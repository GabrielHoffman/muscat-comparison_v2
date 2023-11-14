suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
    library(dreamlet)
    library(edgeR)
    library(matrixStats)
    library(tidyverse)
    library(Matrix)
})

apply_pb <- function(sce, pars, ds_only = TRUE) {
    t1 <- system.time({
        a <- pars$assay
        if (!ds_only) 
            assay(sce, a) <- switch(a, 
                counts = counts(sce),
                cpm = calculateCPM(counts(sce)),
                logcounts = normalizeCounts(computeLibraryFactors(sce)),
                vstresiduals = vst(counts(sce), show_progress = FALSE)$y)
        pb <- aggregateData(sce, a, fun = pars$fun, scale = pars$scale)
    })[[3]]
    t2 <- system.time({
        if( pars$method %in% c("dreamlet", "dreamlet_no_cell_weights") ){

            useCountsWeights = ifelse(pars$method == "dreamlet", TRUE, FALSE)

            library(dreamlet)
            # use dreamlet pseudobulk command here
            pb <- aggregateToPseudoBulk(sce, a, fun = pars$fun, scale = pars$scale, cluster_id = "cluster_id",sample_id = "sample_id")

            if( useCountsWeights ){
                
                # two part correction
                # 1) counts: scale pseudocount so zero counts -> zero var
                # 2) sigSq: shrink sample var to handle small n
                V.list1 = getVarList( sce, "cluster_id", "sample_id", shrink=TRUE, 0.5)

                # W.list = lapply(V.list1, function(x){
                #     x = 1 / ( x + quantile(x, .2))
                #     x / rowMeans(x)
                #     })
            }else{
                W.list = NULL
            }

            # df = merge(t(W.list[[1]][1:4,]), df_pc[df_pc$cellType == "cluster1",], by.x="row.names", by.y="ID")

            vobj <- processAssays(pb, ~ 1, verbose=FALSE, min.count=3, weightsList = W.list)
            fit <- dreamlet(vobj, ~ group_id, verbose=FALSE )
            tab <- topTable(fit, coef='group_idB', number=Inf, sort.by="none")

            tab2 = with(tab, data.frame(gene = ID, cluster_id = assay, logFC, AveExpr, t, p_val=P.Value, B, contrast='B'))

            tab2$p_adj.glb = p.adjust(tab2$p_val, "BH")
            tab2$p_adj.loc = rep(NA, nrow(tab2))

            for( CT in unique(tab2$cluster_id) ){
                idx = which(tab2$cluster_id==CT)
                tab2$p_adj.loc[idx] = p.adjust(tab2$p_val[idx], "BH")
            }
            res = tab2    
        }else{
            res <- tryCatch(
                do.call(pbDS, c(
                    list(pb = pb, filter = "none", verbose = FALSE),
                    pars[names(pars) %in% names(formals(pbDS))])),
                error = function(e) e)
            if (!inherits(res, "error"))
                res <- dplyr::bind_rows(res$table[[1]])
        }
    })[[3]]
    list(rt = c(t1, t2), tbl = res)
}



getBootLCPM = function(sce, cluster_id, sample_id, df_pc, ndraws = NULL){
    # interate thu donors, cell types and bootstrap reps
    df_grid = expand.grid(cellType = unique(sce[[cluster_id]]),
                        ID =  unique(sce[[sample_id]]))

    # bootstrap indeces
    idx = sapply( seq(nrow(df_grid)), function(i){

        # filter
        idx = which(df_grid$cellType[i] == sce[[cluster_id]] & df_grid$ID[i] == sce[[sample_id]])

        # bootstrap cells
        if( is.null(ndraws) ){
            idx2 = idx[sample.int(length(idx), length(idx), replace=TRUE)]
        }else{          
            idx2 = idx[sample.int(length(idx), min(length(idx), ndraws), replace=TRUE)]
        }
        idx2
        })
    idx = sort(unlist(idx))

    # pseudobulk of boostrap
    pb <- aggregateToPseudoBulk(sce[,idx],
      assay = "counts",
      cluster_id = cluster_id,
      sample_id = sample_id,
      verbose = FALSE)

    geneExpr = lapply( assayNames(pb), function(CT){

        df_sub = df_pc %>% filter(cellType == CT)
        pc = df_sub$prior.count
        names(pc) = df_sub$ID
            
        countMatrix = assay(pb, CT)
        pcMat = lapply(colnames(countMatrix), function(id)
                    rpois(nrow(countMatrix), pc[id]))
        pcMat = do.call(cbind, pcMat)
        colnames(pcMat) = colnames(countMatrix)

        dge = DGEList(counts = countMatrix + pcMat, 
                            lib.size = colSums2(countMatrix))
        # dge = calcNormFactors(dge)
        edgeR::cpm(dge, log=TRUE, prior.count=.25)
        })
    names(geneExpr) = assayNames(pb)

    geneExpr
}

summarizeBootstraps = function(geneExprBoot){
    # interate thu donors, cell types and bootstrap reps
    CT.names = names(geneExprBoot[[1]])
    id.names = colnames(geneExprBoot[[1]][[1]])

    df_var = lapply( CT.names, function(CT){

        df_var = lapply(id.names, function(id){

            # create matrix of boostrap samples for cell type and id
            Y = lapply( seq(length(geneExprBoot)), function(j){
                geneExprBoot[[j]][[CT]][,id,drop=FALSE]
            })
            Y = do.call(cbind, Y)

            # y.mean = rowMeans2(Y, useNames=FALSE)

            # sampling variance of mean from boostraps
            y.var = rowVars(Y, useNames=TRUE) / ncol(Y)

            y.var = data.frame(var = y.var)
            colnames(y.var) = id

            y.var
        })
        as.matrix(do.call(cbind, df_var))
    })
    names(df_var) = CT.names
    df_var
}


trimWeightOutliersGene = function(x, zmax){

    # compute z-score
    zscore = scale(x)

    # extract parameters of transform
    # z-score = (x - mu) / s
    mu = attr(zscore,"scaled:center")
    s = attr(zscore,"scaled:scale")
    
    # if x exceeds original value giving z-score of zmax, 
    # replace with that orginal value
    x[x > zmax * s + mu] = zmax * s + mu

    # normalize values to have a mean of 1
    x / mean(x)
}


trimWeightOutliers = function(W, zmax = 5){

    t(apply(W, 1, trimWeightOutliersGene, zmax = zmax))
}




#' @export
getVarFromCounts = function(countMatrix, prior.count = .25){

    count.gene <- rowSums2(countMatrix, useNames=FALSE)
    count.lib <- colSums2(countMatrix, useNames=FALSE)
    ncell <- ncol(countMatrix)
    sclSq <- sum(count.lib^2)

    # add pseudocount for each cell
    count.gene <- count.gene + prior.count
    count.lib <- count.lib + 1

    # normalize counts by library size
    # add pseudocount to counts here
    normCounts <- scale(countMatrix + prior.count, 
                    scale = count.lib, 
                    cente = FALSE)
    # compute variance for each row
    sigmaSq.hat.gene <- rowVars(normCounts, useNames=FALSE)
    sigmaSq.hat.gene[is.na(sigmaSq.hat.gene)] <- 0

    # compute variance
    # vectorize
    # v.hat <- 1 / count.gene + (sigmaSq.hat.gene * sclSq) / (ncell^2 *count.gene^2)

    tibble(Gene = rownames(countMatrix), count.gene, sigmaSq.hat.gene, sclSq, ncell)
}

#' @export
getVarForCellType = function(sce, cluster_id, sample_id, CT, prior.count){

    idx = which(sce[[cluster_id]] == CT)
    lib.size <- colSums2(counts(sce), cols=idx)
    
    # scale prior count so that an observed count of 0
    # gives zero variance across samples
    df_pc = data.frame(ID = sce[[sample_id]][idx], 
        cellType = sce[[cluster_id]][idx], 
        prior.count = prior.count * lib.size/mean(lib.size)) %>%
        group_by(cellType, ID) %>%
        summarize(n=length(ID), prior.count = sum(prior.count) / n)

    # get variance estimates for each ID and gene
    df <- lapply( unique(sce[[sample_id]]), function(ID){

        idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == ID
        countMatrix = counts(sce)[,idx,drop=FALSE]

        pc = df_pc$prior.count[df_pc$ID == ID]
        
        res <- getVarFromCounts( countMatrix, pc)
        res$ID <- ID
        res
        })
    bind_rows(df)
}


# df = bind_rows(df)
# df[df$Gene == gene,] %>% arrange(ID)
# df_pc

# idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == 'sample1.A'
# countMatrix1 = counts(sce)[,idx,drop=FALSE]

# idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == 'sample15.A'
# countMatrix2 = counts(sce)[,idx,drop=FALSE]


# countMatrix1[gene,]
# countMatrix2[gene,]



# hist(df[df$Gene == gene,]$count.gene)

#' @export
getVarList = function(sce, cluster_id, sample_id, shrink, prior.count){

    if( ! cluster_id %in% colnames(colData(sce)) ){
        msg <- paste0("sample_id entry not found in colData(sce): ", cluster_id)
        stop( msg )
    }
    if( ! sample_id %in% colnames(colData(sce)) ){
        msg <- paste0("sample_id entry not found in colData(sce): ", sample_id)
        stop( msg )
    }

    var.list <- lapply( unique(sce[[cluster_id]]), function(CT){
        df <- getVarForCellType( sce, cluster_id, sample_id, CT, prior.count) %>%
                mutate(Gene = factor(Gene, rownames(sce)),
                        ID = factor(ID))

        if( shrink ){           
            res <- limma::squeezeVar( df$sigmaSq.hat.gene, df$ncell-1, robust=FALSE)
            # plot(df$sigmaSq.hat.gene, res$var.post, main=CT, log="xy")
            # abline(0, 1, col="red")
            # browser()
            df$sigmaSq.hat.gene <- res$var.post
        }
        df$vhat <- with(df, 1 / count.gene + (sigmaSq.hat.gene * sclSq) / (ncell^2 *count.gene^2))
        # only count variance
        # df$vhat <- with(df, 1 / count.gene)

        mat <- sparseMatrix(df$Gene, df$ID, 
            x = df$vhat, 
            dims = c(nlevels(df$Gene), nlevels(df$ID)),
            dimnames = list(levels(df$Gene), levels(df$ID)))
        as.matrix(mat)
    })
    names(var.list) <- unique(sce[[cluster_id]])
    var.list
}



