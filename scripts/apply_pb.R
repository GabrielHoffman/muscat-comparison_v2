suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
    library(dreamlet)
    library(edgeR)
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

            if( useCountsWeights ){
                # W.list = getWeightsList(sce, "cluster_id", "sample_id", 10)
                # W.list = lapply( W.list, trimWeightOutliers, zmax=3)

                # Perform Bootstraps
                geneExprBoot = lapply(seq(100), function(i) 
                                    getBootLCPM(sce))

                # Summarize Bootstraps
                bootVars = summarizeBootstraps( geneExprBoot )

                W.list = lapply( names(bootVars), function(x){
                    W = 1 / bootVars[[x]]
                    W / rowMeans(W)
                    })
                names(W.list) = names(bootVars)

            }else{
                W.list = NULL
            }

            library(dreamlet)
            # use dreamlet pseudobulk command here
            pb <- aggregateToPseudoBulk(sce, a, fun = pars$fun, scale = pars$scale, cluster_id = "cluster_id",sample_id = "sample_id")

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



getBootLCPM = function(sce){
    # interate thu donors, cell types and bootstrap reps
    df_grid = expand.grid(cellType = unique(sce$cluster_id),
                        ID =  unique(sce$sample_id))

    # bootstrap indeces
    idx = sapply( seq(nrow(df_grid)), function(i){

        # filter
        idx = which(df_grid$cellType[i] == sce$cluster_id & df_grid$ID[i] == sce$sample_id)

        # bootstrap sample
        idx[sample.int(length(idx), length(idx), replace=TRUE)]
        })
    idx = sort(unlist(idx))

    # pseudobulk of boostrap
    pb <- aggregateToPseudoBulk(sce[,idx],
      assay = "counts",
      cluster_id = "cluster_id",
      sample_id = "sample_id",
      verbose = FALSE)

    geneExpr = lapply( assayNames(pb), function(CT){

        dge = DGEList(counts = assay(pb, CT))
        dge = calcNormFactors(dge)
        edgeR::cpm(dge, log=TRUE)
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

            # variance from boostraps
            y.var = rowVars(Y, useNames=TRUE)

            y.var = data.frame(var = y.var)
            colnames(y.var) = id

            y.var
        })
        as.matrix(do.call(cbind, df_var))
    })
    names(df_var) = CT.names
    df_var
}
