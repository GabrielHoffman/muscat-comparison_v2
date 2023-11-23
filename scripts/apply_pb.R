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
                V.list1 = dreamlet:::getVarList( sce, sample_id = "sample_id", cluster_id = "cluster_id", shrink=TRUE, 0.5)

                W.list = lapply(V.list1, function(x){
                    x = 1 / ( x + quantile(x, .1))
                    x / rowMeans(x)
                    })
            }else{
                W.list = NULL
            }

            vobj <- processAssays(pb, ~ 1, verbose=FALSE, weightsList = W.list)#, min.cells=1, min.prop=.1, min.count=1)
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

            # get gene/cluster pairs that are retained
            library(dreamlet)
            pb.tmp <- dreamlet::aggregateToPseudoBulk(sce, "counts", cluster_id = "cluster_id",sample_id = "sample_id")
            vobj <- dreamlet::processAssays(pb.tmp, ~ 1, verbose=FALSE)
            fit <- dreamlet(vobj, ~ group_id, verbose=FALSE )
            tab <- topTable(fit, coef='group_idB', number=Inf, sort.by="none")
            tab$key = with(tab, paste(assay, ID))

            res <- tryCatch(
                do.call(pbDS, c(
                    list(pb = pb, filter = "none", verbose = FALSE),
                    pars[names(pars) %in% names(formals(pbDS))])),
                error = function(e) e)

            if (!inherits(res, "error")){
                 res <- dplyr::bind_rows(res$table[[1]])
            }

            # retain only gene/cluster pairs from dreamlet
            keep = with(res, paste(cluster_id, gene)) %in% tab$key
            res = res[keep,]
        }
    })[[3]]
    list(rt = c(t1, t2), tbl = res)
}




