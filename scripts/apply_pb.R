suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
    library(dreamlet)
    library(edgeR)
    library(matrixStats)
    library(dreamlet)
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
        if( pars$method %in% c("dreamlet_deltaW", "dreamlet_delta", "dreamlet_ncells", "dreamlet_none") ){

            # use dreamlet pseudobulk command here
            pb <- aggregateToPseudoBulk(sce, a, cluster_id = "cluster_id", sample_id = "sample_id", fun = pars$fun, scale = pars$scale)

            # Gene expressed genes for each cell type
            geneList = getExprGeneNames(pb, min.cells=1)

            # Precision weights
            pc = .5
            W.list = NULL
            W.list <- switch(pars$method, 
                    "dreamlet_delta" = pbWeights( sce, 
                                    sample_id = "sample_id", 
                                    cluster_id = "cluster_id", 
                                    method = "delta", 
                                    geneList = geneList,
                                    details = TRUE,
                                    prior.count = pc), 
                    "dreamlet_ncells" = pbWeights( sce, 
                                    sample_id = "sample_id", 
                                    cluster_id = "cluster_id", 
                                    method = "ncells", 
                                    geneList = geneList), 
                    "dreamlet_none" = {w = pbWeights( sce, 
                                    sample_id = "sample_id", 
                                    cluster_id = "cluster_id", 
                                    method = "ncells", 
                                    geneList = geneList);
                        lapply(w, function(x){x[] = 1; x})})

            priorWeightsAsCounts = ifelse(pars$method == "dreamlet_deltaW", TRUE, FALSE)
            rescaleWeightsAfter = ifelse(pars$method == "dreamlet_deltaW", FALSE, TRUE)

            vobj <- processAssays(pb, ~ group_id, verbose=FALSE, 
                    weightsList = W.list, 
                    prior.count = pc, 
                    priorWeightsAsCounts = priorWeightsAsCounts, 
                    rescaleWeightsAfter = rescaleWeightsAfter,
                    min.cells = 2,
                    min.count = 2,
                    min.samples = 4,
                    min.prop = 0.2)
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

            # In order to keep the same genes for muscat as dreamlet
            # get gene/cluster pairs that are retained
            pb.tmp <- dreamlet::aggregateToPseudoBulk(sce, "counts", cluster_id = "cluster_id",sample_id = "sample_id")
            vobj <- dreamlet::processAssays(pb.tmp, ~ group_id, verbose=FALSE,
                    min.cells = 2,
                    min.count = 2,
                    min.samples = 4,
                    min.prop = 0.2)
            fit <- dreamlet(vobj, ~ group_id, verbose=FALSE )
            tab <- topTable(fit, coef='group_idB', number=Inf, sort.by="none")
            tab$key = with(tab, paste(assay, ID))

            res <- tryCatch(
                do.call(pbDS, c(
                    list(pb = pb, filter = "none", verbose = TRUE, min_cells=1),
                    pars[names(pars) %in% names(formals(pbDS))])),
                error = function(e) e)

            if (!inherits(res, "error")){
                 res <- dplyr::bind_rows(res$table[[1]])
            }

            # retain only gene/cluster pairs from dreamlet
            keep = with(res, paste(cluster_id, gene)) %in% tab$key
            res = res[keep,]
            res
        }
    })[[3]]
    list(rt = c(t1, t2), tbl = res)
}



