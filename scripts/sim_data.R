suppressMessages({
    library(dplyr)
    library(jsonlite)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment) 
})

# load data & simulation parameters
sce <- readRDS(args$sce)
sim_pars <- fromJSON(args$sim_pars)
set.seed(sim_pars$seed + as.numeric(wcs$i))

assignInNamespace( ".check_args_simData", function(u)
    return(list(nk = u$nk, ns = u$ns)), ns="muscat")

# Simulate more cells than needed
# Then downsample later
k_scaling = 10

sim <- simData(sce, 
    paired = FALSE, lfc = 0.5,
    ng = nrow(sce), nc = k_scaling*sim_pars$nc,
    ns = sim_pars$ns, nk = sim_pars$nk,
    p_dd = sim_pars$p_dd, probs = sim_pars$probs,
    force=TRUE)

sim <- sim[rowSums(counts(sim) > 0) >= 10, ]

# don't subsample genes
tab = table(sim$sample_id, sim$cluster_id)

# Downsample cell to get sim_pars$nc total
# using Dirichlet-multinomial
#############################

rdmn = function(counts, alpha){

    stopifnot(identical(ncol(counts), length(alpha)))

    df = lapply(seq(nrow(counts)), function(i){
        prob = dirmult::rdirichlet(1, alpha)
        t(rmultinom(1, counts[i,], prob = prob))
    })
    df = do.call(rbind, df)
    rownames(df) = rownames(counts)
    colnames(df) = colnames(counts)
    df
}

if( k_scaling > 1){
    # overdispersion parameter 
    alpha = 10

    countTarget = rdmn(tab/k_scaling*2, rep(alpha, ncol(tab)))

    df_grid = expand.grid(sid = levels(sim$sample_id), 
                            cid = levels(sim$cluster_id))

    keep = lapply( seq(nrow(df_grid)), function(i){

        keep = which( sim$sample_id == df_grid$sid[i] & sim$cluster_id == df_grid$cid[i])

        ncells = countTarget[df_grid$sid[i],df_grid$cid[i]]
        ncells = max(5, ncells)

        if( ncells < length(keep)){
            keep <- sample(keep, ncells)
        }
        keep
    })
    keep = sort(unlist(keep))

    # Subsample
    sim2 = sim[,keep]

    # rename cells
    colnames(sim2) = paste0("cell", seq(ncol(sim2)))

    # filter genes
    sim2 <- sim2[rowSums(counts(sim2) > 0) >= 10, ]

    # set number of cells
    metadata(sim2)$n_cells = table(sim2$sample_id)
    metadata(sim)$args$nc = sim_pars$nc

}else{
    sim2 = sim
}

# Subsample genes 
sim2 <- sim2[sample(nrow(sim2), min(nrow(sim2), sim_pars$ng)), ]

# back to standard processing
gi <- metadata(sim2)$gene_info 
gi <- dplyr::filter(gi, gene %in% rownames(sim2))
metadata(sim2)$gene_info <- gi

sim2 <- computeLibraryFactors(sim2)
sim2 <- logNormCounts(sim2)
assays(sim2)$cpm <- calculateCPM(sim2)

vst_values <- suppressWarnings(sctransform::vst(counts(sim2))$y)

assays(sim2)$vstresiduals <- vst_values

saveRDS(sim2, args$sim)
