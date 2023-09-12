suppressMessages({
    library(dplyr)
    library(jsonlite)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment) 
    library(HMP) 
})

# load data & simulation parameters
sce <- readRDS(args$sce)
sim_pars <- fromJSON(args$sim_pars)
set.seed(sim_pars$seed + as.numeric(wcs$i))

assignInNamespace( ".check_args_simData", function(u)
    return(list(nk = u$nk, ns = u$ns)), ns="muscat")

sim <- simData(sce, 
    paired = FALSE, lfc = 0.5,
    ng = nrow(sce), nc = 10*sim_pars$nc,
    ns = sim_pars$ns, nk = sim_pars$nk,
    p_dd = sim_pars$p_dd, probs = sim_pars$probs,
    force=TRUE)

sim <- sim[rowSums(counts(sim) > 0) >= 10, ]

# don't subsample genes
# sim <- sim[sample(nrow(sim), min(nrow(sim), sim_pars$ng)), ]

# Downsample cell to get sim_pars$nc total
# using Dirichlet-multinomial
#############################

alpha = 2
tab = table(sim$sample_id, sim$cluster_id)
countTarget = Dirichlet.multinomial(tab/10, c(alpha, alpha))

df_grid = expand.grid(sid = levels(sim$sample_id), 
                        cid = levels(sim$cluster_id))

keep = lapply( seq(nrow(df_grid)), function(i){

    keep = which( sim$sample_id == df_grid$sid[i] & sim$cluster_id == df_grid$cid[i])

    ncells = countTarget[df_grid$sid[i],df_grid$cid[i]]

    sample(keep, ncells)
    })
keep = sort(unlist(keep))

sim = sim[,keep]

# back to standard processing
gi <- metadata(sim)$gene_info 
gi <- dplyr::filter(gi, gene %in% rownames(sim))
metadata(sim)$gene_info <- gi

sim <- computeLibraryFactors(sim)
sim <- logNormCounts(sim)
assays(sim)$cpm <- calculateCPM(sim)
assays(sim)$vstresiduals <- suppressWarnings(
    vst(counts(sim), show_progress = FALSE)$y)

saveRDS(sim, args$sim)

