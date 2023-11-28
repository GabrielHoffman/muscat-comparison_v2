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

# simulate effect size heterogeneity
# with a normal offset added to the logFC 
.nb.replace <- function(cs, d, m, lfc = NULL, f = 1) {
    n_gs <- length(d)
    n_cs <- length(cs)
    if (is.null(lfc)) {
        lfc <- rep(0, n_gs)
    } else {
        lfc[lfc < 0] <- 0
    }

    # effect size heterogeneity
    lfc <- lfc + rnorm(length(lfc), 0, 1) 

    fc <- f * (2 ^ lfc)
    fc <- rep(fc, each = n_cs)
    ds <- rep(1/d, each = n_cs)
    ms <- c(t(m[, cs])) * fc 
    y <- rnbinom(n_gs * n_cs, size = ds, mu = ms)
    y <- matrix(y, byrow = TRUE, 
        nrow = n_gs, ncol = n_cs, 
        dimnames = list(names(d), cs))
    ms <- split(ms, rep(seq_len(nrow(m)), each = n_cs))
    list(counts = y, means = ms)
}

assignInNamespace(".nb", .nb.replace, ns="muscat")


# Simulate more cells than needed
# Then downsample later
k_scaling = 2

sim <- simData(sce, 
    paired = FALSE, lfc = 0.5,
    ng = nrow(sce), nc = k_scaling*sim_pars$nc,
    ns = sim_pars$ns, nk = sim_pars$nk,
    p_dd = sim_pars$p_dd, probs = sim_pars$probs,
    force=TRUE)

# don't subsample genes
tab = table(sim$sample_id, sim$cluster_id)

sim <- sim[rowSums(counts(sim) > 0) >= 10, ]

if( k_scaling > 1){
    df_grid = expand.grid(sid = levels(sim$sample_id), 
                            cid = levels(sim$cluster_id))

    keep = lapply( seq(nrow(df_grid)), function(i){

        keep = which( sim$sample_id == df_grid$sid[i] & sim$cluster_id == df_grid$cid[i])

        target = tab[df_grid$sid[i],df_grid$cid[i]]/k_scaling

        # sample cell counts from Negative Binomial 
        # Poisson if theta = Inf
        # additive overdispersion is mu^2/theta
        # variance is 'a' times the Poisson variance 
        # a = 10
        # theta = target / (a-1)
        theta = 1
        ncells = MASS::rnegbin(1, mu=target, theta=theta)
         
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
    metadata(sim2)$args$nc = sim_pars$nc

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




# # test sampling
# with(colData(sim2), xtabs(~sample_id + cluster_id))


# library(dreamlet)
# # use dreamlet pseudobulk command here
# pb <- aggregateToPseudoBulk(sim2, cluster_id = "cluster_id",sample_id = "sample_id")

# vobj <- processAssays(pb, ~ 1, verbose=FALSE, min.count=3, useCountsWeights=TRUE)
# # plotVoom(vobj)
# fit <- dreamlet(vobj, ~ group_id, verbose=FALSE )

# topTable(fit, coef='group_idB', number=Inf, sort.by="none") %>%
#         as_tibble %>%
#         mutate(cluster_id = assay, gene = ID) %>%
#         left_join(metadata(sim2)$gene_info, by=c("cluster_id", "gene")) %>%
#         filter(category == "ee") %>%
#         ggplot(aes(P.Value)) +
#             geom_histogram() + 
#             theme_classic() +
#             theme(aspect.ratio=1) +
#             facet_wrap(~cluster_id)


# pb = aggregateToPseudoBulk(sce[,sce$subclass=="Astro"],
#     assay = "counts", 
#     cluster_id = "subclass",
#     sample_id = "SubID")

# pb$value = rnorm(ncol(pb))

# res.proc = processAssays( pb, ~ 1, assays="Astro", useCountsWeights=TRUE)
# # plotVoom(res.proc)

# fit = dreamlet(res.proc, ~ value)

# hist(topTable(fit, coef="value", number=Inf)$P.Value)



# assignInNamespace('processOneAssay', value = f, ns="dreamlet")
# checkFormula = dreamlet:::checkFormula
# library(variancePartition)
# library(edgeR)



