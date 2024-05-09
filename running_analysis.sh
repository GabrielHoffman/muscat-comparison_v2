
# Use PsychAD for sims
#######################
cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison_v2/data/raw_data/pad
R
library(zellkonverter)
library(SingleCellExperiment)
file = "PsychAD_r0_Dec_28_2022.h5ad"
sce = readH5AD(file, use_hdf5 = TRUE)
counts(sce) = assay(sce, 'X')
assay(sce, 'X') = NULL
reducedDims(sce) <- NULL

# only keep male controls in the same batch
CTs = c("IN_VIP", "Astro" )
CTs = c("IN_PVALB_CHC", "EN_L5_6_NP" )
# sort(table(sce$subclass))
idx = with(colData(sce), Dx_AD == "Control" & Sex == "Male" & subclass %in% CTs & poolID %in% c("NPSAD-169-A1", "NPSAD-243-A2", "NPSAD-20201106-C1", "NPSAD-169-A1", "NPSAD-243-A1", "NPSAD-217-C2", "NPSAD-228-A1", "NPSAD-119-A1", "NPSAD-141-A2", "NPSAD-235-A2"))
sceSub = sce[,idx]
colData(sceSub) = droplevels(colData(sceSub))
sceSub$lib.size = colSums2(counts(sceSub))

# Test pseudobulk and voom
library(dreamlet)
idx = sample(ncol(sceSub), 1000)
idx = which(sceSub$lib.size < 4000)
idx = 1:ncol(sceSub)
pb = aggregateToPseudoBulk(sceSub[,idx], sample_id="SubID", cluster_id = "subclass")
vobj = processAssays(pb, ~1, 
                    verbose=TRUE, 
                    priorWeightsAsCounts = TRUE, 
                    prior.count.for.weights = .5,
                    rescaleWeightsAfter = FALSE,
                    min.cells = 0,
                    min.count = 0,
                    min.samples = 4,
                    min.prop = 0,
                    min.total.count = 1)
plotVoom(vobj)






tab = with(colData(sceSub), table(subclass, Channel))
cs = colSums(tab)
chs = names(cs)
# chs = names(cs)[cs > 500]
# chs = grep("-1$", chs, value=TRUE)[1:8]

# Subset based on donors
idx = with(colData(sceSub), Channel %in% chs)
sceSub = sceSub[,idx]
colData(sceSub) = droplevels(colData(sceSub))

table(sceSub$SubID, sceSub$subclass)

sceSub$stim = sceSub$Dx_AD
sceSub$ind = sceSub$SubID
sceSub$cell = sceSub$subclass
sceSub$sample_id <- factor(paste0(sceSub$stim, sceSub$ind)) # construct sample IDs
sceSub <- sceSub[, sceSub$stim == "Control"]
colData(sceSub) = droplevels(colData(sceSub))

counts(sceSub) = as.matrix(counts(sceSub))

rs = rowSums(counts(sceSub))
# [rs > 5,]

saveRDS(sceSub, file="kang_sce0.rds")


cd /sc/arion/scratch/hoffmg01/muscat-comparison_v2
Rscript scripts/prep_kang.R --input_sce data/raw_data/sce0_kang.rds --output_sce data/raw_data/ref_kang.rds




#------------------

# cd /Users/gabrielhoffman/workspace/repos/eval_methods/muscat-comparison_v2

# https://hoffmg01.u.hpc.mssm.edu/muscat-comparison_v3/plots/

# cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison_v2/
cd /sc/arion/scratch/hoffmg01/muscat-comparison_v2
ml python gcc/11.2.0
git pull

# clear results
# rm -f plots/* results/* logs/* meta/* meta/*/* data/sim_data/*
rm -f plots/*  data/sim_data/* 
Rscript setup.R

git pull #origin new

snakemake --rerun-incomplete --jobs 500 --cluster 'bsub -q premium -R "rusage[mem=16000]" -R span[hosts=1] -W 16:00 -P acc_CommonMind -n 8' 


# snakemake -R sim_qc --jobs 500 --cluster 'bsub -q premium -R "rusage[mem=24000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 5' 



pars = list(assay = "counts", fun = "sum", scale = "false", method = "limma-voom", treat=FALSE)


pars = list(assay = "counts", fun = "sum", scale = FALSE, method = "DESeq2", treat=FALSE)


snakemake -j12 --rerun-incomplete

# https://hoffmg01.u.hpc.mssm.edu/muscat-comparison_v2/


snakemake -j1 -R plot_perf_by_nx
snakemake -j1 -R sim_qc



## Feb 5 2024
###############

cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison/
ml git python gcc/11.2.0
git pull origin new

rm -f plots/* results/* logs/* meta/sim_pars/*
mkdir -p results
Rscript setup.R

snakemake --rerun-incomplete --jobs 500 --cluster 'bsub -q premium -R "rusage[mem=16000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 5' 




args = list(sce = "data/raw_data/ref_kang.rds", sim_pars = "meta/sim_pars/db10.json")


sim <- simData(sce, 
    paired = FALSE, lfc = .5 ,
    force = TRUE,
    ng = nrow(sce), nc = sim_pars$nc * k_scaling,
    ns = 10, nk = 1,
    p_dd = sim_pars$p_dd, probs = sim_pars$probs)

sim$SubID = gsub("\\..*", "", sim$sample_id)


pb <- aggregateToPseudoBulk(sim, "counts", cluster_id = "cluster_id", sample_id = "sample_id")


snakemake --rerun-incomplete -j1


# how many reads per cell?
# how many cells?

cd /Users/gabrielhoffman/Dropbox/projects/dreamlet/v1/figures/simulations/v2

rsync -avPz sklar1:"/sc/arion/projects/CommonMind/hoffman/muscat-comparison_v2/plots/*.pdf" .




R --args res="results/kang,de*_ns*" wcs="did=kang,x=s" ggp="plots/kang-perf_by_ns.rds" fig="plots/kang-perf_by_ns.pdf" did="kang" x="s"


data/raw_data/sce_kang.rds


library(zellkonverter)

file = "/sc/arion/projects/CommonMind/aging/hui/files/aging_lister_combined_ds1000_rmno_beforeUMAT_full_adata.h5ad"

sce = readH5AD(file, verbose=TRUE, use_hdf5=TRUE)

Error: numpy.core._exceptions._ArrayMemoryError: Unable to allocate 30.5 GiB for an array with shape (4097142408,) and data type float64

files = paste0("results/", list.files("results"))
files = files[grep("s(\\d+).rds", files)]

args = list(res = files)

wcs = list(x = 's')



R --args sim=data/sim_data/kang,de10_ns,7.rds fun=scripts/apply_pb.R wcs=c=x,did=kang,g=x,i=7,j=1,k=x,mid=dreamlet.sum.counts,s=5,sid=de10_ns            meth_pars=meta/meth_pars/dreamlet.sum.counts.json run_pars=meta/run_pars/kang,de10_ns.json res=results/kang,de10_ns,7,dreamlet.sum.counts,1,gx,cx,kx,s5.rds

sce = readRDS("data/sim_data/kang,de10_ns,7.rds")

pb <- aggregateToPseudoBulk(sce, "counts", cluster_id = "cluster_id",sample_id = "sample_id")

vobj1 <- processAssays(pb, ~ 1, verbose=FALSE, min.count=3)
fit1 <- dreamlet(vobj1, ~ group_id, verbose=FALSE )
tab1 <- topTable(fit1, coef='group_idB', number=Inf, sort.by="none")


vobj2 <- processAssays(pb, ~ 1, verbose=FALSE, min.count=3, weightsList = W.list)
fit2 <- dreamlet(vobj2, ~ group_id, verbose=FALSE )
tab2 <- topTable(fit2, coef='group_idB', number=Inf, sort.by="none")

plot(tab1$logFC, tab2$logFC)
abline(0, 1, col="red")


plot(tab1$t, tab2$t)
abline(0, 1, col="red")

hist(tab1$P.Value, 100)


hist(tab2$P.Value, 100)

a = sapply(seq(nrow(assay(vobj1, 1))), function(i){
	cor(assay(vobj1, 1)$weights[i,], assay(vobj2, 1)$weights[i,])
	})




i = which.min(a)

plot(assay(vobj1, 1)$weights[i,], assay(vobj2, 1)$weights[i,])

tab1[i,]
tab2[i,]

R CMD BATCH "--args sim=data/sim_data/kang,de10_nc,4.rds fun=scripts/apply_pb.R wcs=c=25,did=kang,g=x,i=4,j=1,k=x,mid=dreamlet.sum.counts,s=x,sid=de10_nc  meth_pars=meta/meth_pars/dreamlet.sum.counts.json run_pars=meta/run_pars/kang,de10_nc.json res=results/kang,de10_nc,4,dreamlet.sum.counts,1,gx,c25,kx,sx.rds" scripts/run_meth.R logs/run_meth-kang,de10_nc,4,dreamlet.sum.counts,1,gx,c25,kx,sx.Rout


tabSub = tab[tab$assay=="cluster1",]
v = apply(W.list[[1]][tabSub$ID,], 1, function(x) max(x))

plot(v, -log10(tabSub$P.Value))



rs = rowSums(
df[df$Gene == gene,] %>% arrange(ID))

which.min(rs)
gene3435 
      97 
> which.min(rs)
gene = "gene3435"
W.list[[1]][gene,]

V.list1 = getVarList( sce, "cluster_id", "sample_id", shrink=TRUE, 0.01)

V.list1[[1]][gene,]




gene = "gene2341"
idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == 'sample1.A'
countMatrix1 = counts(sce)[,idx,drop=FALSE]

idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == 'sample25.A'
countMatrix2 = counts(sce)[,idx,drop=FALSE]

countMatrix1[gene,]
countMatrix2[gene,]



df_pc %>%
	filter(ID %in% c('sample1.A', 'sample25.A'))


df %>%
	filter(ID %in% c('sample1.A', 'sample25.A')) %>%
	filter(Gene == gene)


V.list1 = getVarList( sce, "cluster_id", "sample_id", shrink=TRUE,5)

W.list = lapply(V.list1, function(x){
    x = 1 / ( x + quantile(x, 0.2))
    x / rowMeans(x)
    })

vobj <- processAssays(pb, ~ 1, verbose=FALSE, min.count=3, weightsList = W.list)
fit <- dreamlet(vobj, ~ group_id, verbose=FALSE )
tab <- topTable(fit, coef='group_id', number=Inf, sort.by="none")
hist(tab$P.Value)


hist(apply(W.list[[1]], 1, sd))


hist(apply(W.list[[1]], 1, max))



7fd3e8ec7fb9b92ea0fa2dcb2f6d343c936fe545


 res <- tryCatch(
                do.call(pbDS, c(
                    list(pb = pb, filter = "none", verbose = TRUE, min_cells=10),
                    pars[names(pars) %in% names(formals(pbDS))])),
                error = function(e) e)

debug(muscat:::.DESeq2)

res = pbDS(pb, method = "limma-voom")

res = pbDS(pb, method = "DESeq2")


Error in rule run_meth:
    jobid: 440
    output: results/kang,dm10,6,DESeq2.sum.counts,1,gx,cx,kx,sx.rds
    log: logs/run_meth-kang,dm10,6,DESeq2.sum.counts,1,gx,cx,kx,sx.Rout (check log file(s) for error message)
    shell:
        /hpc/packages/minerva-centos7/R/4.3.0/lib64/R/bin/R CMD BATCH --no-restore --no-save            "--args sim=data/sim_data/kang,dm10,6.rds fun=scripts/apply_pb.R wcs=c=x,did=kang,g=x,i=6,j=1,k=x,mid=DESeq2.sum.counts,s=x,sid=dm10   meth_pars=meta/meth_pars/DESeq2.sum.counts.json run_pars=meta/run_pars/kang,dm10.json res=results/kang,dm10,6,DESeq2.sum.counts,1,gx,cx,kx,sx.rds"              scripts/run_meth.R logs/run_meth-kang,dm10,6,DESeq2.sum.counts,1,gx,cx,kx,sx.Rout
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)



 sim_pars = list(nc = 1000, ns = 60, nk=1)

library(muscat)
data(example_sce)
library(SingleCellExperiment)

# prep. SCE for simulation
ref <- prepSim(example_sce)

# simulate data
sim <- simData(ref, nc = 200,
	p_dd = c(0.9, 0, 0.1, 0, 0, 0),
	ng = 100, force = TRUE,
	probs = list(NULL, NULL, c(1, 0)))






