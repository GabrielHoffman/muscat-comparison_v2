
# Use PsychAD for sims
#######################
cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison_v2/data/raw_data
R
library(zellkonverter)
library(SingleCellExperiment)
file = "PsychAD_r0_Dec_28_2022.h5ad"
sce = readH5AD(file, use_hdf5 = TRUE)
counts(sce) = assay(sce, 'X')
assay(sce, 'X') = NULL
reducedDims(sce) <- NULL

# only keep controls
CTs = c("Astro", "EN_L2_3_IT")
idx = with(colData(sce), Dx_AD == "Control" & Sex == "Male" & subclass %in% CTs)
sceSub = sce[,idx]
colData(sceSub) = droplevels(colData(sceSub))

tab = with(colData(sceSub), table(subclass, Channel))
cs = colSums(tab)
chs = names(cs)[cs > 500]
chs = grep("-1$", chs, value=TRUE)[1:8]

# Subset based on donors
idx = with(colData(sceSub), Channel %in% chs)
sceSub = sceSub[,idx]
colData(sceSub) = droplevels(colData(sceSub))

table(sceSub$Channel, sceSub$subclass)

sceSub$stim = sceSub$Dx_AD
sceSub$ind = sceSub$SubID
sceSub$cell = sceSub$subclass
sceSub$sample_id <- factor(paste0(sceSub$stim, sceSub$ind)) # construct sample IDs
sceSub <- sceSub[, sceSub$stim == "Control"]
colData(sceSub) = droplevels(colData(sceSub))

counts(sceSub) = as.matrix(counts(sceSub))

rs = rowSums(counts(sceSub))

saveRDS(sceSub[rs > 1000,], file="kang_sce0.rds")
#------------------

# cd /Users/gabrielhoffman/workspace/repos/eval_methods/muscat-comparison_v2

cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison_v2/
ml git python gcc/11.2.0

# clear results
# rm -f plots/* results/* logs/* meta/* meta/*/* data/sim_data/*
rm -f plots/*  data/sim_data/*
Rscript setup.R

git pull


snakemake -R sim_data --jobs 500 --cluster 'bsub -q premium -R "rusage[mem=24000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 5' 



snakemake --rerun-incomplete --jobs 500 --cluster 'bsub -q premium -R "rusage[mem=24000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 5' 



snakemake -j1 --rerun-incomplete

# https://hoffmg01.u.hpc.mssm.edu/muscat-comparison_v2/


snakemake -j1 -R plot_lfc


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






