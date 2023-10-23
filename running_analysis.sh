

# cd /Users/gabrielhoffman/workspace/repos/eval_methods/muscat-comparison_v2

cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison_v2/
ml git python gcc/11.2.0

# clear results
# rm -f plots/* results/* logs/* meta/* meta/*/* data/sim_data/*
Rscript setup.R


git pull

snakemake -j1 --rerun-incomplete

# https://hoffmg01.u.hpc.mssm.edu/muscat-comparison_v2/


snakemake --rerun-incomplete --jobs 500 --cluster 'bsub -q premium -R "rusage[mem=16000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 5' 


snakemake -j1 -R plot_lfc


# how many reads per cell?
# how many cells?

cd /Users/gabrielhoffman/Dropbox/projects/dreamlet/v1/figures/simulations/v2

rsync -avPz sklar1:"/sc/arion/projects/CommonMind/hoffman/muscat-comparison_v2/plots/*.pdf" .




R --args res="results/kang,de*_ns*" wcs="did=kang,x=s" ggp="plots/kang-perf_by_ns.rds" fig="plots/kang-perf_by_ns.pdf" did="kang" x="s"



# Set alpha = 1e5: 2:13 pm
# Set alpha = 100: 3:07 pm
# Set alpha = 20: 4:30 pm
# Set alpha = 1: 10:00 pm


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



