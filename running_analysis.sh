

cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison_v2/
ml git python gcc/11.2.0

# clear results
# rm -f plots/* results/* logs/* meta/* meta/*/* data/sim_data/*
Rscript setup.R


git pull

snakemake -j1 --rerun-incomplete

# https://hoffmg01.u.hpc.mssm.edu/muscat-comparison/


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