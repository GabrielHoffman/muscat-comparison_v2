
cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison

# clear results
# rm -f plots/* results/* logs/* meta/sim_pars/* data/sim_data/*
Rscript setup.R


ml git
git pull

snakemake -j1 --rerun-incomplete

# https://hoffmg01.u.hpc.mssm.edu/muscat-comparison/


snakemake --jobs 500 --cluster 'bsub -q premium -R "rusage[mem=16000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 2' 


snakemake --rerun-incomplete --jobs 500 --cluster 'bsub -q premium -R "rusage[mem=16000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 2' 


snakemake -j1 -R plot_lfc


# how many reads per cell?
# how many cells?

cd /Users/gabrielhoffman/Dropbox/projects/dreamlet/v1/figures/simulations

rsync -avPz sklar1:"/sc/arion/projects/CommonMind/hoffman/muscat-comparison/plots/*.pdf" .

