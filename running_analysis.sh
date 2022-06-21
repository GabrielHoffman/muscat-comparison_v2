
cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison

# clear results
# rm -f plots/* results/* logs/*

ml git
git pull

snakemake -j1 --rerun-incomplete

# https://hoffmg01.u.hpc.mssm.edu/muscat-comparison/


snakemake  --jobs 200 --cluster 'bsub -q premium -R "rusage[mem=6000]" -R span[hosts=1] -W 6

snakemake  --rerun-incomplete --jobs 200 --cluster 'bsub -q premium -R "rusage[mem=6000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 2' 




