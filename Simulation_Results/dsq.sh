partition=scavenge,pi_zhao,general;
module load dSQ
#dSQ --jobfile /gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/corrupt_dsq_jobs.sh -p ${partition} -n 1 -C avx2 --mem-per-cpu=64g -t 24:00:00 --mail-type=ALL --batch-file corrupt.pbs
#sbatch corrupt.pbs


dSQ --jobfile /gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/tissue_number_dsq_jobs.sh -p ${partition} -n 1 -C avx2 --mem-per-cpu=64g -t 24:00:00 --mail-type=ALL --batch-file tissue_number.pbs
sbatch tissue_number.pbs

dSQ --jobfile /gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/gene_number_dsq_jobs.sh -p ${partition} -n 1 -C avx2 --mem-per-cpu=64g -t 24:00:00 --mail-type=ALL --batch-file gene_number.pbs
sbatch gene_number.pbs

#dSQ --jobfile /gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/pi_dsq_jobs.sh -p ${partition} -n 1 -C avx2 --mem-per-cpu=64g -t 24:00:00 --mail-type=ALL --batch-file pi.pbs
#sbatch pi.pbs

#dSQ --jobfile /gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/tau_v_dsq_jobs.sh -p ${partition} -n 1 -C avx2 --mem-per-cpu=64g -t 24:00:00 --mail-type=ALL --batch-file tau_v.pbs
#sbatch tau_v.pbs
