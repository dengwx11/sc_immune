#!/bin/bash

#SBATCH --partition=scavenge,general,bigmem,pi_zhao
#SBATCH --job-name=BAL
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=100000
#SBATCH --mail-type=ALL

#module restore basic-conda
#source activate scrnaseq

module load R/4.0.3-foss-2020b
Rscript --vanilla /gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune/Real_application/BAL/run_transig.R