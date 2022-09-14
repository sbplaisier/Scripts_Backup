#!/bin/bash
#SBATCH --job-name=avgmeth # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL               # notifications 
#SBATCH --mail-user=splaisie@asu.edu # send-to address
#SBATCH -n 2
#SBATCH --mem=28G

cd /home/splaisie/scripts/
python average_methylation_range_gene2.py 
