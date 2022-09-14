#!/bin/bash
#SBATCH --job-name=annotate # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL               # notifications 
#SBATCH --mail-user=splaisie@asu.edu # send-to address
#SBATCH -n 2
#SBATCH --mem=48G

cd /home/splaisie/scripts/
python check_xist_methylation.py
