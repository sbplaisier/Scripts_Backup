#!/bin/bash
#SBATCH --job-name=annotate # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL               # notifications 
#SBATCH --mail-user=splaisie@asu.edu # send-to address
#SBATCH -n 2
#SBATCH --mem=32G
#SBATCH -p phi

source activate /home/splaisie/mambaforge/envs/placenta_RNA_env3/
cd /home/splaisie/scripts/
python annotate_features.py
