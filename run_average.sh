#!/bin/bash
#SBATCH --job-name=avgmeth # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL               # notifications 
#SBATCH --mail-user=splaisie@asu.edu # send-to address

sbatch -n 2 --mem=48G --job-name=MW11 --wrap "python average_methylation.py MW-11"
sbatch -n 2 --mem=48G --job-name=MW21 --wrap "python average_methylation.py MW-21"
sbatch -n 2 --mem=48G --job-name=MW31 --wrap "python average_methylation.py MW-31"
sbatch -n 2 --mem=48G --job-name=OBG55P1 --wrap "python average_methylation.py OBG0055-P1"
