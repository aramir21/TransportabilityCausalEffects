#!/bin/bash

#SBATCH -J STAR_LOSO_metrics
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -p longjobs
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -t 14-00:00:00
#SBATCH --mail-user=aramir21@eafit.edu.co
#SBATCH --mail-type=ALL

module load miniconda3/25.11.1
conda activate ambienteR

cd /home/aramir21/SimulationBHM

Rscript STAR_superpopulation_LOSO_metrics.R
