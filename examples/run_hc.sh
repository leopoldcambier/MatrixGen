#!/bin/bash
#
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH -J hc
#SBATCH -p mc,owners
#SBATCH -o hc.%j.out

using julia

srun julia high_contrast_diffusion.jl
