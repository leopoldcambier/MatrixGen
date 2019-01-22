#!/bin/bash
#
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH -J gen
#SBATCH -p mc
#SBATCH -o gen.%j.out

using julia

srun julia gen.jl
