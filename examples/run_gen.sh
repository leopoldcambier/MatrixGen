#!/bin/bash
#
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH -J gen
#SBATCH -p mc,owners
#SBATCH -o gen.%j.out

using julia

srun julia gen.jl