#!/bin/bash
#
#SBATCH --job-name=pellet_combustion
#SBATCH --output=out/gpu1.out
#SBATCH --error=err/gpu1.err
#
#SBATCH --ntasks=16
#SBATCH --gres=gpu:2
#SBATCH --partition=gpu
#SBATCH -v

scripts/multi-gpu-execute.sh bin/gpu/Pellet-Flame-Propagation 0.68 0.62 0.52