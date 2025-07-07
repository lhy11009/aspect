#!/bin/bash -l
#SBATCH -N 15
#SBATCH -n 720
#SBATCH --threads-per-core=1
#SBATCH --tasks-per-node=48
#SBATCH -o task-%j.stdout
#SBATCH -e task-%j.stderr
#SBATCH -t 48:00:00
#SBATCH --partition=skx-normal
#SBATCH --switches=1
#SBATCH -A TG-EES220051
#SBATCH --job-name=eba3d_width51_c22_AR4

>&2 echo "list of modules:"
>&2 module list
>&2 echo "aspect source: ${ASPECT_SOURCE_DIR}"

ibrun ${ASPECT_SOURCE_DIR}/build_master_TwoD/aspect ./case.prm
