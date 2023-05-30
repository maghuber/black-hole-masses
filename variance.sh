#!/bin/bash
#SBATCH --time=0:90:00			
#SBATCH --mail-type=all
#SBATCH --mail-user=mahu8801@colorado.edu
#SBATCH --qos=normal	      		
#SBATCH --partition=amilan		
#SBATCH --output=~/scratch/alpine/$USER/output/python_%j.out

# Load the python module
module load python
module load anaconda
conda activate blackhole_env

# Run Python Script
cd ~/scratch/alpine/$USER/blackhole_data/
python errorbars.py --ngalaxies=1000 --ndraws=100 --i='~/scratch/alpine/$USER/blackhole_data/errors.dat' --o='~/scratch/alpine/$USER/blackhole_data/bhmasses_variances.dat'