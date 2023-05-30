#!/bin/bash
#SBATCH --time=0:90:00			
#SBATCH --mail-type=all
#SBATCH --mail-user=mahu8801@colorado.edu
#SBATCH --qos=normal	      		
#SBATCH --partition=amilan		
#SBATCH --output=/scratch/alpine/mahu8801/output/python_%j.out

# Load the python module
module load python
module load anaconda
conda activate blackhole_env

# Run Python Script
cd /scratch/alpine/mahu8801/blackhole_data/scripts
python errorbars.py --ngalaxies=1000 --ndraws=100 --i='/scratch/alpine/mahu8801/blackhole_data/data/errors.dat' --o='/scratch/alpine/mahu8801/blackhole_data/data/bhmasses.dat'