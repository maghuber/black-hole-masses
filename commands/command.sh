#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00			
#SBATCH --mail-type=all
#SBATCH --mail-user=mahu8801@colorado.edu
#SBATCH --qos=normal	      			
#SBATCH --output=/scratch/alpine/mahu8801/output/python_%j.out

# Load the python module
module load python
module load anaconda
conda activate blackhole_env

# Run Python Script
cd /scratch/alpine/mahu8801/blackhole_data/scripts
python mbulge_sigma.py --slice1=512000 --slice2=512821 --ndraws=500 --relation='kormendyho' --i='/scratch/alpine/mahu8801/blackhole_data/photometry.dat' --o='/scratch/alpine/mahu8801/blackhole_data/data/512000_512821.dat'