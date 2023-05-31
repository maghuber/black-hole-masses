#!/usr/bin/env python
from astropy.table import Table
import numpy as np
from astropy.io import ascii
from astropy import constants as const
import matplotlib.pyplot as plt
from scipy.stats import rv_continuous, norm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--ngalaxies",type=int,default=1000,
                    help="number of galaxies")
parser.add_argument("--all",action='store_true',
                    help="takes all galaxies in the sample")
parser.add_argument("--ndraws",type=int,
                    help="number of random draws")
parser.add_argument("--i",type=str,
                    help="input file")
parser.add_argument("--o",type=str,
                    help="output file")
args = parser.parse_args()

n_draws = args.ndraws
in_file = args.i
out_file = args.o

# constants
Msun = 1.988435e30
kpctokm = 3.086e16
G = const.G.to('km3 / (kg s2)').value

## Make random distributions from values and uncertainties

def make_distribution_normal(quantity,sigma,s):
    errordist = np.zeros((len(quantity),s))
    for i, q in enumerate(quantity):
        test_draws = np.random.normal(quantity[i], sigma[i], s)
        errordist[i] = (test_draws)
    return errordist

# calculate the black hole mass

def Mbh_Hb(FWHMHb,L5100):
    L = 10**L5100/1e44
    FWHM = FWHMHb/1000
    logMbh = np.log10(FWHM**2*L**0.533)+6.91
    return logMbh
    
errors = Table.read(in_file,
        format="ascii")

if args.all:
    n_galaxies = len(errors)
    indices = np.s_[:n_galaxies]
else:
    n_galaxies = args.ngalaxies
    indices = np.s_[-n_galaxies:]

fwhm = np.asarray(errors['FWHM-BHb'][indices],dtype=float)
e_fwhm = np.asarray(errors['e_FWHM-BHb'][indices],dtype=float)
logL5100 = np.asarray(errors['logL5100'][indices],dtype=float)
e_logL5100 = np.asarray(errors['e_logL5100'][indices],dtype=float)

fwhm_dist = make_distribution_normal(fwhm,e_fwhm,n_draws)
logL5100_dist = make_distribution_normal(logL5100,e_logL5100,n_draws)

# make black hole mass distributions
logMbh_Hb_dist = np.zeros((n_galaxies,n_draws))
for i in np.arange(n_galaxies):
    for j in np.arange(n_draws):
        logMbh_Hb_dist[i][j] = Mbh_Hb(fwhm_dist[i][j],logL5100_dist[i][j])

logMbh_Hb_median = np.zeros(n_galaxies)
logMbh_Hb_err = np.zeros(n_galaxies)

for i in np.arange(n_galaxies):
    logMbh_Hb_median[i] = np.around(np.median(logMbh_Hb_dist[i]),3)
    logMbh_Hb_err[i] = np.around(np.std(logMbh_Hb_dist[i]),3)

arr = {'logMbh_Hb': logMbh_Hb_median,
       'err_logMbh_Hb': logMbh_Hb_err}
table = Table(arr)

ascii.write(table, out_file, overwrite=True) 