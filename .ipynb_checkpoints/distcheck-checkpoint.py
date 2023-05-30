#!/usr/bin/env python
from astropy.table import Table
import numpy as np
from astropy.io import ascii
from astropy import constants as const
import matplotlib.pyplot as plt
from scipy.stats import rv_continuous, norm
import argparse
import pickle

parser = argparse.ArgumentParser()
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

# mcconnell and ma coefficients
# M-sigma for all galaxies
alpha_s = 8.32
e_alpha_s = 0.05
beta_s = 5.64
e_beta_s = 0.32
scatter_s = 0.38
# M-Mbulge for dynamical masses
alpha_b = 8.46
e_alpha_b = 0.08
beta_b = 1.05
e_beta_b = 0.11
scatter_b = 0.34
# constants
Msun = 1.988435e30
kpctokm = 3.086e16
G = const.G.to('km3 / (kg s2)').value

## Make random distributions from values and uncertainties
class twosided_gaussian(rv_continuous):
    "Two-Sided Gaussian distribution"
    def _pdf(self, x, mu, sig1, sig2):
        A = np.sqrt(2 / np.pi) / (sig1 + sig2)
        if x < mu:
            return A * np.exp(-(x-mu)**2 / (2 * sig1**2))
        else:
            return A * np.exp(-(x-mu)**2 / (2 * sig2**2))

def make_distribution(quantity,s1,s2,s):
    errordist = np.zeros((len(quantity),s))
    test_rvs = twosided_gaussian(a=5,b=15)
    for i, q in enumerate(quantity):
        print(s1)
        print(s2)
        test_draws = test_rvs.rvs(mu=quantity[i], sig1=s1[i], sig2=s2[i], size=s)
        errordist[i] = (test_draws)
    return errordist

def make_distribution_normal(quantity,sigma,s):
    errordist = np.zeros((len(quantity),s))
    for i, q in enumerate(quantity):
        test_draws = np.random.normal(quantity[i], sigma[i], s)
        errordist[i] = (test_draws)
    return errordist

# calculate the black hole masses
def Mbh_sigma(Mtotal,n,r_e,a,b):
    Kstar = 0.577*(73.32/(10.465+(n-0.94)**2)+0.954)
    rad = r_e*kpctokm
    M = Mtotal*Msun
    sigma = np.sqrt((G*M)/(Kstar*rad))
    logMbh = a + b*np.log10(sigma/200)
    logMbh_scatter = np.random.normal(logMbh, scatter_s)
    return logMbh_scatter

def Mbh_bulge(Mbulge,a,b):
    M = Mbulge*Msun
    logMbh = a + b*np.log10(Mbulge/(10**11))
    logMbh_scatter = np.random.normal(logMbh, scatter_b)
    return logMbh_scatter

def Mbh_bulge_noscatter(Mbulge,a,b):
    M = Mbulge*Msun
    logMbh = a + b*np.log10(Mbulge/(10**11))
    return logMbh

def Mbh_sigma_noscatter(Mtotal,n,r_e,a,b):
    Kstar = 0.577*(73.32/(10.465+(n-0.94)**2)+0.954)
    rad = r_e*kpctokm
    M = Mtotal*Msun
    sigma = np.sqrt((G*M)/(Kstar*rad))
    logMbh = a + b*np.log10(sigma/200)
    return logMbh

errors = Table.read(in_file,
        format="ascii")

logMb = np.asarray(errors['logMb'],dtype=float)
b_logMb = np.asarray(errors['b_logMb'],dtype=float)
B_logMb = np.asarray(errors['B_logMb'],dtype=float)
sig1_Mb = np.asarray(errors['sig1_Mb'],dtype=float)
sig2_Mb = np.asarray(errors['sig2_Mb'],dtype=float)
sig1_Mt = np.asarray(errors['sig1_Mt'],dtype=float)
sig2_Mt = np.asarray(errors['sig1_Mt'],dtype=float)
logMt = np.asarray(errors['logMt'],dtype=float)
b_logMt = np.asarray(errors['b_logMt'],dtype=float)
B_logMt = np.asarray(errors['B_logMt'],dtype=float)
ng = np.asarray(errors['ng'],dtype=float)
e_ng = np.asarray(errors['e_ng'],dtype=float)
Rhl = np.asarray(errors['Rhlr'],dtype=float)
e = np.asarray(errors['e'],dtype=float)
e_e = np.asarray(errors['e_e'],dtype=float)

logMb_dist = make_distribution(logMb,sig1_Mb,sig2_Mb,n_draws)
logMt_dist = make_distribution(logMt,sig1_Mt,sig2_Mt,n_draws)
ng_dist = make_distribution_normal(ng,e_ng,n_draws)
e_dist = make_distribution_normal(e,e_e,n_draws)

# make coefficient distributions
alpha_s_dist = make_distribution_normal(np.full((n_galaxies,),alpha_s),np.full((n_galaxies,),e_alpha_s),n_draws)
beta_s_dist = make_distribution_normal(np.full((n_galaxies,),beta_s),np.full((n_galaxies,),e_beta_s),n_draws)
alpha_b_dist = make_distribution_normal(np.full((n_galaxies,),alpha_b),np.full((n_galaxies,),e_alpha_b),n_draws)
beta_b_dist = make_distribution_normal(np.full((n_galaxies,),beta_b),np.full((n_galaxies,),e_beta_b),n_draws)

# make black hole mass distributions
logMbh_s_dist = np.zeros((n_galaxies,n_draws))
for i in np.arange(n_galaxies):
    for j in np.arange(n_draws):
        alpha = alpha_s_dist[i][j]
        beta = beta_s_dist[i][j]
        r_e = Rhl[i]*np.sqrt(1-e[i])
        Mtot = 10**(logMt_dist[i][j])
        logMbh_s_dist[i][j] = Mbh_sigma(Mtot,ng[i],r_e,alpha,beta)
        
logMbh_b_dist = np.zeros((n_galaxies,n_draws))
for i in np.arange(n_galaxies):
    for j in np.arange(n_draws):
        alpha = alpha_b_dist[i][j]
        beta = beta_b_dist[i][j]
        Mbulge = 10**logMb_dist[i][j]
        logMbh_b_dist[i][j] = Mbh_bulge(Mbulge,alpha,beta)

logMbh_s_median = np.zeros(n_galaxies)
logMbh_s_err = np.zeros(n_galaxies)
logMbh_b_median = np.zeros(n_galaxies)
logMbh_b_err = np.zeros(n_galaxies)
for i in np.arange(n_galaxies):
    logMbh_s_median[i] = np.around(np.median(logMbh_s_dist[i]),3)
    logMbh_b_median[i] = np.around(np.median(logMbh_b_dist[i]),3)
    logMbh_s_err[i] = np.around(np.std(logMbh_s_dist[i]),3)
    logMbh_b_err[i] = np.around(np.std(logMbh_b_dist[i]),3)

logMbh_s_noerr = np.zeros(n_galaxies)
for i in np.arange(n_galaxies):
    alpha = alpha_s
    beta = beta_s
    r_e = Rhl[i]*np.sqrt(1-e[i])
    Mtot = 10**(logMt[i])
    logMbh_s_noerr[i] = Mbh_sigma_noscatter(Mtot,ng[i],r_e,alpha,beta)
        
logMbh_b_noerr = np.zeros(n_galaxies)
for i in np.arange(n_galaxies):
    alpha = alpha_b
    beta = beta_b
    Mbulge = 10**logMb[i]
    logMbh_b_noerr[i] = Mbh_bulge_noscatter(Mbulge,alpha,beta)

my_data = {'logMbh_s_dist': logMbh_s_dist,
           'logMbh_b_dist': logMbh_b_dist,
           'b_median': logMbh_b_median,
           's_median': logMbh_s_median,
           's_noerr' : logMbh_s_noerr,
           'b_noerr' : logMbh_b_noerr,
           'logMtdist' : logMt_dist,
           'logMbdist' : logMb_dist,
           'logMt' : logMt,
           'logMb' : logMb}

output = open(out_file, 'wb')
pickle.dump(my_data, output)
output.close()