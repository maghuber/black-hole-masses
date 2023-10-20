from astropy.table import Table
import numpy as np
from astropy.io import ascii
from astropy import constants as const
import matplotlib.pyplot as plt
from scipy.stats import rv_continuous, norm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--slice1",type=int,default=0,
                    help="lower slice of galaxy table")
parser.add_argument("--slice2",type=int,default=1000,
                    help="upper slice of galaxy table")
parser.add_argument("--ndraws",type=int,
                    help="number of random draws")
parser.add_argument("--i",type=str,
                    help="input file")
parser.add_argument("--o",type=str,
                    help="output file")
args = parser.parse_args()

n_draws = args.ndraws
n_galaxies = args.slice2 - args.slice1
in_file = args.i
out_file = args.o

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
    test_rvs = twosided_gaussian(a=4,b=13)
    for i, q in enumerate(quantity):
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

def Mbh_bulge(Mbulge,a,b):
    M = Mbulge*Msun
    logMbh = a + b*np.log10(Mbulge/(10**11))
    logMbh_scatter = np.random.normal(logMbh, scatter_b)
    return logMbh_scatter

def Mbulge_inf(logMt,g,r):
    color_cut = 0.06*logMt-0.01
    gr_color = abs(g-r)
    if gr_color > color_cut:
        logMb_inf = 1.*logMt
    else:
        logMb_inf = 0.31*logMt
    return logMb_inf

table = Table.read(in_file,
        format="ascii")

indices = np.s_[args.slice1:args.slice2]
    
objid = np.asarray(table['objID'][indices],dtype=str)
sig1_Mt = np.asarray(table['sig1_Mt'][indices],dtype=float)
sig2_Mt = np.asarray(table['sig2_Mt'][indices],dtype=float)
logMt = np.asarray(table['logMt'][indices],dtype=float)
gg2d = np.asarray(table['gg2d'][indices],dtype=float)
rg2d = np.asarray(table['rg2d'][indices],dtype=float)

logMt_dist = make_distribution(logMt,sig1_Mt,sig2_Mt,n_draws)

# make coefficient distributions
alpha_b_dist = make_distribution_normal(np.full((n_galaxies,),alpha_b),np.full((n_galaxies,),e_alpha_b),n_draws)
beta_b_dist = make_distribution_normal(np.full((n_galaxies,),beta_b),np.full((n_galaxies,),e_beta_b),n_draws)

# make black hole mass distributions
logMb_inf_dist = np.zeros((n_galaxies,n_draws))
for i in np.arange(n_galaxies):
    for j in np.arange(n_draws):
        gg = gg2d[i]
        rg = rg2d[i]
        Mt = logMt_dist[i][j]
        logMb_inf_dist[i][j] = Mbulge_inf(Mt,gg,rg)
        
logMbh_b_dist = np.zeros((n_galaxies,n_draws))
for i in np.arange(n_galaxies):
    for j in np.arange(n_draws):
        alpha = alpha_b_dist[i][j]
        beta = beta_b_dist[i][j]
        Mbulge = 10**logMb_inf_dist[i][j]
        logMbh_b_dist[i][j] = Mbh_bulge(Mbulge,alpha,beta)

logMbh_b_median = np.zeros(n_galaxies)
logMbh_b_err = np.zeros(n_galaxies)
for i in np.arange(n_galaxies):
    logMbh_b_median[i] = np.around(np.median(logMbh_b_dist[i]),3)
    logMbh_b_err[i] = np.around(np.std(logMbh_b_dist[i]),3)

arr = {'objID' : objid,
       'logMbh_b': logMbh_b_median,
       'err_logMbh_b': logMbh_b_err}
table = Table(arr)

ascii.write(table, out_file, overwrite=True) 
