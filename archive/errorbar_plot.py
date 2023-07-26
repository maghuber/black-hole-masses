#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.table import Table
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['axes.labelsize'] = 14
matplotlib.rcParams['legend.fontsize'] = 14
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['savefig.dpi'] = 400
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\mathchardef\mhyphen="2D']
matplotlib.rcParams['legend.fontsize'] = 14
matplotlib.rcParams["axes.edgecolor"] = 'black'
matplotlib.rcParams["legend.edgecolor"] = '0.8'
colors = sns.color_palette("PuBu")


# In[6]:


data = Table.read("tables/1000_galaxies.dat",
        format="ascii")


# In[7]:


x = data['logMbh_b']
y = data['logMbh_s']
yerr = data['err_logMbh_s']
xerr = data['err_logMbh_b']


# In[8]:


plt.figure(dpi=200)
plt.errorbar(x,y,xerr=xerr,yerr=yerr,marker='.',linewidth=0,elinewidth=.3,color='k',alpha=0.2)
plt.ylabel(r'$\log(M_{\bullet,Bulge}) (M_{\odot})$')
plt.xlabel(r'$\log(M_{\bullet,\sigma}) (M_{\odot})$')
plt.title(r'$\mathrm{1000 \ most \ massive } \  M_{Bulge} \ \mathrm{galaxies}$')
plt.savefig('plots/1000_n_galaxies.pdf');

