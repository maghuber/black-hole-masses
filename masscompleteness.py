#!/usr/bin/env python
# coding: utf-8

# In[10]:


from astropy.table import Table, join, vstack, ColumnGroups
import numpy as np
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
from astropy.cosmology import WMAP9 as cosmo
import scipy.integrate as integrate
from scipy.stats import norm, linregress, chisquare
from astropy import constants as const


# In[2]:


t = Table.read("550K Sample",
        format="ascii")
petro = Table.read('petro2_maghuber.csv',
        format="ascii")
petro.rename_column('ObjID','objID')
sample = join(t,petro,keys='objID')
logMbd = np.asarray(sample['LogMbd'])


# In[3]:


mu50 = []
for i, obj in enumerate(sample['petroMag_r']):
    rp = float(sample['petroMag_r'][i])
    theta50 = float(sample['petroR50_r'][i])
    mu50.append(rp+2.5*np.log10(2*np.pi*theta50**2))
mu50 = np.asarray(mu50)


# In[44]:


d = np.digitize(logMbd,np.arange(8.2,11.3,0.2))
xmean = []
ymean = []
ysig = []
for i in range(0,max(d)+1):
    ind = np.where(d == i)
    x = logMbd[ind]
    y = mu50[ind]
    xmu, xstd = norm.fit(x)
    ymu, ystd = norm.fit(y)
    xmean.append(xmu)
    ymean.append(ymu)
    ysig.append(ystd)


# In[53]:


x2 = []
y2 = []
for i,x in enumerate(xmean):
    if x > 9.1 and x < 10.5:
        x2.append(x)
        y2.append(ymean[i])


# In[68]:


result = linregress(x2, y2)
xvals = np.linspace(9.1,10.4,100)
yvals = result.slope*xvals+result.intercept
print(result.slope,chisquare(yvals))


# In[70]:


plt.errorbar(xmean,ymean,yerr = ysig,fmt = '^',linewidth=1, capsize=6, color='r')
xvals = np.linspace(8.0,11.12)
yvals = result.slope*xvals+result.intercept
plt.plot(xvals, yvals,'--',color='k',)
plt.xlabel('Galaxy stellar mass [Log(M/Msun)]')
plt.ylabel(r'$\mathrm{Surface \ brightness}, [\mu_r/mag.arcsec^{-1}]$')
plt.ylim(18,24)
plt.xlim(8,12)
plt.savefig('massincompleteness.png')


# In[72]:


print(sample.colnames)

