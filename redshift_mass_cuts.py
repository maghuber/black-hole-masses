#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.table import Table, join, vstack
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from astroquery.vizier import Vizier
from astropy.cosmology import WMAP9 as cosmo
import scipy.integrate as integrate
from astropy import constants as const
import random


# In[2]:


t4 = Table.read("table4.dat.gz",
        format="ascii")


# In[3]:


t4names = ['objID','z','logMt','b_logMt','B_logMt','logMb','b_logMb','B_logMb','logMd','b_logMd','B_logMd','zmin','zmax','PpS','Type','dBD']
for i,name in enumerate(t4.colnames):
    t4.rename_column(name,t4names[i])


# In[4]:


mask1 = t4['z'] <= 0.2
t41 = t4[mask1]
mask2 = t41['z'] >= 0.02
tredshiftcut = t41[mask2]


# In[5]:


mask3 = tredshiftcut['logMb'] > -99
tmasscut = tredshiftcut[mask3]
print(len(tredshiftcut),len(tmasscut))


# In[6]:


Mbd = []
logMbd = []
for i,obj in enumerate(tmasscut['logMb']):
    if tmasscut['dBD'][i] > 1:
        logMt = tmasscut['logMt'][i]
        Mbd.append(10**logMt)
        logMbd.append(logMt)
    else:
        logMb = tmasscut['logMb'][i]
        Mb = 10**logMb
        logMd = tmasscut['logMd'][i]
        Md = 10**logMd
        Mbd.append(Mb+Md)
        logMbd.append(np.log10(Mb+Md))


# In[9]:


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


# In[17]:


fig,ax = plt.subplots()
ax.hist(logMbd,bins=20,weights=np.ones_like(logMbd) / len(logMbd),edgecolor='white', linewidth=1.2,alpha=0.5,color='tab:cyan')
ax.set_ylabel(r'$\mathrm{Relative \ Distribution \ [n_{bin}/n_{total}]}$')
ax.set_xlabel(r'$\mathrm{Stellar \ Mass \ }[\log(M/M_{\odot})]$')
ax.set_yticks(np.arange(0,0.25,0.05))
ax.set_xticks(np.arange(8,13,1))
ax.set_xlim(8,12)
ax2=ax.twinx()
def div_104(x, *args):
    """
    The function that will you be applied to your y-axis ticks.
    """
    x = float(x)/10**4
    return "{:.1f}".format(x) 
ax2.hist(logMbd,bins=20,edgecolor='black', linewidth=1.2,alpha=0)
ax2.yaxis.set_major_formatter(mtick.FuncFormatter(div_104))
ax2.set_ylabel(r'$\mathrm{Number \ of \ galaxies} [10^4 n_{bin}]$')


# In[13]:


fig,ax = plt.subplots()
ax.hist(tmasscut['z'],bins=20,weights=np.ones_like(tmasscut['z']) / len(tmasscut['z']),edgecolor='white', linewidth=1.2,alpha=0.5,color='tab:pink')
ax.set_ylabel(r'$\mathrm{Relative \ Distribution \ [n_{bin}/n_{total}]}$')
ax.set_xlabel(r'$\mathrm{Redshift, \ z}$')
ax2=ax.twinx()
ax2.hist(tmasscut['z'],bins=20,edgecolor='black', linewidth=1.2,alpha=0)
ax2.yaxis.set_major_formatter(mtick.FuncFormatter(div_104))
ax2.set_ylabel(r'$\mathrm{Number \ of \ galaxies} [10^4 n_{bin}]$')
plt.savefig('redshifthist.png')


# In[9]:


s11 = Table.read("asu.tsv",
        format="ascii")


# In[10]:


joined = join(tmasscut,s11,keys='objID')


# In[11]:


BT = []
for i,obj in enumerate(joined['logMb']):
    if joined['dBD'][i] > 1 and joined['(B/T)r'][i]!='      ':
        BTr = joined['(B/T)r'][i]
        BT.append(float(BTr))
    else:
        logMb = joined['logMb'][i]
        logMd = joined['logMd'][i]
        Mb = 10**logMb
        Md = 10**logMd
        BT.append(float(Mb/(Mb+Md)))
joined.add_column(BT,name='(B/T)')


# In[12]:


for i,obj in enumerate(joined['(B/T)']):
    if joined['PpS'][i] > 0.32 and joined['(B/T)'][i] > 0.7:
        joined['(B/T)'][i] = 1
        logMt = joined['logMt'][i]
        Mbd[i] = 10**logMt


# In[13]:


BT20 = []
BT50 = []
BT80 = []
BT100 = []
for i,obj in enumerate(joined['(B/T)']):
    if obj < 0.2:
        BT20.append(obj)
    if obj >= 0.2 and obj < 0.5:
        BT50.append(obj)
    if obj >= 0.5 and obj < 0.8:
        BT80.append(obj)
    if obj >= 0.8:
        BT100.append(obj)


# In[14]:


print(('BT20:%i  BT50:%i  BT80:%i  BT100:%i')%(len(BT20),len(BT50),len(BT80),len(BT100)))


# In[15]:


def Vmax(galaxy):
   skyfrac = 8200
   angular = skyfrac/41253
   zmin = max(0.02,joined['zmin'][galaxy])
   zmax = min(0.2,joined['zmax'][galaxy])
   zgalaxy = joined['z'][galaxy]
   Vint = integrate.quad(lambda z: (cosmo.angular_diameter_distance(z).value**2*const.c.to('Mpc/s').value)/(cosmo.H(z).to('1/s').value*(1+z)),0.02,zgalaxy)
   Vmaxint = integrate.quad(lambda z: (cosmo.angular_diameter_distance(z).value**2*const.c.to('Mpc/s').value)/(cosmo.H(z).to('1/s').value*(1+z)),zmin,zmax)
   V = angular*Vint[0]
   Vmax = angular*Vmaxint[0]
   return Vmax**-1, V/Vmax


# In[16]:


number = []
for i,obj in enumerate(joined['logMb']):
    number.append(i)
joined.add_column(number,name='Number')


# In[17]:


get_ipython().run_cell_magic('time', '', "V_max = []\nV_Vmax = []\nfor obj in joined['Number']:\n    Vresult = Vmax(obj)\n    V_max.append(Vresult[0])\n    V_Vmax.append(Vresult[1])\n")


# In[18]:


print(np.median(V_Vmax))


# In[19]:


joined.add_column(V_max,name='1/Vmax')


# In[20]:


from astropy.io import ascii
ascii.write(joined, '603K Sample')

