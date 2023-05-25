#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.table import Table, join, vstack, ColumnGroups, Column
import numpy as np
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
from astropy.cosmology import WMAP9 as cosmo
import scipy.integrate as integrate
from astropy import constants as const


# In[8]:


sample = Table.read("602K Sample",
        format="ascii")
sigma_0 = Table.read("veldisp_maghuber.csv",
        format="ascii")
r_e = Table.read("s11radsersic.tsv",
        format="ascii")

samplewre = join(sample,r_e,keys='objID')

sigma_0.rename_column('bestObjID','objID')


# In[9]:


Msun = 1.988435e30
G = const.G.to('km3 / (kg s2)').value
kpctokm = 3.086e16
rap = 1.5
radtoarcsec = 206265
def sigma_inf(Mstar,n,re):
    Kvn = 73.32/(10.465+(n-0.94)**2)+0.954
    Kstar = 0.577*Kvn
    M = Mstar*Msun
    rad = re*kpctokm
    frac = (G*M)/(Kstar*rad)
    return np.sqrt(frac)


# In[10]:


infvd1 = []
for i,obj in enumerate(samplewre['objID']):
    Mstar = float(10**samplewre['LogMbd'][i])
    re = float(samplewre['Rchl_r'][i])
    n = float(samplewre['ng'][i])
    infvd1.append(sigma_inf(Mstar,n,re))


# In[11]:


samplewre.add_column(infvd1,name='sig_inf')


# In[12]:


samplewre = samplewre[samplewre['dBD'] < 1]


# In[13]:


from astropy.io import ascii
ascii.write(samplewre, '550K Sample with siginf')


# In[3]:


vdcut = sigma_0['velDisp'] > 0
validvd = sigma_0[vdcut]
vdsample = join(samplewre,validvd,keys='objID')


# In[4]:


pc = []
for i,obj in enumerate(vdsample['velDisp']):
    pc.append(vdsample['velDispErr'][i]/vdsample['velDisp'][i])
vdsample.add_column(pc,name='error')


# In[5]:


errcut = vdsample['error'] < 0.1
errcorrected = vdsample[errcut]
upper = errcorrected['z'] < 0.07
uppercut = errcorrected[upper]
lower = uppercut['z'] > 0.05
sample = uppercut[lower]


# In[7]:


infvd = []
for i,obj in enumerate(sample['objID']):
    Mstar = float(10**sample['LogMbd'][i])
    re = float(sample['Rchl_r'][i])
    n = float(sample['ng'][i])
    infvd.append(sigma_inf(Mstar,n,re))


# In[8]:


sigma_corr = []
for i,obj in enumerate(sample['velDisp']):
    dA = cosmo.angular_diameter_distance(sample['z'][i]).to('kpc').value
    rkpc = float(sample['Rchl_r'][i])
    r_as = (rkpc/dA)*radtoarcsec 
    corr = (8*rap/r_as)**0.066
    sigma_corr.append(sample['velDisp'][i]*corr)


# In[9]:


plt.scatter(np.log10(infvd),np.log10(sigma_corr),marker='.',s=5,color='k',alpha=0.05)
#plt.hist2d(np.log10(infvd),np.log10(sigma_corr), cmap=plt.cm.Greys)
x = np.linspace(1.25,3,1000)
plt.xlim(1.8,2.8)
plt.ylim(1.8,2.8)
plt.plot(x,x,'--',color='k')
plt.xlabel(r'$\mathrm{log}(\sigma_{\mathrm{inf}}) \ \mathrm{[km/s]}$')
plt.ylabel(r'$\mathrm{log}(\sigma_{0}) \ \mathrm{[km/s]}$')
#plt.savefig('veldisp.png')


# In[10]:


decomp = vdsample['dBD'] < 1
reducedsample = vdsample[decomp]
print(len(reducedsample))


# In[11]:


infvd2 = []
for i,obj in enumerate(reducedsample['objID']):
    Mstar = float(10**reducedsample['LogMbd'][i])
    re = float(reducedsample['Rchl_r'][i])
    n = float(reducedsample['ng'][i])
    infvd2.append(sigma_inf(Mstar,n,re))


# In[12]:


sigma_corr2 = []
for i,obj in enumerate(reducedsample['velDisp']):
    dA = cosmo.angular_diameter_distance(reducedsample['z'][i]).to('kpc').value
    rkpc = float(reducedsample['Rchl_r'][i])
    r_as = (rkpc/dA)*radtoarcsec
    if reducedsample['z'][i] >= 0.14:
        corr = (8*rap/r_as)**0.04
    corr = (8*rap/r_as)**0.066
    sigma_corr2.append(reducedsample['velDisp'][i]*corr)


# In[13]:


fig, axs = plt.subplots(3,2,figsize=(9,11))
redshift = np.asarray(reducedsample['z'])
d = np.digitize(redshift,np.arange(0.02,0.2,0.03))
for i,ax in zip(range(1,max(d)),axs.reshape(-1)):
    ind = np.where(d == i)
    x = np.asarray(np.log10(infvd2))[ind]
    y = np.asarray(np.log10(sigma_corr2))[ind]
    z = redshift[ind]
    ax.plot([0], label = '%.2f to %.2f'%(min(z),max(z)), visible=False)
    ax.legend(loc='upper left', handlelength=0, handletextpad=0)
    ax.scatter(x,y,marker='.',s=5,color='k',alpha=0.05)
    x2 = np.linspace(1.25,3,1000)
    ax.set_xlim(1.5,2.8)
    ax.set_ylim(1.5,2.8)
    ax.plot(x2,x2,'--',color='k')
    ax.set_xlabel(r'$\mathrm{log}(\sigma_{\mathrm{inf}}) \ \mathrm{[km/s]}$')
    ax.set_ylabel(r'$\mathrm{log}(\sigma_{0}) \ \mathrm{[km/s]}$')
#fig.savefig('vdplots.png')


# In[14]:


col1 = Column(name='sigma_0', data=sigma_corr2)
col2 = Column(name='sigma_inf',data=infvd2)
reducedsample.add_columns([col1,col2])
print(reducedsample.colnames)


# In[15]:


sample550k = Table.read("550K Sample",
        format="ascii")
infvd3 = []
for i,obj in enumerate(sample550k['objID']):
    Mstar = float(10**sample550k['LogMbd'][i])
    re = float(sample550k['Rchl_r'][i])
    n = float(sample550k['ng'][i])
    infvd3.append(sigma_inf(Mstar,n,re))
sigma_corr3 = []
for i,obj in enumerate(sample550k['velDisp']):
    dA = cosmo.angular_diameter_distance(sample550k['z'][i]).to('kpc').value
    rkpc = float(sample550k['Rchl_r'][i])
    r_as = (rkpc/dA)*radtoarcsec
    if sample550k['z'][i] >= 0.14:
        corr = (8*rap/r_as)**0.04
    corr = (8*rap/r_as)**0.066
    sigma_corr3.append(sample550k['velDisp'][i]*corr)
col3 = Column(name='sigma_0', data=sigma_corr3)
col4 = Column(name='sigma_inf',data=infvd3)
sample550k.add_columns([col3,col4])


# In[ ]:


# from astropy.io import ascii
# ascii.write(reducedsample, '290K Sample')

