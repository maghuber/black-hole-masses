#!/usr/bin/env python
# coding: utf-8

# In[2]:


from astropy.table import Table, join, vstack, ColumnGroups, Column
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
from astropy.cosmology import WMAP9 as cosmo
import scipy.integrate as integrate
from scipy.stats import norm, linregress, chisquare
from astropy import constants as const
import seaborn as sns


# # M_bulge sigma agreement

# In[3]:


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


# In[3]:


sample = Table.read("550K Sample with MBH",
        format="ascii")


# In[4]:


sample = sample[ sample['logMbd'] >= 9.1 ]


# In[5]:


BT20 = sample[ sample['(B/T)'] < 0.2 ]
BT50l = sample[ sample['(B/T)'] >= 0.2 ]
BT50 = BT50l[ BT50l['(B/T)'] < 0.5 ]
BT80l = sample[ sample['(B/T)'] >= 0.5 ]
BT80 = BT80l[ BT80l['(B/T)'] < 0.8 ]
BT100 = sample[ sample['(B/T)'] >= 0.8 ]


# In[6]:


passive = sample[ sample['color_type'] == 'passive' ]
active = sample[ sample['color_type'] == 'active' ]
tables2 = [passive,active]
tablenames2 = ['passive','active']


# In[7]:


early = join(BT80l,passive)
lowbulge = sample[ sample['(B/T)'] <= 0.5 ]
late = join(lowbulge,active)


# In[14]:


tables = [early,late]
tablenames = ['B/T greater than 0.5 and passive','B/T less than 0.5 and active']
figure = ['early','late']
for table,name,f in zip(tables,tablenames,figure):
    
    t3 = table[table['logMbd'] > 11]
    t22 = table[table['logMbd'] >= 10]
    t2 = t22[t22['logMbd'] <= 11]
    t1 = table[table['logMbd'] < 10]
    
    massrange = ['10^{10} \ < \ M_* \ <  \ 10^{11}','M_* \ > \ 10^{11}']
    
    fig, axs = plt.subplots(1,2,figsize=(10,5),sharey=True,sharex=True,dpi=300)
    fig.suptitle('%s'%(name),fontsize=22)
    x = np.linspace(1,15,10)

    for t, ax, m in zip((t2,t3),axs,massrange):
        inds = np.where(np.isfinite(t['logMbh_siginf']) == True)
        ax.hist2d(t['logMbh_siginf'][inds],t['logMbh_bulge'][inds],bins=[30,30],cmap='bone_r')
        ax.set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ \sigma_{inf} (M_\odot)$')
        ax.set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ M_{bulge} (M_\odot)$')
        ax.plot(x,x,'--',color='k')
        ax.set_xlim(4,12)
        ax.set_ylim(5,10)
        ax.plot([0], label = r'$\mathrm{%s, \ %i \ galaxies}$'%(m,len(t)), visible=False)
        ax.legend(loc='upper left', handlelength=0, handletextpad=0)
        ax.grid(alpha=0.5,which='both')
    fig.savefig('%s.png'%(f))


# In[6]:


tables = [BT20,BT50,BT80,BT100]
tablenames = ['BT20','BT50','BT80','BT100']
for table,name in zip(tables,tablenames):
    
    t3 = table[table['logMbd'] > 11]
    t22 = table[table['logMbd'] >= 10]
    t2 = t22[t22['logMbd'] <= 11]
    t1 = table[table['logMbd'] < 10]
    
    massrange = ['M_* \ < \ 10^{10}','10^{10} \ < \ M_* \ <  \ 10^{11}','M_* \ > \ 10^{11}']
    
    fig, axs = plt.subplots(1,3,figsize=(14,5),sharey=True,sharex=True)
    fig.suptitle(r'$\mathrm{Bulge \ fraction: \ %s}$'%(name),fontsize=16)
    x = np.linspace(1,15,10)

    for t, ax, m in zip((t1,t2,t3),axs,massrange):
        inds = np.where(np.isfinite(t['logMbh_siginf']) == True)
        ax.hist2d(t['logMbh_siginf'][inds],t['logMbh_bulge'][inds],bins=[30,30],cmap='bone_r')
        ax.set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ \sigma_{inf}$')
        ax.set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ M_{bulge}$')
        ax.plot(x,x,'--',color='k')
        ax.set_xlim(4,12)
        ax.set_ylim(5,10)
        ax.plot([0], label = r'$\mathrm{%s, \ %i \ galaxies}$'%(m,len(t)), visible=False)
        ax.legend(loc='upper left', handlelength=0, handletextpad=0)
        ax.grid(alpha=0.5,which='both')
    #fig.savefig('%s.png'%(name))


# # Fundamental Plane

# In[8]:


def FundamentalPlane(Mb,Mbe,sig,sige):
    a = 0.37
    ae = 0.09
    b = 3.19
    be = 0.52
    ZP = -2.93
    ZPe = 0.66
    log10Mbh = a*np.log10(Mb)+b*np.log10(sig)+ZP
    dMdMb = a/(Mb*np.log(10))
    dMds = b/(sig*np.log(10))
    dMda = np.log10(Mb)
    dMdb = np.log10(sig)
    err = np.sqrt(ZPe**2+(dMdMb*Mbe)**2+(dMds*sige)**2+(dMda*ae)**2+(dMdb*be)**2)
    return log10Mbh, err
def Sifoni_MMbulge(Mb,Mbe):
    alpha = 0.962
    ae = 0.066 
    ZPe = 0.716
    ZP = -2.099
    log10Mbh = alpha*np.log10(Mb)+ZP
    dMda = np.log10(Mb)
    dMdMb = alpha/(Mb*np.log(10))
    err = np.sqrt(ae**2+(dMdMb*Mbe)**2+ZPe**2)
    return log10Mbh, err
def Sifoni_Msigma(sig,sige):
    alpha = 5.246
    ae = 0.274 
    ZPe = 0.631
    ZP = -3.77
    log10Mbh = alpha*np.log10(sig)+ZP
    dMda = np.log10(sig)
    dMds = alpha/(sig*np.log(10))
    err = np.sqrt(ae**2+(dMds*sige)**2+ZPe**2)
    return log10Mbh, err


# In[9]:


Mbh_fp = []
Mbh_bulge = []
Mbh_sig = []
e_Mbh_fp = []
e_Mbh_bulge = []
e_Mbh_sig = []
for i in range(len(sample)):
    Mb = 10**sample['logMb'][i]
    Mbe = 10**sample['e_logMb'][i]
    sig = sample['siginf'][i]
    sige = sample['e_siginf'][i]
    fp = FundamentalPlane(Mb,Mbe,sig,sige)
    mbhb = Sifoni_MMbulge(Mb,Mbe)
    mbhs = Sifoni_Msigma(sig,sige)
    Mbh_fp.append(fp[0])
    e_Mbh_fp.append(fp[1])
    Mbh_bulge.append(mbhb[0])
    e_Mbh_bulge.append(mbhb[1])
    Mbh_sig.append(mbhs[0])
    e_Mbh_sig.append(mbhs[1])


# In[10]:


col1 = Column(name='logMbh_fp',data=Mbh_fp)
col2 = Column(name='e_logMbh_fp',data=e_Mbh_fp)
col3 = Column(name='logMbh_sbulge',data=Mbh_bulge)
col4 = Column(name='e_logMbh_sbulge',data=e_Mbh_bulge)
col5 = Column(name='logMbh_ssig',data=Mbh_sig)
col6 = Column(name='e_logMbh_ssig',data=e_Mbh_sig)
sample.add_columns([col1,col2,col3,col4,col5,col6])


# In[24]:


BT20 = sample[ sample['(B/T)'] < 0.2 ]
BT50l = sample[ sample['(B/T)'] >= 0.2 ]
BT50 = BT50l[ BT50l['(B/T)'] < 0.5 ]
BT80l = sample[ sample['(B/T)'] >= 0.5 ]
BT80 = BT80l[ BT80l['(B/T)'] < 0.8 ]
BT100 = sample[ sample['(B/T)'] >= 0.8 ]
tables = [BT20,BT50,BT80,BT100]
tablenames = ['BT20','BT50','BT80','BT100']


# In[25]:


for table,name in zip(tables,tablenames):
    
    t2 = table[table['logMbd'] >= 10]
    
    fig, axs = plt.subplots(1,2,figsize=(10,5),sharey=True,sharex=True)
    fig.suptitle(r'$\mathrm{Bulge \ fraction: \ %s}$'%(name),fontsize=16)
    x = np.linspace(1,15,10)

    inds = np.where(np.isfinite(t2['logMbh_ssig']) == True)
    axs[0].hist2d(t2['logMbh_ssig'][inds],t2['logMbh_fp'][inds],bins=[80,80],cmap='bone_r')
    axs[0].set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ \sigma_{inf}$')
    axs[1].hist2d(t2['logMbh_sbulge'][inds],t2['logMbh_fp'][inds],bins=[30,30],cmap='bone_r')
    axs[1].set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ M_{bulge}$')
    for ax in axs:
        ax.set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ M_{bulge} - \sigma_{inf}$')
        ax.plot(x,x,'--',color='k')
        ax.grid(alpha=0.5,which='both')
    fig.savefig('fundamental_plane_%s.png'%(name))

