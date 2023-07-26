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
from scipy.stats import gaussian_kde


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


# In[4]:


agnsample = Table.read("AGN Sample",
        format="ascii")


# In[5]:


BT20 = agnsample[ agnsample['(B/T)'] < 0.2 ]
BT50l = agnsample[ agnsample['(B/T)'] >= 0.2 ]
BT50 = BT50l[ BT50l['(B/T)'] < 0.5 ]
BT80l = agnsample[ agnsample['(B/T)'] >= 0.5 ]
BT80 = BT80l[ BT80l['(B/T)'] < 0.8 ]
BT100 = agnsample[ agnsample['(B/T)'] >= 0.8 ]


# In[6]:


passive = agnsample[ agnsample['color_type'] == 'passive' ]
active = agnsample[ agnsample['color_type'] == 'active' ]


# In[7]:


early = join(BT80l,passive)
lowbulge = agnsample[ agnsample['(B/T)'] <= 0.5 ]
late = join(lowbulge,active)


# In[9]:


tables= [early,late]
tablenames = [r'$\mathrm{B/T \ greater \ than \ 0.5 \ and \ passive}$',r'$\mathrm{B/T \ less \ than \ 0.5 \ and \ active}$']
figure = ['blearly','bllate']
for t,name,f in zip(tables,tablenames,figure):

    sigma = np.asarray(t['logMBHa']).std(axis=0)
    
    fig, axs = plt.subplots(1,2,figsize=(10,5),sharey=True,sharex=True,dpi=300)
    fig.suptitle('%s'%(name),fontsize=22)
    x = np.linspace(1,15,10)

    xy = np.vstack([t['logMbh_siginf'],t['logMBHa']])
    z = gaussian_kde(xy)(xy)
    
    xy2 = np.vstack([t['logMbh_bulge'],t['logMBHa']])
    z2 = gaussian_kde(xy2)(xy2)
    
    axs[0].scatter(t['logMBHa'],t['logMbh_siginf'],c=z,cmap='BuPu',s=10)
    axs[0].set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ \sigma_{inf} (M_\odot)$')
    axs[1].scatter(t['logMBHa'],t['logMbh_bulge'],c=z2,cmap='BuPu',s=10)
    axs[1].set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ M_{bulge} (M_\odot)$')
    for ax in axs:
        ax.fill_between(x, x+sigma, x-sigma, facecolor='gray', alpha=0.5)
        ax.set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from \ broad \ line \ emission}$')
        ax.plot(x,x,'--',color='k')
        ax.grid(alpha=0.5,which='both')
        ax.plot(x,x,'--',color='k')
        ax.set_xlim(4,12)
        ax.set_ylim(5,10)
        ax.grid(alpha=0.5,which='both')
    fig.savefig('%s.png'%(f))


# In[27]:


fig, axs = plt.subplots(1,2,figsize=(10,5),sharey=True,sharex=True,dpi=300)
fig.suptitle(r'$\mathrm{All \ AGN}$',fontsize=22)
x = np.linspace(1,15,10)

sigma = np.asarray(t['logMBHa']).std(axis=0)

xy = np.vstack([agnsample['logMbh_siginf'],agnsample['logMBHa']])
z = gaussian_kde(xy)(xy)

xy2 = np.vstack([agnsample['logMbh_bulge'],agnsample['logMBHa']])
z2 = gaussian_kde(xy2)(xy2)

axs[0].scatter(agnsample['logMBHa'],agnsample['logMbh_siginf'],c=z,cmap='BuPu',s=10)
axs[0].set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ \sigma_{inf} (M_\odot)$')
axs[1].scatter(agnsample['logMBHa'],agnsample['logMbh_bulge'],c=z2,cmap='BuPu',s=10)
axs[1].set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ M_{bulge} (M_\odot)$')
for ax in axs:
    ax.fill_between(x, x+sigma, x-sigma, facecolor='gray', alpha=0.5)
    ax.set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from \ broad \ line \ emission}$')
    ax.plot(x,x,'--',color='k')
    ax.grid(alpha=0.5,which='both')
    ax.plot(x,x,'--',color='k')
    ax.set_xlim(4,12)
    ax.set_ylim(5,10)
    ax.grid(alpha=0.5,which='both')


# In[43]:


fig, axs = plt.subplots(1,2,figsize=(10,5),sharey=True,sharex=True,dpi=300)
fig.suptitle(r'$\mathrm{All \ AGN}$',fontsize=22)
x = np.linspace(1,15,10)


axs[0].hist2d(agnsample['logMBHa'],agnsample['logMbh_siginf'],bins=[20,70],cmap='BuPu')
axs[0].set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ \sigma_{inf} (M_\odot)$')
axs[1].hist2d(agnsample['logMBHa'],agnsample['logMbh_bulge'],bins=[20,30],cmap='BuPu')
axs[1].set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ M_{bulge} (M_\odot)$')
for ax in axs:
    ax.fill_between(x, x+sigma, x-sigma, facecolor='gray', alpha=0.5)
    ax.set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from \ broad \ line \ emission}$')
    ax.plot(x,x,'--',color='k')
    ax.grid(alpha=0.5,which='both')
    ax.plot(x,x,'--',color='k')
    ax.set_xlim(4,12)
    ax.set_ylim(5,10)
    ax.grid(alpha=0.5,which='both')


# In[37]:


tables = [BT20,BT50,BT80,BT100]
tablenames = ['BT20','BT50','BT80','BT100']
for t,name in zip(tables,tablenames):
    
    
    fig, axs = plt.subplots(1,2,figsize=(14,5),sharey=True,sharex=True)
    fig.suptitle(r'$\mathrm{Bulge \ fraction: \ %s}$'%(name),fontsize=16)
    x = np.linspace(1,15,10)

    xy = np.vstack([t['logMbh_siginf'],t['logMBHa']])
    z = gaussian_kde(xy)(xy)
    
    xy2 = np.vstack([t['logMbh_bulge'],t['logMBHa']])
    z2 = gaussian_kde(xy2)(xy2)
    
    axs[0].scatter(t['logMBHa'],t['logMbh_siginf'],c=z,cmap='BuPu',s=10)
    axs[0].set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ \sigma_{inf} (M_\odot)$')
    axs[1].scatter(t['logMBHa'],t['logMbh_bulge'],c=z2,cmap='BuPu',s=10)
    axs[1].set_ylabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from} \ M_{bulge} (M_\odot)$')
    for ax in axs:
        ax.fill_between(x, x+sigma, x-sigma, facecolor='gray', alpha=0.5)
        ax.set_xlabel(r'$\mathrm{log}(M_{BH}) \ \mathrm{from \ broad \ line \ emission}$')
        ax.plot(x,x,'--',color='k')
        ax.grid(alpha=0.5,which='both')
        ax.plot(x,x,'--',color='k')
        ax.set_xlim(4,12)
        ax.set_ylim(5,10)
        ax.grid(alpha=0.5,which='both')

