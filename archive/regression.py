#!/usr/bin/env python
# coding: utf-8

# In[2]:


from astropy.table import Table
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


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


# In[12]:


def Deming(x,y,e_x,e_y):
    n = len(x)
    delta = np.var(e_y)/np.var(e_x)
    xbar = np.mean(x)
    ybar = np.mean(y)
    sxx =  (1/n)*sum([(xi-xbar)**2 for xi in x]) 
    sxy =  (1/n)*sum([(xi-xbar)*(y[i]-ybar) for i,xi in enumerate(x)])
    syy = (1/n)*sum([(yi-ybar)**2 for yi in y])
    beta1 = (syy-delta*sxx+np.sqrt((syy-delta*sxx)**2+4*delta*sxy**2))/(2*sxy)
    beta0 = ybar-beta1*xbar
    xstar = [xi + beta1/(beta1**2+delta)*(y[i]-beta0-beta1*xi) for i,xi in enumerate(x)]
    return xstar,beta0,beta1


# In[23]:


data = Table.read("tables/1000_galaxies.dat",
        format="ascii")


# In[24]:


xvals = data['logMbh_b']
yvals = data['logMbh_s']
yerr = data['err_logMbh_s']
xerr = data['err_logMbh_b']


# In[59]:


plt.figure(dpi=200)
plt.errorbar(xvals,yvals,xerr=xerr,yerr=yerr,marker='.',linewidth=0,elinewidth=.3,color='k',alpha=0.2)
plt.xlabel(r'$\log(M_{\bullet,Bulge}) (M_{\odot})$')
plt.ylabel(r'$\log(M_{\bullet,\sigma}) (M_{\odot})$')
plt.title(r'$\mathrm{1000 \ most \ massive } \  M_{Bulge} \ \mathrm{galaxies}$');


# In[89]:


(xstar,b0,b1) = Deming(xvals,yvals,xerr,yerr)
xrange = np.linspace(2,13,100)
plt.figure(dpi=200)
plt.errorbar(xvals,yvals,xerr=xerr,yerr=yerr,marker='.',linewidth=0,elinewidth=.3,color='k',alpha=0.2)
plt.xlabel(r'$\log(M_{\bullet,Bulge}) (M_{\odot})$')
plt.ylabel(r'$\log(M_{\bullet,\sigma}) (M_{\odot})$')
plt.title(r'$\mathrm{1000 \ most \ massive } \  M_{Bulge} \ \mathrm{galaxies}, \ \mathrm{with \ errorbars}$')
plt.plot(xrange,[b1*x+b0 for x in xrange],linewidth=0.8,color='xkcd:red',zorder=4,label=r'$y = %.2f + %.2fx$'%(b0,b1))
plt.plot(xrange,xrange,linestyle='--',linewidth=0.8,color='xkcd:blue',zorder=3,label=r'$y = x$')
plt.xlim(3,11)
plt.ylim(0,11)
plt.legend()
plt.grid(alpha=0.5)
plt.savefig('plots/final1000galaxies.pdf');


# In[55]:


sample = Table.read("tables/1000_galaxies_noerrbars.dat",
        format="ascii")


# In[65]:


xvals1 = sample['logMbh_b']
yvals1 = sample['logMbh_s']


# In[70]:


from sklearn.linear_model import LinearRegression
model = LinearRegression().fit(xvals1.reshape((-1, 1)), yvals1.reshape((-1, 1)))


# In[90]:


xrange1 = np.linspace(2,13,100)
plt.figure(dpi=200)
plt.scatter(xvals1,yvals1,marker='.',color='k',alpha=0.2)
plt.xlabel(r'$\log(M_{\bullet,Bulge}) (M_{\odot})$')
plt.ylabel(r'$\log(M_{\bullet,\sigma}) (M_{\odot})$')
plt.title(r'$\mathrm{1000 \ most \ massive } \  M_{Bulge} \ \mathrm{galaxies}, \ \mathrm{without \ errorbars}$')
plt.plot(xrange1,[model.coef_[0][0]*x+model.intercept_[0] for x in xrange1],linewidth=0.8,color='xkcd:red',zorder=4,label=r'$y = %.2f + %.2fx$'%(model.intercept_[0],model.coef_[0][0]))
plt.plot(xrange1,xrange1,linestyle='--',linewidth=0.8,color='xkcd:blue',zorder=3,label=r'$y = x$')
plt.xlim(3,11)
plt.ylim(0,11)
plt.legend()
plt.grid(alpha=0.5)
plt.savefig('plots/1000galaxiesnoerr.pdf');

