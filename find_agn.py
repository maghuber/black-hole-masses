#!/usr/bin/env python
# coding: utf-8

# In[17]:


from astropy.table import Table, join, vstack, ColumnGroups, Column, unique
from astropy.coordinates import SkyCoord, match_coordinates_3d
import numpy as np
import astropy.units as u
from astropy import constants as const
from astropy.io import ascii


# In[63]:


agn1 = ascii.read('agn.txt')
obj1 = ascii.read('550K Sample with MBH')
agnradec = ascii.read('agnradec.txt')
objradec = ascii.read('objradec.txt')
obj = join(obj1,objradec,keys='objID')
agn = join(agn1,agnradec,keys='ID')


# In[61]:


agnskycoord = SkyCoord(ra=agn['RAdeg'], dec=agn['DEdeg'])
objskycoord = SkyCoord(ra=obj['ra']*u.degree, dec=obj['dec']*u.degree)


# In[ ]:


max_sep = 1 * u.arcsec
idx, d2d, d3d = agnskycoord.match_to_catalog_3d(objskycoord)
sep_constraint = d2d < max_sep
obj_matches = agnskycoord[sep_constraint]
agn_matches = objskycoord[idx[sep_constraint]]
broadlines = agn[sep_constraint]
scaling = obj[idx[sep_constraint]]


# In[ ]:


agnsample = Table()
columns = [scaling['objID'],broadlines['ID'],scaling['(B/T)'],scaling['color_type'],scaling['logMbd'],scaling['logMbh_bulge'],scaling['logMbh_siginf'],broadlines['logMBHa'],broadlines['logMBHb'],broadlines['logMBH']]
for column in columns:
    agnsample.add_column(column)


# In[ ]:


mask2 = np.isfinite(agnsample['logMbh_siginf']) == True
agnsample = agnsample[mask2]


# In[ ]:


len(agnsample)


# In[10]:


ascii.write(agnsample, 'AGN Sample', overwrite=True)

