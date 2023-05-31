#!/usr/bin/env python
# coding: utf-8
from astropy.table import Table, join, vstack, ColumnGroups, Column, unique
from astropy.coordinates import SkyCoord, match_coordinates_3d
import numpy as np
import astropy.units as u
from astropy import constants as const
from astropy.io import ascii
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--agnfile",type=str,
                    help="agn input file")
parser.add_argument("--samplefile",type=str,
                    help="sample input file")
parser.add_argument("--o",type=str,
                    help="output file")
args = parser.parse_args()

agn = ascii.read(args.agnfile)
obj = ascii.read(args.samplefile)

agnskycoord = SkyCoord(ra=agn['RAdeg'], dec=agn['DEdeg'])
objskycoord = SkyCoord(ra=np.asarray(obj['_RA'],dtype=float)*u.degree, dec=np.asarray(obj['_DE'],dtype=float)*u.degree)

max_sep = 1 * u.arcsec
idx, d2d, d3d = agnskycoord.match_to_catalog_3d(objskycoord)
sep_constraint = d2d < max_sep
obj_matches = agnskycoord[sep_constraint]
agn_matches = objskycoord[idx[sep_constraint]]
broadlines = agn[sep_constraint]
scaling = obj[idx[sep_constraint]]

columns = [broadlines['ID'],broadlines['FWHM-BHb'],broadlines['e_FWHM-BHb'],broadlines['logL5100'],broadlines['e_logL5100']]
for column in columns:
    scaling.add_column(column)
    
ascii.write(scaling, args.o, overwrite=True)