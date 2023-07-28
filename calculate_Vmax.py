from astropy.cosmology import WMAP9 as cosmo
import scipy.integrate as integrate
from astropy import constants as const
from astropy.table import Table
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--i",type=str,
                    help="input file")
parser.add_argument("--o",type=str,
                    help="output file")
args = parser.parse_args()
input = args.i
output = args.o

data = Table.read(input,
        format="ascii")

def Vmax(galaxy,data):
    skyfrac = 8200
    angular = skyfrac/41253
    zmin = max(0.02,data['zmin'][galaxy])
    zmax = min(0.2,data['zmax'][galaxy])
    zgalaxy = data['z'][galaxy]
    Vint = integrate.quad(lambda z: (cosmo.angular_diameter_distance(z).value**2*const.c.to('Mpc/s').value)/(cosmo.H(z).to('1/s').value*(1+z)),0.02,zgalaxy)
    Vmaxint = integrate.quad(lambda z: (cosmo.angular_diameter_distance(z).value**2*const.c.to('Mpc/s').value)/(cosmo.H(z).to('1/s').value*(1+z)),zmin,zmax)
    V = angular*Vint[0]
    Vmax = angular*Vmaxint[0]
    return Vmax**-1, V/Vmax

V_max = [Vmax(g,data)[0] for g in range(len(data))]
V_Vmax = [Vmax(g,data)[1] for g in range(len(data))]

arr = {'objID' : data['objID'],
       '1/Vmax' : V_max,
       'V/Vmax' : V_Vmax
      }
table = Table(arr)

ascii.write(table, output, overwrite=True) 