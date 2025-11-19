#!/home/silvia/PhD/bin/python3
#Reduce axes in Racs header.

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS  
import numpy as np
import sys
    
def add_beam(fitsfile):
	hdu = fits.open(fitsfile)
	bpa = 8.0
	bmaj = 2.2222e-03
	bmin = 2.2222e-03
	hdu[0].header.set("BMAJ", bmaj)
	hdu[0].header.set("BMIN", bmin)
	hdu[0].header.set("BPA", bpa)
	hdu.writeto(fitsfile, overwrite=True)

inputfile=sys.argv[1]
output=sys.argv[2]

hdu=fits.open(inputfile)

w=WCS(hdu[0].header, naxis=2)
data=hdu[0].data

header=w.to_header()
hdu=fits.PrimaryHDU(header=header)
hdu.data=data

hdu.writeto(output, overwrite=True)

add_beam(output)