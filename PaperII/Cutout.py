#!/home/silvia/PhD/bin/python3
#Create a cutout of a source.

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS  
import numpy as np
import sys
    

file1=sys.argv[1]
output=sys.argv[2]

hdu=fits.open(file1)
hdu2=fits.open('/home/silvia/DS9_file/Summary_Plot/Candidates/G321.3-3.9/SpectralIndex/FIts_file/SNR_G321.3-3.9_SUMSS.fits')

image=hdu[0].data
w=WCS(hdu[0].header)
boxsize=2*u.deg
coords=SkyCoord('15:32:42.4509', '-60:51:49.735', unit=(u.hourangle, u.deg), frame='fk5')
cutout=Cutout2D(image, coords, boxsize, wcs=w)

hdu[0].data=cutout.data
hdu[0].header.update(cutout.wcs.to_header())

if 'PC1_1' in dict(hdu[0].header).keys():
	hdu[0].header['CD1_1'] = hdu[0].header['PC1_1']
	hdu[0].header['CD2_2'] = hdu[0].header['PC2_2']
	hdu[0].header['CDELT1'] = hdu[0].header['CD1_1']
	hdu[0].header['CDELT2'] = hdu[0].header['CD2_2']
	hdu[0].header.pop('PC1_1')
	hdu[0].header.pop('PC2_2')
	hdu[0].header.pop('CDELT1')
	hdu[0].header.pop('CDELT2')
        
hdu[0].writeto(output, overwrite=True)