#!/home/silvia/anaconda3/envs/racs/bin/python3
#Create a cutout of a source for Meerkat cube.

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS  
import numpy as np
import sys
    

inputfile=sys.argv[1]
output=sys.argv[2]

hdu=fits.open(inputfile)

w=WCS(hdu[0].header, naxis=2)
#image=np.squeeze(hdu[0].data)
image=hdu[0].data[0,0,:,:]

boxsize=(1*u.deg, 1*u.deg)
coords=SkyCoord('18:17:09.7573', '-14:39:13.826', unit=(u.hourangle, u.deg), frame='fk5')
cutout=Cutout2D(image, position=coords, size=boxsize, wcs=w)

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