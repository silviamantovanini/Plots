#!/home/silvia/anaconda3/envs/racs/bin/python3

#Plot GPM before and after the subtraction of the sources with WISE contours.

import os
import sys
import shutil
import glob

import matplotlib
matplotlib.use('Agg') # So does not use display -- only good if just making plots
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Circle, Rectangle

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import wcs
from astropy.visualization import make_lupton_rgb
from astropy.visualization import PercentileInterval
import numpy as np

from matplotlib.colors import LogNorm

plt.rcParams.update({
    "font.size": 8,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54  # centimeters in inches

def normalize(arr, vmin, vmax):
    nor = (arr - vmin) / (vmax - vmin)
    nor[np.where(nor<0.0)] = 0.0
    nor[np.where(nor>1.0)] = 1.0
    return nor
    
vmin_gpm, vmax_gpm = -0.02*1e3, 1e3*0.07
#vmin_radio, vmax_radio = -0.002*1e3, 1e3*0.008
vmin_radio, vmax_radio = -0.0002*1e3, 1e3*0.00045

#Load images and define intervals.
gpm=sys.argv[1]
wise12=sys.argv[2]
wise22=sys.argv[3]
wise34=sys.argv[4]
radio=sys.argv[5]

hdu_gpm = fits.open(gpm)
hdu_wise12 = fits.open(wise12)
hdu_wise22 = fits.open(wise22)
hdu_wise34 = fits.open(wise34)
hdu_radio = fits.open(radio)

pct = 95.0
interval = PercentileInterval(pct)

w_gpm = wcs.WCS(hdu_gpm[0].header)
w_wise = wcs.WCS(hdu_wise22[0].header)
w_radio = wcs.WCS(hdu_radio[0].header, naxis=2)

image_gpm=hdu_gpm[0].data
image_radio=hdu_radio[0].data
data_radio=np.squeeze(image_radio)

i = interval.get_limits(hdu_wise22[0].data)
r = normalize(hdu_wise22[0].data, *i)
i = interval.get_limits(hdu_wise12[0].data.data)
g = normalize(hdu_wise12[0].data, *i)
i = interval.get_limits(hdu_wise34[0].data.data)
b = normalize(hdu_wise34[0].data, *i)

wise_rgb = np.dstack([r,g,b])

#Now print them.
fig = plt.figure(figsize=(20*cm,5*cm))

ax_gpm = fig.add_axes([0.1,0.11,0.25,1.0], projection = w_gpm)
ax_radio = fig.add_axes([0.42,0.11,0.25,1.0], projection = w_radio)
ax_wise = fig.add_axes([0.77,0.11,0.25,1.0], projection = w_wise)

cbaxes_gpm = fig.add_axes([0.36, 0.11, 0.01, 1.0])
cbaxes_radio = fig.add_axes([0.68, 0.11, 0.01, 1.0])

img_gpm = ax_gpm.imshow(1e3*image_gpm, origin="lower", vmin=vmin_gpm, vmax=vmax_gpm, cmap="cubehelix")
img_wise = ax_wise.imshow(wise_rgb, origin="lower")
img_radio = ax_radio.imshow(1e3*data_radio, origin="lower", vmin=vmin_radio, vmax=vmax_radio, cmap="cubehelix")

cb_gpm = plt.colorbar(img_gpm, cax = cbaxes_gpm, orientation="vertical")
plt.setp(plt.getp(cb_gpm.ax.axes, 'yticklabels'))
cb_radio = plt.colorbar(img_radio, cax = cbaxes_radio, orientation="vertical")
plt.setp(plt.getp(cb_radio.ax.axes, 'yticklabels'))

cb_gpm.set_label(' ')
cb_radio.set_label('Flux density / mJy beam$^{-1}$')

for ax in ax_gpm, ax_wise, ax_radio:
    lon = ax.coords['ra']
    lon.set_major_formatter('hh:mm')
    lon.set_axislabel("Right Ascension")

lat = ax_gpm.coords['dec']
lat.set_axislabel("Declination")
lat.set_major_formatter('dd:mm')
    
for ax in ax_wise, ax_radio:
    lat = ax.coords['dec']
    lat.set_ticklabel_visible(False)
    lat.set_ticks_visible(False)

#Add ellipse around the shell just in gpm image - coordinates are fk5 in degrees
for ax in ax_gpm, ax_wise, ax_radio:
    ax.add_patch(Ellipse((274.2906556, -14.6538405),
                         width=0.2, height=0.2, angle=360.0,
                         facecolor='none',
                         edgecolor='white',
                         ls = "dashed",
                         alpha = 0.9,
                         transform=ax.get_transform('fk5')))
       
for ax in ax_gpm, ax_radio:                       
    boxwidth = 92
    ax.add_patch(Rectangle((22, 18), width=boxwidth, height=boxwidth, facecolor = 'white', edgecolor='black'))
    
ax_gpm.add_patch(Ellipse((274.52, -14.88),
                     width=2.16263402253E-02, height=1.65142789483E-02, angle=90.0,
                     facecolor='black',
                     edgecolor='black',
                     transform=ax_gpm.get_transform('fk5')))
                      
ax_radio.add_patch(Ellipse((274.52, -14.88),
                     width=2.2222e-03, height=2.2222e-03, angle=90.0,
                     facecolor='black',
                     edgecolor='black',
                     transform=ax_radio.get_transform('fk5')))
                     
ax_gpm.annotate(
    '',
    xy=(145, 110), xycoords='figure points',
    xytext=(165, 120), textcoords='figure points',
    arrowprops=dict(arrowstyle="->", color='cyan'))

ax_radio.annotate(
    '',
    xy=(325, 110), xycoords='figure points',
    
    xytext=(345, 120), textcoords='figure points',
    arrowprops=dict(arrowstyle="->", color='cyan'))    
                     
fig.savefig("SNR_G16.0+0.7.pdf", bbox_inches="tight", dpi=150)