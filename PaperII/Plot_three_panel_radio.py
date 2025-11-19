#!/home/silvia/PhD/bin/python3

__author__ = "Natasha Hurley-Walker"
__date__ = "30/08/2021"

import os
import sys
import shutil
import glob

import matplotlib
matplotlib.use('Agg') # So does not use display -- only good if just making plots
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Ellipse, Rectangle
from astropy.visualization import PercentileInterval, AsinhStretch

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import wcs
import numpy as np

plt.rcParams.update({
    "font.size": 10,
    "font.sans-serif": ["Helvetica"]})
#    "font.sans-serif": ["Arial"]})
#    "font.family": "sans-serif",

cm = 1/2.54  # centimeters in inches

def normalize(arr, vmin, vmax):
    nor = (arr - vmin) / (vmax - vmin)
    nor[np.where(nor<0.0)] = 0.0
    nor[np.where(nor>1.0)] = 1.0
    return nor
    
bkg = 0.9

vmin_chipass, vmax_chipass = 13.9, 17.4
vmin_spass, vmax_spass = 0.9, 1.e3*(1.9 - bkg)

fig = plt.figure(figsize=(25*cm,25*cm))

hdu_chipass = fits.open("CHIPASS_rg_subtracted_nb.fits")
hdu_spass = fits.open("SPASS_model_subtracted_nb.fits")

d_chipass = hdu_chipass[0].data
d_spass = hdu_spass[0].data

w_chipass = wcs.WCS(hdu_chipass[0].header)
w_spass = wcs.WCS(hdu_spass[0].header)

ax_chipass = fig.add_axes([0.1,0.11,0.37,0.35], projection = w_chipass)
ax_spass = fig.add_axes([0.52,0.11,0.37,0.35], projection = w_spass)

cbaxes_chipass = fig.add_axes([0.47, 0.11, 0.01, 0.35])
cbaxes_spass = fig.add_axes([0.89, 0.11, 0.01, 0.35])

img_chipass = ax_chipass.imshow(d_chipass, origin="lower", vmin=vmin_chipass, vmax=vmax_chipass, cmap="inferno")
img_spass = ax_spass.imshow(1.e3*d_spass - 1.e3*bkg, origin="lower", vmin=vmin_spass, vmax=vmax_spass, cmap="inferno")

cb_chipass = plt.colorbar(img_chipass, cax = cbaxes_chipass, orientation="vertical")
#cb_chipass.ax.tick_params(labelsize=8)
cb_spass = plt.colorbar(img_spass, cax = cbaxes_spass, orientation="vertical")
#cb_spass.ax.tick_params(labelsize=8)

cb_chipass.set_label(' ')
cb_spass.set_label('Flux density / mJy beam$^{-1}$')

for ax in ax_spass, ax_chipass:
    lon = ax.coords['ra']
    lon.set_major_formatter('hh:mm')
    lon.set_axislabel("Right Ascension (J2000)")

lat = ax_chipass.coords['dec']
lat.set_axislabel("Declination (J2000)")
lat.set_major_formatter('dd')

lat = ax_spass.coords['dec']
lat.set_ticklabel_visible(False)
lat.set_ticks_visible(False)
    
for ax in ax_chipass, ax_spass:
    ax.add_patch(Ellipse((233.0323196, -60.8420575),
                         width=1.9, height=3.1, angle=103.0,
                         facecolor='none',
                         edgecolor='white',
                         ls = "dashed",
                         linewidth=1.5,
                         alpha = 0.9,
                         transform=ax.get_transform('fk5')))
                       
    boxwidth = 75
    ax.add_patch(Rectangle((22, 18), width=boxwidth, height=boxwidth, facecolor = 'white', edgecolor='black'))
                  
ax_chipass.add_patch(Ellipse((236.4797049, -62.3716549),
                     width=(14.4/60.), height=(8.9/60.), angle=0.0,
                     facecolor='black',
                     edgecolor='black',
                     transform=ax_chipass.get_transform('fk5')))
                      
ax_spass.add_patch(Ellipse((236.4797049, -62.3716549),
                     width=(8.9/60.), height=(8.9/60.), angle=0.0,
                     facecolor='black',
                     edgecolor='black',
                     transform=ax_spass.get_transform('fk5')))
                       
fig.savefig("CHIPASS_SPASS.pdf", bbox_inches="tight")
fig.savefig("CHIPASS_SPASS.png", bbox_inches="tight")



###############################################################################################################################
###############################################################################################################################

#Print GPM, MGPS & WISE
hdu_mwa=fits.open("SNR_G321.3-3.9_MWA_rg.fits")
hdu_mgps=fits.open("SNR_G321.3-3.9_SUMSS.fits")
hdu_wise12 = fits.open("SNR_G321.3-3.9_WISE12.fits")
hdu_wise22 = fits.open("SNR_G321.3-3.9_WISE22.fits")
hdu_wise34 = fits.open("SNR_G321.3-3.9_WISE34.fits")

pct = 95.0
interval = PercentileInterval(pct)

w_wise = wcs.WCS(hdu_wise22[0].header)
i = interval.get_limits(hdu_wise22[0].data)
r = normalize(hdu_wise22[0].data, *i)
i = interval.get_limits(hdu_wise12[0].data)
g = normalize(hdu_wise12[0].data, *i)
i = interval.get_limits(hdu_wise34[0].data)
b = normalize(hdu_wise34[0].data, *i)

wise_rgb = np.dstack([r,g,b])

image_mwa=hdu_mwa[0].data
image_mgps=hdu_mgps[0].data

w_mwa=wcs.WCS(hdu_mwa[0].header)
w_mgps=wcs.WCS(hdu_mgps[0].header)

vmin_mwa, vmax_mwa = -0.017e3, 0.03e3
vmin_mgps, vmax_mgps = -2, 4.

fig_t = plt.figure(figsize=(23*cm,5*cm))

ax_mwa = fig_t.add_axes([0.1,0.11,0.25,1.0], projection = w_mwa)
ax_mgps = fig_t.add_axes([0.41,0.11,0.25,1.0], projection = w_mgps)
ax_wise = fig_t.add_axes([0.77,0.11,0.25,1.0], projection = w_wise)

cbaxes_mwa = fig_t.add_axes([0.35, 0.11, 0.01, 1.0])
cbaxes_mgps = fig_t.add_axes([0.66, 0.11, 0.01, 1.0])

img_mwa = ax_mwa.imshow(1e3*image_mwa, origin="lower", vmin=vmin_mwa, vmax=vmax_mwa, cmap="inferno")
img_mgps = ax_mgps.imshow(1e3*image_mgps, origin="lower", vmin=vmin_mgps, vmax=vmax_mgps, cmap="inferno")
img_wise = ax_wise.imshow(wise_rgb, origin="lower")

cb_mwa = plt.colorbar(img_mwa, cax = cbaxes_mwa, orientation="vertical")
cb_mgps = plt.colorbar(img_mgps, cax = cbaxes_mgps, orientation="vertical")

cb_mwa.set_label(' ')
cb_mgps.set_label('Flux density / mJy beam$^{-1}$')

for ax in ax_mwa, ax_mgps, ax_wise:
    lon = ax.coords['ra']
    lon.set_major_formatter('hh:mm')
    lon.set_axislabel("Right Ascension (J2000)")

lat = ax_mwa.coords['dec']
lat.set_axislabel("Declination (J2000)")
lat.set_major_formatter('dd')

lat = ax_mgps.coords['dec']
lat.set_ticklabel_visible(False)
lat.set_ticks_visible(False)

lat = ax_wise.coords['dec']
lat.set_ticklabel_visible(False)
lat.set_ticks_visible(False)
    
for ax in ax_mwa, ax_mgps, ax_wise:
    ax.add_patch(Ellipse((233.0323196, -60.8420575),
                         width=1.9, height=2.6, angle=105.0,
                         facecolor='none',
                         edgecolor='white',
                         ls = "dashed",
                         linewidth=1.5,
                         alpha = 0.9,
                         transform=ax.get_transform('fk5')))

fig_t.savefig("GPM-MGPS.pdf", bbox_inches="tight")
fig_t.savefig("GPM-MGPS.png", bbox_inches="tight")


#width=1.6, height=2.3, angle=105.0,


###############################################################################################################################
###############################################################################################################################

#Print RACS, EMU
hdu_racs=fits.open("SNR_G321.3-3.9_RACS_rg_nb_cut.fits")
hdu_emu=fits.open("SNR_G321.3-3.9_EMU_cut_rg_nb.fits")

w_racs=wcs.WCS(hdu_racs[0].header, naxis=2)
w_emu=wcs.WCS(hdu_emu[0].header, naxis=2)

data_racs=np.squeeze(hdu_racs[0].data)
interval=PercentileInterval(97.).get_limits(data_racs)
image_racs=normalize(data_racs, *interval)

data_emu=np.squeeze(hdu_emu[0].data)
interval_emu=PercentileInterval(97.).get_limits(data_emu)
image_emu=normalize(data_emu, *interval_emu)

fig_r = plt.figure(figsize=(25*cm,25*cm))

ax_racs = fig_r.add_axes([0.1,0.11,0.37,0.35], projection = w_racs)
ax_emu = fig_r.add_axes([0.52,0.11,0.37,0.35], projection = w_emu)
cbaxes_racs = fig_r.add_axes([0.47, 0.11, 0.01, 0.35])
cbaxes_emu = fig_r.add_axes([0.89, 0.11, 0.01, 0.35])

img_racs = ax_racs.imshow(image_racs, origin="lower", cmap="inferno")
img_emu = ax_emu.imshow(image_emu, origin="lower", cmap="inferno")

cb_racs = plt.colorbar(img_racs, cax = cbaxes_racs, orientation="vertical")
cb_emu = plt.colorbar(img_emu, cax = cbaxes_emu, orientation="vertical")

cb_racs.set_label(' ')
cb_emu.set_label('Flux density / mJy beam$^{-1}$')

for ax in ax_racs, ax_emu:
    lon = ax.coords['ra']
    lon.set_major_formatter('hh:mm')
    lon.set_axislabel("Right Ascension (J2000)")

lat = ax_racs.coords['dec']
lat.set_axislabel("Declination (J2000)")
lat.set_major_formatter('dd')

lat = ax_emu.coords['dec']
lat.set_ticklabel_visible(False)
lat.set_ticks_visible(False)
    
for ax in ax_racs, ax_emu:
    ax.add_patch(Ellipse((233.0323196, -60.8420575),
                         width=1.9, height=3.1, angle=103.0,
                         facecolor='none',
                         edgecolor='white',
                         ls = "dashed",
                         linewidth=1.5,
                         alpha = 0.9,
                         transform=ax.get_transform('fk5')))

fig_r.savefig("RACS-EMU.pdf", bbox_inches="tight")
fig_r.savefig("RACS-EMU.png", bbox_inches="tight")


###############################################################################################################################
###############################################################################################################################

#Print HALPHA
bkg_h=15
vmin_h, vmax_h = 10, 190

hdu_h=fits.open("Halpha_rg_nb.fits")
w_h=wcs.WCS(hdu_h[0].header)
d_h = hdu_h[0].data

fig_h = plt.figure(figsize=(14,14))

ax_h = fig_h.add_axes([0.1,0.11,0.45, 0.45], projection = w_h)
cbaxes_h = fig_h.add_axes([0.56, 0.11, 0.01, 0.45])

img_h = ax_h.imshow(d_h, vmin=vmin_h, vmax=vmax_h, origin="lower", cmap="gnuplot2_r")	#twilight

cb_h = plt.colorbar(img_h, cax = cbaxes_h, orientation="vertical")
cb_h.ax.tick_params(labelsize=5)
cb_h.set_label('Flux density / mJy beam$^{-1}$')

lon = ax_h.coords['ra']
lon.set_major_formatter('hh:mm')
lon.set_axislabel("Right Ascension (J2000)")

lat = ax_h.coords['dec']
lat.set_axislabel("Declination (J2000)")
lat.set_major_formatter('dd')
    
ax_h.add_patch(Ellipse((233.0323196, -60.8420575),
                     width=1.9, height=3.1, angle=103.0,
                     facecolor='none',
                     edgecolor='cyan',
                     ls = "dashed",
                     linewidth=1.5,
                     alpha = 0.9,
                     transform=ax_h.get_transform('fk5')))

fig_h.savefig("HALPHA.pdf", bbox_inches="tight")
fig_h.savefig("HALPHA.png", bbox_inches="tight")