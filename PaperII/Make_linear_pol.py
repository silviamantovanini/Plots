#!/home/silvia/PhD/bin/python3
# Make a total linear polarisation spass image

from astropy.io import fits
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from astropy.visualization import PercentileInterval

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import wcs
import numpy as np
from matplotlib.ticker import ScalarFormatter

plt.rcParams.update({
    "font.size": 11,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54 

def normalize(arr, vmin, vmax):
    nor = (arr - vmin) / (vmax - vmin)
    nor[np.where(nor<0.0)] = 0.0
    nor[np.where(nor>1.0)] = 1.0
    return nor
    
hq = fits.open("SPASS_q_rg_mod_nb.fits")
q = hq[0].data
hu = fits.open("SPASS_u_rg_mod_nb.fits")
u = hu[0].data
hi = fits.open("SPASS_model_subtracted_nb.fits")
i = hi[0].data

q = q - 0.01
u = u - 0.01

l = np.sqrt(q**2 + u**2)
p = l / i
p[i<0.02] = np.nan

hq[0].data = l
hq.writeto("SPASS_l_rg_mod.fits", overwrite=True)
hq[0].data = p
hq.writeto("SPASS_p_rg_mod.fits", overwrite=True)


#Print total linear polarization
hdu=fits.open("SPASS_l_rg_mod.fits")

w=wcs.WCS(hdu[0].header)
interval=PercentileInterval(99.).get_limits(hdu[0].data)
image=normalize(hdu[0].data, *interval)

fig = plt.figure(figsize=(25*cm,25*cm))

ax = fig.add_axes([0.1,0.11,0.45, 0.45], projection = w)
cbaxes = fig.add_axes([0.56, 0.11, 0.01, 0.45])

img = ax.imshow(image, origin="lower", cmap="inferno")

cb = plt.colorbar(img, cax = cbaxes, orientation="vertical")
#cb.ax.tick_params(labelsize=5)
cb.set_label('Flux density / mJy beam$^{-1}$')

lon = ax.coords['ra']
lon.set_major_formatter('hh:mm')
lon.set_axislabel("Right Ascension (J2000)")

lat = ax.coords['dec']
lat.set_axislabel("Declination (J2000)")
lat.set_major_formatter('dd')
    
ax.add_patch(Ellipse((233.0323196, -60.8420575),
                     width=1.9, height=3.1, angle=103.0,
                     facecolor='none',
                     edgecolor='white',
                     ls = "dashed",
                     linewidth=1.5,
                     alpha = 0.9,
                     transform=ax.get_transform('fk5')))
                     
ax.add_patch(Ellipse((233.4770904, -61.2639339),
                     width=1.0, height=1.7, angle=120.0,
                     facecolor='none',
                     edgecolor='green',
                     ls = "dashed",
                     linewidth=1.5,
                     alpha = 0.9,
                     transform=ax.get_transform('fk5')))

fig.savefig("TotLinearPol.pdf", bbox_inches="tight")
fig.savefig("TotLinearPol.png", bbox_inches="tight")


################################################################################################################################
################################################################################################################################
#Print Q and U stokes images
vmin_pol, vmax_pol = -0.05e3, 0.09e3
fig2 = plt.figure(figsize=(10,10))

hdu_q = fits.open("SPASS_q_rg_mod_nb.fits")
hdu_u = fits.open("SPASS_u_rg_mod_nb.fits")

wq = wcs.WCS(hdu_q[0].header)
wu = wcs.WCS(hdu_q[0].header)

ax_q = fig2.add_axes([0.1,0.11,0.37,0.35], projection = wq)
ax_u = fig2.add_axes([0.52,0.11,0.37,0.35], projection = wu)
cbaxes_pol = fig2.add_axes([0.90, 0.11, 0.01, 0.35])

img_q = ax_q.imshow(1.e3*hdu_q[0].data, origin="lower", vmin=vmin_pol, vmax=vmax_pol, cmap="inferno")
img_u = ax_u.imshow(1.e3*hdu_u[0].data, origin="lower", vmin=vmin_pol, vmax=vmax_pol, cmap="inferno")

cb_pol = plt.colorbar(img_q, cax = cbaxes_pol, orientation="vertical")
cb_pol.set_label('Flux density / mJy beam$^{-1}$')

for ax in ax_q, ax_u:
    lon = ax.coords['ra']
    lon.set_axislabel("Right Ascension (J2000)")
    lon.set_major_formatter('hh:mm')

lat = ax_q.coords['dec']
lat.set_axislabel("Declination (J2000)")
lat.set_major_formatter('dd')

lat = ax_u.coords['dec']
lat.set_ticklabel_visible(False)
lat.set_ticks_visible(False)

for ax in ax_q, ax_u:
    ax.add_patch(Ellipse((233.0323196, -60.8420575),
                         width=1.9, height=3.1, angle=103.0,
                         facecolor='none',
                         edgecolor='white',
                         ls = "dashed",
                         linewidth=1.5,
                         alpha = 0.9,
                         transform=ax.get_transform('fk5')))

fig2.savefig("Q&U.pdf", bbox_inches="tight")
fig2.savefig("Q&U.png", bbox_inches="tight")