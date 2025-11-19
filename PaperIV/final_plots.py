#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12
})

cm = 1 / 2.54

FULL_CATALOGUE_PATH = "/home/silvia/Desktop/Papers/Project_4/Pulsars/pulsar2_6_1.txt"
MATCHED_NAMES_PATH = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/FINAL/final_matches.txt"

matches_catalogue = '/home/silvia/Desktop/Papers/Project_4/Pulsars/matched_4FGL_out_FLAGS_artefacts_blocking.fits'
gleamx_catalogue = '/home/silvia/Desktop/Papers/Project_4/Pulsars/merged_catalogue_lon_SEDs_compact.fits'
assoc_to_pulsars = '/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/compact_multiple.fits'	#cross-matched_pulsars_no_clusters_artefacts_multiple_compact.fits
fermi_catalogue = '/home/silvia/Desktop/Papers/Project_4/Pulsars/fermi_compact_unassoc.fits'

PARAMETERS = ["P0", "DM", "RM", "AGE", "DIST", "S1400"]
USE_KS_TEST = True
NAME_COLUMN = "PSRJ"

full_data = {"PSRJ": [], "RAJ": [], "DECJ": [], "P0": [], "P1": []}

def unwrap(a):
    if a > 180:
        return a - 360
    else:
        return a

vunwrap = np.vectorize(unwrap)
		
def parse_float(x):
    return float(x) if x not in ['', '*'] else None

def parse_pp(x):
    return float(x) if x not in ['', '*'] else np.nan
    
with open(FULL_CATALOGUE_PATH, 'r') as f:
    for i, line in enumerate(f):
        if i < 2:
            continue
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(';')
        try:
            name = parts[1].strip()
            ra_sexa = parts[3].strip()
            dec_sexa = parts[6].strip()
            
            coord = SkyCoord(ra=ra_sexa, dec=dec_sexa, unit=(u.hourangle, u.deg))
            ra_deg = coord.ra.deg
            dec_deg = coord.dec.deg
            ra_hours = coord.ra.hour 
            
            period = parts[9].strip()
            pdev = parts[12].strip()


            full_data["PSRJ"].append(name)
            full_data["RAJ"].append(parse_float(ra_hours))
            full_data["DECJ"].append(parse_float(dec_deg))
            full_data["P0"].append(parse_pp(period))
            full_data["P1"].append(parse_pp(pdev))
        except IndexError:
            continue

df_full = pd.DataFrame(full_data)

with open(MATCHED_NAMES_PATH, 'r') as f:
    matched_names = [line.split(';')[0].strip() for line in f if not line.startswith("#")]

df_subset = df_full[df_full[NAME_COLUMN].isin(matched_names)]

mask_replace = df_full["P1"].isna() & (df_full["P0"] < 0.01)
df_full.loc[mask_replace, "P1"] = 1e-21

mask_replace_subset = df_subset["P1"].isna() & (df_subset["P0"] < 0.01)
df_subset.loc[mask_replace_subset, "P1"] = 1e-21

converted_names = df_subset.loc[mask_replace_subset, "PSRJ"]

##################################################################
# PLOTTING THE PULSAR DISTRIBUTION WITH ON TOP THE MATCHED PULSARS
##################################################################
fig = plt.figure(figsize=(16*cm, 10*cm))
ax = fig.add_subplot(1, 1, 1)

coord = SkyCoord(ra=df_full["RAJ"], dec=df_full["DECJ"], unit=(u.hourangle, u.deg))
gal = coord.galactic

l_full = gal.l.deg
b_full = gal.b.deg

b_mask = (b_full >= -11) & (b_full <= 11)
l_mask = ((l_full >= 233) & (l_full <= 360)) | ((l_full <= 44) & (l_full >= 0))

plane_mask = b_mask & l_mask
l_plane = l_full[plane_mask]
b_plane = b_full[plane_mask]

coord_subset = SkyCoord(ra=df_subset["RAJ"], dec=df_subset["DECJ"], unit=(u.hourangle, u.deg))
gal_subset = coord_subset.galactic

l_subset = gal_subset.l.deg
b_subset = gal_subset.b.deg

ax.scatter(vunwrap(l_plane), b_plane, color='lightsteelblue', s=30, label='ATNF v2.6.2')
ax.scatter(vunwrap(l_subset), b_subset, color='forestgreen', edgecolor='black', s=30, label='Matched pulsars')

msp_mask = (df_subset["P0"] <= 0.03) & (df_subset["P1"] < 1e-16)
coord_msp = SkyCoord(ra=df_subset.loc[msp_mask, "RAJ"], dec=df_subset.loc[msp_mask, "DECJ"], unit=(u.hourangle, u.deg))
gal_msp = coord_msp.galactic

ax.scatter(vunwrap(gal_msp.l.deg), gal_msp.b.deg, color='firebrick', edgecolor='black', s=30, marker='s', label='Matched MSPs')

ax.set_xlabel('Galactic longitude (°)')
ax.set_ylabel('Galactic latitude (°)')
ax.set_xlim([-128, 46])
ax.legend(loc='best')
ax.invert_xaxis()

fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/plane_matches.pdf", dpi=300)

################################################
# PLOTTING THE P-PDOT DIAGRAM (KNOWN VS MATCHED)
################################################  
fig = plt.figure(figsize=(13*cm, 10*cm))
ax = fig.add_subplot(1, 1, 1)

P_vals = np.logspace(-3, 2, 500)
B_vals = [1e8, 1e10, 1e12, 1e14]

for B in B_vals:
    Pdot_vals = (B / 3.2e19)**2 / P_vals
    ax.plot(P_vals, Pdot_vals, '--', color="black", alpha=0.6, lw=1, zorder=-1)

B_plot = [1e10, 1e12, 1e14]   
for B in B_plot:
    label_x = P_vals[int(0.88 * len(P_vals))]  # 80% along the line
    label_y = (B / 3.2e19)**2 / label_x

    B_exp = int(np.log10(B))
    label = rf'$10^{{{B_exp}}}$ G'

    ax.text(label_x, label_y, label, fontsize=9, verticalalignment='bottom', horizontalalignment='left', rotation=340, color='black', alpha=0.6)

ax.scatter(df_full["P0"], df_full["P1"], color='lightsteelblue', s=30, label='ATNF v2.6.2', zorder=1)
ax.scatter(df_subset["P0"], df_subset["P1"], color='forestgreen', edgecolor='black', s=30, label='Matched pulsars', zorder=2)

ax.set_xlabel('Period (s)')
ax.set_xscale('log')
ax.set_ylabel('Period derivative (s/s)')
ax.set_yscale('log')

ax.set_xlim([1e-3, 1e2])
ax.set_ylim([6e-23, 1e-9])

msp_mask = (df_subset["P0"] <= 0.03) & (df_subset["P1"] < 1e-16)

ax.scatter(df_subset.loc[msp_mask, "P0"], df_subset.loc[msp_mask, "P1"], color='firebrick', edgecolor='black', s=30, marker='s', label='Matched MSPs', zorder=2)
ax.legend(loc="lower right")

fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/p-pdot.pdf", bbox_inches="tight", dpi=300)
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/p-pdot.png", bbox_inches="tight", dpi=300)

##################################################################
# PLOTTING THE GAMMA DISTRIBUTION WITH ON TOP THE MATCHED SOURCES
##################################################################

# Loading the data:
with fits.open(matches_catalogue) as hdul:
    data = Table(hdul[1].data)
    names = [name for name in data.colnames if len(data[name].shape) <= 1]
    df = data[names].to_pandas()

with fits.open(assoc_to_pulsars) as hdu:
    data = hdu[1].data
    dff = pd.DataFrame(data)
    
coords = SkyCoord(ra=dff['ref_ra']*u.deg, dec=dff['ref_dec']*u.deg)
formatted_names = [f'GP J{ra.to_string(unit=u.hourangle, sep="", precision=1, pad=True)} {dec.to_string(sep="", precision=0, alwayssign=True, pad=True)}' for ra, dec in zip(coords.ra, coords.dec)]
            
filtered_df = df[~df['Source_Name'].isin(formatted_names)]
        
with fits.open(fermi_catalogue) as h:
    tutte = Table(h[1].data)
    names = [name for name in tutte.colnames if len(tutte[name].shape) <= 1]
    gamma = tutte[names].to_pandas()

###################################################
# Plotting the distribution of the sources in gamma  
###################################################
fig = plt.figure(figsize=(16*cm, 10*cm))
ax = fig.add_subplot(1, 1, 1)

gamma_coord = SkyCoord(ra=gamma['RAJ2000']*u.deg, dec=gamma['DEJ2000']*u.deg)
gal = gamma_coord.galactic

l_full = gal.l.deg
b_full = gal.b.deg

b_mask = (b_full >= -11) & (b_full <= 11)
l_mask = ((l_full >= 233) & (l_full <= 360)) | ((l_full <= 44) & (l_full >= 0))

plane_mask = b_mask & l_mask
l_plane = l_full[plane_mask]
b_plane = b_full[plane_mask]

coords = SkyCoord(ra=filtered_df['ref_ra']*u.deg, dec=filtered_df['ref_dec']*u.deg)
gal_subset = coords.galactic

l_subset = gal_subset.l.deg
b_subset = gal_subset.b.deg        

ax.scatter(vunwrap(l_plane), b_plane, color='grey', s=30, label='4FGL-DR4 unassociated')
ax.scatter(vunwrap(l_subset), b_subset, color='#084081', edgecolor='black', s=30, label='Matched sources')

ax.set_xlabel('Galactic longitude (°)')
ax.set_ylabel('Galactic latitude (°)')
ax.set_xlim([-128, 46])
ax.legend(loc='best')
ax.invert_xaxis()

fig.tight_layout()
fig.savefig(f"/home/silvia/Desktop/Papers/Paper_IV/plane_gamma.pdf", dpi=300)

######################################################
# Separation in arcminutes between pulsars and matched  
######################################################
fig = plt.figure(figsize=(13*cm, 10*cm))
ax = fig.add_subplot(1, 1, 1)

sep = []
with open('/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/FINAL/separation.txt', 'r') as f:
    for line in f:
        sep.append(float(line))

sep = np.array(sep)

bins = np.logspace(np.log10(sep.min()), np.log10(sep.max()), num=16)
ax.hist(sep, color='#084081', edgecolor='black', bins=bins, alpha=0.7)

ax.set_ylabel('Count')
ax.set_xlabel("Separation (')")
ax.set_xscale('log')

fig.tight_layout()
fig.savefig(f"/home/silvia/Desktop/Papers/Paper_IV/separation.pdf", dpi=300)