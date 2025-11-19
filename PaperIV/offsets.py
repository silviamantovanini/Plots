#!/usr/bin/env python3

##################################################################
# SCRIPT TO PLOT THE OFFSETS IN RA AND DEC FOR THE MATCHED SOURCES
##################################################################

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp, mannwhitneyu
import seaborn as sns
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.coordinates import search_around_sky
import astropy.units as u
from collections import Counter

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12
})

cm = 1 / 2.54

FULL_CATALOGUE_PATH = "/home/silvia/Desktop/Papers/Project_4/Pulsars/pulsar2_6_1.txt"
MATCHED_NAMES_PATH = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/final_matches.txt"	#"matched_pulsars.txt"
SOURCE_CATALOGUE = '/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/cross-matched_pulsars_no_clusters_artefacts_multiple_compact.fits'
#"/home/silvia/Desktop/Papers/Project_4/Pulsars/matched_sources.fits"

NAME_COLUMN = "PSRJ"
NAME_COLUMN_SOURCE = 'Pulsar_Name'

with open(MATCHED_NAMES_PATH, 'r') as f:
    matched_names = [line.split(';')[0].strip() for line in f if not line.startswith("#")]
    
full_data = {"PSRJ": [], "RAJ": [], "DECJ": [], "err_RA": [], "err_DEC": []}
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
            ra_err = float(parts[4].strip())*15
            dec_err = float(parts[7].strip())*15
            
            coord = SkyCoord(ra=ra_sexa, dec=dec_sexa, unit=(u.hourangle, u.deg))
            ra_deg = coord.ra.deg
            dec_deg = coord.dec.deg
            ra_hours = coord.ra.hour 

            full_data["PSRJ"].append(name)
            full_data["RAJ"].append(ra_deg)
            full_data["DECJ"].append(dec_deg)
            full_data["err_RA"].append(ra_err)
            full_data["err_DEC"].append(dec_err)

        except IndexError:
            continue

df_full = pd.DataFrame(full_data)
df_subset = df_full[df_full[NAME_COLUMN].isin(matched_names)]

source_data = {"Pulsar_Name": [], "ra_str": [], "dec_str": []}
with fits.open(SOURCE_CATALOGUE) as hdul:
    data = hdul[1].data
    name = data['PSRJ']
    ra = data['ra_str']
    dec = data['dec_str']
    
    name_counts = Counter(name)
    unique_indices = [i for i, n in enumerate(name) if name_counts[n] == 1]
    unique_names = name[unique_indices]
    unique_ra = ra[unique_indices]
    unique_dec = dec[unique_indices]
    
    coords = SkyCoord(ra=unique_ra, dec=unique_dec, unit=(u.hourangle, u.deg))
    ra_deg = coords.ra.deg
    dec_deg = coords.dec.deg
    
    source_data["Pulsar_Name"].extend(unique_names)
    source_data["ra_str"].extend(ra_deg)
    source_data["dec_str"].extend(dec_deg)
    
df_source = pd.DataFrame(source_data)
df_source_subset = df_source[df_source[NAME_COLUMN_SOURCE].isin(matched_names)]

unique_names_clean = [str(n).strip().upper() for n in unique_names]
df_subset = df_subset.copy()
df_subset[NAME_COLUMN] = df_subset[NAME_COLUMN].str.strip().str.upper()
df_subset_unique = df_subset[df_subset[NAME_COLUMN].isin(unique_names_clean)]

print(len(df_subset_unique))
print(len(df_source_subset))

removed_pulsars = set(name) - set(name[unique_indices])
#print(f"Removed pulsars due to duplicates: {removed_pulsars}")
print(len(unique_names))
print(len(removed_pulsars))

df_merged = pd.merge(df_subset_unique, df_source_subset, left_on="PSRJ", right_on="Pulsar_Name", suffixes=('_atnf', '_mwa'))

df_merged["delta_ra"] = (df_merged["RAJ"] - df_merged["ra_str"]) * 3600
df_merged["delta_dec"] = (df_merged["DECJ"] - df_merged["dec_str"]) * 3600

#df_merged["total_offset_arcsec"] = np.sqrt(df_merged["delta_ra_arcsec"]**2 + df_merged["delta_dec_arcsec"]**2)

print(np.mean(df_merged["delta_ra"]*np.cos(df_merged["dec_str"])), np.mean(df_subset_unique["err_RA"]))
print(np.mean(df_merged["delta_dec"]), np.mean(df_subset_unique["err_DEC"]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(32*cm, 16*cm))

ax1.hist(df_merged["delta_ra"], bins=20, color='#084081', edgecolor='black', alpha=0.7)
ax1.set_title(f'Mean RA offset = {np.mean(df_merged["delta_ra"]):.2f}')
ax1.set_xlabel("RA offset ('')")
ax1.set_ylabel('Count')
ax1.grid(True)

ax2.hist(df_merged["delta_dec"], bins=20, color='#084081', edgecolor='black', alpha=0.7)
ax2.set_title(f'Mean Dec offset = {np.mean(df_merged["delta_dec"]):.2f}')
ax2.set_xlabel("Dec offset ('')")
ax2.set_ylabel('Count')
ax2.grid(True)

fig.tight_layout()
fig.savefig("offsets.pdf", dpi=200)