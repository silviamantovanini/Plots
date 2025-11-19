#! /usr/bin/env python3

############################################
# GOALS OF THE SCRIPT
# 1 - COMPARISON WITH KNOWN PULSARS IN GAMMA
############################################

import re
import pandas
import corner
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12})

cm = 1/2.54

############################################################################
# Link to all the files and directories I am gonna use:
matches_catalogue = '/home/silvia/Desktop/Papers/Project_4/Pulsars/matched_4FGL_out_FLAGS_artefacts.fits'
pulsar_catalogue = '/home/silvia/Desktop/Papers/Project_4/Pulsars/pulsar2_6_1.txt'
psr_gamma = '/home/silvia/Desktop/Papers/Project_4/Pulsars/gamma_pulsars.fits'
msp_gamma = '/home/silvia/Desktop/Papers/Project_4/Pulsars/gamma_msps.fits'

paper_dir = '/home/silvia/Desktop/Papers/Paper_IV'
plots_dir = '/home/silvia/Desktop/Papers/Project_4/Pulsars/paper_plots'

############################################################################
# Loading the data:
with fits.open(psr_gamma) as hdul:
    data = Table(hdul[1].data)
    names = [name for name in data.colnames if len(data[name].shape) <= 1]
    df_psr = data[names].to_pandas()

with fits.open(msp_gamma) as hdu:
    datam = Table(hdu[1].data)
    names = [name for name in datam.colnames if len(datam[name].shape) <= 1]
    df_msp = datam[names].to_pandas()
 
with fits.open(matches_catalogue) as h:
    datac = Table(h[1].data)
    names = [name for name in datac.colnames if len(datac[name].shape) <= 1]
    df_cand = datac[names].to_pandas()

df_cand = df_cand.drop_duplicates(subset='Source_Name')
     
coords_psr = SkyCoord(l=df_psr['GLON']*u.deg, b=df_psr['GLAT']*u.deg, frame='galactic')
coords_msp = SkyCoord(l=df_msp['GLON']*u.deg, b=df_msp['GLAT']*u.deg, frame='galactic')
coords_cand = SkyCoord(l=df_cand['GLON']*u.deg, b=df_cand['GLAT']*u.deg, frame='galactic')

sm_psr = df_psr['Conf_95_SemiMajor']
sm_msp = df_msp['Conf_95_SemiMajor']
sm_cand = df_cand['Conf_95_SemiMajor']

############################################################################
# Galactic latitude distribution as a function of the semimajor axis

fig = plt.figure(figsize=(18*cm, 15*cm))
ax = fig.add_subplot(1, 1, 1)

ax.scatter(coords_psr.b, sm_psr, color='grey', s=20, label='Pulsars')
ax.scatter(coords_msp.b, sm_msp, color='darkmagenta', s=20, label='MSPs')
ax.scatter(coords_cand.b, sm_cand, color='#084081', s=20, label='Pulsar-like sources')

ax.set_xlabel('Galactic Latitude (°)')
ax.set_ylabel('Semi major axis (°)')
ax.legend(loc='best')

fig.tight_layout()
fig.savefig(f"{plots_dir}/lat_distribution.pdf", dpi=300)

############################################################################
# Distribution of the 3 parameters used for filtering
pl_psr = df_psr['PL_Index']
pl_msp = df_msp['PL_Index']
pl_cand = df_cand['PL_Index']

lp_psr = df_psr['LP_SigCurv']
lp_msp = df_msp['LP_SigCurv']
lp_cand = df_cand['LP_SigCurv']

var_psr = df_psr['Variability_Index']
var_msp = df_msp['Variability_Index']
var_cand = df_cand['Variability_Index']

fig, axs = plt.subplots(3, 3, figsize=(40*cm, 25*cm))
axs = axs.flat

bins = np.histogram_bin_edges(pandas.concat([pl_psr, pl_msp, pl_cand]), bins=20)

axs[0].hist(pl_psr, bins=bins, color='grey', label='Pulsars')
axs[3].hist(pl_msp, bins=bins, color='darkmagenta', label='MSPs')
axs[6].hist(pl_cand, bins=bins, color='#084081', label='Pulsar-like sources')
axs[6].set_xlabel('Photon index')

for ax in (axs[0], axs[3], axs[6]):
    ax.legend(loc='upper left')
    ax.set_xlim([1.9, 3.0])
    ax.set_ylabel('Counts')
    
bins = np.histogram_bin_edges(pandas.concat([np.log10(lp_psr), np.log10(lp_msp), np.log10(lp_cand)]), bins=20)

axs[1].hist(np.log10(lp_psr), bins=bins, color='grey', label='Pulsars')
axs[4].hist(np.log10(lp_msp), bins=bins, color='darkmagenta', label='MSPs')
axs[7].hist(np.log10(lp_cand), bins=bins, color='#084081', label='Pulsar-like sources')
axs[7].set_xlabel('Significance curve')

for ax in (axs[1], axs[4], axs[7]):
    ax.set_xlim([0.0, 2.5])
    ax.legend()

bins = np.histogram_bin_edges(pandas.concat([np.log10(var_psr), np.log10(var_msp), np.log10(var_cand)]), bins=20)

axs[2].hist(np.log10(var_psr), bins=bins, color='grey', label='Pulsars')
axs[5].hist(np.log10(var_msp), bins=bins, color='darkmagenta', label='MSPs')
axs[8].hist(np.log10(var_cand), bins=bins, color='#084081', label='Pulsar-like sources')
axs[8].set_xlabel('Variability index')

for ax in (axs[2], axs[5], axs[8]):
    ax.set_xlim([0.5, 2.5])
    ax.legend()
    
fig.tight_layout()
fig.savefig(f"{plots_dir}/params.pdf", dpi=300)

############################################################################
# Scatter plots of the 3 parameters used for filtering
fig, axs = plt.subplots(3, 1, figsize=(12*cm, 30*cm))
axs = axs.flat

axs[0].scatter(pl_psr, np.log10(lp_psr), color='grey', marker='s', label='Pulsars', s=50)
axs[0].scatter(pl_msp, np.log10(lp_msp), color='darkmagenta', marker='^', label='MSPs', s=50)
axs[0].scatter(pl_cand, np.log10(lp_cand), color='#084081', edgecolor='black', label='Pulsar-like sources', s=50)
axs[0].set_xlabel('Photon index')
axs[0].set_ylabel('Significance curve')
axs[0].legend(loc='upper right')

axs[1].scatter(np.log10(lp_psr), np.log10(var_psr), color='grey', marker='s', label='Pulsars', s=50)
axs[1].scatter(np.log10(lp_msp), np.log10(var_msp), color='darkmagenta', marker='^', label='MSPs', s=50)
axs[1].scatter(np.log10(lp_cand), np.log10(var_cand), color='#084081', edgecolor='black', label='Pulsar-like sources', s=50)
axs[1].set_xlabel('Significance curve')
axs[1].set_ylabel('Variability index')
axs[1].legend(loc='upper right')

axs[2].scatter(np.log10(var_psr), pl_psr, color='grey', marker='s', label='Pulsars', s=50)
axs[2].scatter(np.log10(var_msp), pl_msp, color='darkmagenta', marker='^', label='MSPs', s=50)
axs[2].scatter(np.log10(var_cand), pl_cand, color='#084081', edgecolor='black', label='Pulsar-like sources', s=50)
axs[2].set_xlabel('Variability index')
axs[2].set_ylabel('Photon index')
axs[2].legend(loc='upper right')
    
fig.tight_layout()
fig.savefig(f"{plots_dir}/params_scatter.pdf", dpi=300)

#
df_psr = pandas.DataFrame({"Photon index": pl_psr, "Significance curve": np.log10(lp_psr), "Variability index": np.log10(var_psr), "Class": "Pulsars"})
df_msp = pandas.DataFrame({"Photon index": pl_msp, "Significance curve": np.log10(lp_msp), "Variability index": np.log10(var_msp), "Class": "MSPs"})
df_cand = pandas.DataFrame({"Photon index": pl_cand, "Significance curve": np.log10(lp_cand), "Variability index": np.log10(var_cand), "Class": "Pulsar-like sources"})

df_all = pandas.concat([df_psr, df_msp, df_cand], ignore_index=True)

sns.set(style="ticks")
plot = sns.pairplot(df_all, hue="Class", palette={"Pulsars": "grey", "MSPs": "darkmagenta", "Pulsar-like sources": "#084081"}, plot_kws={"s": 50, "alpha": 0.7}, diag_kind="hist", diag_kws={"bins": 15}, corner=True)

#diag_kws=dict(fill=False),    
if plot._legend is not None:
    plot._legend.set_title("")  
    plot._legend.set_bbox_to_anchor((0.85, 0.85))
    plot._legend._loc = 1 
    plot._legend.set_frame_on(True)
    #plot._legend.get_frame().set_edgecolor('black')
    #plot._legend.get_frame().set_linewidth(1)

#plot._legend.remove()
plot.fig.set_size_inches(18*cm, 18*cm)
plot.fig.tight_layout()
plot.fig.savefig(f"{plots_dir}/params_triangle.pdf", bbox_inches="tight", dpi=300)

############################################################################
# Plotting gamma energy flux vs radio flux density
energy_cand = df_cand['Energy_Flux100']
index = df_cand['pl_alpha']
flux200 = df_cand['S_200']

flux = flux200 * (pow(7, -index))

def normalize_name(name):
    name = name.replace("4FGL", "").strip()
    match = re.match(r'J(\d+)\.\d+([+-]\d+)', name)
    if match:
        return f"J{match.group(1)}{match.group(2)}"
    
    return name
    
fluxes_psr = []
energy_psr = []

df_psr['Normalized_Name'] = df_psr['Source_Name'].apply(normalize_name)
df_msp['Normalized_Name'] = df_msp['Source_Name'].apply(normalize_name)

normalized_names_psr = df_psr['Normalized_Name']
with open("/home/silvia/Desktop/Papers/Project_4/Pulsars/pulsar2_6_1.txt", 'r') as f:
    for line in f:
        line = line.strip()
        parts = line.split(';')
        name = parts[1].strip()
        
        if name in normalized_names_psr.tolist():
            fluxes_psr.append(parts[24].strip())
            energy_value = df_psr.loc[df_psr['Normalized_Name'] == parts[1].strip(), 'Energy_Flux100'].values[0]
            energy_psr.append(energy_value)

fluxes_msp = []
energy_msp = []
normalized_names_msp = df_msp['Normalized_Name']
with open("/home/silvia/Desktop/Papers/Project_4/Pulsars/pulsar2_6_1.txt", 'r') as f:
    for line in f:
        line = line.strip()
        parts = line.split(';')
        name = parts[1].strip()
        
        if name in normalized_names_msp.tolist():
            fluxes_msp.append(parts[24].strip())
            energy_value = df_msp.loc[df_msp['Normalized_Name'] == parts[1].strip(), 'Energy_Flux100'].values[0]
            energy_msp.append(energy_value)
                                                         
fig = plt.figure(figsize=(18*cm, 15*cm))
ax = fig.add_subplot(1, 1, 1)

ax.scatter(energy_cand, flux, color='#084081', edgecolor='black', label='Pulsar-like sources', s=50)
ax.scatter(energy_psr, fluxes_psr, color='grey', marker='s', label='Pulsars', s=50)
ax.scatter(energy_msp, fluxes_msp, color='darkmagenta', marker='^', label='MSPs', s=50)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([0.01, 100])
ax.set_xlabel('$\gamma$-ray energy flux (erg/cm$^2$/s)')
ax.set_ylabel('1400 MHz radio flux density (mJy)')

fig.tight_layout()
fig.savefig(f"{plots_dir}/fluxes.pdf", dpi=300)