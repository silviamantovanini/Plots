#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp, mannwhitneyu
import seaborn as sns
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
import astropy.units as u

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 18 #12
})

cm = 1 / 2.54

FULL_CATALOGUE_PATH = "/home/silvia/Desktop/Papers/Project_4/Pulsars/pulsar2_6_2.txt"
MATCHED_NAMES_PATH = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/FINAL/final_matches.txt"	#matched_pulsars.txt

PARAMETERS = ["P0", "DM", "RM", "AGE", "DIST", "S1400"]
USE_KS_TEST = True
NAME_COLUMN = "PSRJ"

full_data = {"PSRJ": [], "RAJ": [], "DECJ": [], "P0": [], "DM": [], "RM": [], "AGE": [], "DIST": [], "S1400": [], "P1": [], "W10": []}

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
            dm = parts[15].strip()
            rm = parts[18].strip()	
            dist = parts[36].strip()	#30
            age = parts[41].strip()	#35
            s1400 = parts[30].strip()	#24
            w10 = parts[24].strip()	#w50=21

            full_data["PSRJ"].append(name)
            full_data["RAJ"].append(parse_float(ra_hours))
            full_data["DECJ"].append(parse_float(dec_deg))
            full_data["P0"].append(parse_pp(period))
            full_data["DM"].append(parse_float(dm))
            full_data["RM"].append(parse_float(rm))
            full_data["AGE"].append(parse_float(age))
            full_data["DIST"].append(parse_float(dist))
            full_data["S1400"].append(parse_float(s1400))
            full_data["P1"].append(parse_pp(pdev))
            full_data["W10"].append(parse_float(w10))
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
print("Pulsars with missing P1 replaced with 1e-21 (P0 < 0.01):")
print(converted_names.to_string(index=False))
print(df_subset.loc[converted_names.index, "P0"])

print(f"Loaded full catalogue: {df_full.shape[0]} sources")
print(f"Extracted matched subset: {df_subset.shape[0]} sources using {len(matched_names)} names")

# PLOTTING COMPARISON OF DISTRIBUTION OF VARIOUS PARAMETERS (KNOWN VS MATCHED)
fig, axes = plt.subplots(2, 3, figsize=(48*cm, 25*cm))
axes = axes.flatten()

plot_idx = 0
transform_log = {"AGE", "RM", "P0", "DM", "DIST", "S1400"}

for param in PARAMETERS:
    if param in df_full.columns and param in df_subset.columns:
        x_full = df_full[param].dropna()
        x_subset = df_subset[param].dropna()

        if len(x_full) == 0 or len(x_subset) == 0:
            print(f"Skipping {param}: missing or empty values.")
            continue

        if param in transform_log:
            x_full = np.log10(x_full[x_full > 0])
            x_subset = np.log10(x_subset[x_subset > 0])
            label_suffix = " (log)"
        else:
            label_suffix = ""

        if USE_KS_TEST:
            stat, p_value = ks_2samp(x_full, x_subset)
            test_name = "Kolmogorov-Smirnov"
        else:
            stat, p_value = mannwhitneyu(x_full, x_subset, alternative='two-sided')
            test_name = "Mann-Whitney U"

        #print(f"{test_name} Test for '{param}':")
        #print(f"Statistic = {stat:.4f}, P-value = {p_value:.4g}")
        #print("  " + ("Significant" if p_value < 0.05 else "Not significant"))

        ax = axes[plot_idx]
        sns.kdeplot(x_full, label='ATNF v2.6.1 pulsars', lw=2, ax=ax)
        sns.kdeplot(x_subset, label='Matched pulsars', lw=2, linestyle='--', ax=ax)
        ax.set_xlabel(f"{param}{label_suffix}")
        ax.set_ylabel("Density")
        ax.legend()
        ax.grid(True)
        plot_idx += 1

# Remove empty subplot if fewer than 6 plots
if plot_idx < len(axes):
    for i in range(plot_idx, len(axes)):
        fig.delaxes(axes[i])

fig.tight_layout()
fig.savefig("all_distributions.pdf", dpi=300)
print("Comparison complete and plot saved as 'all_distributions.pdf'.")


####################################################
# SAME AS ABOVE BUT DIFFERENT STYLE AND ONLY 3 PLOTS
####################################################
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(54 * cm, 15 * cm))

if USE_KS_TEST:
    stat, p_value = ks_2samp(x_full, x_subset)
    test_name = "Kolmogorov-Smirnov"
else:
    stat, p_value = mannwhitneyu(x_full, x_subset, alternative='two-sided')
    test_name = "Mann-Whitney U"

sns.kdeplot(np.log10(df_full["P0"]), label='ATNF v2.6.1 pulsars', lw=2, ax=ax1, fill=True, color='grey', alpha=0.9)
sns.kdeplot(np.log10(df_subset["P0"]), label='Matched pulsars', lw=2, linestyle='--', ax=ax1, fill=True, color='#084081', alpha=0.75)

ax1.set_xlabel(f"P [s]")
ax1.set_ylabel("Density")
ax1.legend()

sns.kdeplot(np.log10(df_full["DM"]), label='ATNF v2.6.1 pulsars', lw=2, ax=ax2, fill=True, color='grey', alpha=0.9)
sns.kdeplot(np.log10(df_subset["DM"]), label='Matched pulsars', lw=2, linestyle='--', ax=ax2, fill=True, color='#084081', alpha=0.75)

ax2.set_xlabel(r"DM [pc/cm$^3$]")
ax2.set_ylabel('')
ax2.legend()

sns.kdeplot(np.log10(df_full["S1400"]), label='ATNF v2.6.1 pulsars', lw=2, ax=ax3, fill=True, color='grey', alpha=0.9)
sns.kdeplot(np.log10(df_subset["S1400"]), label='Matched pulsars', lw=2, linestyle='--', ax=ax3, fill=True, color='#084081', alpha=0.75)

ax3.set_xlabel(f"Flux density [mJy]")
ax3.set_ylabel('')
ax3.legend()

fig.tight_layout()
fig.savefig("P-DM-FLUX_distributions.pdf", dpi=300)

# and histograms
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(54 * cm, 15 * cm))

sns.histplot(np.log10(df_full["P0"]), label='ATNF v2.6.1 pulsars', lw=2, ax=ax1, fill=True, color='grey', alpha=0.9, stat='density')
sns.histplot(np.log10(df_subset["P0"]), label='Matched pulsars', lw=2, linestyle='--', ax=ax1, fill=True, color='#084081', alpha=0.75, stat='density')

ax1.set_xlabel(f"P [s]")
ax1.set_ylabel("Density")
ax1.legend()

sns.histplot(np.log10(df_full["DM"]), label='ATNF v2.6.1 pulsars', lw=2, ax=ax2, fill=True, color='grey', alpha=0.9, stat='density')
sns.histplot(np.log10(df_subset["DM"]), label='Matched pulsars', lw=2, linestyle='--', ax=ax2, fill=True, color='#084081', alpha=0.75, stat='density')

ax2.set_xlabel(r"DM [pc/cm$^3$]")
ax2.set_ylabel('')
ax2.legend()

sns.histplot(np.log10(df_full["S1400"]), label='ATNF v2.6.1 pulsars', lw=2, ax=ax3, fill=True, color='grey', alpha=0.9, stat='density')
sns.histplot(np.log10(df_subset["S1400"]), label='Matched pulsars', lw=2, linestyle='--', ax=ax3, fill=True, color='#084081', alpha=0.75, stat='density')

ax3.set_xlabel(f"Flux density [mJy]")
ax3.set_ylabel('')
ax3.legend()

fig.tight_layout()
fig.savefig("P-DM-FLUX_histograms.pdf", dpi=300)

####################################################################
# and only two
####################################################################
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(54 * cm, 15 * cm)) #54

valid_flux_full = df_full["S1400"].replace([np.inf, -np.inf], np.nan).dropna()
valid_flux_full = valid_flux_full[valid_flux_full > 0]

valid_flux_subset = df_subset["S1400"].replace([np.inf, -np.inf], np.nan).dropna()
valid_flux_subset = valid_flux_subset[valid_flux_subset > 0]

log_flux_full = np.log10(valid_flux_full)
log_flux_subset = np.log10(valid_flux_subset)

combined_log_flux = np.concatenate([log_flux_full, log_flux_subset])
min_val = combined_log_flux.min()
max_val = combined_log_flux.max()
bins = np.linspace(min_val, max_val, num=20)
#bin_edges = np.histogram_bin_edges(combined_log_flux, bins='scott')

sns.histplot(np.log10(df_full["S1400"]), label='ATNF v2.6.1', ax=ax1, fill=True, color='lightsteelblue', alpha=0.9, stat='density', bins=bins)
sns.histplot(np.log10(df_subset["S1400"]), label='Matched pulsars', ax=ax1, fill=True, color='forestgreen', alpha=0.75, stat='density', bins=bins)

ax1.set_xlabel(f"Flux density [mJy]")
ax1.set_ylabel("Density")
ax1.legend()

valid_dm_full = df_full["DM"].replace([np.inf, -np.inf], np.nan).dropna()
valid_dm_full = valid_dm_full[valid_dm_full > 0]

valid_dm_subset = df_subset["DM"].replace([np.inf, -np.inf], np.nan).dropna()
valid_dm_subset = valid_dm_subset[valid_dm_subset > 0]

combined_dm = np.concatenate([valid_dm_full, valid_dm_subset])
min_val = combined_dm.min()
max_val = combined_dm.max()
bins = np.linspace(min_val, max_val, num=30)
#bin_edges = np.histogram_bin_edges(combined_dm, bins='scott')

sns.histplot(df_full["DM"], label='ATNF v2.6.1', ax=ax2, fill=True, color='lightsteelblue', alpha=0.9, stat='density', bins=bins)
sns.histplot(df_subset["DM"], label='Matched sources', ax=ax2, fill=True, color='forestgreen', alpha=0.75, stat='density', bins=bins)

ax2.set_xlabel(r"DM [pc/cm$^3$]")
ax2.set_ylabel('')
ax2.set_xlim([0, 1350])
#ax2.set_xscale('log')
ax2.legend()

valid_full = df_full[["W10", "P0"]].replace([np.inf, -np.inf], np.nan).dropna()
valid_full = valid_full[(valid_full["W10"] > 0) & (valid_full["P0"] > 0)]

valid_subset = df_subset[["W10", "P0"]].replace([np.inf, -np.inf], np.nan).dropna()
valid_subset = valid_subset[(valid_subset["W10"] > 0) & (valid_subset["P0"] > 0)]

frac_width_full = valid_full["W10"] / (valid_full["P0"]*1000)
frac_width_subset = valid_subset["W10"] / (valid_subset["P0"]*1000)

combined_frac_width = np.concatenate([frac_width_full, frac_width_subset])

min_val = combined_frac_width.min()
max_val = combined_frac_width.max()
log_bins = np.logspace(np.log10(min_val), np.log10(max_val), num=20)

sns.histplot(frac_width_full, label='ATNF v2.6.1', ax=ax3, fill=True, color='lightsteelblue', alpha=0.9, stat='density', bins=log_bins)
sns.histplot(frac_width_subset, label='Matched pulsars', ax=ax3, fill=True, color='forestgreen', alpha=0.75, stat='density', bins=log_bins)

ax3.set_xlabel(r"W$_{10}$ / P")
ax3.set_ylabel('')
ax3.set_xscale('log')
ax3.legend()

for a in (ax1, ax2, ax3):
    print(a.get_position())
    
fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/DM-FLUX.pdf", dpi=300)

##################################################################
# PLOTTING THE PULSAR DISTRIBUTION WITH ON TOP THE MATCHED PULSARS
##################################################################
fig = plt.figure(figsize=(20*cm, 13*cm))
ax = fig.add_subplot(1, 1, 1)

ax.scatter(df_full["RAJ"], df_full["DECJ"], color='lightsteelblue', s=30, label='ATNF v2.6.1')
ax.scatter(df_subset["RAJ"], df_subset["DECJ"], color='forestgreen', edgecolor='black', s=30, label='Matched pulsars')

msp_mask = (df_subset["P0"] < 0.03) & (df_subset["P1"] < 1e-16)
ax.scatter(df_subset.loc[msp_mask, "RAJ"], df_subset.loc[msp_mask, "DECJ"], color='forestgreen', edgecolor='black', s=30, marker='s', label='Matched MSPs', zorder=2)

l_vals_full = np.linspace(0, 360, 1000)
b_vals = np.zeros_like(l_vals_full)

gal_up_full = SkyCoord(l=l_vals_full * u.deg, b=11 * u.deg, frame='galactic')
gal_down_full = SkyCoord(l=l_vals_full * u.deg, b=-11 * u.deg, frame='galactic')
eq_up_full = gal_up_full.transform_to('fk5')
eq_down_full = gal_down_full.transform_to('fk5')

def split_ra_discontinuity(ra_array, dec_array, threshold=12):
    ra = np.array(ra_array)
    dec = np.array(dec_array)
    diff = np.abs(np.diff(ra))
    split_indices = np.where(diff > threshold)[0] + 1
    segments = np.split(np.column_stack((ra, dec)), split_indices)
    return segments
    
ra_up = eq_up_full.ra.hour
dec_up = eq_up_full.dec.deg
ra_down = eq_down_full.ra.hour
dec_down = eq_down_full.dec.deg

for segment in split_ra_discontinuity(ra_up, dec_up):
    ax.plot(segment[:, 0], segment[:, 1], color='black', lw=1.5, linestyle='--', zorder=1)
for segment in split_ra_discontinuity(ra_down, dec_down):
    ax.plot(segment[:, 0], segment[:, 1], color='black', lw=1.5, linestyle='--', zorder=1)

l_vals_inrange = np.concatenate([np.linspace(233, 360, 500, endpoint=False), np.linspace(0, 44, 500)])
b_vals_inrange = np.zeros_like(l_vals_inrange)

gal_up_in = SkyCoord(l=l_vals_inrange * u.deg, b=11 * u.deg, frame='galactic')
gal_down_in = SkyCoord(l=l_vals_inrange * u.deg, b=-11 * u.deg, frame='galactic')
eq_up_in = gal_up_in.transform_to('fk5')
eq_down_in = gal_down_in.transform_to('fk5')

ra_up_in = eq_up_in.ra.hour
dec_up_in = eq_up_in.dec.deg
ra_down_in = eq_down_in.ra.hour
dec_down_in = eq_down_in.dec.deg

for segment in split_ra_discontinuity(ra_up_in, dec_up_in):
    ax.plot(segment[:, 0], segment[:, 1], color='firebrick', lw=2, linestyle='--', zorder=1)
for segment in split_ra_discontinuity(ra_down_in, dec_down_in):
    ax.plot(segment[:, 0], segment[:, 1], color='firebrick', lw=2, linestyle='--', zorder=1)

b_range = np.linspace(-11, 11, 300)
l_bounds = [233, 44]

for l_val in l_bounds:
    gal_line = SkyCoord(l=np.full_like(b_range, l_val)*u.deg, b=b_range*u.deg, frame='galactic')
    eq_line = gal_line.transform_to('fk5')
    
    ra_line = eq_line.ra.hour
    dec_line = eq_line.dec.deg
    
    ax.plot(ra_line, dec_line, color='firebrick', linestyle='--', linewidth=2, zorder=1)

ax.set_xlabel('Right Ascension (h)')
ax.set_ylabel('Declination (°)')
ax.set_xlim([0,24])
ax.legend()
ax.invert_xaxis()

fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/matched_pulsar.pdf", dpi=300)

################################################
# PLOTTING THE P-PDOT DIAGRAM (KNOWN VS MATCHED)
################################################
mag_coords = []

with open("/home/silvia/Desktop/Papers/Project_4/Pulsars/magnetars.txt", 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(',')
        try:
            ra_sexa = parts[0].strip()
            dec_sexa = parts[1].strip()
            coord = SkyCoord(ra=ra_sexa, dec=dec_sexa, unit=(u.hourangle, u.deg))
            mag_coords.append(coord)
        except Exception as e:
            print(f"Skipping malformed line: {line} | Error: {e}")

mag_coords = SkyCoord(mag_coords)

pulsar_coords = SkyCoord(ra=df_subset["RAJ"].values, dec=df_subset["DECJ"].values, unit=(u.hourangle, u.deg))
idx_mag, idx_pulsar, sep2d, _ = search_around_sky(mag_coords, pulsar_coords, 0.5 * u.arcmin)

print(f"Matched magnetars with matched pulsars: {len(idx_mag)}")
matched_mag_pulsars = df_subset.iloc[idx_pulsar]
          

fig = plt.figure(figsize=(18*cm, 15*cm))
ax = fig.add_subplot(1, 1, 1)

P_vals = np.logspace(-3, 2, 500)
B_vals = [1e8, 1e10, 1e12, 1e14]

for B in B_vals:
    Pdot_vals = (B / 3.2e19)**2 / P_vals
    ax.plot(P_vals, Pdot_vals, '--', color="black", alpha=0.6, lw=1, zorder=-1)

B_plot = [1e10, 1e12, 1e14]   
for B in B_plot:
    label_x = P_vals[int(0.9 * len(P_vals))]  # 80% along the line
    label_y = (B / 3.2e19)**2 / label_x

    B_exp = int(np.log10(B))
    label = rf'$10^{{{B_exp}}}$ G'

    ax.text(label_x, label_y, label, fontsize=9, verticalalignment='bottom', horizontalalignment='left', rotation=340, color='black', alpha=0.6)
#f'$B={B:.0e}$ G'

ax.scatter(df_full["P0"], df_full["P1"], color='lightsteelblue', s=50, label='ATNF v2.6.1', zorder=1)
ax.scatter(df_subset["P0"], df_subset["P1"], color='forestgreen', edgecolor='black', s=50, label='Matched pulsars', zorder=2)

print(f'lenght of the subset: {len(df_subset)}')
ax.set_xlabel('Period (s)')
ax.set_xscale('log')
ax.set_ylabel('Period derivative (s/s)')
ax.set_yscale('log')

ax.set_xlim([1e-3, 1e2])
ax.set_ylim([6e-23, 1e-9])

# Count MSPs in the full dataset
msp_mask = (df_subset["P0"] < 0.03) & (df_subset["P1"] < 1e-16)
num_msp = msp_mask.sum()
print(f"Number of MSPs: {num_msp}")

ax.scatter(df_subset.loc[msp_mask, "P0"], df_subset.loc[msp_mask, "P1"], color='forestgreen', edgecolor='black', s=50, marker='s', label='Matched MSPs', zorder=2)
ax.legend(loc="lower right")

fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/p-pdot.pdf", bbox_inches="tight", dpi=300)

##################################################################################################
# scatter plot + dm 
##################################################################################################
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(23 * cm, 44 * cm))

valid_full = df_full[["W10", "P0"]].replace([np.inf, -np.inf], np.nan).dropna()
valid_full = valid_full[(valid_full["W10"] > 0) & (valid_full["P0"] > 0)]

valid_subset = df_subset[["W10", "P0"]].replace([np.inf, -np.inf], np.nan).dropna()
valid_subset = valid_subset[(valid_subset["W10"] > 0) & (valid_subset["P0"] > 0)]

#valid_full = df_full[["W10", "P0", "S1400"]].replace([np.inf, -np.inf], np.nan)
#valid_full = valid_full.dropna(subset=["W10", "P0"])
#valid_full = valid_full[(valid_full["W10"] > 0) & (valid_full["P0"] > 0)]

#valid_subset = df_subset[["W10", "P0", "S1400"]].replace([np.inf, -np.inf], np.nan)
#valid_subset = valid_subset.dropna(subset=["W10", "P0"])
#valid_subset = valid_subset[(valid_subset["W10"] > 0) & (valid_subset["P0"] > 0)]

#valid_full["S1400"] = valid_full["S1400"].fillna(1)
#valid_full = valid_full[valid_full["S1400"] > 0]

#valid_subset["S1400"] = valid_subset["S1400"].fillna(1)
#valid_subset = valid_subset[valid_subset["S1400"] > 0]

#valid_flux_full = df_full["S1400"].replace([np.inf, -np.inf], np.nan).dropna()
#valid_flux_full = valid_flux_full[(valid_flux_full > 0)]

#valid_flux_subset = df_subset["S1400"].replace([np.inf, -np.inf], np.nan).dropna()
#valid_flux_subset = valid_flux_subset[(valid_flux_subset > 0)]

frac_width_full = valid_full["W10"] / (valid_full["P0"]*1000)
frac_width_subset = valid_subset["W10"] / (valid_subset["P0"]*1000)

combined_frac = np.concatenate([frac_width_full, frac_width_subset])
min_val = combined_frac.min()
max_val = combined_frac.max()
bins = np.logspace(np.log10(min_val), np.log10(max_val), num=18)

ax1.hist(frac_width_full, label='ATNF v2.6.2', color='lightsteelblue', alpha=0.9, bins=bins, density=True)
ax1.hist(frac_width_subset, label='Matched pulsars', color='forestgreen', alpha=0.6, bins=bins, density=True)

#ax1.scatter(valid_full["S1400"], frac_width_full, label='ATNF v2.6.1', color='lightsteelblue', alpha=0.9)
#ax1.scatter(valid_subset["S1400"], frac_width_subset, label='Matched pulsars', color='forestgreen', alpha=0.8)

#ax1.axvline(np.median(valid_flux_full), color='darkblue', linestyle='-.', linewidth=1.5)
#ax1.axvline(np.median(valid_flux_subset), color='darkgreen', linestyle='--', linewidth=1.5)

#ax1.axhline(np.median(frac_width_full), color='darkblue', linestyle='-.', linewidth=1.5)
#ax1.axhline(np.median(frac_width_subset), color='darkgreen', linestyle='--', linewidth=1.5)

ax1.set_xlabel(r"W$_{10}$ / P")
ax1.set_ylabel("Fraction")
ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.legend()

valid_flux_full = df_full["S1400"].replace([np.inf, -np.inf], np.nan).dropna()
valid_flux_full = valid_flux_full[valid_flux_full > 0]

valid_flux_subset = df_subset["S1400"].replace([np.inf, -np.inf], np.nan).dropna()
valid_flux_subset = valid_flux_subset[valid_flux_subset > 0]

print(valid_flux_full)

combined_flux = np.concatenate([valid_flux_full, valid_flux_subset])
min_val = combined_flux.min()
max_val = combined_flux.max()
bins = np.logspace(np.log10(min_val), np.log10(max_val), num=18)

ax2.hist(valid_flux_full, label='ATNF v2.6.2', color='lightsteelblue', alpha=0.9, bins=bins, density=True)
ax2.hist(valid_flux_subset, label='Matched pulsars', color='forestgreen', alpha=0.6, bins=bins, density=True)

ax2.set_xlabel(r"Flux density [mJy]")
ax2.set_ylabel("Fraction")
ax2.set_xscale('log')
ax2.legend()

print('Done')

valid_dm_full = df_full["DM"].replace([np.inf, -np.inf], np.nan).dropna()
valid_dm_full = valid_dm_full[valid_dm_full > 0]

valid_dm_subset = df_subset["DM"].replace([np.inf, -np.inf], np.nan).dropna()
valid_dm_subset = valid_dm_subset[valid_dm_subset > 0]

combined_dm = np.concatenate([valid_dm_full, valid_dm_subset])
min_val = combined_dm.min()
max_val = combined_dm.max()
bins = np.logspace(np.log10(min_val), np.log10(max_val), num=16)

ax3.hist(valid_dm_full, label='ATNF v2.6.2', color='lightsteelblue', alpha=0.9, bins=bins, density=True)
ax3.hist(valid_dm_subset, label='Matched pulsars', color='forestgreen', alpha=0.6, bins=bins, density=True)

ax3.set_xlabel(r"DM [pc/cm$^3$]")
ax3.set_ylabel("Fraction")
ax3.set_xscale('log')
#ax3.set_xlim([0, 1500])
ax3.legend()
    
fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/DM-FLUX-FRAC.pdf", dpi=300)