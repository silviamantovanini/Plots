#! /usr/bin/env python3

import os
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 15
})

cm = 1 / 2.54

def extract_alpha(param_str):
    try:
        params = ast.literal_eval(param_str)
        return params.get('a')
    except (ValueError, SyntaxError):
        return None

def extract_params(param_str):
    try:
        return ast.literal_eval(param_str)
    except:
        return {}

# Consider only the pulsars that have not been discarded
with open("/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/FINAL/final_matches.txt", "r") as f:
    pulsar_name_list = [line.split(';')[0].strip() for line in f if not line.startswith("#")]
    #[line.strip() for line in f if line.strip() and not line.startswith("#")]

pulsar_name_set = set(pulsar_name_list)

input_file = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/FINAL/fit_results.csv"
df = pd.read_csv(input_file)

df["Pulsar_Name"] = df["Pulsar_Name"].apply(lambda x: x.encode().decode("unicode_escape").strip("b'\""))
df = df[df["Pulsar_Name"].isin(pulsar_name_set)].copy()

plot_dir = "paper_plots"

# === PLOT 1: Histogram for all the model used ===
def func(pct, allvals):
    absolute = int(pct/100.*np.sum(allvals))
    return "{:.1f}%".format(pct, absolute)
    
fig, ax = plt.subplots(figsize=(18*cm, 15*cm))
model_counts = df['Best_Model'].value_counts()
models = model_counts.index
counts = model_counts.values

print(model_counts)

#cmap = plt.get_cmap("GnBu")
#colors = [cmap((i+1.2) / len(models)) for i in range(len(models))]
colors = ['#08519c', '#980043', '#2171b5', '#4292c6', '#6baed6']
x = np.arange(len(models))
ax.pie(counts, colors=colors, autopct=lambda pct: func(pct, counts), radius=1.5)

#ax.bar(x, counts, color=colors)

#ax.set_ylabel("Count")
#ax.set_xlabel("Model")
#ax.set_xticks(x)
#ax.set_xticklabels([])
#ax.tick_params(axis='x', length=0)

#handles = [mpatches.Patch(color=colors[i], label=model) for i, model in enumerate(models)]
#ax.legend(handles=handles, loc='best')

ax.legend(['pl hard cut-off', 'simple pl', 'pl low turn-over', 'broken pl', 'double turn-over spectrum'], loc='upper left', bbox_to_anchor=(-0.35, 1.2), frameon=True)

fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/model_histogram.pdf", dpi=300)
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/model_histogram.png", dpi=300)

# === PLOT 2: Alpha Histogram for Simple Power Law ===
powerlaw_df = df[df["Best_Model"] == "simple_power_law"].copy()
powerlaw_df["Spectral_Index"] = powerlaw_df["Parameters"].apply(extract_alpha)
powerlaw_df = powerlaw_df.dropna(subset=["Spectral_Index"])

print(f"Mean pl spectral index: {np.mean(powerlaw_df['Spectral_Index'])} $\pm$ {np.std(powerlaw_df['Spectral_Index'])}")
fig, ax = plt.subplots(figsize=(18 * cm, 16 * cm))
ax.hist(powerlaw_df["Spectral_Index"], bins=20, color="#2b8cbe", edgecolor="#08589e", alpha=0.7)
ax.set_xlabel(r"Spectral Index ($\alpha$)")
ax.set_ylabel("Count")
ax.grid(True)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "simple_powerlaw.png"))

# === PLOT 3: High-Frequency Cutoff ===
high_cutoff_df = df[df['Best_Model'] == 'high_frequency_cut_off_power_law'].copy()
high_cutoff_df['params'] = high_cutoff_df['Parameters'].apply(extract_params)
high_cutoff_df['alpha'] = high_cutoff_df['params'].apply(lambda x: x.get('a'))
high_cutoff_df['nu_c'] = high_cutoff_df['params'].apply(lambda x: x.get('vc'))
high_cutoff_df = high_cutoff_df.dropna(subset=['alpha', 'nu_c'])

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(54 * cm, 16 * cm))
ax1.scatter(high_cutoff_df['nu_c'] / 1e6, high_cutoff_df['alpha'], color='#2b8cbe', s=60)
ax1.set_xlabel("Cutoff Frequency (nu) [MHz]")
ax1.set_ylabel(r"Spectral Index ($\alpha$)")
ax1.grid(True)

ax2.hist(high_cutoff_df['nu_c'] / 1e6, bins=20, color="#2b8cbe", edgecolor="#08589e", alpha=0.7)
ax2.set_xlabel("Cutoff Frequency (nu) [MHz]")
ax2.set_ylabel("Count")
ax2.grid(True)

ax3.hist(high_cutoff_df['alpha'], bins=20, color="#2b8cbe", edgecolor="#08589e", alpha=0.7)
ax3.set_xlabel(r"Spectral Index ($\alpha$)")
ax3.set_ylabel("Count")
ax3.grid(True)

fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "high_freq_cutoff.png"))

# === PLOT 4: Low-Frequency Turnover ===
turnover_df = df[df['Best_Model'] == 'low_frequency_turn_over_power_law'].copy()
turnover_df['params'] = turnover_df['Parameters'].apply(extract_params)
turnover_df['alpha'] = turnover_df['params'].apply(lambda x: x.get('a'))
turnover_df['nu_turnover'] = turnover_df['params'].apply(lambda x: x.get('vpeak'))
turnover_df['curv'] = turnover_df['params'].apply(lambda x: x.get('beta'))
turnover_df = turnover_df.dropna(subset=['alpha', 'nu_turnover', 'curv'])

fig, axs = plt.subplots(2, 2, figsize=(36 * cm, 36 * cm))
ax1 = axs[0, 0]
ax2 = axs[0, 1]
ax3 = axs[1, 0]
ax4 = axs[1, 1]

ax1.scatter(turnover_df['nu_turnover'] / 1e6, turnover_df['alpha'], color='#2b8cbe', s=60)
ax1.set_xlabel(r"Turnover Frequency ($\nu$) [MHz]")
ax1.set_ylabel(r"Spectral Index ($\alpha$)")
ax1.grid(True)

ax2.hist(turnover_df['nu_turnover'] / 1e6, bins=20, color="#2b8cbe", edgecolor="#08589e", alpha=0.7)
ax2.set_xlabel(r"Turnover Frequency ($\nu$) [MHz]")
ax2.set_ylabel("Count")
ax2.grid(True)

ax3.hist(turnover_df['alpha'], bins=20, color="#2b8cbe", edgecolor="#08589e", alpha=0.7)
ax3.set_xlabel(r"Spectral Index ($\alpha$)")
ax3.set_ylabel("Count")
ax3.grid(True)

ax4.scatter(turnover_df['curv'], turnover_df['alpha'], color='#2b8cbe', s=60)
ax4.set_xlabel(r"Curvature ($\beta$)")
ax4.set_ylabel(r"Spectral Index ($\alpha$)")
ax4.grid(True)

fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "turnover_lowfreq.png"))

# === PLOT 5: Broken Power Law ===
broken_df = df[df['Best_Model'] == 'broken_power_law'].copy()
broken_df['params'] = broken_df['Parameters'].apply(extract_params)
broken_df['alpha_1'] = broken_df['params'].apply(lambda x: x.get('a1'))
broken_df['alpha_2'] = broken_df['params'].apply(lambda x: x.get('a2'))
broken_df['broken_freq'] = broken_df['params'].apply(lambda x: x.get('vb'))
broken_df = broken_df.dropna(subset=['alpha_1', 'alpha_2', 'broken_freq'])

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(54 * cm, 16 * cm))

ax1.hist(broken_df['alpha_1'], bins=20, color="#2b8cbe", edgecolor="#08589e", label=r"$\alpha_1$", alpha=0.6)
ax1.hist(broken_df['alpha_2'], bins=20, color="#7bccc4", edgecolor="#4eb3d3", label=r"$\alpha_2$", alpha=0.6)
ax1.set_ylabel("Count")
ax1.legend()
ax1.grid(True)

ax2.scatter(broken_df['alpha_1'], broken_df['alpha_2'], color="#2b8cbe", s=60)
ax2.set_xlabel(r"$\alpha_1$")
ax2.set_ylabel(r"$\alpha_2$")
ax2.grid(True)

bins = np.logspace(np.log10(min(broken_df['broken_freq'] / 1e6)), np.log10(max(broken_df['broken_freq'] / 1e6)), 15)
ax3.hist(broken_df['broken_freq'] / 1e6, bins=bins, color='#2b8cbe', edgecolor="#08589e", alpha=0.6)
ax3.set_xlabel(r"$\nu_b$ [MHz]")
ax3.set_ylabel("Count")
ax3.set_xscale("log")
ax3.grid(True)

fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "broken_power_law.png"))

# === PLOT 6: Double Turn Over ===
double_turnover_df = df[df['Best_Model'] == 'double_turn_over_spectrum'].copy()
double_turnover_df['params'] = double_turnover_df['Parameters'].apply(extract_params)
double_turnover_df['turnover_1'] = double_turnover_df['params'].apply(lambda x: x.get('vc'))
double_turnover_df['turnover_2'] = double_turnover_df['params'].apply(lambda x: x.get('vpeak'))
double_turnover_df = double_turnover_df.dropna(subset=['turnover_1', 'turnover_2'])

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(54 * cm, 16 * cm))

ax1.hist(double_turnover_df['turnover_1'] / 1e6, bins=20, color='#2b8cbe', edgecolor="#08589e", label=r'High-Freq Turnover ($\nu_1$)', alpha=0.7)
ax1.set_xlabel('Frequency [MHz]')
ax1.set_ylabel('Count')
ax1.legend()
ax1.grid(True)

ax2.scatter(double_turnover_df['turnover_1'] / 1e6, double_turnover_df['turnover_2'] / 1e6, color="#2b8cbe", s=60)
ax2.set_xlabel(r"High-Freq Turnover ($\nu_1$)")
ax2.set_ylabel(r"Low-Freq Turnover ($\nu_2$)")
ax2.set_ylim([0, 350])
ax2.grid(True)

ax3.hist(double_turnover_df['turnover_2'] / 1e6, bins=20, color='#7bccc4', edgecolor="#4eb3d3", label=r'Low-Freq Turnover ($\nu_2$)', alpha=0.7)
ax3.set_xlabel('Frequency [MHz]')
ax3.set_ylabel('Count')
ax3.legend()
ax3.grid(True)

fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "double_turn_over.png"))

# === PLOT 7: DM (dispersion measure) and period distributions of the known pulsar detected ===
name_list = []
dm_list = []
period_list = []
dm_err_list = []
period_err_list = []
flux_list = []

with open("/home/silvia/Desktop/Papers/Project_4/Pulsars/pulsar2_6_1.txt", 'r') as f:
    for i, line in enumerate(f):
        if i < 2:
            continue
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            parts = line.split(';')
            name = parts[1].strip()
            dm_str = parts[15].strip()
            period_str = parts[9].strip()
            dm_err_str = parts[16].strip()
            period_err_str = parts[10].strip()
            flux_str = parts[24].strip()	# at 1400
            
            if not dm_str or not period_str or dm_str == '*' or period_str == '*' or not flux_str or flux_str == '*':
                continue
                
            dm = float(dm_str)
            dm_err = float(dm_err_str)
            period = float(period_str)
            period_err = float(period_err_str)
            flux = float(flux_str)
            
            name_list.append(name)
            dm_list.append(dm)
            dm_err_list.append(dm_err)
            period_list.append(period)
            period_err_list.append(period_err)
            flux_list.append(flux)
        except IndexError:
            print(f"Skipping malformed line: {line}")
            continue

df_names = set(df["Pulsar_Name"].dropna().str.strip())
matched_dm = [dm for name, dm in zip(name_list, dm_list) if name in df_names]
matched_dm_err = [dm_err for name, dm_err in zip(name_list, dm_err_list) if name in df_names]
matched_period = [p for name, p in zip(name_list, period_list) if name in df_names]
matched_period_err = [p_err for name, p_err in zip(name_list, period_err_list) if name in df_names]
matched_flux = [fl for name, fl in zip(name_list, flux_list) if name in df_names]

# --- Histogram of Dispersion Measure ---
fig1 = plt.figure(figsize=(18*cm, 15*cm))
ax1 = fig1.add_subplot(1, 1, 1)

ax1.hist(dm_list, bins=40, color='grey', label='ATNF v2.6.1 pulsars', alpha=0.7, density=False)
ax1.hist(matched_dm, bins=20, color="#084081", label='Matched pulsars', alpha=0.7, density=False)
ax1.set_xlabel('Dispersion Measure')
ax1.set_ylabel('Count')
ax1.set_xlim([0, 1250])
ax1.legend()
fig1.tight_layout()
fig1.savefig(f"{plot_dir}/dm.pdf", dpi=300)

# --- Errorbar plot ---
fig2 = plt.figure(figsize=(18*cm, 15*cm))
ax2 = fig2.add_subplot(1, 1, 1)

ax2.errorbar(dm_list, period_list, xerr=dm_err_list, yerr=period_err_list, fmt='o', linestyle='none', color="grey", label="ATNF v2.6.1 pulsars")
ax2.errorbar(matched_dm, matched_period, xerr=matched_dm_err, yerr=matched_period_err, fmt='o', linestyle='none', markerfacecolor="#084081", markeredgecolor='black', label="Matched pulsars")
ax2.set_xlabel(r"Dispersion Measure")
ax2.set_ylabel(r"Period")
ax2.set_yscale('log')
ax2.legend()
fig2.tight_layout()
fig2.savefig(f"{plot_dir}/fig_errorbar_dm_period.png")

# --- DM-FLUX ---
fig = plt.figure(figsize=(18*cm, 15*cm))
ax = fig.add_subplot(1, 1, 1)

ax.scatter(flux_list, dm_list, color="grey", s=20, label="ATNF v2.6.1 pulsars")
ax.scatter(matched_flux, matched_dm, color='#084081', edgecolor='black', s=20, label="Matched pulsars")
ax.set_xlabel("Flux density (mJy)")
ax.set_xscale('log')
ax.set_ylabel(r"Dispersion measure pc cm^(-3)")
ax.legend()
fig.tight_layout()
fig.savefig(f"{plot_dir}/dm_flux.pdf", dpi=300)

# --- FLUX / DM / PERIOD ---
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(54 * cm, 15 * cm))

ax1.hist(dm_list, bins=20, color="grey", label='ATNF v2.6.1 pulsars', alpha=0.7)
ax1.hist(matched_dm, bins=20, color="#084081", label='Matched pulsars', alpha=0.7)
ax1.set_xlabel(r'DM [pc/cm$^3$')
ax1.set_ylabel('Count')
ax1.legend()

bins_all = np.logspace(np.log10(min(flux_list)), np.log10(max(flux_list)), 20)
bins = np.logspace(np.log10(min(matched_flux)), np.log10(max(matched_flux)), 20)

ax2.hist(flux_list, bins=bins_all, color="grey", label='ATNF v2.6.1 pulsars', alpha=0.7, density=True)
ax2.hist(matched_flux, bins=bins, color="#084081", label='Matched pulsars', alpha=0.7, density=True)
ax2.set_xlabel('Flux [mJy]')
ax2.set_xscale('log')
ax2.set_ylabel('Count')
ax2.legend()

bins_all = np.logspace(np.log10(min(period_list)), np.log10(max(period_list)), 20)
bins = np.logspace(np.log10(min(matched_period)), np.log10(max(matched_period)), 20)

ax3.hist(period_list, bins=bins_all, color="grey", label='ATNF v2.6.1 pulsars', alpha=0.7)
ax3.hist(matched_period, bins=bins, color="#084081", label='Matched pulsars', alpha=0.7)
ax3.set_xlabel('P [s]')
ax3.set_xscale('log')
ax3.set_ylabel('Count')
ax3.legend()

fig.tight_layout()
fig.savefig(f"{plot_dir}/dm-p-flux.pdf", dpi=300)

################################################################################
fig = plt.figure(figsize=(20*cm, 20*cm))
gs = gridspec.GridSpec(4, 4)
ax_main = fig.add_subplot(gs[1:4, 0:3])
ax_xhist = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
ax_yhist = fig.add_subplot(gs[1:4, 3], sharey=ax_main) 

plt.setp(ax_xhist.get_xticklabels(), visible=False)
plt.setp(ax_yhist.get_yticklabels(), visible=False)

ax_main.grid(True)
ax_main.errorbar(dm_list, period_list, xerr=dm_err_list, yerr=period_err_list, fmt='.', color="#2b8cbe", label="ATNF pulsars v2.6.1", alpha=0.7)
ax_main.errorbar(matched_dm, matched_period, xerr=matched_dm_err, yerr=matched_period_err, fmt='.', color="#7bccc4", label="Matched pulsars", alpha=0.7)
ax_main.set_xlabel("Dispersion Measure")
ax_main.set_ylabel("Period [s]")
ax_main.set_yscale('log')
ax_main.legend(loc='best')

ax_xhist.grid(True)
ax_xhist.hist(dm_list, bins=20, color='#2b8cbe', edgecolor="#08589e", alpha=0.7)
ax_xhist.hist(matched_dm, bins=20, color='#7bccc4', edgecolor="#4eb3d3", alpha=0.7)

log_bins = np.logspace(np.log10(min(period_list)), np.log10(max(period_list)), num=20)
                     
ax_yhist.grid(True)
ax_yhist.hist(period_list, bins=log_bins, orientation='horizontal', color='#2b8cbe', edgecolor="#08589e", alpha=0.7)
ax_yhist.hist(matched_period, bins=log_bins, orientation='horizontal', color='#7bccc4', edgecolor="#4eb3d3", alpha=0.7)
ax_yhist.set_yscale('log')

fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "dm_period_joint.png"))