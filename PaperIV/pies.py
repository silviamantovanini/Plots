#!/usr/bin/env python3

########
# Pies #
########

import os
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 15
})

cm = 1 / 2.54

with open("/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/FINAL/final_matches.txt", "r") as f:
    pulsar_name_list = [line.split(';')[0].strip() for line in f if not line.startswith("#")]

pulsar_name_set = set(pulsar_name_list)

input_file_mine = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/FINAL/fit_results.csv"
df = pd.read_csv(input_file_mine)

df["Pulsar_Name"] = df["Pulsar_Name"].apply(lambda x: x.encode().decode("unicode_escape").strip("b'\""))
df = df[df["Pulsar_Name"].isin(pulsar_name_set)].copy()

########################################################################################################
# === PLOT 1: Pie of the models for the sources identified by me with the addition of the radio data ===
def func(pct, allvals):
    absolute = int(pct/100.*np.sum(allvals))
    return "{:.1f}%".format(pct, absolute)
      
fig, ax = plt.subplots(figsize=(18*cm, 15*cm))
model_counts = df['Best_Model'].value_counts()
models = model_counts.index
counts = model_counts.values
explode = (0, 0.1, 0, 0, 0)

colors = ['#081d58', '#980043', '#225ea8', '#016c59', '#014636']	#['#08519c', '#980043', '#2171b5', '#4292c6', '#6baed6']
x = np.arange(len(models))
pies = ax.pie(counts, colors=colors, autopct=lambda pct: func(pct, counts), radius=1.5, explode=explode, shadow=True, textprops={'color': 'white', 'fontsize': 18, 'weight': 'bold'})
ax.legend(['pl hard cut-off', 'simple pl', 'pl low turn-over', 'broken pl', 'double turn-over spectrum'], loc='upper left', bbox_to_anchor=(-0.35, 1.2), frameon=True)

for w in pies[0]:
    w.set_edgecolor("black")
    
fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/pie_my_low.pdf", dpi=300)
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/pie_my_low.png", dpi=300)

########################################################################################################
# === PLOT 2: Pie of the models for the sources identified by me with no low freq ===
input_file_all = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/FINAL/fit_results_nolow.csv"
df_no = pd.read_csv(input_file_all)

df_no["Pulsar_Name"] = df_no["Pulsar_Name"].apply(lambda x: x.encode().decode("unicode_escape").strip("b'\""))
    
fig, ax = plt.subplots(figsize=(18*cm, 15*cm))

model_counts = df_no['Best_Model'].value_counts()
model_counts = model_counts.reindex(index=['high_frequency_cut_off_power_law', 'simple_power_law', 'low_frequency_turn_over_power_law', 'broken_power_law'])
counts = model_counts.values
models = model_counts.index
explode = (0, 0.1, 0, 0)

colors = ['#081d58', '#980043', '#225ea8', '#016c59']
x = np.arange(len(models))
pies = ax.pie(counts, colors=colors, autopct=lambda pct: func(pct, counts), radius=1.5, explode=explode, shadow=True, textprops={'color': 'white', 'fontsize': 18, 'weight': 'bold'})
ax.legend(['pl hard cut-off', 'simple pl', 'pl low turn-over', 'broken pl'], loc='upper left', bbox_to_anchor=(-0.35, 1.2), frameon=True)	

for w in pies[0]:
    w.set_edgecolor("black")
    
fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/pie_my_before.pdf", dpi=300)
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/pie_my_before.png", dpi=300)

########################################################################################################
# === PLOT 3: Pie of the models for the sources in the ATNF catalogue (with my ones, (low freq)) ===
input_file_all = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/known_fit_results.csv"
df_all = pd.read_csv(input_file_all)

df_all["Pulsar_Name"] = df_all["Pulsar_Name"].apply(lambda x: x.encode().decode("unicode_escape").strip("b'\""))

df_combined = pd.concat([df, df_all], ignore_index=True)
    
fig, ax = plt.subplots(figsize=(18*cm, 15*cm))

model_counts = df_combined['Best_Model'].value_counts()
model_counts = model_counts.reindex(index=['high_frequency_cut_off_power_law', 'simple_power_law', 'low_frequency_turn_over_power_law', 'broken_power_law', 'double_turn_over_spectrum'])
counts = model_counts.values
models = model_counts.index
explode = (0, 0.1, 0, 0, 0)

colors = ['#081d58', '#980043', '#225ea8', '#016c59', '#014636']
x = np.arange(len(models))
ax.pie(counts, colors=colors, autopct=lambda pct: func(pct, counts), radius=1.5, explode=explode, shadow=True, textprops={'color': 'white', 'fontsize': 18, 'weight': 'bold'})
ax.legend(['pl hard cut-off', 'simple pl', 'pl low turn-over', 'broken pl', 'double turn-over spectrum'], loc='upper left', bbox_to_anchor=(-0.35, 1.2), frameon=True)	

fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/pie_all_low.pdf", dpi=300)
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/pie_all_low.png", dpi=300)

########################################################################################################
# === PLOT 4: Pie of the models for the sources in the ATNF catalogue (with my ones, (no, low freq)) ===
#input_file_before = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/known_fit_results_before.csv"
#df_before = pd.read_csv(input_file_before)

#df_before["Pulsar_Name"] = df_before["Pulsar_Name"].apply(lambda x: x.encode().decode("unicode_escape").strip("b'\""))

df_before = pd.concat([df_no, df_all], ignore_index=True)
    
fig, ax = plt.subplots(figsize=(18*cm, 15*cm))

model_counts = df_before['Best_Model'].value_counts()
model_counts = model_counts.reindex(index=['high_frequency_cut_off_power_law', 'simple_power_law', 'low_frequency_turn_over_power_law', 'broken_power_law', 'double_turn_over_spectrum'])
counts = model_counts.values
models = model_counts.index
explode = (0, 0.1, 0, 0, 0)

colors = ['#081d58', '#980043', '#225ea8', '#016c59', '#014636']
x = np.arange(len(models))
ax.pie(counts, colors=colors, autopct=lambda pct: func(pct, counts), radius=1.5, explode=explode, shadow=True, textprops={'color': 'white', 'fontsize': 18, 'weight': 'bold'})
ax.legend(['pl hard cut-off', 'simple pl', 'pl low turn-over', 'broken pl', 'double turn-over spectrum'], loc='upper left', bbox_to_anchor=(-0.35, 1.2), frameon=True)	

fig.tight_layout()
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/pie_all_before.pdf", dpi=300)
fig.savefig("/home/silvia/Desktop/Papers/Paper_IV/pie_all_before.png", dpi=300)