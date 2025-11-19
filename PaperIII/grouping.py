#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from matplotlib.lines import Line2D

# === Plot config ===
plt.rcParams.update({
    "font.size": 11,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54
fig, ax = plt.subplots(figsize=(35*cm, 15*cm))

color = plt.cm.GnBu

# === Group drawing helper ===
group_height = 5
ra_span = 5

def draw_group(ax, ra_start, dec, num_groups, color_edge, color_fill):
    for i in range(num_groups):
        ra = ra_start + i * ra_span
        # 30° grey circle overlay
        circle = patches.Circle((ra + ra_span / 2, dec), 15, facecolor='grey', edgecolor='none', alpha=0.2, zorder=1)
        # Rectangle
        rect = patches.Rectangle((ra, dec - group_height / 2), ra_span, group_height, edgecolor=color_edge, facecolor=color_fill, zorder=2)
        ax.add_patch(rect)
        ax.text(ra + ra_span / 2, dec, " ", ha='center', va='center', zorder=3)
        ax.add_patch(circle)

# === Draw all groups ===
draw_group(ax, 185, -72, 5, color(0.95), color(0.85))
draw_group(ax, 140, -55, 22, color(0.85), color(0.75))
draw_group(ax, 124, -40.5, 4, color(0.75), color(0.65))
draw_group(ax, 253, -40.5, 4, color(0.75), color(0.65))
draw_group(ax, 115, -26.7, 3, color(0.65), color(0.55))
draw_group(ax, 266, -26.7, 3, color(0.65), color(0.55))
draw_group(ax, 274, -13, 3, color(0.55), color(0.45))
draw_group(ax, 284, +1.6, 2, color(0.45), color(0.35))

# === Galactic Plane (b=0°) ===
l_vals = np.linspace(0, 360, 1000)
b_vals = np.zeros_like(l_vals)
gal_plane = SkyCoord(l=l_vals * u.deg, b=b_vals * u.deg, frame='galactic')
equatorial = gal_plane.transform_to('fk5')
ra = equatorial.ra.deg
dec = equatorial.dec.deg
ax.plot(ra, dec, color='black', lw=1.5, linestyle='--', label='Galactic Plane (b=0°)')

# === Plot formatting ===
ax.set_xlabel("Right Ascension (°)")
ax.set_ylabel("Declination (°)")
ax.set_xlim(110, 300)
ax.set_ylim(-77, +5)
ax.invert_xaxis()

# === Legend ===
legend_elements = [
    Line2D([0], [0], color=color(0.35), lw=8, label='Dec = +1.6°'),
    Line2D([0], [0], color=color(0.45), lw=8, label='Dec = -13°'),
    Line2D([0], [0], color=color(0.55), lw=8, label='Dec = -26.7°'),
    Line2D([0], [0], color=color(0.65), lw=8, label='Dec = -40.5°'),
    Line2D([0], [0], color=color(0.75), lw=8, label='Dec = -55°'),
    Line2D([0], [0], color=color(0.85), lw=8, label='Dec = -72°'),
    Line2D([0], [0], color='black', lw=1.5, linestyle='--', label='Galactic Plane (b=0°)')
]
ax.legend(handles=legend_elements, loc='best')

# Convert RA axis from degrees to hours
xticks_deg = np.arange(300, 110-1, -30)  # Original deg ticks, inverted axis
xticks_hour = xticks_deg / 15  # Convert to hours
xtick_labels = [f"{int(h)}h" for h in xticks_hour]

ax.set_xticks(xticks_deg)
ax.set_xticklabels(xtick_labels)
ax.set_xlabel("Right Ascension (h)")

# === Output ===
plt.tight_layout()
plt.savefig("grouping.png", dpi=300, bbox_inches='tight')
plt.savefig("grouping.pdf", dpi=300, bbox_inches='tight')