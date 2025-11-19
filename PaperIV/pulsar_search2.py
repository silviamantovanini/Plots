#! /usr/bin/env python3

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from regions import EllipseSkyRegion, PixCoord
from astropy.coordinates import search_around_sky
import astropy.units as u
from collections import defaultdict
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.spatial import cKDTree
import os
import subprocess
from pathlib import Path
from astroquery.skyview import SkyView
 
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12})

cm = 1/2.54

def load_txt_catalog(txt_path):
    # Loading the pulsar catalogue (saved in txt format)
    ra_list = []
    dec_list = []
    name_list = []
    ra_err_list = []
    dec_err_list = []
    
    with open(txt_path, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                parts = line.split(';')
                name = parts[1].strip()
                ra_str = parts[3].strip()
                ra_err = float(parts[4].strip())*15
                dec_str = parts[6].strip()
                dec_err = float(parts[7].strip())*15

                #total_err = (ra_err**2 + dec_err**2)**0.5

                if ra_err > 60 or dec_err > 60:
                    continue
                    
                name_list.append(name)
                ra_list.append(ra_str)
                dec_list.append(dec_str)
                ra_err_list.append(ra_err)
                dec_err_list.append(dec_err)
                
            except IndexError:
                print(f"Skipping malformed line: {line}")
                continue
    print(len(ra_list))
    coords = SkyCoord(ra=ra_list, dec=dec_list, unit=(u.hourangle, u.deg), frame='fk5')
    return coords, np.array(name_list)

def load_fits_catalog(fits_path):
    # Loading my catalogue (saved in fits format)
    # All the sources have a S/N > 5, and are limited in the GP area 
    with fits.open(fits_path) as hdul:
        data = hdul[1].data
        ra = data['ref_ra']
        dec = data['ref_dec']
        spix = data['pl_alpha']
        snr = data['peak_flux'] / data['local_rms']
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    return coords, spix, snr, data

def match_catalogs(cat1_coords, cat1_names, cat2_coords, cat2_data, max_sep_arcmin=2.0):
    # Cross-matching the catalogues - return match if within 2 arcminutes 
    idx, sep2d, _ = cat1_coords.match_to_catalog_sky(cat2_coords)
    matched = sep2d < (max_sep_arcmin * u.arcmin)

    matched_indices_in_cat1 = np.where(matched)[0]
    matched_cat1_names = cat1_names[matched]
    matched_cat2_indices = idx[matched]

    cat1_matched = cat1_coords[matched]
    cat1_unmatched = cat1_coords[~matched]

    cat2_table = Table(cat2_data)
    matched_cat2_data = cat2_table[matched_cat2_indices]

    return cat1_matched, cat1_unmatched, matched_cat2_data, matched_cat1_names, matched_cat2_indices

def plot_cutout_with_source(mosaic_path, pulsar_coords, names, matched_catalog_coords, output_dir="cutout_plots", size_deg=0.2):
    # Plotting the matched sources for visual inspection - I'll have to discard all sources that are probably artefacts
    os.makedirs(output_dir, exist_ok=True)
    
    with fits.open(mosaic_path) as hdul:
        wcs = WCS(hdul[0].header)
        data = hdul[0].data

        for i, coord_eq in enumerate(pulsar_coords):
            name = names[i]
            matched_eq = matched_catalog_coords[i]

            coord_gal = coord_eq.galactic
            matched_gal = matched_eq.galactic

            framesize = u.Quantity((0.2, 0.2), u.deg)
            cutout = Cutout2D(data, coord_gal, framesize, wcs=wcs)
            
            ra_ha = coord_eq.ra.to(u.hourangle)
            dec_deg = coord_eq.dec.to(u.deg)

            coord_eq_ha = SkyCoord(ra=ra_ha.value, dec=dec_deg.value, unit=(u.hourangle, u.deg))
            position_str = coord_eq_ha.to_string("hmsdms", sep=" ")

            wise_fits_path = f"wise_images/{name}_WISE34_fk5.fits"
            wise_file = Path(wise_fits_path)

            if not wise_file.exists():
                print(f"{wise_fits_path} not found. Downloading...")

                try:
                    images = SkyView.get_images(position=coord_eq_ha, survey=["WISE 3.4"], coordinates="J2000", projection="Sin", pixels=1309, height=0.5*u.deg, width=0.5*u.deg)
                    images[0].writeto(wise_fits_path, overwrite=True)
                    #subprocess.run(["perl", "/home/silvia/Desktop/Papers/Paper_I/skvbatch_wget.pl", f"file={wise_fits_path}", f"position={position_str}", f"Survey=WISE3.4", f"Coordinates=J2000", f"Projection=SIN", f"Size=0.5", f"Pixels=1309"], check=True)
                except Exception as e:
                    print(f"Download failed:\n{e}")
                    continue
                
            try:
                with fits.open(wise_fits_path) as wise_hdul:
                    wise_wcs = WCS(wise_hdul[0].header)
                    wise_data = wise_hdul[0].data
                    wise_cutout = Cutout2D(wise_data, coord_eq, framesize, wcs=wise_wcs)
            except Exception as e:
                print(f"Error reading WISE image for {name}: {e}")
                continue

            fig = plt.figure(figsize=(30*cm, 15*cm))
            ax1 = fig.add_subplot(1, 2, 1, projection=cutout.wcs)
            ax2 = fig.add_subplot(1, 2, 2, projection=wise_cutout.wcs)
            
            im1 = ax1.imshow(cutout.data, origin='lower', cmap='inferno', vmin=np.nanpercentile(cutout.data, 5), vmax=np.nanpercentile(cutout.data, 99))   
            #cb = fig.colorbar(im, cax=ax, orientation='vertical')
            ax1.scatter([coord_gal.l.deg], [coord_gal.b.deg], transform=ax1.get_transform('world'), s=50, edgecolor='cyan', facecolor='none', linewidth=1.5, label='Pulsar')
            ax1.scatter([matched_gal.l.deg], [matched_gal.b.deg], transform=ax1.get_transform('world'), s=50, color='red', marker='x', linewidth=1.5, label='Matched Source')

            ax1.set_xlabel("Galactic Longitude (deg)")
            ax1.set_ylabel("Galactic Latitude (deg)")
            ax1.grid(color='white', ls='dotted')
            ax1.invert_xaxis()
            ax1.legend(loc='upper right')

            im2 = ax2.imshow(wise_cutout.data, origin='lower', cmap='gray', vmin=np.nanpercentile(wise_cutout.data, 5), vmax=np.nanpercentile(wise_cutout.data, 99))
            ax2.scatter([coord_eq.ra.deg], [coord_eq.dec.deg], transform=ax2.get_transform('world'), s=50, edgecolor='cyan', facecolor='none', linewidth=1.5, label='Pulsar')
            ax2.scatter([matched_eq.ra.deg], [matched_eq.dec.deg], transform=ax2.get_transform('world'), s=50, color='red', marker='x', linewidth=1.5, label='Matched Source')
                    
            ax2.set_xlabel("Right Ascension")
            ax2.set_ylabel("Declination")
            ax2.grid(color='white', ls='dotted')
            ax2.invert_xaxis()
            ax2.legend(loc='upper right')
                
            fig.tight_layout()
            fig.savefig(f"{output_dir}/{name.replace(' ', '_')}_cutout.png", dpi=200)
            plt.close(fig)
        
def save_region_file(coords, filename="matched_sources.reg"):
    # Save sources in a region file 
    with open(filename, 'w') as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=cyan dashlist=8 3 width=1 font='helvetica 10 normal' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
        f.write("fk5\n")
        for c in coords:
            f.write(f"circle({c.ra.deg},{c.dec.deg},1')\n")

def save_matched_fits(data, names=None, filename="matched_sources.fits"):
    # Save sources in a fits file - adding the name of the pulsar to which they corrrespond
    table = Table(data)
    if names is not None:
        table['Pulsar_Name'] = names
    table.write(filename, format='fits', overwrite=True)

def plot_sources(matched, unmatched, extra=None, filename="matched_plot.pdf"):
    # Plotting the whole pulsar population, the matched sources and all the compact sources left in my catalogue
    fig = plt.figure(figsize=(20*cm, 13*cm))
    ax = fig.add_subplot(1, 1, 1)

    ax.scatter(unmatched.ra.deg, unmatched.dec.deg, color='grey', s=15, label='Unmatched pulsars')
    ax.scatter(matched.ra.deg, matched.dec.deg, color='darkblue', s=15, marker='+', label='Matched pulsars')

    if extra is not None and len(extra) > 0:
        ax.scatter(extra.ra.deg, extra.dec.deg, color='magenta', s=15, marker='x', label='Steep spectrum')

    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.legend()
    ax.set_title('Pulsar Matching and Steep Spectrum Search')
    ax.grid(True)
    ax.invert_xaxis()
    fig.tight_layout()
    fig.savefig(filename)

def is_within_ellipse(source_coord, fermi_coord, semi_major, semi_minor, pos_angle_deg):
    # Check if coordinates fall within ellipse
    
    # Offset in degrees using small-angle approximation
    dra = (source_coord.ra - fermi_coord.ra).to(u.deg).value * np.cos(np.deg2rad(fermi_coord.dec.deg))
    ddec = (source_coord.dec - fermi_coord.dec).to(u.deg).value

    # Rotate offsets by -theta
    theta = np.deg2rad(pos_angle_deg)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    x_rot = dra * cos_t + ddec * sin_t
    y_rot = -dra * sin_t + ddec * cos_t

    # Elliptical distance check
    return (x_rot / semi_major)**2 + (y_rot / semi_minor)**2 <= 1.0  

def match_steep_to_fermi(steep_coords, fermi_coords, fermi_table, max_search_radius_deg=2):
    matches = []

    idx_steep, idx_fermi, sep2d, _ = search_around_sky(steep_coords, fermi_coords, max_search_radius_deg * u.deg)

    for i, j in zip(idx_steep, idx_fermi):
        row = fermi_table[j]
        semi_major = row['Conf_95_SemiMajor']
        semi_minor = row['Conf_95_SemiMinor']
        pos_angle = row['Conf_95_PosAng']

        if np.isnan(semi_major) or np.isnan(semi_minor) or np.isnan(pos_angle):
            continue

        if is_within_ellipse(steep_coords[i], fermi_coords[j], semi_major, semi_minor, pos_angle):
            matches.append((i, j))

    return matches

def plot_fermi_ellipse_with_match(fermi_row, fermi_coord, steep_coords, index, output_dir="fermi_match_plots"):
    os.makedirs(output_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(18*cm, 15*cm))

    t = np.linspace(0, 2 * np.pi, 200)

    for conf, color, alpha in [(68, 'blue', 0.5), (95, 'red', 0.3)]:
        semi_major = fermi_row[f'Conf_{conf}_SemiMajor']
        semi_minor = fermi_row[f'Conf_{conf}_SemiMinor']
        pos_angle = fermi_row[f'Conf_{conf}_PosAng']

        if np.isnan(semi_major) or np.isnan(semi_minor) or np.isnan(pos_angle):
            print(f"Skipping source {fermi_row['Source_Name']} because of NaN values")
            continue  

        if isinstance(semi_major, (list, np.ndarray)):
            print(f"Invalid type for semi_major, skipping source {fermi_row['Source_Name']}")
            continue

        theta = np.deg2rad(pos_angle)

        x = semi_major * np.cos(t)
        y = semi_minor * np.sin(t)

        x_rot = x * np.cos(theta) - y * np.sin(theta)
        y_rot = x * np.sin(theta) + y * np.cos(theta)

        ra_ellipse = fermi_coord.ra.deg + x_rot / np.cos(np.deg2rad(fermi_coord.dec.deg))
        dec_ellipse = fermi_coord.dec.deg + y_rot

        ax.plot(ra_ellipse, dec_ellipse, color=color, alpha=alpha, label=f'{conf}% Confidence Ellipse')

    ax.scatter(fermi_coord.ra.deg, fermi_coord.dec.deg, c='black', marker='x', label='Fermi Source')

    matched_coords = []
    for steep_src in steep_coords:
        if is_within_ellipse(steep_src, fermi_coord, semi_major, semi_minor, pos_angle):
            matched_coords.append(steep_src)

    if matched_coords:
        ax.scatter([s.ra.deg for s in matched_coords], [s.dec.deg for s in matched_coords], c='magenta', marker='o', label='Matched Sources')

    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")
    ax.set_title(f"{fermi_row['Source_Name']}")
    ax.invert_xaxis()
    ax.legend(loc='upper right')
    ax.grid(True)

    fig.tight_layout()
    fig.savefig(f"{output_dir}/{fermi_row['Source_Name']}.png")
    plt.close(fig)

def main(txt_path, fits_path):
    print("Loading catalogs...")
    cat1_coords, cat1_names = load_txt_catalog(txt_path)
    cat2_coords, spix2, snr2, cat2_data = load_fits_catalog(fits_path)

    print("Matching catalogs...")
    matched, unmatched, matched_cat2_data, matched_names, matched_cat2_indices = match_catalogs(cat1_coords, cat1_names, cat2_coords, cat2_data)

    print(f"Matched sources: {len(matched)}")

    print("Creating cutout plots around matched pulsars...")
    plot_cutout_with_source('/home/silvia/GLEAM-X/Joint_Deconvolution/GALACTICPLANE/Mosaics/Round7/release/JDGP_mosaic_170-231MHz_MFS.fits', matched, matched_names, cat2_coords[matched_cat2_indices])

    print("Saving matched region file and FITS...")
    save_region_file(matched)
    save_matched_fits(matched_cat2_data, names=matched_names)

    print("Creating unmatched catalog for steep-spectrum search...")
    cat2_table = Table(cat2_data)
    all_indices = np.arange(len(cat2_coords))
    unmatched_mask = np.ones(len(cat2_coords), dtype=bool)
    unmatched_mask[matched_cat2_indices] = False

    cat2_unmatched_table = cat2_table[unmatched_mask]

    steep_mask = (cat2_unmatched_table['int_flux'] / cat2_unmatched_table['peak_flux'] < 1.5)	#(cat2_unmatched_table['pl_alpha'] < -1.8) &
    steep_sources_table = cat2_unmatched_table[steep_mask]

    print(f"New steep-spectrum candidates: {len(steep_sources_table)}")

    print("Saving steep-spectrum region and FITS file...")
    extra_coords = SkyCoord(ra=steep_sources_table['ref_ra'], dec=steep_sources_table['ref_dec'], unit='deg')
    save_region_file(extra_coords, filename="steep_sources.reg")
    save_matched_fits(steep_sources_table, filename="steep_sources.fits")

    print("Plotting results...")
    plot_sources(matched, unmatched, extra=extra_coords)

    fermi_data = fits.getdata('fermi_compact_unassoc_PL_var_LP.fits', 1)
    fermi_table = Table(fermi_data)
    #non_variable_mask = (fermi_table['Variability_Index'] <= 18.48) & (fermi_table['PL_Index'] <= 2.2)
    #fermi_table = fermi_table[non_variable_mask]
    
    fermi_coords = SkyCoord(ra=fermi_table['RAJ2000'] * u.deg, dec=fermi_table['DEJ2000'] * u.deg)

    print("Cross-matching steep-spectrum sources with Fermi...")
    matches = match_steep_to_fermi(extra_coords, fermi_coords, fermi_table)
    
    print("Plotting ellipses for matched Fermi sources...")

    fermi_to_steep = defaultdict(list)
    for i, j in matches:
        fermi_to_steep[j].append(i)

    for j, steep_idxs in fermi_to_steep.items():
        steep_coords_subset = extra_coords[steep_idxs]
        fermi_row = fermi_table[j]
        fermi_c = fermi_coords[j]
        plot_fermi_ellipse_with_match(fermi_row, fermi_c, steep_coords_subset, j)

    # Save results
    matched_steep_indices = [i for i, _ in matches]
    matched_fermi_indices = [j for _, j in matches]

    steep_matches = steep_sources_table[matched_steep_indices]
    fermi_matches = fermi_table[matched_fermi_indices]

    steep_coords_matched = extra_coords[matched_steep_indices]
    
    save_region_file(steep_coords_matched, filename="steep_fermi_matches.reg")
    save_matched_fits(steep_matches, filename="steep_fermi_matched_steep.fits")
    save_matched_fits(fermi_matches, filename="steep_fermi_matched_fermi.fits")

    print(f"Steep sources matched to Fermi: {len(matches)}")

    print("Plotting results...")
    plot_sources(matched, unmatched, extra=steep_coords_matched, filename='pulsar&steep.pdf')
    
main('pulsar2_6_0.txt', 'merged_catalogue_lon_SEDs.fits')