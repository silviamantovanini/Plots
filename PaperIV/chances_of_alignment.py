#!/usr/bin/env python3

# Script to check what are the chances of alignment of pulsars with stars

import pandas as pd
import re
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.gaia import Gaia
import time
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12
})

cm = 1 / 2.54

def parse_atnf_txt(file_path):
    full_data = {"PSRJ": [], "RAJ": [], "DECJ": []}

    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(';')

            try:
                name = parts[1].strip()
                ra = parts[3].strip()
                dec = parts[6].strip()

                full_data["PSRJ"].append(name)
                full_data["RAJ"].append(ra)
                full_data["DECJ"].append(dec)

            except IndexError:
                continue

    df_full = pd.DataFrame(full_data)

    with open("matched_pulsars.txt", 'r') as f:
        matched_names = [line.strip() for line in f.readlines()]

    df_subset = df_full[df_full["PSRJ"].isin(matched_names)]

    return df_subset
            
def create_skycoords_from_df(df):
    coords = SkyCoord(ra=df['RAJ'].values, dec=df['DECJ'].values, unit=(u.hourangle, u.deg))
    return coords

def query_gaia_around(coord, radius):
    query = f"""SELECT COUNT(*) FROM gaiadr3.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {coord.ra.deg}, {coord.dec.deg}, {radius.to(u.deg).value})) = 1"""
    job = Gaia.launch_job(query)
    return job.results[0][0]

def main():
    atnf_file = "pulsar2_6_1.txt"
    search_radius = 30 * u.arcsec

    print("Parsing ATNF data...")
    atnf_df = parse_atnf_txt(atnf_file)
    coords_pulsars = create_skycoords_from_df(atnf_df)

    print(f"Querying Gaia for {len(coords_pulsars)} pulsars...")
    alignment_count = 0

    counts = []
    for i, coord in enumerate(coords_pulsars):
        try:
            count = query_gaia_around(coord, search_radius)
            if count > 0:
                alignment_count += 1
            #print(f"[{i+1}/{len(coords_pulsars)}] Pulsar: {atnf_df.iloc[i]['PSRJ']} - Matches: {count}")
            counts.append(count)
        except Exception as e:
            print(f"Query failed at {coord.to_string('hmsdms')}: {e}")
            counts.append(0)
            time.sleep(1)

    print(f"Total pulsars with at least one Gaia star within {search_radius}: {alignment_count}")

    # After collecting counts for each pulsar in a list called `counts`:
    fig, ax1 = plt.subplots(1, 1, figsize=(25 * cm, 16 * cm))
    
    ax1.hist(counts, bins=40, align='left', color='mediumblue', edgecolor="darkblue", alpha=0.7)
    ax1.set_xlabel('Number of Gaia Stars within 30 arcsec')
    ax1.set_ylabel('Number of Pulsars')

    fig.tight_layout()
    fig.savefig("alignment.png")

if __name__ == "__main__":
    main()