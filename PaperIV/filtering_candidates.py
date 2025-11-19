#! /usr/bin/env python3

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
import pandas as pd

fits_file = "steep_fermi_matched_steep.fits"
with fits.open(fits_file) as hdul:
    data = hdul[1].data
    df = pd.DataFrame(data)

Simbad.add_votable_fields('otype')

# Object types to exclude
non_pulsar_types = {'Galaxy', 'AGN', 'QSO', 'BLLac', 'LINER', 'Seyfert_1', 'Seyfert_2', 'StarburstG', 'IRAS', 'RadioG', 'XrayG', 'OptG', 'GClstr', 'GravLensSystem', 'Star', 'HighPM*', 'LowPM*', 'WD*', 'BYDra', 'Flare*', 'V*', 'Eclipsing_binary', 'EB*', 'EB*Algol', 'EB*WUMa', 'EB*BetaLyr', 'RotV*', 'PulsV*', 'Pre-main sequence*', 'RSCVn', 'TTau*', 'S*', 'BlueStr*'}

# OTYPEs seen so far: {'RR*', 'AB?', 'X', 'HS?', 'G', 'cor', 'PN', '*', 'SB*', 'DNe', 'AGN', 'LP*', 'MoC', 'Rad', 'No*', 'OpC', 'PM*', 'pA?', 'BLL', 'Mi*', 'EB*', 'cC*', 'mm', 'RG*', 'Y*?', 'OH*', 'gam', 'cm', 'FIR', 'gB', 'dS*', 'Mas', '**', 'HII', 'gD*', 'El*', 'Y*O', 'WD?', 'IR', 'V*', 'Em*', 'bub', 'LP?', 'smm', 'Cl*', 'C*?'}

clean_candidates = []

#otype_seen = set()
with open("simbad_query_errors.log", "w") as error_log:
    for idx, row in df.iterrows():
        try:
            coord = SkyCoord(ra=row['ref_ra'] * u.deg, dec=row['ref_dec'] * u.deg)
            result = Simbad.query_region(coord, radius=1 * u.arcmin)
            #print(f"SIMBAD result: {result}")
            
            if result is None:
                print(f"No SIMBAD match at RA={row['ref_ra']}, Dec={row['ref_dec']}")
                clean_candidates.append(row)
            else:
                matched = False
                for obj in result:                    
                    obj_type = obj['otype'].decode() if isinstance(obj['otype'], bytes) else str(obj['otype'])
                    obj_name = obj['main_id'].decode() if isinstance(obj['main_id'], bytes) else str(obj['main_id'])

                    print(f"Found SIMBAD match: {obj_name} - {obj_type}")

                    if obj_type in non_pulsar_types:
                        matched = True
                        print(f"Probably an AGN or a star : Index {idx} | RA={row.get('ref_ra')} DEC={row.get('ref_dec')}")
                        break

                if not matched:
                    clean_candidates.append(row)

        except Exception as e:
            error_msg = f"Index {idx} | RA={row.get('ref_ra')} DEC={row.get('ref_dec')} | Error: {e}\n"
            print(error_msg)
            error_log.write(error_msg)

filtered_df = pd.DataFrame(clean_candidates)
filtered_table = Table.from_pandas(filtered_df)
filtered_table.write("filtered_fermi_candidates.fits", format="fits", overwrite=True)

print(f"Done. Saved {len(filtered_table)} likely pulsar candidates to FITS.")