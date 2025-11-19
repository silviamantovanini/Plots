#! /usr/bin/env python3

from astropy.io import fits
from astropy.table import Table

import pandas as pd

###############################################
# GOALS OF THE SCRIPT
# 1 - FILTER OUT SOURCES WITH FLAGS: 1, 2, 5, 7
###############################################

matches_catalogue = '/home/silvia/Desktop/Papers/Project_4/Pulsars/matched_fermi_out.fits'

def get_set_flags(flag_value):
    # Returns the list of bit positions (flags) set in the given integer.
    return [bit for bit in range(flag_value.bit_length()) if flag_value & (1 << bit)]

with fits.open(matches_catalogue) as hdul:
    data = Table(hdul[1].data)
    names = [name for name in data.colnames if len(data[name].shape) <= 1]
    df = data[names].to_pandas()
 
disallowed_flag_plus1 = {1, 2, 5, 6}  
df_filtered = df[~((df['Flags'] + 1).isin(disallowed_flag_plus1))]
 
output_table = Table.from_pandas(df_filtered)
output_table.write('filtered_catalogue.fits', overwrite=True)