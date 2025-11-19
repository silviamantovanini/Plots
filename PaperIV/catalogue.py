#!/usr/bin/env python3
#Creating catalogue

from astropy.io import fits
from astropy.table import Table

infile = "/home/silvia/Desktop/Papers/Project_4/Pulsars/CROSS-MATCHED/compact_multiple.fits"
outfile = "/home/silvia/Desktop/Papers/Paper_IV/testo/matched_pulsars.fits"


table = Table.read(infile, hdu=1)

selected_columns = ["PSRJ", "RAJ", "DECJ", "err_RA_1", "err_DEC_1", "S1400", "src_name", "ref_ra", "err_ra_2", "ref_dec", "err_dec_2", "peak_flux_69_0000", "err_peak_flux_69_0000", "peak_flux_69_0001", "err_peak_flux_69_0001", "peak_flux_69_0002", "err_peak_flux_69_0002", "peak_flux_69_0003", "err_peak_flux_69_0003", "peak_flux_93_0000", "err_peak_flux_93_0000", "peak_flux_93_0001", "err_peak_flux_93_0001", "peak_flux_93_0002", "err_peak_flux_93_0002", "peak_flux_93_0003", "err_peak_flux_93_0003", "peak_flux_121_0000", "err_peak_flux_121_0000", "peak_flux_121_0001", "err_peak_flux_121_0001", "peak_flux_121_0002", "err_peak_flux_121_0002", "peak_flux_121_0003", "err_peak_flux_121_0003", "peak_flux_145_0000", "err_peak_flux_145_0000", "peak_flux_145_0001", "err_peak_flux_145_0001", "peak_flux_145_0002", "err_peak_flux_145_0002", "peak_flux_145_0003", "err_peak_flux_145_0003", "peak_flux_169_0000", "err_peak_flux_169_0000", "peak_flux_169_0001", "err_peak_flux_169_0001", "peak_flux_169_0002", "err_peak_flux_169_0002", "peak_flux_169_0003", "err_peak_flux_169_0003"]

new_table = table[selected_columns].copy()

# Define renaming dictionary
rename_map = {
    "peak_flux_69_0000": "peak_flux_076MHz",
    "err_peak_flux_69_0000": "err_peak_flux_076MHz",
    "peak_flux_69_0001": "peak_flux_084MHz",
    "err_peak_flux_69_0001": "err_peak_flux_084MHz",
    "peak_flux_69_0002": "peak_flux_092MHz",
    "err_peak_flux_69_0002": "err_peak_flux_092MHz",
    "peak_flux_69_0003": "peak_flux_099MHz",
    "err_peak_flux_69_0003": "err_peak_flux_099MHz",
    "peak_flux_93_0000": "peak_flux_107MHz",
    "err_peak_flux_93_0000": "err_peak_flux_107MHz",
    "peak_flux_93_0001": "peak_flux_115MHz",
    "err_peak_flux_93_0001": "err_peak_flux_115MHz",
    "peak_flux_93_0002": "peak_flux_122MHz",
    "err_peak_flux_93_0002": "err_peak_flux_122MHz",
    "peak_flux_93_0003": "peak_flux_130MHz",
    "err_peak_flux_93_0003": "err_peak_flux_130MHz",
    "peak_flux_121_0000": "peak_flux_143MHz",
    "err_peak_flux_121_0000": "err_peak_flux_143MHz",
    "peak_flux_121_0001": "peak_flux_151MHz",
    "err_peak_flux_121_0001": "err_peak_flux_151MHz",
    "peak_flux_121_0002": "peak_flux_158MHz",
    "err_peak_flux_121_0002": "err_peak_flux_158MHz",
    "peak_flux_121_0003": "peak_flux_166MHz",
    "err_peak_flux_121_0003": "err_peak_flux_166MHz",
    "peak_flux_145_0000": "peak_flux_174MHz",
    "err_peak_flux_145_0000": "err_peak_flux_174MHz",
    "peak_flux_145_0001": "peak_flux_181MHz",
    "err_peak_flux_145_0001": "err_peak_flux_181MHz",
    "peak_flux_145_0002": "peak_flux_189MHz",
    "err_peak_flux_145_0002": "err_peak_flux_189MHz",
    "peak_flux_145_0003": "peak_flux_197MHz",
    "err_peak_flux_145_0003": "err_peak_flux_197MHz",
    "peak_flux_169_0000": "peak_flux_204MHz",
    "err_peak_flux_169_0000": "err_peak_flux_204MHz",
    "peak_flux_169_0001": "peak_flux_212MHz",
    "err_peak_flux_169_0001": "err_peak_flux_212MHz",
    "peak_flux_169_0002": "peak_flux_220MHz",
    "err_peak_flux_169_0002": "err_peak_flux_220MHz",
    "peak_flux_169_0003": "peak_flux_227MHz",
    "err_peak_flux_169_0003": "err_peak_flux_227MHz",
}

for old, new in rename_map.items():
    if old in new_table.colnames:
        new_table.rename_column(old, new)
        
new_table.write(outfile, overwrite=True)