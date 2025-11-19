#!/home/silvia/PhD/bin/python3

from astropy.io import fits
from astropy.table import Table

input_fits = 'merged_catalogue_lon_SEDs.fits'
output_fits = 'merged_catalogue_sample.fits'

match_string = ['_69_0000', '_69_0001', '_69_0002', '_69_0003', '_69_MFS', '93_0000', '93_0001', '93_0002', '93_0003', '93_MFS', '121_0000', '121_0001', '121_0002', '121_0003', '121_MFS', '145_0000', '145_0001', '145_0002', '145_0003', '145_MFS', '169_0000', '169_0001', '169_0002', '169_0003', '169_MFS']
replacement = ['_072-080MHz', '_080-088MHz', '_088-095MHz', '_095-103MHz', '_072-103MHz', '103-111MHz', '111-118MHz', '118-126MHz', '126-134MHz', '103-134MHz', '139-147MHz', '147-154MHz', '154-162MHz', '162-170MHz', '139-170MHz', '170-177MHz', '177-185MHz', '185-193MHz', '193-200MHz', '170-200MHz', '200-208MHz', '208-216MHz', '216-223MHz', '223-231MHz', '200-231MHz']

matching = ['src_name', 'background', 'local_rms', 'ref_ra', 'ref_dec', 'err_ra', 'err_dec', 'peak_flux', 'err_peak_flux', 'int_flux', 'err_int_flux', 'a', 'err_a', 'b', 'err_b', 'pa', 'err_pa', 'residual_mean', 'residual_std', 'psf_a', 'psf_b', 'psf_pa', ]
replace = ['Name', 'background_wide', 'local_rms_wide', 'RAJ2000', 'DEJ2000', 'err_RAJ2000', 'err_DEJ2000', 'peak_flux_wide', 'err_peak_flux_wide', 'int_flux_wide', 'err_int_flux_wide', 'a_wide', 'err_a_wide', 'b_wide', 'err_b_wide', 'pa_wide', 'err_pa_wide', 'residual_mean_wide', 'residual_std_wide', 'psf_a_wide', 'psf_b_wide', 'psf_pa_wide', ]

table = Table.read(input_fits)
table = table[:1000]

replace_map = dict(zip(match_string, replacement))

new_colnames = []
for name in table.colnames:
    new_name = name
    for match_str, repl in replace_map.items():
        if match_str in new_name:
            new_name = new_name.replace(match_str, repl)
    if new_name != name:
        table.rename_column(name, new_name)

if 'src_name' in table.colnames:
    table['src_name'] = [name.replace('GP ', 'GLEAM-X ') for name in table['src_name']]
else:
    print("Warning: 'src_name' column not found.")
    
# Rename columns using exact match
exact_replace_map = dict(zip(matching, replace))
for old_name, new_name in exact_replace_map.items():
    if old_name in table.colnames:
        table.rename_column(old_name, new_name)

table.write(output_fits, overwrite=True)
print(f"Saved first 1000 rows to '{output_fits}' with updated column names.")