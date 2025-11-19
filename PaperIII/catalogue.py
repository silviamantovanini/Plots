#!/usr/bin/env python3
#Creating catalogue

from astropy.io import fits
from astropy.table import Table
import os

infile = "/home/silvia/Desktop/Papers/Paper_III/merged_catalogue_lon_SEDs.fits"
outfile = "/home/silvia/Desktop/Papers/Paper_III/Paper/testo/catalogue_gleam-x-gp.fits"
final = "/home/silvia/Desktop/Papers/Paper_III/Paper/testo/GLEAM-X-GP_catalogue.fits"

table = Table.read(infile, hdu=1)

for col in table.colnames:
    if "_69_0000" in col:
        new_name = col.replace("_69_0000", "_072-080MHz")
        table.rename_column(col, new_name)
    if "_69_0001" in col:
        new_name = col.replace("_69_0001", "_080-088MHz")
        table.rename_column(col, new_name)
    if "_69_0002" in col:
        new_name = col.replace("_69_0002", "_088-095MHz")
        table.rename_column(col, new_name)
    if "_69_0003" in col:
        new_name = col.replace("_69_0003", "_095-103MHz")
        table.rename_column(col, new_name)
        
    if "_93_0000" in col:
        new_name = col.replace("_93_0000", "_103-111MHz")
        table.rename_column(col, new_name)
    if "_93_0001" in col:
        new_name = col.replace("_93_0001", "_111-118MHz")
        table.rename_column(col, new_name)
    if "_93_0002" in col:
        new_name = col.replace("_93_0002", "_118-126MHz")
        table.rename_column(col, new_name)
    if "_93_0003" in col:
        new_name = col.replace("_93_0003", "_126-134MHz")
        table.rename_column(col, new_name)
        
    if "_121_0000" in col:
        new_name = col.replace("_121_0000", "_139-147MHz")
        table.rename_column(col, new_name)
    if "_121_0001" in col:
        new_name = col.replace("_121_0001", "_147-154MHz")
        table.rename_column(col, new_name)
    if "_121_0002" in col:
        new_name = col.replace("_121_0002", "_154-162MHz")
        table.rename_column(col, new_name)
    if "_121_0003" in col:
        new_name = col.replace("_121_0003", "_162-170MHz")
        table.rename_column(col, new_name)
        
    if "_145_0000" in col:
        new_name = col.replace("_145_0000", "_170-177MHz")
        table.rename_column(col, new_name)
    if "_145_0001" in col:
        new_name = col.replace("_145_0001", "_177-185MHz")
        table.rename_column(col, new_name)
    if "_145_0002" in col:
        new_name = col.replace("_145_0002", "_185-193MHz")
        table.rename_column(col, new_name)
    if "_145_0003" in col:
        new_name = col.replace("_145_0003", "_193-200MHz")
        table.rename_column(col, new_name)
        
    if "_169_0000" in col:
        new_name = col.replace("_169_0000", "_200-208MHz")
        table.rename_column(col, new_name)
    if "_169_0001" in col:
        new_name = col.replace("_169_0001", "_208-216MHz")
        table.rename_column(col, new_name)
    if "_169_0002" in col:
        new_name = col.replace("_169_0002", "_216-223MHz")
        table.rename_column(col, new_name)
    if "_169_0003" in col:
        new_name = col.replace("_169_0003", "_223-231MHz")
        table.rename_column(col, new_name)
    
    if "_69_MFS" in col:
        new_name = col.replace("_69_MFS", "_072-103MHz")
        table.rename_column(col, new_name)
    if "_93_MFS" in col:
        new_name = col.replace("_93_MFS", "_103-134MHz")
        table.rename_column(col, new_name)
    if "_121_MFS" in col:
        new_name = col.replace("_121_MFS", "_139-170MHz")
        table.rename_column(col, new_name)
    if "_145_MFS" in col:
        new_name = col.replace("_145_MFS", "_170-200MHz")
        table.rename_column(col, new_name)
    if "_169_MFS" in col:
        new_name = col.replace("_169_MFS", "_200-231MHz")
        table.rename_column(col, new_name)

table.rename_column("ref_ra", "RAJ2000")
table.rename_column("ref_dec", "DEJ2000")
table.rename_column("err_ra", "err_RAJ2000")
table.rename_column("err_dec", "err_DEJ2000")

table.write(outfile, overwrite=True)
table = Table.read(outfile, hdu=1)
      
freqs = """072-080MHz
    080-088MHz
    088-095MHz
    095-103MHz
    103-111MHz
    111-118MHz
    118-126MHz
    126-134MHz
    139-147MHz
    147-154MHz
    154-162MHz
    162-170MHz
    170-177MHz
    177-185MHz
    185-193MHz
    193-200MHz
    200-208MHz
    208-216MHz
    216-223MHz
    223-231MHz
    072-103MHz
    103-134MHz
    139-170MHz
    170-200MHz
    200-231MHz
    """.split()
    
metas = {}
metas['src_name']     = {'ucd':'meta.id;src',                                     'units':''       ,'description':'IAU Name'}
metas['background']   = {'ucd':'stat.median;phot.flux.density',                   'units':'Jy/beam','description':'Background level in wide image'}
metas['local_rms']    = {'ucd':'stat.variance;phot.flux.density',                 'units':'Jy/beam','description':'Local noise leve in wide image'}
metas['ra_str']       = {'ucd':'pos.eq.ra;meta.main',                             'units':'h:m:s',  'description':'Right Ascension J2000'}
metas['dec_str']      = {'ucd':'pos.eq.dec;meta.main',                            'units':'d:m:s',  'description':"Declination J2000"}
metas['RAJ2000']      = {'ucd':'pos.eq.ra',                                       'units':'deg',    'description':'Right Ascension J2000'}
metas['DEJ2000']      = {'ucd':'pos.eq.dec',                                      'units':'deg',    'description':"Declination J2000"}
metas['err_RAJ2000']  = {'ucd':'stat.error;pos.eq.ra',                            'units':'deg',    'description':"Uncertainty in Right Ascension"}
metas['err_DEJ2000']  = {'ucd':'stat.error;pos.eq.dec',                           'units':'deg',    'description':"Uncertainty in Declination"}
metas['peak_flux']    = {'ucd':'phot.flux.density;em.radio.170-231MHz',           'units':'Jy/beam','description':'Peak flux in wide image'}
metas['err_peak_flux']= {'ucd':'stat.error;phot.flux.density;em.radio.170-231MHz','units':'Jy/beam','description':'Uncertainty in fit for peak flux in wide image'}
metas['int_flux']     = {'ucd':'phot.flux;em.radio.170-231MHz',                   'units':'Jy',     'description':'Integrated flux in wide image'}
metas['err_int_flux'] = {'ucd':'stat.error;phot.flux;em.radio.170-231MHz',        'units':'Jy',     'description':'Uncertainty in fit for integrated flux in wide image'}
metas['a']            = {'ucd':'phys.angSize.smajAxis;meta.modelled',             'units':'arcsec', 'description':'Fitted semi-major axis in wide image'}
metas['err_a']        = {'ucd':'stat.error;phys.angSize.smajAxis',                'units':'arcsec', 'description':'Uncertainty in fitted semi-major axis in wide image'}
metas['b']            = {'ucd':'phys.angSize.sminAxis;meta.modelled',             'units':'arcsec', 'description':'Fitted semi-minor axis in wide image'}
metas['err_b']        = {'ucd':'stat.error;phys.angSize.sminAxis',                'units':'arcsec', 'description':'Uncertainty in fitted semi-minor axis in wide image'}
metas['pa']           = {'ucd':'pos.posAng;meta.modelled',                        'units':'deg',    'description':'Fitted position angle in wide image'}
metas['err_pa']       = {'ucd':'stat.error;pos.posAng',                           'units':'deg',    'description':'Uncertainty in fitted position angle in wide image'}
metas['residual_mean']= {'ucd':'stat.fit.residual',                               'units':'Jy/beam','description':'Mean value of data-model in wide image'}
metas['residual_std'] = {'ucd':'stat.fit.residual',                               'units':'Jy/beam','description':'Standard deviation of data-model in wide image'}
metas['psf_a']        = {'ucd':'phys.angSize.smajAxis;instr.det.psf',             'units':'arcsec', 'description':'Semi-major axis of the point spread function in wide image'}
metas['psf_b']        = {'ucd':'phys.angSize.sminAxis;instr.det.psf',             'units':'arcsec', 'description':'Semi-minor axis of the point spread function in wide image'}
metas['psf_pa']       = {'ucd':'pos.posAng;instr.det.psf',                        'units':'deg',    'description':'Position angle of the point spread function in wide image'}
metas['S_200']        = {'ucd':'phot.flux;em.radio.200MHz',                       'units':'Jy',     'description':'Power-law fitted flux density at 200 MHz'}
metas['pl_alpha']     = {'ucd':'spec.index;em.radio',                             'units':'',       'description':'Fitted spectral index assuming power-law SED'}
metas['pl_alpha_err'] = {'ucd':'stat.error;spec.index;em.radio',                  'units':'',       'description':'Error on fitted spectral index assuming power-law SED'}
metas['pl_rchi2']     = {'ucd':'stat.fit.goodness',                               'units':'',       'description':'Reduced chi-squared statistic for power-law SED fit'}
metas['cpl_alpha']    = {'ucd':'spec.index;em.radio',                             'units':'',       'description':'Fitted spectral index assuming curved SED'}
metas['cpl_alpha_err']= {'ucd':'stat.error;spec.index;em.radio',                  'units':'',       'description':'Error on fitted spectral index assuming curved SED'}
metas['cpl_rchi2']    = {'ucd':'stat.fit.goodness',                               'units':'',       'description':'Reduced chi-squared statistic for curved SED fit'}

for frange  in freqs:
    suffix = '_{0}'.format(frange)
    metas['background'+suffix]   = {'ucd':'stat.median;src;phot.flux.density',               'units':'Jy/beam','description':'Background level in '+frange+' image'}
    metas['local_rms'+suffix]    = {'ucd':'stat.variance;src;phot.flux.density',             'units':'Jy/beam','description':'Local noise level in '+frange+' image'}
    metas['peak_flux'+suffix]    = {'ucd':'phot.flux.density;em.radio.'+frange,              'units':'Jy/beam','description':'Peak flux in '+frange+' image'}
    metas['err_peak_flux'+suffix]= {'ucd':'stat.error;phot.flux.density;em.radio.'+frange,   'units':'Jy/beam','description':'Uncertainty in fit for peak flux in '+frange+' image'}
    metas['int_flux'+suffix]     = {'ucd':'phot.flux;em.radio.'+frange,                      'units':'Jy',     'description':'Integrated flux in '+frange+' image'}
    metas['err_int_flux'+suffix] = {'ucd':'stat.error;phot.flux;em.radio.'+frange,           'units':'Jy',     'description':'Uncertainty in fit for integrated flux in '+frange+' image'}
    metas['a'+suffix]            = {'ucd':'phys.angSize.smajAxis;src;meta.modelled',         'units':'arcsec', 'description':'Fitted semi-major axis in '+frange+' image'}
    metas['b'+suffix]            = {'ucd':'phys.angSize.sminAxis;src;meta.modelled',         'units':'arcsec', 'description':'Fitted semi-minor axis in '+frange+' image'}
    metas['pa'+suffix]           = {'ucd':'pos.posAng;meta.modelled',                        'units':'degree', 'description':'Fitted position angle in '+frange+' image' }
    metas['residual_mean'+suffix]= {'ucd':'stat.fit.residual',                               'units':'Jy/beam','description':'Mean value of data-model in '+frange+' image'}
    metas['residual_std'+suffix] = {'ucd':'stat.fit.residual',                               'units':'Jy/beam','description':'Standard deviation of data-model in '+frange+' image'}
    metas['psf_a'+suffix]        = {'ucd':'phys.angSize.smajAxis;instr.det.psf',             'units':'arcsec', 'description':'Semi-major axis of the point spread function in '+frange+' image'}
    metas['psf_b'+suffix]        = {'ucd':'phys.angSize.sminAxis;instr.det.psf',             'units':'arcsec', 'description':'Semi-minor axis of the point spread function in '+frange+' image'}
    metas['psf_pa'+suffix]       = {'ucd':'pos.posAng;instr.det.psf',                        'units':'deg',    'description':'Position angle of the point spread function in '+frange+' image'}
       
cmd = f'stilts tpipe in={outfile} out={final} ofmt=fits-plus cmd="'
for k, val in metas.items():
    cmd += f"colmeta -units '{val['units']}' -ucd '{val['ucd']}' -desc '{val['description']}' '{k}'; \\\n"
cmd += '"'
os.system(cmd)