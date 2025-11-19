#!/home/silvia/anaconda3/envs/racs/bin/python3
#Download a RACS cutout from Casda, centered in the coordinates I'm giving in input.

from astropy import units as u, wcs
from astropy.coordinates import SkyCoord
from astroquery.casda import Casda
import sys, numpy as np, os

RaJ2000=sys.argv[1]
DecJ2000=sys.argv[2]

centre = SkyCoord(RaJ2000, DecJ2000, unit=(u.hourangle, u.deg), frame='fk5')

casda = Casda()
casda.login(username='silvia.mantovanini@postgrad.curtin.edu.au')

result = Casda.query_region(centre, radius=150*u.arcmin)
public_data = Casda.filter_out_unreleased(result)

#subset = public_data[((public_data['obs_collection'] == 'The Rapid ASKAP Continuum Survey') & (np.char.startswith(public_data['filename'], 'RACS-DR1_')) & (np.char.endswith(public_data['filename'], 'A.fits')))]
subset = public_data[((public_data['obs_collection'] == 'EMU') & (np.char.startswith(public_data['filename'], 'image.i.EMU_')) & (np.char.endswith(public_data['filename'], '0.restored.conv.fits')))]

url_list = casda.cutout(subset[:1], coordinates=centre, radius=90*u.arcmin)
filelist = casda.download_files(url_list, savedir='/home/silvia/DS9_file/Summary_Plot/Candidates/G321.3-3.9/SpectralIndex/FIts_file/Plot-cutted')

#os.rename(filelist[1], '/home/silvia/Desktop/Papers/Paper_I/G341.4-0.2/SNR_G341.4-0.2_RACS.fits')
os.rename(filelist[1], '/home/silvia/DS9_file/Summary_Plot/Candidates/G321.3-3.9/SpectralIndex/FIts_file/Plot-cutted/SNR_G321.3-3.9_EMU.fits')
os.remove(filelist[0])