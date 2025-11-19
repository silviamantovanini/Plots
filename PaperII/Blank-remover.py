#!/home/silvia/PhD/bin/python3

#Remove 'BLANK' keyword from headers.

import sys
from astropy.io import fits
from numpy import sin, radians
import argparse
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

def remove_blank(fitsfile, outname):
    hdu = fits.open(fitsfile)
    
    hdr = hdu[0].header
    del hdr['BLANK']

    if outname is None:
        outname = fitsfile.replace(".fits", "_nb.fits")
    hdu.writeto(outname, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Input/output files")
    group1.add_argument("--infits", dest='infits', type=str, default=None,
                        help="The fits image to add a beam to.")
    group1.add_argument("--outfits", dest='outfits', type=str, default=None,
                        help="The output filename (default = _wpsf).")
    results = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    remove_blank(results.infits, results.outfits)