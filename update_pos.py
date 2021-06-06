#import fitsio
from astropy.io.fits import update
from astropy.io import fits as pyfits
import os, sys

filename = sys.argv[1]
ra = sys.argv[2]
dec = sys.argv[3]

hdulist = pyfits.open(filename, mode='update')
hd0 = hdulist[0]
hd0.header['RA'] = ra
hd0.header['DEC'] = dec
print(hd0.header['RA'])
print(hd0.header['DEC'])

hdulist.close()
