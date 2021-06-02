import fitsio
import os, sys

filename = sys.argv[1]
ra = sys.argv[2]
dec = sys.argv[3]

hdulist = fitsio.FITS(filename, 'rw')
hdulist[0].write_key('RA',  ra)
hdulist[0].write_key('DEC', dec)

hdulist.close()
