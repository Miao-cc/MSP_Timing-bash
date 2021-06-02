#!/bin/usr/evn python
import numpy as np
import matplotlib.pyplot as plt 
import os,sys
from astropy.io import fits

datapath=sys.argv[1]
print('='*80)
print(datapath)
print('='*80)
hdu=fits.open(datapath)
data=hdu[1].data['DATA']

if len(data.shape)>3:
    print ('you checking is PRS backend FITS!')
    print('='*80)
    tsamp=hdu[1].header['TBIN']
    NSUB=hdu['SUBINT'].header['NAXIS2']
    NBIN=hdu['SUBINT'].header['NSBLK']
    timescale=NSUB*NBIN
    NCHAN=hdu['SUBINT'].header['NCHAN']
    cutnchan=np.int32(NCHAN*0.1)
    POLA=data[:,:,0,:]
    POLB=data[:,:,1,:]
    AA=POLA.reshape(timescale,NCHAN)
    BB=POLB.reshape(timescale,NCHAN)	
    #cutAA=AA[:,cutnchan:-cutnchan]
    #cutBB=BB[:,cutnchan:-cutnchan]
    cutAA=AA[:,:]
    cutBB=BB[:,:]
    obstime=timescale*tsamp
else:
    print ('you checking is spectrum backend FITS!')
    print('='*80)
    NSPEC=len(hdu['SINGLE DISH'].data['OBSNUM'])
    EXPOSURE=hdu['SINGLE DISH'].data['EXPOSURE'][0]
    CHAN_BW=hdu['SINGLE DISH'].data['CHAN_BW'][0]
    NCHAN=hdu['SINGLE DISH'].data['NCHAN'][0]
    cutNCHAN=np.int32(NCHAN*0.1)
    AA=data[:,:,0]
    BB=data[:,:,1]
    if (sys.argv[1][-11])=='W':
        cutAA=AA[:,cutNCHAN:-cutNCHAN]
        cutBB=BB[:,cutNCHAN:-cutNCHAN]
    	
    else:
        cutAA=AA
        cutBB=BB
    obstime=NSPEC*EXPOSURE

plt.clf()
fig = plt.figure(figsize=(19,12))
font1={'weight': 'normal','size': 30}
font2={'weight': 'normal','size': 25}
time=np.linspace(0,obstime,(len(cutAA.mean(1))))
plt.subplot(211)
print(cutAA.shape)
plt.plot(time,cutAA.mean(1),'k')
plt.title('POLA Noise_diode (output power spectrum)',font1)
plt.ylabel('Intensity [Arbitrary Units]',font2)
plt.xticks([])
plt.tick_params(labelsize=18)
plt.subplot(212)
plt.plot(time,cutBB.mean(1),'k')
plt.title('POLB Noise_diode (output power spectrum)',font1)
plt.xlabel('duration time(seconds)',font1)
plt.ylabel('Intensity [Arbitrary Units]',font2)
plt.subplots_adjust(wspace=0.0, hspace=0.25)
plt.tick_params(labelsize=18)
#plt.show()
figname = datapath[:-5] + 'checkTime.png'
plt.savefig(figname, dpi=300)
