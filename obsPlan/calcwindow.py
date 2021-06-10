import numpy as np
import sys
import ephem
import datetime
import argparse


def calcwindow(date, ra, dec, zamax=38.0):

    ZAangle = 90.0 - zamax
    Dawodang=ephem.Observer()
    Dawodang.lat='25.652939'
    Dawodang.lon='106.856594'
    Dawodang.horizon= str(ZAangle) # =(90deg-Maximum Za)
    Dawodang.date= date
    src = ephem.FixedBody()
    src._ra = ra
    src._dec = dec

    src.compute()
    trise  = Dawodang.next_rising(src)
    ttrans = Dawodang.next_transit(src)
    tset   = Dawodang.next_setting(src)
    #trise  = trise.strftime('%Y-%m-%dT%H:%M:%S')
    #ttrans = ttrans.strftime('%Y-%m-%dT%H:%M:%S')
    #tset   = tset.strftime('%Y-%m-%dT%H:%M:%S')

    trise  = str(trise).replace(' ','T').replace('/','-')
    ttrans = str(ttrans).replace(' ','T').replace('/','-')
    tset   = str(tset).replace(' ','T').replace('/','-')
    return trise, ttrans, tset
