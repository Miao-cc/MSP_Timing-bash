#!/usr/bin/env python
import sys
import numpy as num
from builtins import range
from presto import parfile
from presto import bestprof
import matplotlib.pyplot as plt
from presto import psr_utils as pu
from presto import psr_constants as pc
from calcwindow import *
from astropy.time import Time

#import astropy.units as u
#import pint.models as models
#from astropy.visualization import quantity_support

#quantity_support()

time = num.asarray([])
period = num.asarray([])

def getobsTime(firstMJD, lastMJD, RA, DEC):
    obsTimeList = []
    for mjd in range(int(firstMJD), int(lastMJD)):
        Data = Time(mjd, format='mjd').iso
        riseTime, transitTime, setTime = calcwindow(Data, RA, DEC, 36)
        riseTimeMJD = Time(str(riseTime), format='isot', scale='utc').mjd
        transitTimeMJD = Time(str(transitTime), format='isot', scale='utc').mjd
        setTimeMJD = Time(str(setTime), format='isot', scale='utc').mjd
        obsTimeList.append([riseTimeMJD, transitTimeMJD, setTimeMJD])
    return obsTimeList

def orbeqn(Ppxt, times):
    # P = Ppsr, p = Porb, x = a*sin(i)/s, t = T_o
    phi = pc.TWOPI*(times - Ppxt[3])*86400.0/Ppxt[1]
    return Ppxt[0]*(1.0+pc.TWOPI*Ppxt[2]/Ppxt[1]*num.cos(phi))*1e3

def funct(Ppxt, times, measured):
    return orbeqn(Ppxt, times) - measured

if __name__ == '__main__':
    if len(sys.argv)==1:
        print("\nusage: fit_circular_orbit.py parfiles MJD_list")
        print("Or")
        print("\nusage: fit_circular_orbit.py parfiles MJD_list Period_list")

        exit(0)

    
    if len(sys.argv[1:]) == 2:
        print("Input file number: 2")
        Ppsr_bestprof = None
    elif len(sys.argv[1:]) == 3:
        print("Input file number: 3")
        Ppsr_bestprof = sys.argv[3]
    else:
        print("Too many input files")

    parFile = sys.argv[1]
    MJDList = sys.argv[2]


    MJDs = num.loadtxt("MJDList.txt")
    startMJD = min(MJDs)
    endMJD = max(MJDs)

    x = parfile.psr_par(parFile)

    Ppsr = 1./x.F0
    Porb = x.PB *86400.0
    Torb = x.T0
    Xorb = x.A1
    P_RAJ = x.RAJ 
    P_DECJ = x.DECJ
    
    print("PSR period: %s ms orbit period: %s d  A1: %s  Tepoch: %s" %(Ppsr, Porb/86400.0, Xorb, Torb))
    print("MJD range: %s to %s" %(startMJD, endMJD))
    newts = num.arange(startMJD-0.1, endMJD+0.2,0.01)
    time = orbeqn([Ppsr, Porb, Xorb, Torb], newts)

    plt.figure(figsize=(8,6))
    plt.plot(newts, time,'r')
    newts = MJDs
    times = orbeqn([Ppsr, Porb, Xorb, Torb], newts)
    plt.scatter(newts, times,c='b')

    obsWindoList = getobsTime(startMJD, endMJD, P_RAJ, P_DECJ)
    for obsRange in obsWindoList:
        plt.axvline(x=obsRange[0], c="gray", ls="-", lw=1)
        plt.axvline(x=obsRange[1], c="gray", ls="--", lw=1)
        plt.axvline(x=obsRange[2], c="gray", ls="-", lw=1)

    if Ppsr_bestprof is not None:
        print("Open file: ", Ppsr_bestprof)
        data = num.loadtxt(Ppsr_bestprof)
        print(data)
        plt.scatter(data[:,0], data[:,1],c='k')

    plt.ylabel("Period (ms)")
    plt.xlabel("MJD")
    plt.xlim(startMJD-0.15, endMJD+0.25)
    plt.savefig("result.eps", format='eps')
    plt.show()
