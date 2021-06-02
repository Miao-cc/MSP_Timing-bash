import math 
import time
import sys, os
#import psrchive
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from astropy import constants as const

def loadData(filename):
    """
    load the data
    """
    # set the global parmaters
    setParams()
    filename = sys.argv[1]

    # load data from ASCII text
    with open(filename, 'r') as f:
        first_line = f.readline()
    #File: puppi_59055_19C64_fold.zap.pF Src: J2016+0758 Nsub: 79 Nch: 1 Npol: 1 Nbin: 128 RMS: 0.000468897
    # MJD(mid): 59289.108033144711 Tsub: 3540.000 Freq: 1255.467 BW: 500
    # get the file info
    archiveFileName = first_line.split(' ')[1]
    first_line = first_line.replace(':','').split(' ')
    print(first_line)
    for num, item in enumerate(first_line):
        if item == "Nsub":
            subintNum = int(first_line[num+1])
        if item == "Nch":
            chanNum = int(first_line[num+1]) 
        if item == "Npol":
            polNum = int(first_line[num+1])
        if item == "Nbin":
            profBin = int(first_line[num+1])

    print("No found: %s, set the freqence low: 1000 MHz, freqence high: 1500 MHz" %(archiveFileName))
    # get frequence information

    # get frequence in hz and wavelength in meter


    data_pdv = np.loadtxt(filename, comments='#')
    data_pdv = data_pdv[:,3:7]
    print("data_pdv shape: ", data_pdv.shape)
    data = np.zeros((subintNum, polNum, chanNum, profBin))
    data[:,0,:,:] = data_pdv[:,0].reshape(subintNum, chanNum, profBin)
    data[:,1,:,:] = data_pdv[:,1].reshape(subintNum, chanNum, profBin)
    data[:,2,:,:] = data_pdv[:,2].reshape(subintNum, chanNum, profBin)
    data[:,3,:,:] = data_pdv[:,3].reshape(subintNum, chanNum, profBin)




    pulsarData['subintNum'] = subintNum
    pulsarData['polNum'] = polNum
    pulsarData['chanNum'] = chanNum
    pulsarData['profBin'] = profBin
    pulsarData['data'] = data
    #print data.shape
    #plt.plot(data[:,0,:,:].sum(0).sum(0))
    #plt.show()


    #return pulsarData

def plot_profile(data, showFile=True):
    matplotlib.rcParams.update({'font.size': 16,})
    start_subint = pulsarData['start_subint']
    end_endint = pulsarData['end_endint']
    profBin = pulsarData['profBin']

    I_prof = data[start_subint:end_endint,0,:,:].sum(0).sum(0)
    Q_prof = data[start_subint:end_endint,1,:,:].sum(0).sum(0)
    U_prof = data[start_subint:end_endint,2,:,:].sum(0).sum(0)
    V_prof = data[start_subint:end_endint,3,:,:].sum(0).sum(0)
    print("Profile shape: \nI %10s \nQ %10s \nU %10s \nV %10s" %(I_prof.shape,  Q_prof.shape, U_prof.shape, V_prof.shape))
    Liner_prof = np.sqrt(Q_prof**2 + U_prof**2)


    #PA = 0.5*np.arctan2(U_prof, Q_prof) / np.pi * 180 % 180
    PA = 0.5*np.arctan2(U_prof, Q_prof) * 180 / np.pi
    off_pulse_rms = get_off_pulse_rms(I_prof)
    index, PAsigma = get_Ltrue(Q_prof, U_prof, off_pulse_rms)/off_pulse_rms
    print(index)
    PAsigma = PAsigma[index]
    PAerr = np.array([PA[i]*(1 - stats.norm.cdf(PAsigma[i]*1/2.)) for i in range(profBin)])
    PAindex = PAsigma > 1
    #PAindex = PAerr < 5
    
    #fig = plt.figure(figsize=(9, 16))
    fig = plt.figure()
    gs = gridspec.GridSpec(nrows=4, ncols=1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1:, 0])
    ax1.set_xlim(0,profBin)
    ax2.set_xlim(0,profBin)
    ax1.set_ylim(-100, 100)
    ax1.set_yticks(np.arange(-90, 90, 30))
    ax1.set_ylabel('P.A. (deg.)')
    ax2.set_xlabel("Phase bin")
    ax2.set_ylabel("Profile Value")
    
    ax1.errorbar(np.arange(profBin)[PAindex], PA[PAindex], yerr=PAerr[PAindex], fmt="bo:",color="blue",ecolor='grey', elinewidth=2,capsize=4 )
    
    ax2.plot(I_prof, 'k')
    ax2.plot(Liner_prof, 'r')
    ax2.plot(V_prof, 'b')
    if showFile is True:
        plt.show()

def setParams():
    global pulsarData
    pulsarData = {}
    ################
    # physical constant
    pulsarData['c'] = const.c.to('m/s').value
    pulsarData['start_subint'] = None
    pulsarData['end_endint'] = None
    pulsarData['startPhase'] = None
    pulsarData['endPhase'] = None

def main(filename):
    """Each function below is independent.
    """

    loadData(filename)

    # get the zapped channel list

    # plot the profile
    data = pulsarData['data']

    print(data.shape)
    plt.imshow(data[:,0,0,:], origin='lower', aspect='auto', cmap='binary')
    plt.show()
    plt.plot(data[0,0,:,:].sum(0))
    plt.plot(np.sqrt(data[0,1,:,:]**2 + data[0,2,:,:]**2).sum(0))
    plt.plot(data[0,3,:,:].sum(0))
    plt.show()

if __name__ == "__main__":
    filename = sys.argv[1]
    main(filename)
