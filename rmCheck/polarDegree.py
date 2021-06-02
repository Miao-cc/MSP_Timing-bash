#! /usr/bin/env python
import psrchive
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def get_off_pulse_rms(allPulse):
    all_rms = np.std(allPulse)
    offPulseIndex = allPulse/all_rms < 1
    offPulse = allPulse[offPulseIndex]
    off_pulse_rms = np.std(offPulse)
    print "Profile RMS; %f, Off Pulse RMS: %s" %(all_rms, off_pulse_rms)
    return off_pulse_rms, offPulseIndex

def RMScal(allPulse, index):
    all_rms = np.std(allPulse[index])
    return all_rms

try:
    filename = sys.argv[1]
except:
    print "usage: python plotArch.py filename.calibP"

arch = psrchive.Archive_load(filename)

# archive file dedisperse
arch.dedisperse()
# archive file remove baseline
arch.remove_baseline()

# seunch in nbin, nchan ,nsubint
#arch.bscrunch_to_nbin(512)
#arch.fscrunch_to_nchan(128)
arch.convert_state('Stokes')
arch.centre_max_bin()

# get the data 
data = arch.get_data()

# get frequence information
freq_lo = arch.get_centre_frequency() - arch.get_bandwidth()/2.0
freq_hi = arch.get_centre_frequency() + arch.get_bandwidth()/2.0
# get data shape
subintNum, polNum, chanNum, profBin = data.shape
obsLength  = arch.integration_length()
print obsLength

print "Load file: %s" %(filename)
print "Low Frequence: %s, High Frequence: %s" %(freq_lo, freq_hi)
print "Subint Number: %s, Polarization Number: %s, Channel Number: %s, Profile bin: %s" %(subintNum, polNum, chanNum, profBin)

freq_hz = np.linspace(freq_lo*1e6, freq_hi*1e6, num=chanNum, endpoint=True)

# get I, Q, U, V
# subint, frequence, nbin
I = data[0,0,:,:]
Q = data[0,1,:,:]
U = data[0,2,:,:]

I_phase = I.sum(0)

RMS, offPulseIndex = get_off_pulse_rms(I_phase)

onPulseIndex = (I_phase/RMS > 3)
peakIndex = I_phase.argmax()

for num in range(peakIndex,0, -1):
    #print num, onPulseIndex[num]
    if onPulseIndex[num] == False:
        onPulseIndex[0:num] = False
        break

for num in range(peakIndex,len(onPulseIndex)):
    #print num, onPulseIndex[num]
    if onPulseIndex[num] == False:
        onPulseIndex[num:] = False
        break

phase = np.linspace(0,360,profBin)


downSample_I = I[:, onPulseIndex].sum(1)
downSample_Q = Q[:, onPulseIndex].sum(1)
downSample_U = U[:, onPulseIndex].sum(1)

#print downSample_I_err.shape, downSample_Q_err.shape, downSample_U_err.shape


I_RMS = np.array([RMScal(I[i,:], offPulseIndex) for i in range(chanNum)])
Q_RMS = np.array([RMScal(Q[i,:], offPulseIndex) for i in range(chanNum)])
U_RMS = np.array([RMScal(U[i,:], offPulseIndex) for i in range(chanNum)])

noneZeroIndex = (downSample_I != 0)
outArray = np.zeros((len(downSample_I[noneZeroIndex]), 7))
outArray[:,0] = freq_hz[noneZeroIndex]
outArray[:,1] = downSample_I[noneZeroIndex]
outArray[:,2] = downSample_Q[noneZeroIndex]
outArray[:,3] = downSample_U[noneZeroIndex]
outArray[:,4] = np.sqrt(I_RMS[noneZeroIndex])
outArray[:,5] = np.sqrt(Q_RMS[noneZeroIndex])
outArray[:,6] = np.sqrt(U_RMS[noneZeroIndex])

#outArray[:,4] = np.sqrt(RMS / downSample_I[noneZeroIndex].sum(0)**2)
#outArray[:,5] = np.sqrt(RMS / downSample_I[noneZeroIndex].sum(0)**2)
#outArray[:,6] = np.sqrt(RMS / downSample_I[noneZeroIndex].sum(0)**2)

plotPic = False
#plotPic = True 
if plotPic:
    plt.scatter(phase, I_phase)
    plt.scatter(phase[onPulseIndex], I_phase[onPulseIndex], color='r')
    plt.show()

filename = filename + ".txt"
np.savetxt(filename,outArray,fmt='%s',newline='\n')
