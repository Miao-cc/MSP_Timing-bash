#! /usr/bin/env python
import psrchive
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def get_Ltrue(Q, U, off_pulse_rms):
    L_meas = np.sqrt(U**2 + Q**2)
    L_true = np.zeros(len(L_meas))
    for i in range(len(L_meas)):
        if L_meas[i] / off_pulse_rms >= 1.57:
            L_true[i] = np.sqrt((L_meas[i]/off_pulse_rms)**2 -off_pulse_rms)
            #L_true[i] = np.sqrt((L_meas[i]/off_pulse_rms)**2)
    return L_true

def get_PA_sigma(L_true, off_pulse_rms):
    PA_sigma = np.zeros(len(L_true))
    for i in range(len(L_true)):
        if L_true[i] != 0:
            PA_sigma[i] = 28.65 * off_pulse_rms / L_true[i]
    return PA_sigma

def get_off_pulse_rms(allPulse):
    all_rms = np.std(allPulse)
    offPulse = []
    for i in range(len(allPulse)):
        if allPulse[i]/all_rms < 1:
            offPulse.append(allPulse[i])
    off_pulse_rms = np.std(offPulse)
    print "Profile RMS; %f, Off Pulse RMS: %s" %(all_rms, off_pulse_rms)
    return off_pulse_rms


def getRatio(RM):
    PA = 2.*RM*wave*wave
    #print PA
    Q_ = np.cos(PA)*Q + np.sin(PA)*U
    U_ = -np.sin(PA)*Q + np.cos(PA)*U
    return Q_, U_

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
arch.bscrunch_to_nbin(512)
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

# get I, Q, U, V
# subint, frequence, nbin
I = data[:,0,:]
Q = data[:,1,:]
U = data[:,2,:]
V = data[:,3,:] 


I_phase = data[:,0,:].sum(1)
Q_phase = data[:,1,:].sum(1)
U_phase = data[:,2,:].sum(1)
V_phase = data[:,3,:].sum(1)

#I_phase = np.freeze(data[:,0,:])
#Q_phase = np.freeze(data[:,1,:])
#U_phase = np.freeze(data[:,2,:])
#V_phase = np.freeze(data[:,3,:])

print I_phase.shape, I_phase.max(), I_phase.min()

L_phase = np.sqrt(Q_phase**2 + Q_phase**2)


fig, ax = plt.subplots(1, 3, figsize=(16, 9))
ax = ax.flatten()
ax1 = ax[0]
im1 = ax1.imshow(I_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
ax1.set_title("Total Intensity", fontsize=20)
ax1.set_xlabel("Pulse Phase", fontsize=20)
ax1.set_ylabel("Time (seconds)", fontsize=20)
ax1.set_yticks(np.linspace(0, subintNum, num=10, endpoint=True))
ax1.set_yticklabels([str(int(obsLength*i/10.)) for i in range(10)], fontsize = 20)
ax1.set_xticks(np.linspace(0.1, 1, num=9, endpoint=False))
ax1.set_xticklabels([str(i/10.) for i in range(1, 10)], fontsize = 13)

ax2 = ax[1]
#ax2 = plt.subplot(1,3,2, sharey=ax1)
im2 = ax2.imshow(L_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
#im2 = ax2.imshow(L_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
ax2.set_xlabel("Pulse Phase", fontsize=20)
ax2.set_title("Linear Polarization", fontsize=20)
#ax1.set_xlabel("Linear Polarization")
ax2.set_yticks([])
ax2.set_xticks(np.linspace(0.1, 1, num=9, endpoint=False))
ax2.set_xticklabels([str(i/10.) for i in range(1, 10)], fontsize = 13)


I_max_index = I_phase.sum(1).argmax()
ax3 = ax[2]
#ax3 = plt.subplot(1,3,3, sharey=ax2)
#im3 = ax3.imshow(V_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
#im3 = ax3.imshow(V_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
im3 = ax3.imshow(L_phase/I_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum),vmin=0,vmax=1)
#im3 = ax3.imshow(V_phase/I_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum),vmin=-1,vmax=1)
print (L_phase[:, I_max_index:I_max_index+2]/I_phase[:, I_max_index:I_max_index+2]).max(), (L_phase[:, I_max_index:I_max_index+2]/I_phase[:, I_max_index:I_max_index+2]).min()
print (L_phase[:, I_max_index:I_max_index+2]/I_phase[:, I_max_index:I_max_index+2]).mean(), (L_phase[:, I_max_index:I_max_index+2]/I_phase[:, I_max_index:I_max_index+2]).std()
ax3.set_xlabel("Pulse Phase",fontsize=20)
ax3.set_title("Circular Polarization",fontsize=20)
ax3.set_yticks([])
ax3.set_xticks(np.linspace(0.1, 1, num=9, endpoint=False))
ax3.set_xticklabels([str(i/10.) for i in range(1, 10)], fontsize = 13)
plt.subplots_adjust(wspace=0)

ax1.set_xlim(0.4,0.8)
ax2.set_xlim(0.4,0.8)
ax3.set_xlim(0.4,0.8)

#fig.colorbar(im3, ax=[ax1, ax2, ax3])
fig.colorbar(im3, ax=[ax1, ax2, ax3], fraction=0.03, pad=0.01)
#plt.tight_layout()
#plt.savefig(filename + "_Time2phase.eps", format='eps', dpi=300)
plt.show()
