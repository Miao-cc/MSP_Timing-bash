#! /usr/bin/env python
import psrchive
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

try:
    filename = sys.argv[1]
except:
    print "usage: python plotArch.py filename.calibP"

arch = psrchive.Archive_load(filename)

# archive file dedisperse
#arch.dedisperse()
arch.defaraday()
arch.remove_baseline()
# archive file remove baseline

# seunch in nbin, nchan ,nsubint
#arch.bscrunch_to_nbin(512)
#arch.fscrunch_to_nchan(128)
#arch.tscrunch_to_nsub(32)
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

#data = np.roll(data, 700, axis=3)

# get I, Q, U, V
# subint, frequence, nbin
I = data[:,0,0,:]
Q = data[:,1,0,:]
U = data[:,2,0,:]
V = data[:,3,0,:]

I_phase = data[:,0,0,:].sum(0)
Q_phase = data[:,1,0,:].sum(0)
U_phase = data[:,2,0,:].sum(0)
V_phase = data[:,3,0,:].sum(0)

print I_phase.shape, I_phase.max(), I_phase.min()

L = np.sqrt(Q**2 + Q**2)
L_phase = np.sqrt(Q_phase**2 + Q_phase**2)

plt.plot(I_phase, 'k')
plt.plot(L_phase, 'r')
plt.plot(V_phase, 'b')
plt.show()

plt.plot(I[0,:], 'k')
plt.plot(L[0,:], 'r')
plt.plot(V[0,:], 'b')
plt.show()


fig, ax = plt.subplots(1, 3, figsize=(16, 9))
ax = ax.flatten()
ax1 = ax[0]
im1 = ax1.imshow(I, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
ax1.set_title("Total Intensity", fontsize=20)
ax1.set_xlabel("Pulse Phase", fontsize=20)
ax1.set_ylabel("Time (seconds)", fontsize=20)
ax1.set_yticks(np.linspace(0, subintNum, num=10, endpoint=True))
ax1.set_yticklabels([str(int(obsLength*i/10.)) for i in range(10)], fontsize = 20)
ax1.set_xticks(np.linspace(0.1, 1, num=9, endpoint=False))
ax1.set_xticklabels([str(i/10.) for i in range(1, 10)], fontsize = 13)

ax2 = ax[1]
#ax2 = plt.subplot(1,3,2, sharey=ax1)
im2 = ax2.imshow(L, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
#im2 = ax2.imshow(L_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
ax2.set_xlabel("Pulse Phase", fontsize=20)
ax2.set_title("Linear Polarization", fontsize=20)
#ax1.set_xlabel("Linear Polarization")
ax2.set_yticks([])
ax2.set_xticks(np.linspace(0.1, 1, num=9, endpoint=False))
ax2.set_xticklabels([str(i/10.) for i in range(1, 10)], fontsize = 13)

ax3 = ax[2]
#ax3 = plt.subplot(1,3,3, sharey=ax2)
im3 = ax3.imshow(V, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
#im3 = ax3.imshow(V_phase, origin='lower', aspect='auto', cmap = 'binary', extent=(0, 1, 0, subintNum))
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
