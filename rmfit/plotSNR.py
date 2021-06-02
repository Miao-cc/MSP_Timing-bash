from __future__ import division
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
font1 = {'family' : 'Ubuntu Mono',
'weight' : 'normal',
'size'   : 16,
}

matchName = sys.argv[1]
fileList = glob.glob(matchName+"*.npz")
print matchName
#fileList = glob.glob("19C13/19C13-RM_-500_100_subint_1_2*.npz")
#print fileList
# sort the filename
RMstrs = [filename.split("__RM_")[-1].split(".npz")[0] for filename in fileList]
rms = list(map(float, RMstrs))
rms.sort()
fileBase = fileList[0].split("__RM_")[0] + '__RM_'

dataImage = []
snrList = np.zeros(len(fileList))
rmList = np.zeros(len(fileList))

for num, RM in enumerate(rms):
    filename = fileBase+str(RM)+'.npz'
    npzfile = np.load(filename)
    data = npzfile['arr_0']
    rmValue, snr = npzfile['header']
    dataImage.append(data)
    snrList[num] = snr
    rmList[num] = rmValue
    #print filename, len(data), snrList[num]
dataImage = np.array(dataImage)
snrList[np.isnan(snrList)] = 0
snrList[np.isinf(snrList)] = 0
snrList[snrList < 0.0] = 0
print dataImage.shape
#plt.show()
#plt.plot(rmList, snrList)
#plt.xlabel("RM")
#plt.ylabel("SNR")
#plt.show()

maxLinPol = dataImage[snrList.argmax(),:]
print rms[snrList.argmax()]
print snrList.argmax()
print snrList.max()
print maxLinPol

rows=cols = 4
extent = [0, len(maxLinPol), rmList[0], rmList[-1]]
figsize = (10, 8)
gs = gridspec.GridSpec(rows, cols)
fig = plt.figure(num=1, figsize=figsize)
ax1 = fig.add_subplot(gs[1:,0])
ax2 = fig.add_subplot(gs[0,1:])
ax3 = fig.add_subplot(gs[1:,1:])
ax1.plot(snrList, rmList, "+-k")
ax1.set_ylim(rmList[0], rmList[-1])
ax1.invert_xaxis()
ax2.plot(maxLinPol,"k")
ax2.set_xlim(0, len(maxLinPol))
ax3.imshow(dataImage, origin='lower', aspect='auto', cmap='hot', extent=extent)

#ax1.set_axis_off()
#ax2.set_axis_off()
#ax3.set_axis_off()

ax1.get_shared_y_axes().join(ax1, ax3)
ax3.get_shared_x_axes().join(ax2, ax3)
#ax1.set_yticks([])
ax2.set_xticks([])
ax3.set_yticks([])

ax1.set_xlabel("SNR", font1)
ax1.set_ylabel(r"Rotation Measure [rad $m^{-2}$]", font1)
#ax2.set_ylabel("linear polarization\n profile", font1)
ax2.set_title("Linear polarization profile", font1)
ax3.set_xlabel("Phase [bin]", font1)


ax4 = fig.add_subplot(gs[0,0])
maxSNRindex = snrList.argmax()
ax4.text(0.0,0.4,"Best SNR: %s\nBest RM: %s" %(int(snrList[maxSNRindex]), int(rms[maxSNRindex])), font1)
ax4.set_axis_off()

#ax[-1].set_title('markevery=%s' % str(case))
#ax[-1].plot(x, y, 'o', ls='-', ms=4, markevery=case)

#fig.subplots_adjust(left=0.05,right=0.9, bottom=0.05, top=0.9, wspace=0, hspace=0)


fig.subplots_adjust(hspace=0, wspace=0)
#fig.tight_layout()
outPNGname = matchName.split('/')[-1] + 'result.png'
#plt.savefig(outPNGname, dpi=300)

#plt.savefig("19C13_RMsearch_1.png", dpi=300)
plt.show()

