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
bandList = np.zeros(len(fileList))
rmList = np.zeros(len(fileList))

for num, RM in enumerate(rms):
    filename = fileBase+str(RM)+'.npz'
    npzfile = np.load(filename)
    data = npzfile['arr_0']
    rmValue, snr = npzfile['header']
    dataImage.append(data)
    bandList[num] = data.sum()
    rmList[num] = rmValue
dataImage = np.array(dataImage)
print dataImage.shape

maxLinPol = dataImage[bandList.argmax(),:]
print rms[bandList.argmax()]
print bandList.argmax()
print bandList.max()
print maxLinPol

rows=cols = 4
extent = [0, len(maxLinPol), rmList[0], rmList[-1]]
figsize = (10, 8)
gs = gridspec.GridSpec(rows, cols)
fig = plt.figure(num=1, figsize=figsize)
ax1 = fig.add_subplot(gs[1:,0])
ax2 = fig.add_subplot(gs[0,1:])
ax3 = fig.add_subplot(gs[1:,1:])
ax1.plot(bandList, rmList, "+-k")
ax1.set_ylim(rmList[0], rmList[-1])
ax1.invert_xaxis()
ax2.plot(maxLinPol,"k")
ax2.set_xlim(0, len(maxLinPol))
ax3.imshow(dataImage, origin='lower', aspect='auto', cmap='hot', extent=extent)
ax1.axhline(y=rms[bandList.argmax()], c = 'r')

#ax1.set_axis_off()
#ax2.set_axis_off()
#ax3.set_axis_off()

ax1.get_shared_y_axes().join(ax1, ax3)
ax3.get_shared_x_axes().join(ax2, ax3)
#ax1.set_yticks([])
ax2.set_xticks([])
ax3.set_yticks([])

ax1.set_xlabel("band Sum", font1)
ax1.set_ylabel(r"Rotation Measure [rad $m^{-2}$]", font1)
#ax2.set_ylabel("linear polarization\n profile", font1)
ax2.set_title("Linear polarization profile", font1)
ax3.set_xlabel("Phase [bin]", font1)


ax4 = fig.add_subplot(gs[0,0])
maxBandindex = bandList.argmax()
ax4.text(0.0,0.4,"Best RM: %s" %(int(rms[maxBandindex])), font1)
ax4.set_axis_off()


fig.subplots_adjust(hspace=0, wspace=0)
outPNGname = matchName.split('/')[-1] + 'result.png'
#plt.savefig(outPNGname, dpi=300)

plt.show()

