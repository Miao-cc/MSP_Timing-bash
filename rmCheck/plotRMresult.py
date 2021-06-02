#! /usr/bin/env python
import sys
import math
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import leastsq
from scipy.optimize import curve_fit


def gaussFunc(x,mean,sd,a):
    norm = a*np.exp(-(x - mean)**2/(2*sd**2))
    return np.array(norm)

def gauss(x,*p):
    return gaussFunc(x,p[0],p[1],p[2])

def gaussFit(rmResult):
    #def gaussFunc(x,mean,sd,a):
        #norm = a*np.exp(-(x - mean)**2/(2*sd**2))
        #return np.array(norm)
    
    #def gauss(x,*p):
        #return gaussFunc(x,p[0],p[1],p[2])
    
    RM = rmResult[:,0]
    Flux_norm = rmResult[:,1]
    Flux_norm_err = rmResult[:,2]
    
    
    # mean, sd, a
    p0 = [RM[Flux_norm.argmax()], 5 ,Flux_norm.max()]
    y_5init = gauss(RM,*p0)
    popt,pcov = curve_fit(gauss,RM[1900:2300],Flux_norm[1900:2300],p0=p0)
    perr = np.sqrt(np.diag(pcov))
    #print "Perr", perr
    #print "popt", popt
    #print Flux_norm.max()
    #print 'fit: mean=%5.3f, sd=%5.3f, A=%5.3f' % tuple(popt)
    #print 'fit: mean_err=%5.3f, sd_err=%5.3f, A_err=%5.3f' % tuple(perr)
    return popt, perr
    
    #plt.figure(figsize=[15,9])
    #plt.errorbar(RM, Flux_norm,yerr=Flux_norm_err,label='data')
    #plt.plot(RM, gauss(RM, *popt), 'r-',label='fit: mean=%5.3f, sd=%5.3f, A=%5.3f' % tuple(popt))
    #
    #plt.legend(fontsize=20)
    #plt.xlabel('RM',fontsize=17)
    #plt.ylabel('Normrized Flux',fontsize=17)
    #plt.xticks(fontsize=15)
    #plt.yticks(fontsize=15)
    #plt.show()



def getData(filename):
    data = filename.split('Timing')[0]
    data = data.split('_')[-1]
    return int(data)


filelist = glob.glob("*.RMresult")
filelist = sorted(filelist, key=getData)
print "File Number: ", len(filelist)

rmResult = {}
keys = []
axises = []
legends = []

maxRM = []

for result in filelist:
    print result
    filename = result.split('.')[0]
    rmResult[filename] = np.loadtxt(result)
    keys.append(filename)


fitResult = []
cmap = plt.cm.get_cmap('tab20')
for num, key in enumerate(keys):
    data = rmResult[key]
    fitResult.append(gaussFit(data))
    print key, data.shape, fitResult[num][0]
    #print fitResult[num]
    ax = plt.errorbar(data[:, 0],data[:, 1],yerr=data[:, 2], color = cmap.colors[num])
    plt.plot(data[:, 0], gauss(data[:, 0],*fitResult[num][0]), color = 'r')
    maxRMIndex = data[:,1].argmax()
    maxRM.append(data[:,0][maxRMIndex])
    axises.append(ax)
    legends.append(str(num).rjust(2) + ': ' + key)


    

plt.legend(axises, legends)
plt.show()



for num, i in enumerate(fitResult):
    ax1 = plt.errorbar(num,i[0][0],yerr=i[1][0], label='mean',fmt='o',ecolor='r',color='b',elinewidth=2,capsize=4)
    ax2 = plt.errorbar(num,i[0][1],yerr=i[1][1], label='variance',fmt='o',ecolor='r',color='r',elinewidth=2,capsize=4)
#    ax3 = plt.errorbar(num,i[0][2],yerr=i[1][2], label='A',fmt='o',ecolor='r',color='k',elinewidth=2,capsize=4)
    #print i[0], i[1]
#plt.legend([ax1, ax2, ax3], ['mean', 'variance', 'A'])
plt.legend([ax1, ax2], ['mean', 'variance'])
plt.axhline(y=23,c='b',ls='--',lw=2)
plt.show()

