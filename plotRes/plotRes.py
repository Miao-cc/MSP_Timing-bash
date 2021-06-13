#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import math
from builtins import str, zip
from optparse import OptionParser
import bisect, os, sys, glob
import scipy, scipy.signal, scipy.stats
#from presto.presto import rfft, next2_to_n
#from presto import infodata
#from presto.Pgplot import *
from Pgplot import *


def getJUMP(file_name):
    jumpPoint = []
    upperLine = None
    firstLine = None
    with open(file_name,'r') as f:
        for line in f.readlines():
            if(len(line)>25):
                if (line[0]!='C'):
                    firstLine = line
                    break
    with open(file_name,'r') as f:
        for line in f.readlines():
            if(len(line)>25):
                if (line[0]!='C'):
                    upperLine = line
            elif(line.startswith("JUMP")):
                jumpPoint.append(upperLine)
    if len(jumpPoint)>0:
        if jumpPoint[0] is None:
            jumpPoint[0] = firstLine
        for num, line in enumerate(jumpPoint):
            line = line.split(' ')[2]
            jumpPoint[num] = float(line)
        #print(jumpPoint) 
        return np.array(jumpPoint)
    else:
        return np.array([])


def backUp(file_name):
#def addJUMP(file_name):
    f=open(file_name,'rb')
    new_file=open('.%s.bak' % file_name,'wb')
    for line in f.readlines():
        new_file.write(line)
    f.close()
    new_file.close()


def plotJUMP(jumpPoint,yLimit):
    #jumpPoint
    jumpPoint = np.array(jumpPoint)
    ppgplot.pgsls(2)
    ppgplot.pgslw(3)
    ppgplot.pgsci(2)
    if len(jumpPoint)>0:
        for mjd in jumpPoint:
            print("JUMP Point:", mjd)
            mjd = math.floor(mjd) -0.5
            ppgplot.pgline(np.array([mjd, mjd]), np.array(yLimit))
    ppgplot.pgsls(1)
    ppgplot.pgsci(1)


# plot the residual
def loadResidual():
    data = np.loadtxt('res.dat')
    print("Load the res.dat. Data shape: ", data.shape)
    mjd = data[:,0]
    res = data[:,1]
    err = data[:,3]
    return mjd, res, err


def devPlot(mjd, res, err, xlimit, ylimit):
    ppgplot.pgupdt()
    #print(xlimit, ylimit)
    xLimit = np.array([min(xlimit), max(xlimit)])
    yLimit = np.array([min(ylimit), max(ylimit)])
    ppgplot.pgeras()
    ppgplot.pgsch(1.5)
    ppgplot.pgswin(*xLimit, *yLimit)
    ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
    ppgplot.pgsch(1.3)
    #ppgplot.pgswin(*xLimit, *yLimit)
    ppgplot.pgpt(mjd, res, 17)
    ppgplot.pgsls(1)
    ppgplot.pgerry(mjd, res+err, res-err, 1.0)
    ppgplot.pgsch(1.5)
    ppgplot.pgmtxt('L', 2.5, 0.5, 0.5, r"Residual \gms")
    ppgplot.pgmtxt('B', 2.0, 0.5, 0.5, "MJD + %s" %(xLimit[0]+minMJD))
    #ppgplot.pgline(xLimit, yLimit)
    #ppgplot.pgline(np.array([10, 10]), np.array(yLimit)) 

    if toafile is not None:
        jumpPoint = getJUMP(toafile) - minMJD
        #jumpPoint = getJUMP(toafile)
        plotJUMP(jumpPoint, yLimit)
    

usage = "usage: %prog -x [-f temp.tim]"
def main():
    parser = OptionParser(usage)
    parser.add_option("-x", "--xwin", action="store_true", dest="xwin",
                      default=False, help="Open an X-window and plot")
    parser.add_option("-f", "--file", type="string", dest="toafile", default=None,
                      help="load a TOA file")
    (opts, args) = parser.parse_args()
    if opts.xwin:
        pgplot_device = "/XWIN"
    else:
        pgplot_device = ""

    if pgplot_device:
        print("Open the pgplot device")
        ppgplot.pgopen(pgplot_device)
    global toafile
    toafile = opts.toafile

    ppgplot.pgpage() # open a new page
    ppgplot.pgpap(8, 0.5)  # Width in inches, aspect
    ppgplot.pgsvp(0.08, 0.95, 0.08, 0.95)
    ppgplot.pgupdt()
    ppgplot.pgeras()
    startPlot()


def startPlot():
    mjd, res, err = loadResidual()
    global minMJD
    minMJD = min(mjd)
    mjd = mjd - minMJD
    mjdStart = min(mjd) - 10
    mjdEnd = max(mjd) + 10
    
    xlimit = [mjdStart, mjdEnd]
    ylimit = [min(res)*0.9, max(res)*1.1]
    devPlot(mjd, res, err, xlimit, ylimit)
    while(True):
        mousePos = ppgplot.pgband(7)
        key = str(mousePos[2],encoding='utf-8')
        mousePos = [mousePos[0], mousePos[1]]
        print("MJD: ", mousePos[0]+minMJD)
        #print("MJD: ", mousePos[0]+minMJD, key)
        #ppgplot.pgtext(mousePos[0], mousePos[1], "We clicked here.")
        ppgplot.pgupdt()
        # quit the plot
        if key=='q':
            break

        # zoom in
        if key=='z':
            replacePos = ppgplot.pgband(2,0, mousePos[0], mousePos[1])
            replacePos = [replacePos[0], replacePos[1]]
            xlimit = np.array([mousePos[0], replacePos[0]])
            ylimit = np.array([min(res)*0.9, max(res)*1.1])
            devPlot(mjd, res, err, xlimit, ylimit)

        # replot
        if key=='u':
            xlimit = np.array([mjdStart, mjdEnd])
            ylimit = np.array([min(res)*0.9, max(res)*1.1])
            devPlot(mjd, res, err, xlimit, ylimit)

        # reload the residual
        if key=='R':
            mjd, res, err = loadResidual()
            minMJD = min(mjd)
            mjd = mjd - minMJD
            mjdStart = min(mjd) - 10
            mjdEnd = max(mjd) + 10
            
            xlimit = [mjdStart, mjdEnd]
            ylimit = [min(res)*0.9, max(res)*1.1]
            devPlot(mjd, res, err, xlimit, ylimit)

        ## add jump
        #if key=='j':
        #    #print("add JMUP to TOA %s" %(mousePos[0]))
        #    jumpMJD = mousePos[0]+minMJD
        #    print("add JMUP to TOA %s" %(jumpMJD))

        # add jump
        if key=='h':
            print(
"""
command:
  z: zoom in
  u: replot
  R: reload 'the res.dat'
""")


if __name__ == '__main__':
    main()
