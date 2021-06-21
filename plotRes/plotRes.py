#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import math
from builtins import str, zip
from optparse import OptionParser
import bisect, os, sys, glob
import scipy, scipy.signal, scipy.stats
from readTOA import get_TOAs
from readRes import Res
#from presto.presto import rfft, next2_to_n
#from presto import infodata
#from presto.Pgplot import *
from Pgplot import *

"""
read res.dat
read res.dat
copy toa.tim to .toa.tim.bak
read .toa.tim.bak

get toa array
get jumps 

"""
def plotUsage():
    """
    print the plot command
    """
    print(
         "command:\n"
       + "  z: zoom in\n"
       +"  u: replot\n"
       +"  R: reload 'the res.dat'\n")

# copy the file and edit it
def backUp(file_name):
    """
    back up the file
    """
    f=open(file_name,'rb')
    new_file=open('.%s.bak' % file_name,'wb')
    for line in f.readlines():
        new_file.write(line)
    f.close()
    new_file.close()


def getJUMP(file_name):
    """
    get the jump in toas
    """
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

def plotSet():
    """
    set the plot window
    """
    ppgplot.pgpage() # open a new page
    ppgplot.pgsch(1.5)
    ppgplot.pgslw(2)
    ppgplot.pgsls(1)
    ppgplot.pgsci(1)
    ppgplot.pgpap(8, 0.5)  # Width in inches, aspect
    ppgplot.pgsvp(0.08, 0.95, 0.08, 0.95)
    ppgplot.pgupdt()
    ppgplot.pgeras()

def plotDots(mjd, res, err):
    ppgplot.pgsls(1)
    ppgplot.pgsci(1)
    ppgplot.pgpt(mjd, res, 17)
    
def plotDotsErr(mjd, res, err):
    ppgplot.pgsls(1)
    ppgplot.pgsci(1)
    ppgplot.pgerry(mjd, res+err, res-err, 1.0)

def plotJUMP(jumpPoint,yLimit):
    """
    plot the jumps
    input: jumpPoint, yLimit
    """
    #jumpPoint
    ppgplot.pgsch(1.5)
    jumpPoint = np.array(jumpPoint)
    ppgplot.pgsls(2)
    ppgplot.pgslw(3)
    ppgplot.pgsci(2)
    if len(jumpPoint)>0:
        for mjd in jumpPoint:
            print("JUMP Point:", mjd)
            mjd = math.floor(mjd)
            ppgplot.pgline(np.array([mjd, mjd]), np.array(yLimit))
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)


#def devPlot(residual, toas, xlimit, ylimit):
def devPlot(mjd, res, err, xlimit, ylimit, minMJD=None):
    """
    get the residual and toas
    get the xlimit and ylimit
    """
    ppgplot.pgupdt()
    # set the x,y limit
    xLimit = np.array([min(xlimit), max(xlimit)])
    yLimit = np.array([min(ylimit), max(ylimit)])
    ppgplot.pgeras()
    ppgplot.pgsch(1.5)
    ppgplot.pgswin(*xLimit, *yLimit)
    ppgplot.pgbox("BCNST", 0, 0, "BCNST", 0, 0)
    ppgplot.pgsch(1.3)
    ppgplot.pgmtxt('L', 2.5, 0.5, 0.5, r"Residual \gms")
    #ppgplot.pgmtxt('B', 2.0, 0.5, 0.5, "MJD + %s" %(xLimit[0]+minMJD))
    print("DEV Min MJD", minMJD)
    ppgplot.pgmtxt('B', 2.0, 0.5, 0.5, "MJD + %s" %(minMJD))

    # plot the dots
    plotDots(mjd, res, err)
    # plot the error bar
    plotDotsErr(mjd, res, err)

#    if toafile is not None:
#        jumpPoint = getJUMP(toafile) - minMJD
#        #jumpPoint = getJUMP(toafile)
#        plotJUMP(jumpPoint, yLimit)
    

def startPlot():
    """
    load the residual
    load the toa file
    set the window start and end
    cal devPlot 
    devPlot: 

    """
    # load the residual and toa
    global toas
    global residual
    if toafile is not None:
        toas = get_TOAs(toafile)
    residual = Res(resfile)

    # set the minimual MJD
    global minMJD

    res_MJD = residual.table['mjd_float'].data
    res = residual.table['res'].data
    res_Error = residual.table['error'].data

    minMJD = res_MJD.min()
    res_MJD = res_MJD - res_MJD.min()
    mjdStart = res_MJD.min()
    mjdEnd = res_MJD.max()

    minRes = res.min()
    maxRes = res.max()
    
    xlimit = [mjdStart, mjdEnd]
    ylimit = [minRes*0.9, maxRes*1.1]

    devPlot(res_MJD, res, res_Error, xlimit, ylimit, minMJD=minMJD)

def plotloop():
    residual = Res(resfile)
    # load the data for the next plot
    res_MJD = residual.table['mjd_float'].data
    res = residual.table['res'].data
    res_Error = residual.table['error'].data

    minMJD = res_MJD.min()
    res_MJD = res_MJD - res_MJD.min()
    mjdStart = res_MJD.min()
    mjdEnd = res_MJD.max()

    minRes = res.min()
    maxRes = res.max()
    
    xlimit = [mjdStart, mjdEnd]
    ylimit = [minRes*0.9, maxRes*1.1]

    while(True):
        key = None
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
            devPlot(res_MJD, res, res_Error, xlimit, ylimit, minMJD=minMJD)

        # replot
        if key=='u':
            xlimit = np.array([mjdStart, mjdEnd])
            ylimit = np.array([min(res)*0.9, max(res)*1.1])
            devPlot(res_MJD, res, res_Error, xlimit, ylimit, minMJD=minMJD)

        # reload the residual
        if key=='R':

            residual = Res(resfile)
            res_MJD = residual.table['mjd_float'].data
            res = residual.table['res'].data
            res_Error = residual.table['error'].data

            minMJD = res_MJD.min()
            print("Minn MJD", minMJD)
            res_MJD = res_MJD - res_MJD.min()
            mjdStart = res_MJD.min()
            mjdEnd = res_MJD.max()

            minRes = res.min()
            maxRes = res.max()
            
            xlimit = [mjdStart, mjdEnd]
            ylimit = [minRes*0.9, maxRes*1.1]
            ppgplot.pgeras()
            devPlot(res_MJD, res, res_Error, xlimit, ylimit, minMJD=minMJD)

        # add jump
        if key=='j':
            #print("add JMUP to TOA %s" %(mousePos[0]))
            jumpMJD = mousePos[0]+minMJD
            print("add JMUP to TOA %s" %(jumpMJD))

        # replot
        if key=='u':
            xlimit = np.array([mjdStart, mjdEnd])
            ylimit = np.array([min(res)*0.9, max(res)*1.1])
            devPlot(res_MJD, res, res_Error, xlimit, ylimit, minMJD=minMJD)


        # add jump
        if key=='h':
            plotUsage()


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

    # set the toa filename
    global toafile
    global resfile 
    toafile = opts.toafile
    if os.path.exists("res.dat"):
        resfile = "res.dat"
    else:
        print("No found residual file")
        sys.exit(0)

    # set the plot window
    plotSet()
    # start plot
    startPlot()
    # plot loop
    plotloop()

if __name__ == '__main__':
    main()
