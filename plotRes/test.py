from __future__ import print_function
import numpy as np
import math
from builtins import str, zip
from optparse import OptionParser
import bisect, os, sys, glob
import scipy, scipy.signal, scipy.stats
from readTOA import get_TOAs
from readRes import Res

def readRes(file_name):
    res = Res(file_name)
    return res


def readTOAs(file_name):
    toas = get_TOAs(file_name)
    return toas


resfile = 'res.dat'
timFile = '19C63-all_clean.tim'


residual = readRes(resfile)
toas = readTOAs(timFile)
#print(resfile)
#print(timFile)
print(toas)
print(residual)
