# uncompyle6 version 3.7.4
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.17 (default, Feb 27 2021, 15:10:58) 
# [GCC 7.5.0]
# Embedded file name: /data/miaocc/FAST-Timing/19C13/RMtoolkit/paERR.py
# Compiled at: 2020-12-22 12:37:16
import numpy as np, matplotlib.pyplot as plt
from masuiFunc import derotate
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

class getLtrue:

    def __init__(self, rm, pa0):
        self.rm = rm
        self.pa0 = pa0

    def getdata(self, I, Q, U):
        self.I = I
        self.Q = Q
        self.U = U
        offPulseRMS = get_off_pulse_rms(I)
        self.Lindex, self.Ltrue = get_Ltrue(self.Q, self.U, offPulseRMS)


def get_off_pulse_rms(allPulse):
    all_rms = np.std(allPulse - np.mean(allPulse))
    offPulse_rms = np.std(allPulse[(allPulse / all_rms < 1)])
    return offPulse_rms


def get_SNR(profile):
    off_pulse_rms = get_off_pulse_rms(I_prof)
    I_sigma_value = Q_sigma = U_sigma = V_sigma = I_prof / off_pulse_rms
    index = I_sigma_value > THRESHOLD
    W_eq = max(profile)


def get_Ltrue(Q, U, off_pulse_rms):
    L_meas = np.sqrt(U ** 2 + Q ** 2)
    Width_eq, pulseIndex, offPulseIndex = get_equivalent_width(L_meas, fitWidth=0.05)
    off_pulse_rms = np.std(L_meas[offPulseIndex])
    L_true = np.zeros(len(L_meas))
    index = L_meas / off_pulse_rms >= 1.57
    L_true = L_true
    L_true[index] = off_pulse_rms * np.sqrt((L_meas[index] / off_pulse_rms) ** 2 - 1.0)
    return (
     index, L_true)


def get_PA_sigma(L_true, off_pulse_rms):
    PA_sigma = np.zeros(len(L_true))
    index = L_true != 0
    PA_sigma[index] = 28.65 * off_pulse_rms / L_true[index]
    return PA_sigma


def I_sigma(I_prof, THRESHOLD=3):
    off_pulse_rms = get_off_pulse_rms(I_prof)
    I_sigma_value = Q_sigma = U_sigma = V_sigma = I_prof / off_pulse_rms
    index = I_sigma_value > THRESHOLD
    I_sigma_value_tmp = np.zeros(len(I_sigma_value))
    I_sigma_value_tmp[index] = I_sigma_value[index]
    return I_sigma_value_tmp


def getSNR(profile, fitWidth=0.05):
    profile = np.array(profile)
    Width_eq, pulseIndex, offPulseIndex = get_equivalent_width(profile, fitWidth=fitWidth)
    meanProf = np.mean(profile[offPulseIndex])
    sigma_offPulse = np.std(profile[offPulseIndex])
    print 'Mean of profile: %s, off pulse RMS: %s, pulse Width: %s' % (meanProf, sigma_offPulse, Width_eq)
    snr = np.sum(profile - meanProf) / np.sqrt(Width_eq) / sigma_offPulse
    return snr


def get_equivalent_width(profile, fitWidth=0.05):
    phs = np.linspace(0, 1.0, len(profile), endpoint=True)
    pulseWidth, pulseIndex, offPulseIndex = get_pulse_index(profile, fitWidth=0.05)
    area = np.trapz(y=profile[pulseIndex], x=phs[pulseIndex])
    return (
     len(profile) * area / max(profile), pulseIndex, offPulseIndex)


def get_pulse_index(profile, fitWidth=0.05):
    pulseProfile = np.array(profile)
    phase = np.linspace(0, 1.0, len(pulseProfile), endpoint=True)
    maxValue = pulseProfile.max()
    normProf = pulseProfile / maxValue
    P0 = [
     np.max(normProf), phase[normProf.argmax()], 1.0 * fitWidth, np.mean(normProf)]
    result_fit1 = leastsq(residuals, P0, args=(normProf, phase))
    pulseCenter = result_fit1[0][1]
    pulseWidth = result_fit1[0][2]
    pulseIndex = np.where(np.logical_and(phase <= pulseCenter + 4 * pulseWidth, phase >= pulseCenter - 4 * pulseWidth))[0]
    index = np.logical_not(np.ones_like(phase))
    index[pulseIndex] = True
    pulseIndex = index
    offPulseIndex = np.logical_not(pulseIndex)
    return (
     pulseWidth, pulseIndex, offPulseIndex)


def funGaussion(x, p):
    a, b, c, d = p
    return a * np.exp(-(x - b) ** 2 / (2 * c ** 2)) + d


def residuals(p, y, x):
    return y - funGaussion(x, p)


def qfunc(wavelength_meter, RM, pa0):
    pa0_rad = pa0 * np.pi / 180.0
    Q = np.cos(2 * (pa0_rad + RM * wavelength_meter ** 2))
    return Q


def ufunc(wavelength_meter, RM, pa0):
    pa0_rad = pa0 * np.pi / 180.0
    U = np.sin(2 * (pa0_rad + RM * wavelength_meter ** 2))
    return U


def chiSquare(Q_norm, U_norm, Q_fun, U_fun, sigma_Q, sigma_U):
    chiSquare_value = (Q_norm - Q_fun) ** 2 / sigma_Q + (U_norm - U_fun) ** 2 / sigma_U
    return chiSquare_value
