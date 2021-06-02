import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const

def rebin_freq(spectra, mask_chans, rebin_fact):
    # get the spectra shape, frequence and polarization
    nfreq ,naxis1 = spectra.shape
    print "spectra shape: ", nfreq ,naxis1

    # deep copy spectra
    spectra = spectra.copy()
    # mask the bad channels
    spectra[mask_chans] = 0
    # rebin the spectra
    spectra_rebin = np.reshape(spectra, (nfreq // rebin_fact, rebin_fact, naxis1))
    print "spectra_rebin", spectra_rebin.shape
    spectra_rebin = np.sum(spectra_rebin, 1)

    # get the mask channels
    norm = np.reshape(np.logical_not(mask_chans), (nfreq // rebin_fact, rebin_fact))
    print "reshaped norm: ", norm.shape
    norm = np.sum(norm, 1).astype(float)
    print "summed norm: ", norm.shape, "if channel hadn't zapped, set the number +1"
    bad = norm == 0
    norm[bad] = 1
    norm = 1./norm
    norm[bad] = 0
    print "spectra_rebin", spectra_rebin.shape
    print "norm.shape", norm.shape
    spectra_rebin *= norm[:,None]
    return spectra_rebin

def integrated_pulse_spectrum(data, freq, time, pars, matched=True):
    t0 = pars[0]
    dm = pars[1]
    amp_800 = pars[2]
    alpha = pars[3]
    width = pars[4]
    scatter_800 = pars[5]

    scatter = scatter_800 * (freq / 800.)**-4
    delay = delay_from_dm(freq, dm, t0)

    delta_t = np.median(np.diff(time))
    delta_f = np.median(np.abs(np.diff(freq)))

    out = np.empty(data.shape[:-1], dtype=float)

    time_selector = RangeSelector(time)

    for ii in range(len(freq)):
        if matched:
            nw = 10.
        else:
            nw = 3.
        window_low = nw * width
        window_high = nw * width + nw**2 / 2 * scatter[ii]
        start_ind, stop_ind = time_selector(delay[ii] - window_low, delay[ii] + window_high)
        near_times = time[start_ind:stop_ind].copy()
        if matched:
            pulse = windowed_pulse(near_times, freq[[ii]], delta_t, delta_f, dm, t0, width, scatter[ii])
            # Matched filter normalization.
            pulse /= np.sum(pulse**2) * delta_t
            # Compensated? Probably don't need to if noise is white.
            #pulse -= np.mean(pulse)
        else:
            pulse = np.ones_like(near_times)
        for jj in range(data.shape[1]):
            out[ii, jj] = np.sum(data[ii,jj,start_ind:stop_ind] * pulse)
    out *= delta_t
    return out



def derotate(data, freq, RM, phi0):
    """
    derotate the data 
    if fit_phi == phi0, then phi == pol_angle
    Q_ = cos(0) = cos(phi-pol_angle) = cos(phi)cos(pol_angle) + sin(phi)sin(pol_angle)
    U_ = sin(0) = sin(phi-pol_angle) = sin(phi)cos(pol_angle) - cos(phi)sin(pol_angle)
    set frequence refer to inf
    """
    c = const.c.to('m/s').value
    pol_angle = 2 * RM *  (c / freq / 1e6)**2 + phi0 * np.pi/180.

    nsubint, npol, nfreq, nbin = data.shape
    #print "data shape: subint: %s pol: %s frequence: %s bin: %s" %(nsubint, npol, nfreq, nbin)

    # deep copy the Q and U and data
    # repeat the data in a new axis
    Q = copy.deepcopy(data[:, 1, :, :])
    U = copy.deepcopy(data[:, 2, :, :])
    data_derot = copy.deepcopy(data)

    #Q = data[:, 1, :, :]
    #U = data[:, 2, :, :]
    #data_derot = data

    data_derot[:, 1, :, :] =  Q * np.cos(pol_angle)[None, :, None] + U * np.sin(pol_angle)[None, :, None]
    data_derot[:, 2, :, :] = -Q * np.sin(pol_angle)[None, :, None] + U * np.cos(pol_angle)[None, :, None]
    return data_derot



def derotate_pa0(data, freq, RM, phi0):
    """
    derotate the data 
    if fit_phi == phi0, then phi == pol_angle
    Q_ = cos(0) = cos(phi-pol_angle) = cos(phi)cos(pol_angle) + sin(phi)sin(pol_angle)
    U_ = sin(0) = sin(phi-pol_angle) = sin(phi)cos(pol_angle) - cos(phi)sin(pol_angle)
    set frequence refer to inf
    """
    c = const.c.to('m/s').value
    pol_angle = np.zeros((len(freq), len(phi0)))
    pol_angle = 2 * RM *  (c / freq[:,None] / 1e6)**2 + phi0[None,:] * np.pi/180.

    nsubint, npol, nfreq, nbin = data.shape
    #print "data shape: subint: %s pol: %s frequence: %s bin: %s" %(nsubint, npol, nfreq, nbin)

    # deep copy the Q and U and data
    # repeat the data in a new axis
    Q = copy.deepcopy(data[:, 1, :, :])
    U = copy.deepcopy(data[:, 2, :, :])
    data_derot = copy.deepcopy(data)

    #Q = data[:, 1, :, :]
    #U = data[:, 2, :, :]
    #data_derot = data

    data_derot[:, 1, :, :] =  Q * np.cos(pol_angle)[None, :, None] + U * np.sin(pol_angle)[None, :, None]
    data_derot[:, 2, :, :] = -Q * np.sin(pol_angle)[None, :, None] + U * np.cos(pol_angle)[None, :, None]
    return data_derot
