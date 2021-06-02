import sys
import math 
import time
import psrchive
import matplotlib
import numpy as np
from paERR import *
from masuiFunc import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.optimize import leastsq
#from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from astropy import constants as const
import multiprocessing as MP

# multi processs func
def multiProcesss(func, params, ncpus):
    num_cpus = max(1, MP.cpu_count() - 1)
    if ncpus>num_cpus:
        print "The max cpu number is: ", num_cpus+1
        ncpus=num_cpus

    p = MP.Pool(processes=ncpus)
    #print params
    for param1 in params[0]:
        new_args =  [param1] + params[1:]
        p.apply_async(func, args=new_args)
    p.close()
    p.join()

def outResult(outFileName, data, header = None):
    if header is None:
        np.savez(outFileName, data)
    else:
        np.savez(outFileName, data, header=header)

    print "Writing file: ", outFileName


def loadData(filename):
    """
    load the data
    """

    # set the global parmaters
    setParams()
    
    filename = sys.argv[1]
    arch = psrchive.Archive_load(filename)
    # archive file dedisperse
    # archive file remove baseline
    # archive file dedisperse
    arch.dedisperse()
    arch.remove_baseline()
    arch.convert_state('Stokes')
    arch.centre_max_bin()
    # get the data 
    data = arch.get_data()
    #arch.convert_state('Stokes')

    subintNum, polNum, chanNum, profBin = data.shape
    # get frequence information
    freq_lo = arch.get_centre_frequency() - arch.get_bandwidth()/2.0
    freq_hi = arch.get_centre_frequency() + arch.get_bandwidth()/2.0
    # get frequence in hz and wavelength in meter
    freq_hz = np.linspace(freq_lo*1e6, freq_hi*1e6, num=chanNum, endpoint=True)
    wavelength_meter = pulsarData['c'] / freq_hz

    pulsarData['freq_lo'] = freq_lo
    pulsarData['freq_hi'] = freq_hi
    pulsarData['wavelength_meter'] = wavelength_meter
    pulsarData['freq_hz'] = freq_hz
    pulsarData['subintNum'] = subintNum
    pulsarData['polNum'] = polNum
    pulsarData['chanNum'] = chanNum
    pulsarData['profBin'] = profBin
    pulsarData['data'] = data
    pulsarData['startPhase'] = 0
    pulsarData['endPhase'] = profBin
    #print data.shape
    #plt.plot(data[:,0,:,:].sum(0).sum(0))
    #plt.show()

    print "Frequence: %s MHz, %s MHz" %(freq_lo, freq_lo)
    print "wavelength: %s m, %s m" %(wavelength_meter[0], wavelength_meter[-1])

    #return pulsarData

def plot_profile(data, showFile=True):
    matplotlib.rcParams.update({'font.size': 16,})
    start_subint = pulsarData['start_subint']
    end_endint = pulsarData['end_endint']
    profBin = pulsarData['profBin']

    I_prof = data[start_subint:end_endint,0,:,:].sum(0).sum(0)
    Q_prof = data[start_subint:end_endint,1,:,:].sum(0).sum(0)
    U_prof = data[start_subint:end_endint,2,:,:].sum(0).sum(0)
    V_prof = data[start_subint:end_endint,3,:,:].sum(0).sum(0)
    print "Profile shape: \nI %10s \nQ %10s \nU %10s \nV %10s" %(I_prof.shape,  Q_prof.shape, U_prof.shape, V_prof.shape)
    Liner_prof = np.sqrt(Q_prof**2 + U_prof**2)


    #PA = 0.5*np.arctan2(U_prof, Q_prof) / np.pi * 180 % 180
    PA = 0.5*np.arctan2(U_prof, Q_prof) * 180 / np.pi
    off_pulse_rms = get_off_pulse_rms(I_prof)
    index, PAsigma = get_Ltrue(Q_prof, U_prof, off_pulse_rms)/off_pulse_rms
    print index
    PAsigma = PAsigma[index]
    PAerr = np.array([PA[i]*(1 - stats.norm.cdf(PAsigma[i]*1/2.)) for i in range(profBin)])
    PAindex = PAsigma > 1
    #PAindex = PAerr < 5
    
    #fig = plt.figure(figsize=(9, 16))
    fig = plt.figure()
    gs = gridspec.GridSpec(nrows=4, ncols=1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1:, 0])
    ax1.set_xlim(0,profBin)
    ax2.set_xlim(0,profBin)
    ax1.set_ylim(-100, 100)
    ax1.set_yticks(np.arange(-90, 90, 30))
    ax1.set_ylabel('P.A. (deg.)')
    ax2.set_xlabel("Phase bin")
    ax2.set_ylabel("Profile Value")
    
    ax1.errorbar(np.arange(profBin)[PAindex], PA[PAindex], yerr=PAerr[PAindex], fmt="bo:",color="blue",ecolor='grey', elinewidth=2,capsize=4 )
    
    ax2.plot(I_prof, 'k')
    ax2.plot(Liner_prof, 'r')
    ax2.plot(V_prof, 'b')
    if showFile is True:
        plt.show()

def plot_FaradayRot():
    matplotlib.rcParams.update({'font.size': 16,})
    wavelength_meter = pulsarData['wavelength_meter']
    PA0 = pulsarData['PA0']
    RM = pulsarData['RM']
    freq_hz = pulsarData['freq_hz']

    #print PA0, RM
    Q_fun = Qfunc(wavelength_meter, PA0, RM)
    U_fun = Ufunc(wavelength_meter, PA0, RM)
    fig = plt.figure()
    ax1 = plt.subplot(311)
    plt.plot(freq_hz/1.e6, Q_fun, '.k', linewidth = 15)
    plt.xticks([])
    plt.ylim(-1.05, 1.05)
    plt.ylabel('Q')
    ax2 = plt.subplot(312, sharex=ax1)
    plt.plot(freq_hz/1.e6, U_fun, '.k', linewidth = 15)
    plt.xlim(freq_hz[0]/1.e6, freq_hz[-1]/1.e6)
    plt.ylim(-1.05, 1.05)
    plt.ylabel("U")
    #ax3 = plt.subplot(313, sharex=ax1)
    ax3 = plt.subplot(313)
    plt.plot(freq_hz/1.e6, np.sqrt(Q_fun**2 + U_fun**2), '.k', linewidth = 15)
    plt.xlim(freq_hz[0]/1.e6, freq_hz[-1]/1.e6)
    plt.ylim(-1.05, 1.05)
    plt.ylabel("Linear")
    plt.xlabel("Frequence (MHz)")
    plt.subplots_adjust(wspace = 0, hspace = 0)
    plt.show()


def get_zap_Chan(data):
    print "***************************Zap channels***************************"

    start_subint = pulsarData['start_subint']
    end_endint = pulsarData['end_endint']

    I_freq = data[start_subint:end_endint,0,:,:].sum(0)
    zap_Chan_Index = I_freq.sum(axis=1) == 0
    Chan_Index = (1-zap_Chan_Index).astype(np.bool)

    pulsarData['non_Chan_Index'] = Chan_Index
    pulsarData['zap_Chan_Index'] = zap_Chan_Index
    print "zapped channel number: %s" %(np.sum(zap_Chan_Index == True))



def RMgrid(rmA):
    t_Start = time.time()
    RM = pulsarData['RM']
    PA0 = pulsarData['PA0']
    data = pulsarData['data']
    profBin = pulsarData['profBin']
    start_subint = pulsarData['start_subint']
    end_endint = pulsarData['end_endint']
    freq = pulsarData['freq_hz']/1e6
    snrs = np.zeros_like(rmA)
    
    data_rotated = data[start_subint:end_endint,:,:,:]
    print data_rotated.shape

    for i, rmValue in enumerate(rmA):
        print data_rotated.shape
        #data_derot = derotate(data, freq, rmValue, PA0)
        #I = data_derot[start_subint:end_endint,0,:,:].sum(0).sum(0)
        #Q = data_derot[start_subint:end_endint,1,:,:].sum(0).sum(0)
        #U = data_derot[start_subint:end_endint,2,:,:].sum(0).sum(0)
###
        data_derot = derotate(data_rotated, freq, rmValue, PA0)
        print data_derot.shape
        I = data_derot[:,0,:,:].sum(0).sum(0)
        Q = data_derot[:,1,:,:].sum(0).sum(0)
        U = data_derot[:,2,:,:].sum(0).sum(0)
###
        #L_true = getLtrue(rmValue, PA0)
        #L_true.getdata(I, Q, U)
        #snr = getSNR(L_true.Ltrue, fitWidth = 0.5)
###

        L = np.sqrt(Q**2 + U **2)
        snr = getSNR(L, fitWidth = 0.5)


        t_End = time.time()
        #print "Using time: %s sec." %(t_End - t_Start)
        #print L_true.Ltrue
        
        outFileName = outRootName+"_RM_"+str(round(rmValue,3))
        print "Out file: ", outFileName
        #headerInfo = {'rmA':rmA, 'rmValue':rmValue, 'SNR':snr}
        headerInfo = np.array([rmValue, snr])
        #outResult(outFileName, L_true.Ltrue, header = headerInfo)
        outResult(outFileName, L, header = headerInfo)

        print "RM %s   SRN %s" %(rmValue, snr)
        snrs[i] = snr

    pulsarData['data'] = None
    #plt.plot(rmA, snrs)
    #plt.axvline(x=rmA[snrs.argmax()])
    #plt.xlabel("RM")
    #plt.ylabel("SNR")
    #plt.show()


"""
def RMgrid(rmA):
    t_Start = time.time()
    RM = pulsarData['RM']
    PA0 = pulsarData['PA0']
    data = pulsarData['data']
    profBin = pulsarData['profBin']
    start_subint = pulsarData['start_subint']
    end_endint = pulsarData['end_endint']
    freq = pulsarData['freq_hz']/1e6
    snrs = np.zeros_like(rmA)

    ncpus = 16
    multiProcesss(get_SNR_from_derotate, [rmA, PA0, data, freq], ncpus)

def get_SNR_from_derotate(rmValue, PA0, data, freq):
    data_derot = derotate(data, freq, rmValue, PA0)
    I = data_derot[start_subint:end_endint,0,:,:].sum(0).sum(0)
    Q = data_derot[start_subint:end_endint,1,:,:].sum(0).sum(0)
    U = data_derot[start_subint:end_endint,2,:,:].sum(0).sum(0)
#    L_true = getLtrue(rmValue, PA0)
#    L_true.getdata(I, Q, U)
#    snr = getSNR(L_true.Ltrue, fitWidth = 0.5)
    L = np.sqrt(Q**2 + U **2)
    snr = getSNR(L, fitWidth = 0.5)
    print "RM %s   SRN %s" %(rmValue, snr)

"""

#    Liner_prof = np.sqrt(Q_prof**2 + U_prof**2)
#    print "Q %10s \nU %10s \nL %10s" %(Q_prof.shape, U_prof.shape, Liner_prof.shape)
#
#    #PA = 0.5*np.arctan2(U_prof, Q_prof) / np.pi * 180 % 180
#    PA = 0.5*np.arctan2(U_prof, Q_prof) * 180 / np.pi
#    off_pulse_rms = get_off_pulse_rms(I_prof)
#    PAsigma = get_Ltrue(Q_prof, U_prof, off_pulse_rms)/off_pulse_rms
#    PAerr = np.array([PA[i]*(1 - stats.norm.cdf(PAsigma[i]*1/2.)) for i in range(profBin)])
#    PAindex = PAsigma > 1
#    #PAindex = PAerr < 5




def plot_spectra(data, zap_Chan_Index):
    matplotlib.rcParams.update({'font.size': 16,})
    spectra = data[start_subint:end_endint, :, :, startPhase:endPhase].sum(axis=0).sum(axis=2).T
    spec_rebin = rebin_freq(spectra, zap_Chan_Index, rebin_fact)
    print spec_rebin.shape

    freq = freq_hz/1e6
    f_resample = freq[rebin_fact//2::rebin_fact]
    wave_resample = wavelength_meter[rebin_fact//2::rebin_fact]
    # get the index of non 0 value
    none_zero_index = (spec_rebin!=0)
    print none_zero_index.shape, wave_resample[0], wave_resample[-1]

    
    Liner_freq = np.sqrt(spec_rebin[:, 1]**2 + spec_rebin[:, 2]**2)
    spec_rebin[none_zero_index[:,1], 1] = spec_rebin[none_zero_index[:,1], 1]/Liner_freq[none_zero_index[:,1]]
    spec_rebin[none_zero_index[:,2], 2] = spec_rebin[none_zero_index[:,1], 2]/Liner_freq[none_zero_index[:,1]]


    #popt, pcov = curve_fit(Qfunc, wave_resample[none_zero_index[:,1]], spec_rebin[:,1][none_zero_index[:,1]], maxfev=10000)
    #popt, pcov = curve_fit(Q_normfunc, wave_resample[none_zero_index[:,1]], spec_rebin[:,1][none_zero_index[:,1]], maxfev=10000)
    #print popt, pcov
    #Liner_freq_func = np.sqrt(Qfunc(wavelength_meter, *popt)**2 + Qfunc(wavelength_meter, *popt)**2)
    print wave_resample[none_zero_index[:,1]].shape, spec_rebin[:,1][none_zero_index[:,1]].shape
    para = leastsq(error, [PA0, RM], args=(wave_resample[none_zero_index[:,1]], spec_rebin[:,1][none_zero_index[:,1]]))
    print para


"""
    #########
    print "spec_rebin shape: ", spec_rebin.shape
    Q_, U_ = derotateRM(spec_rebin[:,:][none_zero_index[:,1]], f_resample[none_zero_index[:,1]], RM, 140)

    # Plot polarized spectra.
    fig = plt.figure(figsize=(12,9))
    ax1 = plt.subplot(211)
    plt.title(filename)
    #plt.plot(f_resample[none_zero_index[:,1]], spec_rebin[:,1][none_zero_index[:,1]], '.k', linewidth = 15)
    plt.plot(f_resample[none_zero_index[:,1]], spec_rebin[:,1][none_zero_index[:,1]], '.k', linewidth = 10)
    plt.plot(f_resample, Q_normfunc(wave_resample, *para[0]), '.r', linewidth = 10)
    plt.plot(f_resample[none_zero_index[:,1]], spec_rebin[:,1][none_zero_index[:,1]] - Q_normfunc(wave_resample[none_zero_index[:,1]], *para[0]), '.b', linewidth = 15)
    plt.legend(["data", "func", "data - func"], loc = 'lower right')
    #plt.plot(f_resample, Q_normfunc(wave_resample, PA0, RM), '.r', linewidth = 15)
    plt.ylabel('Normalized Q')
    ax2 = plt.subplot(212, sharex=ax1)
    #plt.plot(f_resample[none_zero_index[:,2]], spec_rebin[:,2][none_zero_index[:,2]], '.k', linewidth = 15)
    plt.plot(f_resample[none_zero_index[:,2]], spec_rebin[:,2][none_zero_index[:,2]], '.k', linewidth = 10)
    plt.plot(f_resample, U_normfunc(wave_resample, *para[0]), '.r', linewidth = 10)
    plt.plot(f_resample[none_zero_index[:,2]], spec_rebin[:,2][none_zero_index[:,2]] - U_normfunc(wave_resample[none_zero_index[:,2]], *para[0]), '.b', linewidth = 15)
    plt.legend(["data", "func", "data - func"], loc = 'lower right')
    #plt.plot(f_resample, U_normfunc(wave_resample, PA0, RM), '.r', linewidth = 15)
    plt.ylabel("Normalized U")
    plt.xlabel("Frequence (MHz)")
    plt.xlim(freq[0], freq[-1])
    plt.show()
"""


def error(param, wavelength, Q):
    return Q_normfunc(wavelength, *param) - Q


# model for Q, a function of wavelegth 
def Qfunc(wavelength_meter, pa0, rm):
    # change the pa0 to radians, rm: rad*m^-2, wavelength: m^-2
    pa0_rad = pa0*np.pi/180.
    Q_func = np.cos(2*(pa0_rad+rm*wavelength_meter**2)) * wavelength_meter
    return Q_func


# model for U, a function of wavelegth
def Ufunc(wavelength_meter, pa0, rm):
    # change the pa0 to radians, rm: rad*m^-2, wavelength: m^-2
    pa0_rad = pa0*np.pi/180.
    U_func = np.sin(2*(pa0_rad+rm*wavelength_meter**2)) * wavelength_meter
    return U_func


# model for normalized Q, a function of wavelegth 
def Q_normfunc(wavelength_meter, pa0, rm):
    # change the degree pa0 to radians, rm: rad*m^-2, wavelength: m^-2
    pa0_rad = pa0*np.pi/180.
    Q_norm = np.cos(2*(pa0_rad+rm*wavelength_meter**2))
    return Q_norm


# model for normalized U, a function of wavelegth
def U_normfunc(wavelength_meter, pa0, rm):
    # change the degree pa0 to radians, rm: rad*m^-2, wavelength: m^-2
    pa0_rad = pa0*np.pi/180.
    U_norm = np.sin(2*(pa0_rad+rm*wavelength_meter**2))
    return U_norm


## get he chi square of the result
#def chiSquare(Q_norm, U_norm, Q_fun, U_fun, sigma_Q, sigma_U):
#    #chiSquare_array = (Q_norm - Q_fun)**2 + (U_norm - U_fun)**2
#    chiSquare_array = ((Q_norm - Q_fun) / sigma_Q)**2  + ((U_norm - U_fun) / sigma_U)**2  chiSquare_value = chiSquare_array.sum(axis=0).sum(axis=0)
#    return chiSquare_value


"""

# get data shape
subintNum, polNum, chanNum, profBin = data.shape
print "load data: %20s, data shape: %20s" %(filename, data.shape)

I_all = data[start_subint:end_endint,0,:,:].sum(0).sum(0)
Q_all = data[start_subint:end_endint,1,:,:].sum(0).sum(0)
U_all = data[start_subint:end_endint,2,:,:].sum(0).sum(0)
print "Profile shape: \nI %10s \nQ %10s \nU %10s" %(I_all.shape,  Q_all.shape, U_all.shape)




data = data[start_subint:end_endint,:,:,startPhasea:endPhase]

plt.show()
sys.exit(0)


phawin = endPhase - startPhasea

##

# sum in time
I_freq = data[:,0,:,:].sum(axis=0)
Q_freq = data[:,1,:,:].sum(axis=0)
U_freq = data[:,2,:,:].sum(axis=0)
V_freq = data[:,3,:,:].sum(axis=0)



sys.exit(0)



print "freq shape: ", I_freq.shape, Q_freq.shape, U_freq.shape, V_freq.shape
print "freq sum: ", np.sum(Q_freq), np.sum(U_freq)
print "freq min: ", Q_freq.min(), U_freq.min()
print "freq max: ", Q_freq.max(), U_freq.max()

# sum in freq and time
I_prof = I_freq.sum(axis=0)
Q_prof = Q_freq.sum(axis=0)
U_prof = U_freq.sum(axis=0)
V_prof = V_freq.sum(axis=0)
print "prof shape: ", I_prof.shape, Q_prof.shape, U_prof.shape, V_prof.shape

# get the normalized Q and U, no summed in frequence
LinerPol_freq = np.sqrt(U_freq**2 + Q_freq**2)
Q_freq_norm = Q_freq / LinerPol_freq
U_freq_norm = U_freq / LinerPol_freq
print Q_freq[0,0], U_freq[0,0], LinerPol_freq[0,0]
print "norm freq shape: ", Q_freq_norm.shape, U_freq_norm.shape

plt.figure()
#fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
fig, axes = plt.subplots(1, 3)
#fig.subplots_adjust(wspace = .1,hspace = 0)
#im1 = axes[0].imshow(Q_freq_norm, origin='lower', aspect='auto', cmap='binary')
#im2 = axes[1].imshow(U_freq_norm, origin='lower', aspect='auto', cmap='binary')


tmp = np.zeros((chanNum, phawin))
tmp[tmp==0] = float('nan')
tmp[Chan_Index,:] = Q_freq
im1 = axes[0].imshow(tmp, origin='lower', aspect='auto', cmap='binary')
tmp[Chan_Index,:] = U_freq
im2 = axes[1].imshow(tmp, origin='lower', aspect='auto', cmap='binary')
tmp[Chan_Index,:] = LinerPol_freq
im3 = axes[2].imshow(tmp, origin='lower', aspect='auto', cmap='binary')
 
cb = fig.colorbar(im2, ax=axes[2])
cb.set_label("colorbar")
plt.tight_layout()
plt.show()


#####################
# get the variance, no off pulse
Q_freq_sigma = np.std(Q_freq, axis=1)
U_freq_sigma = np.std(U_freq, axis=1)
print "Q_freq_sigma shape: ", Q_freq_sigma.shape, Q_freq_sigma.shape

# calculate normalized L sigma
L_freq_norm_sigma = np.sqrt((np.dot(Q_freq_sigma, U_freq)**2 + np.dot(U_freq_sigma, Q_freq)**2) / np.power(LinerPol_freq, 4))
print "L_freq_norm_sigma shape: ", L_freq_norm_sigma.shape

sys.exit(0)

#################################################
print "#"*50
# get normalized Q and normalized U sigma
Q_freq_norm_sigma = np.sqrt(np.multiply(L_freq_norm_sigma[:,100], U_freq_norm[:,100])**2)
U_freq_norm_sigma = np.sqrt(np.multiply(L_freq_norm_sigma[:,100], Q_freq_norm[:,100])**2)
print Q_freq_norm.shape, L_freq_norm_sigma.shape
print "Q_freq_norm_sigma shape: ", Q_freq_norm_sigma.shape, U_freq_norm_sigma.shape
###################

chiSquare_result = np.zeros((360,500))
for p in range(0,360):
    for rm in range(-500,0):
        pa = p * np.pi/180.
        #res = np.sum((Q_freq_norm[:,100] - Q_normfunc(wavelength_meter, pa, rm))**2 + (U_freq_norm[:,100] - U_normfunc(wavelength_meter, pa, rm))**2)
        res = chiSquare(Q_freq_norm[:,100], U_freq_norm[:,100], Q_normfunc(wavelength_meter, pa, rm), U_normfunc(wavelength_meter, pa, rm), Q_freq_norm_sigma, U_freq_norm_sigma)
        chiSquare_result[p,rm+500] = res
        print "PA: ", p, "rm: " ,rm, res

plt.figure()
plt.imshow(chiSquare_result, origin='lower', aspect='auto', cmap='binary')
plt.colorbar()
minChiSquare = np.unravel_index(chiSquare_result.argmin(), chiSquare_result.shape)
print minChiSquare
#print "PA: %s, RM: %s" %(round(np.arange(180)[minChiSquare[0]],2), round(rmA[minChiSquare[1]], 1))
print "PA: %s, RM: %s" %(round(range(0,360)[minChiSquare[0]],2), round(range(-500,0)[minChiSquare[1]], 1))
plt.figure()

Q_normfunc_array = Q_normfunc(wavelength_meter, 75, -5)
U_normfunc_array = U_normfunc(wavelength_meter, 75, -5)
#plt.plot(wavelength_meter, Q_normfunc_array,'r')
#plt.plot(wavelength_meter, U_normfunc_array,'b')
plt.scatter(wavelength_meter, Q_freq_norm[:,100],'r')
plt.scatter(wavelength_meter, U_freq_norm[:,100],'b')
#plt.plot(wavelength_meter, Q_freq_norm[:,100].reshape(4345/5,5).mean(axis=1)
       #- U_freq_norm[:,100].reshape(4345/5,5).mean(axis=1))

plt.show()


sys.exit(0)
print "freq_norm shape: ", Q_freq_norm.shape, U_freq_norm.shape, LinerPol_freq.shape
print "freq_norm sum: ", np.sum(Q_freq_norm), np.sum(U_freq_norm), np.sum(LinerPol_freq)
print "freq_norm min: ", Q_freq_norm.min(), U_freq_norm.min(), LinerPol_freq.min()
print "freq_norm max: ", Q_freq_norm.max(), U_freq_norm.max(), LinerPol_freq.max()

# get the variance, no off pulse
Q_freq_sigma = np.std(Q_freq, axis=1)
U_freq_sigma = np.std(U_freq, axis=1)
print "Q_freq_sigma shape: ", Q_freq_sigma.shape, Q_freq_sigma.shape

# calculate normalized L sigma
L_freq_norm_sigma = np.sqrt((np.dot(Q_freq_sigma, U_freq)**2 + np.dot(U_freq_sigma, Q_freq)**2) / np.power(LinerPol_freq, 4))
print "L_freq_norm_sigma shape: ", L_freq_norm_sigma.shape

#################################################
print "#"*50
# get normalized Q and normalized U sigma
Q_freq_norm_sigma = np.sqrt(np.multiply(L_freq_norm_sigma, U_freq_norm)**2)
U_freq_norm_sigma = np.sqrt(np.multiply(L_freq_norm_sigma, Q_freq_norm)**2)
print Q_freq_norm.shape, L_freq_norm_sigma.shape
print "Q_freq_norm_sigma shape: ", Q_freq_norm_sigma.shape, U_freq_norm_sigma.shape


#pa0 = np.ones(profBin)*100
#pa0 = np.ones(profBin)*0
#rm = -1000

chiSquare_result = np.zeros((180, rm_len))
count = 0
#'''
for i in range(180):
    for j in range(rm_len):
        #rm = i-50 - 22
        #rm = i-200/2.
        pa0 = np.ones(phawin)*i
        rm = rmA[j]
        Q_normfunc_array = np.array([Q_normfunc(wavelength_meter, pa, rm) for pa in pa0]).T
        U_normfunc_array = np.array([U_normfunc(wavelength_meter, pa, rm) for pa in pa0]).T
        #print "U_normfunc_array shape: ", U_normfunc_array.shape, Q_normfunc_array.shape
        
        
        #print "Q_freq_norm, U_freq_norm, Q_normfunc_array, U_normfunc_array, Q_freq_norm_sigma, U_freq_norm_sigma"
        #print Q_freq_norm.shape, U_freq_norm.shape, Q_normfunc_array.shape, U_normfunc_array.shape, Q_freq_norm_sigma.shape, U_freq_norm_sigma.shape
        #print np.sum(Q_freq_norm[Chan_Index, :]), np.sum(U_freq_norm[Chan_Index, :]), np.sum(Q_normfunc_array[Chan_Index, :]), np.sum(U_normfunc_array[Chan_Index, :]), np.sum(Q_freq_norm_sigma[Chan_Index, :]), np.sum(U_freq_norm_sigma[Chan_Index, :])
        
        chiSquare_value = chiSquare(Q_freq_norm, U_freq_norm, Q_normfunc_array, U_normfunc_array, Q_freq_norm_sigma, U_freq_norm_sigma)
        chiSquare_result[i,j] = (chiSquare_value)
        count += 1
        print count / (180.*rm_len)
        #print chiSquare_value
#'''


plt.figure()
plt.imshow(I_freq, origin='lower', aspect='auto')

plt.figure()
plt.plot(np.arange(phawin)+startPhasea, I_prof, linewidth =5)
plt.plot(I_all, linewidth =2.5)
plt.plot(np.sqrt(Q_all**2 + U_all**2),'r')
plt.axvline(x=startPhasea, color='k')
plt.axvline(x=endPhase, color='k')

#Q_ = np.cos(2*-22*wavelength_meter)*Q_freq + np.sin(2*-22*wavelength_meter)*U_freq
#U_ = -np.sin(2*-22*wavelength_meter)*Q_freq + np.cos(2*-22*wavelength_meter)*U_freq
#Lsum = np.sqrt(Q_.sum(axis=0)**2 + U_.sum(axis=0)**2)
#plt.plot(Lsum, linewidth=3.5)
#plt.plot(np.sqrt(Q_prof**2 + U_prof**2), linewidth=2)

plt.figure()
#plt.plot(np.arange(360), chiSquare_result)
plt.imshow(chiSquare_result, origin='lower', aspect='auto', cmap='binary')
minChiSquare = np.unravel_index(chiSquare_result.argmin(), chiSquare_result.shape)
print minChiSquare
plt.scatter(minChiSquare[1],minChiSquare[0] ,color='b')
plt.text(minChiSquare[1]+3,minChiSquare[0], "PA: %s, RM: %s" %(round(np.arange(180)[minChiSquare[0]],2), round(rmA[minChiSquare[1]], 1)), color = 'b')
xlabel_ticks = np.arange(0,rm_len,5)
plt.xticks(xlabel_ticks, [str(round(rmA[i],1)) for i in xlabel_ticks])
#plt.yticks(  [rm_s, rm_e, 0, 180])
plt.xlabel("RM")
plt.ylabel("PA(degree)")

plt.figure()
plt.plot(wavelength_meter, Q_freq_norm.mean(axis=1),'r')
plt.plot(wavelength_meter, U_freq_norm.mean(axis=1),'b')
Q_normfunc_array = Q_normfunc(wavelength_meter, np.arange(180)[minChiSquare[0]], rmA[minChiSquare[1]])
U_normfunc_array = U_normfunc(wavelength_meter, np.arange(180)[minChiSquare[0]], rmA[minChiSquare[1]])
plt.plot(wavelength_meter, Q_normfunc_array,'r')
plt.plot(wavelength_meter, U_normfunc_array,'b')

#Q_normfunc_array = Q_normfunc(wavelength_meter, 0.0, -22)
#U_normfunc_array = U_normfunc(wavelength_meter, 0.0, -22)
#plt.plot(wavelength_meter, Q_normfunc_array,'r')
#plt.plot(wavelength_meter, U_normfunc_array,'b')


plt.figure()
plt.plot(wavelength_meter, Q_freq.mean(axis=1),'r')
plt.plot(wavelength_meter, U_freq.mean(axis=1),'b')

#from scipy import optimize as so
#so.leastsq()

#plt.plot(wavelength_meter, Q_normfunc_array,'r')
#plt.plot(wavelength_meter, U_normfunc_array,'b')

plt.show()

'''
off_pulse_rms = get_off_pulse_rms(I_prof)
I_sigma = Q_sigma = U_sigma = V_sigma = I_prof / off_pulse_rms
phase = np.arange(len(I_sigma))

#plt.plot(I_sigma[I_sigma<3])
THRESHOLD  = 3
index = I_sigma>THRESHOLD
plt.plot(phase, I_sigma)
plt.plot(phase[index], I_sigma[index])
#plt.plot(Q_sigma)
#plt.plot(U_sigma)
#plt.plot(V_sigma)
plt.show()
'''
"""

def setParams():

    global pulsarData

    pulsarData = {}
    ################
    # physical constant
    pulsarData['c'] = const.c.to('m/s').value
    ################
    # set the subint start and subint end
    pulsarData['start_subint'] = 0
    pulsarData['end_endint'] = 1
    ################
    # get the pulse profile
    # startPhase, endPhase = 484, 502 # for J1917+1353-b2048.calibP.nsub10
    # startPhase, endPhase = 950, 1100 # for J1917+1353-b2048.calibP.nsub10
    pulsarData['startPhase'] = 0
    pulsarData['endPhase'] = 512
    #pulsarData['startPhase'] = 950 
    #pulsarData['endPhase'] = 1100
    ################
    # set the rebin_fact
    #rebin_fact = 4
    #rebin_fact = 128
    pulsarData['rebin_fact'] = 1

    ################
    # set the RM and PA
    pulsarData['PA0'] = 0
    pulsarData['RM'] = -66
    #pulsarData['PA0'] = 150.72
    #pulsarData['RM'] = 230.21

    print "startPhaesa: %s endPhase: %s" %(pulsarData['startPhase'], pulsarData['endPhase'])
    print "PA0 : %s  RM: %s" %(pulsarData['PA0'], pulsarData['RM'])



def main(filename):
    """Each function below is independent.
    """

    loadData(filename)

    # get the zapped channel list

    # plot the profile
    data = pulsarData['data']
    RM = pulsarData['RM']
    PA0 = pulsarData['PA0']
    freq = pulsarData['freq_hz']/1e6

    print data.shape
    #plt.plot(data[0,0,:,:].sum(0))
    #plt.imshow(data[0,0,:,:], origin='lower', aspect='auto', cmap='binary')
    #plt.show()
    #plt.plot(data[0,0,:,:].sum(0))
    #sys.exit(0)

    #plot_profile(data)

    ####
    # check derotate function
    #data_derot = derotate(data, freq, RM, PA0)
    #plot_profile(data_derot)

    #sys.exit(0)

    # RM search range 
    rm_s = -10000
    rm_e = 10000
    
    global outRootName
    outRootName = 'rmResult/%s-RM_%s_%s_subint_%s_%s_' % (filename, rm_s, rm_e, pulsarData['start_subint'], pulsarData['end_endint'])
    rm_len = rm_e - rm_s
    rmA = np.linspace(rm_s,rm_e,rm_len+1, endpoint=True)
    RMgrid(rmA)

    #data_derot = derotate(data, freq, RM, PA0)
    #plot_profile(data_derot)



    # plot the Faraday rotation measure
    #plot_FaradayRot()
    
    # rebin the data
    #plot_spectra(data, zap_Chan_Index)

    """
    #for i in range(14, 16):
    for i in range(1, 2):
        global start_subint
        global end_endint
        start_subint = i
        end_endint = i + 5
        #plot_profile(data, showFile=False)
        plot_profile(data)
        non_Chan_Index, zap_Chan_Index = get_zap_Chan(data)
        #derotate(data)
        plot_spectra(data, zap_Chan_Index)
    """

if __name__ == "__main__":
    filename = sys.argv[1]
    main(filename)
