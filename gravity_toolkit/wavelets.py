#!/usr/bin/env python
u"""
wavelets.py
Written by Hugo Lecomte (02/2021)

Function to apply a wavelets analysis, code based on (Torrence and Compo, 1998)
"""
import numpy as np
import scipy.special

def wave_bases(mother, k, scale, param=-1):
    """Computes the wavelet function as a function of Fourier frequency
       used for the CWT in Fourier space (Torrence and Compo, 1998)

    Arguments
    ---------
        mother:  str equal to 'MORLET' or 'DOG' to choose the wavelet type
        k: vector of the Fourier frequencies
        scale:  wavelet scales
        param: nondimensional parameter for the wavelet function

    Returns
    -------
        daughter: the wavelet function
        fourier_factor: the ratio of Fourier period to scale
        coi: cone-of-influence size at the scale
        dofmin: degrees of freedom for each point in the wavelet power (Morlet = 2)
    """
    mother = mother.upper()
    n = len(k)  # length of Fourier frequencies
    k = np.array(k) # turn k to array

    if mother == 'MORLET':  # choose the wavelet function, in this case Morlet
        if param == -1:
            param = 6 # For Morlet this is k0 (wavenumber), default is 6

        expnt = -(scale*k - param)**2/2*(k > 0) # table 1 Torrence and Compo (1998)
        norm = np.sqrt(scale*k[1])*(np.pi** -0.25)*np.sqrt(len(k))

        daughter = [] # define daughter as a list
        for ex in expnt:  # for each value scale (equal to next pow of 2)
            daughter.append(norm*np.exp(ex))
        daughter = np.array(daughter) # transform in array

        daughter = daughter*(k > 0) # Heaviside step function
        fourier_factor = (4*np.pi)/(param + np.sqrt(2 + param * param)) # scale --> Fourier period
        coi = fourier_factor/np.sqrt(2) # cone-of-influence
        dofmin = 2 # degrees of freedom

    elif mother == 'DOG': # DOG Wavelet
        if param == -1:
            param = 2  # For DOG this is m (wavenumber), default is 2
        m = param

        expnt = -(scale*k)**2/2
        pws = np.array((scale*k)**m)
        # gamma(m+0.5) = 1.3293
        norm = np.sqrt(scale*k[1]/1.3293*np.sqrt(n))

        daughter = []
        for ex in expnt:
            daughter.append(-norm* 1j**m * np.exp(ex))
        daughter = np.array(daughter)
        daughter = daughter[:]*pws

        fourier_factor = 2*np.pi/np.sqrt(m + .5)
        coi = fourier_factor/np.sqrt(2)
        dofmin = 1

    elif mother == 'PAUL':  # Paul Wavelet
        if param == -1:
            param = 4
        m = param

        expnt = -(scale*k)*(k > 0)
        norm = np.sqrt(scale*k[1]) *(2**m /np.sqrt(m*(np.math.factorial(2*m - 1))))*np.sqrt(n)
        pws = np.array((scale*k)**m)

        daughter = []
        for ex in expnt:
            daughter.append(norm*np.exp(ex))
        daughter = np.array(daughter)
        daughter = daughter[:]*pws

        daughter = daughter*(k > 0)  # Heaviside step function
        fourier_factor = 4*np.pi/(2*m + 1)
        coi = fourier_factor*np.sqrt(2)
        dofmin = 2

    return daughter, fourier_factor, coi, dofmin


def wavelet(Y, dt, pad=1, dj=.25, s0=-1, J1=-1, mother='MORLET', param=-1):
    """Computes the wavelet continuous transform of the vector Y,
           by definition:
        W(a,b) = sum(f(t)*psi[a,b](t) dt)        a dilate/contract
        psi[a,b](t) = 1/sqrt(a) psi(t-b/a)       b displace
        The wavelet basis is normalized to have total energy = 1 at all scales

    Arguments
    ---------
        Y: time series
        dt: sampling rate
        pad: bool for zero padding or not
        dj: spacing between discrete scales
        s0: smallest scale of the wavelet
        J1: total number of scales
        mother: the mother wavelet function
        param: the mother wavelet parameter

    Returns
    -------
        wave: wavelet transform of Y
        period: the vector of "Fourier" periods (in time units) that correspond to the scales
        scale: vector of scale indices, given by S0*2(j*DJ), j =0 ...J1
        coi: cone of influence
        """
    n1 = len(Y)  # time series length

    if s0 == -1: # define s0 as 2 times dt (Shannon criteria) if s0 is not given
        s0 = 2 * dt
    if J1 == -1: # define J1 if not provide
        J1 = int((np.log(n1*dt/s0) / np.log(2))/dj)

    x = Y - np.mean(Y) # remove mean of the time serie

    if pad: # if zero padding, add zeros to x
        base2 = int(np.log(n1)/np.log(2) + 0.4999)
        x = np.concatenate((x, np.zeros(2**(base2 + 1) - n1)))

    n = len(x) #update length of x

    k = np.arange(0, int(n/2))
    k = k*(2*np.pi) / (n*dt)
    k = np.concatenate((k, -k[int((n - 1)/2)::-1]))  # be careful for parity

    f = np.fft.fft(x) # fft on the padded time series

    scale = s0 * 2**(np.arange(0, J1 + 1, 1)*dj)
    # define wavelet array
    wave = np.zeros((int(J1 + 1), n))
    wave = wave + 1j * wave # make it complex

    for a1 in range(0, int(J1 + 1)):
        daughter, fourier_factor, coi, dofmin = wave_bases(mother, k, scale[a1], param)
        wave[a1, :] = np.fft.ifft(f * daughter)

    period = fourier_factor * scale

    # cone-of-influence, differ for uneven len of timeseries:
    if n1%2: # uneven
        coi = coi * dt * np.concatenate((np.arange(0, n1/2 - 1), np.arange(0, n1/2)[::-1]))
    else: # even
        coi = coi * dt * np.concatenate((np.arange(0, n1/2), np.arange(0, n1/2)[::-1]))

    # cut zero padding
    wave = wave[:, :n1]

    return wave, period, scale, coi

def wave_signif(Y, dt, scale, dof=-1, lag1=0, siglvl=0.95, mother='MORLET', param=-1):
    """Computes the wavelet significance test at a level of confidence siglvl%

        Arguments
        ---------
            Y: time series
            dt: sampling rate
            scale: scales of the wavelet decomposition
            dof: degrees of freedom
            lag1: assuming lag-1 autocorrelation of the serie (0 for white noise RECOMMENDED, 0.72 for red noise)
            siglvl: percentage of the confidence level
            mother: the mother wavelet function
            param: the mother wavelet parameter

        Returns
        -------
            wave: wavelet transform of Y
            period: the vector of "Fourier" periods (in time units) that correspond to the scales
            scale: vector of scale indices, given by S0*2(j*DJ), j =0 ...J1
            coi: cone of influence
    """
    mother = mother.upper()
    variance = np.var(Y)

    # define default param and fourier factor for the wavelet
    if mother == 'MORLET':
        if param == -1:
            param = 6 # For Morlet this is k0 (wavenumber), default is 6
        if dof == -1:
            dof = 2

        fourier_factor = float(4 * np.pi) / (param + np.sqrt(2 + param**2))

    if mother == 'DOG':
        if param == -1:
            param = 2  # For DOG, default param is 2
        if dof == -1:
            dof = 1

        fourier_factor = float(2 * np.pi / (np.sqrt(param + 0.5)))

    if mother == 'PAUL':
        if param == -1:
            param = 4  # For PAUL, default param is 4
        if dof == -1:
            dof = 2

        fourier_factor = float(4 * np.pi / (2 * param + 1))

    # compute period from scale
    period = [e * fourier_factor for e in scale]

    # compute theoretical fft associated to the theoretical noise of the data given by lag1
    freq = [dt / p for p in period]
    fft_theor = [variance*((1 - lag1**2) / (1 - 2*lag1*np.cos(f * 2 * np.pi) + lag1**2)) for f in freq]

    chisquare = scipy.special.gammaincinv(dof/2.0, siglvl)*2.0/dof
    signif = [ft * chisquare for ft in fft_theor]
    return signif