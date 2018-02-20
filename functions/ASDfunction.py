#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 18:09:18 2018

@author: keithstaylor
"""

def ASDGenerator(ny, d):
    
    
    import matplotlib.mlab as mlab
    import numpy as np
    
     
    t1 = ny[:, 0]
     
    t2 = d[:, 0]
    
    Ts1 = t1[1]-t1[0]
    Fs1 = 1./Ts1
    
    Ts2 = t2[1]-t2[0]
    Fs2 = 1./Ts2

    NFFT1 = 4 * Fs1
    
    NFFT2 = 4 * Fs2
    
    #gaussian noise
    Pxx_gaussian, freqs = mlab.psd(ny, NFFT = int(NFFT1), Fs = Fs1, window=mlab.window_hanning)
    PSD_gaussian = interpolate.interp1d(freqs, Pxx_gaussian)

    #coloured noise
    Pxx_coloured, freqs = mlab.psd(d, NFFT = int(NFFT2), Fs = Fs2, window=mlab.window_hanning)
    PSD_coloured = interpolate.interp1d(freqs, Pxx_coloured)
    
    
    np.savetxt('ASDgaus.txt', np.array([freqs, np.sqrt(Pxx_gaussian)]).T)
    np.savetxt('ASDcol.txt', np.array([freqs, np.sqrt(Pxx_coloured)]).T)
    
    return np.array([freqs, np.sqrt(Pxx_gaussian)]).T, np.array([freqs, np.sqrt(Pxx_coloured)]).T
