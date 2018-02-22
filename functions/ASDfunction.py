#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 18:09:18 2018
@author: keithstaylor
"""

def ASDGenerator(data):
    
    
    import matplotlib.mlab as mlab
    import numpy as np
    from scipy import interpolate
    
     
    t = data[:, 0]
     
    Ts = t[1]-t[0]
    Fs = 1./Ts

    NFFT = 4 * Fs
    

    

    Pxx, freqs = mlab.psd(data, NFFT = int(NFFT), Fs = Fs, window=mlab.window_hanning)
    PSD = interpolate.interp1d(freqs, Pxx)

    
    
    
    np.savetxt('ASDgaus.txt', np.array([freqs, np.sqrt(Pxx)]).T)
    np.savetxt('ASDcol.txt', np.array([freqs, np.sqrt(Pxx)]).T)
    
    return Pxx, PSD
