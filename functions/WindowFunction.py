#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 18:09:18 2018

@author: keithstaylor
"""

def Window(gwtd, cwtd):
    
    import numpy as np
    from scipy import fftpack as fftp
    import scipy as sp
    
    
    t = gwtd[:, 0]
    
    Ts = t[1]-t[0]
    Fs = 1./Ts
    
    NFFT = 4 * Fs
    
    window1 = sp.signal.tukey(NFFT)
    window1fft = fftp.fft(window1) / (len(window1)/2.0)
    freq = np.linspace(-0.5, 0.5, len(window1fft))
    response = 20 * np.log10(np.abs(fftp.fftshift(window1fft / abs(window1fft).max())))

    data1 = gwtd * window1
    data2 = cwtd * window1
    
    return freq, response, data1, data2 
