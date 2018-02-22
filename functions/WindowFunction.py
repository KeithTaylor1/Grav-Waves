#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 18:09:18 2018
@author: keithstaylor
"""

def Window(data):
    
    import numpy as np
    from scipy import fftpack as fftp
    import scipy as sp
    
    
    t = data[:, 0]
    
    Ts = t[1]-t[0]
    Fs = 1./Ts
    
    NFFT = 4 * Fs
    
    window = sp.signal.tukey(NFFT)
    windowfft = fftp.fft(window) / (len(window)/2.0)
    freq = np.linspace(-0.5, 0.5, len(windowfft))
    response = 20 * np.log10(np.abs(fftp.fftshift(windowfft / abs(windowfft).max())))

    windowed = data * window
    
    
    return freq, response, windowed 
