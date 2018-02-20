# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:53:05 2018

@author: njg573
"""

def WhitenFunc(ny, N, Pxx_gaussian, d, Pxx_coloured):
    
    import numpy as np
    from scipy import fftpack as fftp
    
    '''Whitening  - transform to freq domain, divide by asd, transform back
   also plotting the whitened data in the time domain'''

    #gaussian
    #convert to frequency domain
    gfd2 = fftp.fft(ny) #2 sided fft
    gfd1 = gfd2[0:(N//2)+1] #1 sided fft
    
    #divide by sqrt of psd to whiten
    gw = gfd1/ np.sqrt(Pxx_gaussian)
    
    #convert back to time domain
    gwtd = np.fft.irfft(gw)
    
    
    #coloured
    cfd2 = fftp.fft(d) #2 sided fft
    cfd1 = cfd2[0:(N//2)+1] #1 sided fft
    
    #divide by sqrt of psd to whiten
    cw = cfd1/ np.sqrt(Pxx_coloured)
    #convert back to time domain
    cwtd = np.fft.irfft(cw)
    
    return gfd2, gfd1, gw, gwtd, cfd2, cfd1, cw, cwtd
