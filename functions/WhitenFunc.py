# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:53:05 2018
@author: njg573
"""

def WhitenFunc(data, N, Pxx):
    
    import numpy as np
    from scipy import fftpack as fftp
    
    '''Whitening  - transform to freq domain, divide by asd, transform back
   also plotting the whitened data in the time domain'''

    
    #convert to frequency domain
    fft2 = fftp.fft(data) #2 sided fft
    fft1 = fft2[0:(N//2)+1] #1 sided fft
    
    #divide by sqrt of psd to whiten
    whitened_freq= fft1/ np.sqrt(Pxx)
    
    #convert back to time domain
    whitened_time = np.fft.irfft(whitened_freq)
    
    
   
    
    return whitened_freq, whitened_time
