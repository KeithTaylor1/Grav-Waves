def WindowFunc(txtout, NFFT, gwtd, cwtd):

    import numpy as np
    import scipy as sp
    from scipy import signal
    
    '''Windowing of whitened signals and templates'''
    
    window1 = sp.signal.tukey(len(NFFT))
    '''sp.signal.tukey cannot be applied to any 
    data with a different # of points than itself'''
    data1 = gwtd * window1
    window1fft = fftp.fft(window1) / (len(window1)/2.0)
    freq = np.linspace(-0.5, 0.5, len(window1fft))
    response = 20 * np.log10(np.abs(fftp.fftshift(window1fft / abs(window1fft).max())))
    
    
    data2 = cwtd * window1
    



if __name__ == '__main__':
    WindowFunc(txtout=1)  

    
