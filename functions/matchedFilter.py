def matchedFilter(data, template=None, makePlots=False):
    """
    Created on Thu Feb 15 15:20:01 2018
    @author: Michael Stephens
    
    Program performs matched filtering using the methods described in 
    Allen et al., 2012. 
    
    
    Parameters
    ----------
    data: float64 numpy array
        The noisy data to be filtered against.
        
    template: float64 numpy array
        Template to be matched to noisy data, if not supplied as an argument 
        can be input as a text file in time domain, RH column time, LH column 
        strain.
        
    makePlots: Bool
        Toggles whether module returns maximum signal to noise ratio (SNR) or 
        creates plot of SNR in time domain.
    
    
    Returns
    -------
    SNR: tuple (float64,float64,float64)
        The maximum value and index of SNR for a given template, also d_eff
    
    
    Notes
    -----
    The template and sampling frequency MUST be the same. If template is 
    shorter than noisy data template will be zero padded.
    """
    
    import numpy as np
    from scipy import signal
    from numpy import fft, sqrt
    import matplotlib.pyplot as pl
    import matplotlib.mlab as ml
    
    # get template time and strain and additional parameters
    try:  
        t_T = template[:, 0]
        h = template[:, 1]
    except:    
        # load template from *.txt file, checking user input
        filename = input('Enter template file name (including \".txt\"): ')
        while 1:
            try:
                template = np.loadtxt(filename)
                break
            except FileNotFoundError:
                filename = input('File not found, please try again... ')
        
        t_T = template[:, 0]
        h = template[:, 1]
        
    # get data time and strain and parameters
    t_D = data[:, 0]
    d = data[:, 1]

    
    # pad template with zeros
    if len(h) < len(d):
        h = np.pad(h, (0,len(d)-len(h)), 'constant', constant_values=(0,0))
    # check sampling rates equal, else raise exception.
    if (t_T[1]-t_T[0]) != (t_D[1]-t_D[0]):
        raise Exception('template and data sampling rate not equal')
    
    
    t = t_D
    N = len(t)  # number sample points
    Ts = (t[1]-t[0])  # sample spacing
    fs = 1./Ts  # sampling frequency
    # get Template frequency bins
    xf = np.fft.rfftfreq(N, Ts)
    Ts_f = abs(xf[1] - xf[0])
    
    #PSD parameters chosen to minimise spectral leakage, need to test and implement windowing
    NFFT = int(4*fs)
    NOVL = NFFT/2
    
    psd_window = np.blackman(NFFT)
    dwindow = np.hanning(N)
    
    # PSD of the data
    Sn, freqs = ml.psd(d, Fs = fs, NFFT = NFFT, window=psd_window, noverlap=NOVL)
    #Interpolate to get the PSD values at the needed frequencies
    Sn_interp = np.interp(xf, freqs, Sn)

    # FFT of the data and the template
    df = fft.rfft(d*dwindow)/fs
    hf = fft.rfft(h*dwindow)/fs
    
    # Calculate the matched filter and output in the tmerge domain:
    inner_dh = df*np.conjugate(hf)/Sn_interp
    inner_dh_time = 2*fft.irfft(inner_dh)*fs

    # Normalize the matched filter output:
    inner_hh = sum((hf*np.conjugate(hf)/Sn_interp))*Ts_f
    inner_hh = sqrt(abs(inner_hh))
    SNR = inner_dh_time / inner_hh

    # shift the SNR vector by the template length so that the peak is at the END of the template
    peaksample = np.argmax(abs(h))  # location of peak in the template
    SNR = np.roll(SNR, peaksample)
    SNR = abs(SNR)
    
    # effective "distance"
    d_eff = inner_hh/max(SNR)
    
    if makePlots:
        pl.figure()
        pl.plot(t, SNR, label='SNR')
        pl.title(f'Mean SNR: {np.mean(SNR)}')
        pl.grid()
        pl.legend(loc='best')
        pl.xlabel('time')
        pl.ylabel('strain')    
    
    return (np.argmax(SNR), max(SNR), d_eff)
    
    
    
    
    



import numpy as np

if __name__ == '__main__':
    
    filename = input('Enter data file name (including \".txt\"): ')
    while 1:
        try:
            data = np.loadtxt(filename)
            break
        except FileNotFoundError:
            filename = input('File not found, please try again... ')
            
    matchedFilter(data, makePlots=1)