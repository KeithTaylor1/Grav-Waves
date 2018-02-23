def GWExtract(data=None, makePlots = 0):    
    """
    Created on Wed Feb 14 10:28:01 2018
    
    @author: wxw587, MNS543, Keith
    
    Cycles through BH masses and outputs max SNR for each. Then calculates time
    of detection and extracts signal.
    
    
    Parameters
    ----------
    
    data: float64 numpy array
        noisy LIGO data in two columns, LH time RH strain.

    makePlots: bool
        Toggles output of graphs of SNR for matched template, 
    
    
    Notes
    -----
    Maximum and minimum masses to be taken from merger rate subgroup.
    
    Requires matchedFilter and teplateGen and WhitenFunc modules
    """
    
    import numpy as np
    from scipy import signal
    import matplotlib.pyplot as pl
    import matplotlib.mlab as ml
    from matchedFilter import matchedFilter
    from TemplateSignalGenerator import templateGen
    from WhitenFunc import WhitenFunc
    
    
    if data == None:
        # load data from *.txt file, checking user input
        filename = input('Enter data file name (including \".txt\"): ')
        while 1:
            try:
                data = np.loadtxt(filename)
                break
            except FileNotFoundError:
                filename = input('File not found, please try again... ')
        
    # Defining constants
    Ts = data[1,0]-data[0,0]  # Sample spacing
    M_min = 45  # minimum mass in solar masses
    M_max = 56  # maximum mass in solar masses
    
    out = []  # mA, mB, time of max and max SNR
    out_list = np.array([[0,0,0,0]])  # list of mA, mB and no. of iterations
    
    
    #For-loop of mB inside mA
    massesA = range(M_min, M_max)  # range of mass of BHa in solar masses
    for mA in massesA:
        massesB = range(M_min, mA+1)  # range of mass of BHb in solar masses, avoid double counting to reduce runtime
        for mB in massesB:
            
            #calculate Template            
            Template = templateGen(mA, mB, Ts)
            
            # Call Matched filter Function for the Signal-to-Noise Ratio
            SNR = matchedFilter(data,Template) 
            out = [mA, mB, SNR[0], SNR[1]]
            print(out)
            out_list = np.append(out_list, [out],axis=0)
        

    
    Imax = np.argmax(np.array(out_list)[:, 3])  # itration number of the maximum value SNR
    matchedTemplate = templateGen(out_list[Imax, 0], out_list[Imax, 1], Ts) 
    
    print('\nClosest template match:')
    print(f'mA: {out_list[Imax, 0]}')
    print(f'mB: {out_list[Imax, 1]}')    
    print(f'match time: {out_list[Imax, 2]}')
    print(f'SNR peak: {out_list[Imax, 3]}\n')
            
    
    if makePlots:
        
        SNR = matchedFilter(data, matchedTemplate, makePlots=1)
        templateGen(out_list[Imax, 0], out_list[Imax, 1], Ts, makePlots=1)
        
        # PSD of the data
        Sn, freqs = ml.psd(data[:,1], Fs = 1./Ts, NFFT = int(4./Ts), window=np.blackman(int(4./Ts)), noverlap=int(2./Ts))
        #Interpolate to get the PSD values at the needed frequencies
        Sn_interp = np.interp(np.fft.rfftfreq(len(data[:,0]),Ts), freqs, Sn)
        # LIGO frequency range 43-800Hz
        # bb and ab are array_like's
        # bb is the numerator coefficient vector of the filter
        # ab is the denominator coefficient vector of the filter
        bb, ab = signal.butter(4, [43*2.*Ts, 800*2.*Ts], btype='band')
        whitened = WhitenFunc(data[:,1], len(data[:,0]), Sn_interp)
        # filter to get signal
        extract= signal.filtfilt(bb, ab, whitened)
        
        pl.figure()
        pl.plot(data[:, 0],extract, label='extracted signal')
        pl.legend(loc='best')
        pl.grid()
    
  

  
if __name__ == '__main__':
    GWExtract(makePlots=1)    
    
