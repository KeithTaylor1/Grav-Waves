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
        Toggles output of graphs of SNR for matched template
        
        
    Returns
    -------
    
    list of: mA, mB, merge time and SNR peak value
    
    
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
    M_max = 58  # maximum mass in solar masses
    
    out = []  # mA, mB, time of max and max SNR
    out_list = np.array([[0,0,0,0,0]])  # list of mA, mB and no. of iterations
    
    
    #For-loop of mB inside mA
    massesA = range(M_min, M_max+1)  # range of mass of BHa in solar masses
    for mA in massesA:
        massesB = range(M_min, mA+1)  # range of mass of BHb in solar masses, avoid double counting to reduce runtime
        for mB in massesB:
            
            #calculate Template            
            Template = templateGen(mA, mB, Ts)
            
            # Call Matched filter Function for the Signal-to-Noise Ratio
            SNR = matchedFilter(data,Template) 
            out = [mA, mB, SNR[0], SNR[1], SNR[2]]
            print(out)
            out_list = np.append(out_list, [out],axis=0)
        

    
    Imax = np.argmax(np.array(out_list)[:, 3])  # itration number of the maximum value SNR
    matchedTemplate = templateGen(out_list[Imax, 0], out_list[Imax, 1], Ts) 
    
    outs = [out_list[Imax, 0], out_list[Imax, 1], data[int(out_list[Imax, 2]), 0], out_list[Imax, 3], out_list[Imax, 4]]
    print('\nClosest template match:')
    print(f'mA: {outs[0]}')
    print(f'mB: {outs[1]}')    
    print(f'time of merge: {outs[2]}')
    print(f'SNR peak: {outs[3]}')
    print(f'd_eff: {outs[4]}\n')
         
    
    if makePlots:
        
        SNR = matchedFilter(data, matchedTemplate, makePlots=1)
        
        # plot whitened signal and template
        #Interpolated PSD 
        Sn, freqs = ml.psd(data[:,1], Fs = 1./Ts, NFFT = int(4./Ts), window=np.blackman(int(4./Ts)), noverlap=int(2./Ts))
        Sn_interp = np.interp(np.fft.rfftfreq(len(data[:,0]),Ts), freqs, Sn)
        # LIGO frequency range 43-400Hz 
        # bb and ab are array_like
        # bb is the numerator coefficient vector of the filter
        # ab is the denominator coefficient vector of the filter
        bb, ab = signal.butter(4, [43*2.*Ts, 400*2.*Ts], btype='band')
        whitened = WhitenFunc(data[:,1], Sn_interp, Ts)
        # filter to get signal and normalise
        extract= signal.filtfilt(bb, ab, whitened)/np.sqrt((400.-43.)*Ts*2)
        
        matchedTemplate = np.pad(matchedTemplate[:,1], (0,len(data[:,1])-len(matchedTemplate[:,1])), 'constant', constant_values=(0,0))
        matchedTemplate = np.roll(matchedTemplate,int(out_list[Imax, 2]-np.argmax(abs(matchedTemplate))))/outs[4]
        
        pl.figure()
        pl.plot(data[:, 0],extract, label='extracted signal')
        pl.plot(data[:,0],matchedTemplate, label='matched template')
        pl.xlim([outs[2]-1, outs[2]+1])
        pl.ylim([-10, 10])
        pl.legend(loc='best')
        pl.grid()
        pl.show()
    
    
    return outs

  
if __name__ == '__main__':
    GWExtract(makePlots=1)    
    
