def GWExtract(data=None, plt_graphs = 0):    
    """
    Created on Wed Feb 14 10:28:01 2018
    
    @author: wxw587, MNS543
    
    Cycles through BH masses and outputs max SNR for each. Then calculates time
    of detection.
    
    Parameters
    ----------
    
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matchedFilter import matchedFilter
    from TemplateSignalGenerator import templateGen
    
    
    if data == None:
        # load data from *.txt file, checking user input
        filename = input('Enter data file name (including \".txt\"): ')
        while 1:
            try:
                data = np.loadtxt(filename)
                break
            except FileNotFoundError:
                filename = input('File not found, please try again... ')
        
    #Defining constants
    Ts = data[1,0]-data[0,0]  # Sample spacing,
    n = 0
    mA = 50 #Initial definition of mass of BHa in solar masses
    mB = 50 #Initial definition of mass of BHb in solar masses
    
    out = []  # mA, mB, time of max and max SNR
    out_list = []  # list of mA, mB and no. of iterations
    
    
    #For-loop of mB inside mA
    while mA <= 50:
        while mB <= 50:
            
            #calculate Template            
            Template = templateGen(mA, mB, Ts)
            
            # Call Matched filter Function for the Signal-to-Noise Ratio
            SNR = matchedFilter(data,Template) 
            out = [mA, mB, SNR[0], SNR[1]] 
            print(out)
            out_list.append(out) 

            
            mB += 1
            n += 1
        mA += 1
        mB = 20
        
    
    Imax = np.argmax(np.array(out_list)[:,3]) #itration number of the maximum value SNR
    matchedTemplate = templateGen(out_list[Imax, 0], out_list[Imax, 1], Ts) 
    
    print('Closest template match:')
    print(f'mA: {out_list[Imax, 0]}')
    print(f'mB: {out_list[Imax, 1]}')    
    print(f'match time: {out_list[Imax, 2]}')
    print(f'SNR peak: {out_list[Imax, 3]}')
            
    
    if plt_graphs:
        templateGen(out_list[Imax, 0], out_list[Imax, 1], Ts, makePlots=1)
        matchedFilter(data, matchedTemplate, makePlots=1)
        
        
        
bb, ab = signal.butter(4, fband[0]*2./fs, btype='band')
#fband is the range of frequencies
#second parameter is 'array_like' 'A scalar or length-2 sequence giving the critical frequencies.'


whitened  = WhitenFunc(data, N, Pxx)


extract= signal.filtfilt(bb, ab, whitened)

# bb and ab are array_like's
# bb is the numerator coefficient vector of the filter
# ab is the denominator coefficient vector of the filter


    
  

  
if __name__ == '__main__':
    GWExtract()    
    