def soundGeneration(signal=None, outputwav = None, outputshiftwav = None): #consider other parameters as inputs e.g. fs and fshift
    ''' Sound File Generation
    
    Program takes input strain signal and initially produces a wav file 
    representation of the data. The program then proceeds to shift the 
    frequency of the data such that the wav file is in a more easily audible
    range.
    
    Parameters
    ----------
    signal: np.array 
        Either the Gaussian or coloured noisy signal, formatted as a two column
        array with time in RH column and strain in LH column. If not present 
        program asks for a text file.
        
    outputwav: bool
        toggles outputwav, i.e. toggles whether to produce a wav file of the
        signal
        
    outputshiftwav: bool
        toggles outputshiftwav, i.e. toggles whether to produce a wav file of 
        the frequency shifted signal
        
    Output
    ----------
    audio: np.array written as wav file 
        input array is re-written as a wav file, frequency of the signal can be
        shifted and the corresponding array re-written as a wav file 
    '''
   
    '''Converting Numpy Array to Wav file '''

    import numpy as np
    from scipy.io import wavfile
    
    try:
        t = signal[:, 0]
        h = signal[:, 1]
    except:    
        # load signal from *.txt file, checking user input
        filename = input('Enter file name (including \".txt\"): ')
        while 1:
            try:
                signal = np.loadtxt(filename)
                break
            except FileNotFoundError:
                filename = input('File not found, please try again... ')
        
        t = signal[:, 0]
        h = signal[:, 1]
        
    
    N = len(t)               #no. of sample points
    Ts = (t[N-1]-t[0])/N     #sampling interval
    Fs = 1./Ts               #sampling frequency
    
    
    # function to keep the data within integer limits, and write to wavfile:
    def write_wavfile(filename,fs,data):
        d = np.int16(data/np.max(np.abs(data)) * 32767 * 0.9) 
        wavfile.write(filename,int(fs), d)
    
    if outputwav:
        write_wavfile("TimeDomain.wav",int(Fs), h)
    

    '''Frequency Shifting the Wav file to Audible Chirp'''
    
    # function that shifts frequency of a band-passed signal
    def reqshift(data,fshift=100,sample_rate=4096):
        """Frequency shift the signal by constant
        """
        x = np.fft.rfft(data)
        T = len(data)/float(sample_rate)
        df = 1.0/T
        nbins = int(fshift/df)
        # print T,df,nbins,x.real.shape
        y = np.roll(x.real,nbins) + 1j*np.roll(x.imag,nbins)
        y[0:nbins]=0.
        z = np.fft.irfft(y)
        return z
    
    # parameters for frequency shift
    fs = 4096 #why the re-definition?
    fshift = 400.
    #speedup = 1. - not used
    #fss = int(float(fs)*float(speedup)) - not used
    
    # shift frequency of the data
    hshift = reqshift(h,fshift=fshift,sample_rate=fs)
    
    #write the files
    if outputshiftwav:
        write_wavfile("TimeDomainShift.wav",int(fs), hshift)
    
   
    
    
    
    
    
if __name__ == '__main__':
    soundGeneration()
    
