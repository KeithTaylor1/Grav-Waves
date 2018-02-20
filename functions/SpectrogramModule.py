def spectrogram(signal=None, makeplots=True):
    ''' Spectrogram Generator
    
    Program takes input strain signal and outputs figure of a spectrogram of 
    the data i.e. allowing for a visual depiction of the frequency content of 
    the data.
    
    Parameters
    ----------
    signal: np.array 
        Either the Gaussian or coloured noisy signal, formatted as a two column
        array with time in RH column and strain in LH column. If not present 
        program asks for a text file.
    
    makeplots: bool
        Toggles makeplots, default set to true i.e. to produce plots of data.
        
    Output
    ----------
    Spectrogram: im (perhaps use all the outputs of the function)
        The matplotlib.pyplot function spectrogram outputs: 
        a 2D array whereby the columns are the periodgrams of successive 
        segaments, 1D array of the frequencies corresponding to the rows in the
        spectrum, 1D array of the times corresponding to the midpoints of the 
        segments and the image created by imshow containing the spectrogram.
        This spectrogram is the output we are most concerned with but
        program allows option for the other outputs to be viewed / analysed.
        
    '''
    '''Spectrograms - spectrograms of the gaussian and coloured noisy data'''
    
  
    
    import numpy as np
    import matplotlib.pyplot as pl
    
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
     
    NFFT1 = int(Fs/8)
    NOVL = int(NFFT1*15./16)
    window = np.blackman(NFFT1)
    spec_cmap='ocean'
    
    if makeplots:
        pl.figure()
        Pxx, freqs, bins, im = pl.specgram(h, NFFT=NFFT1, Fs=Fs, window=window, noverlap=NOVL, cmap=spec_cmap)
        pl.colorbar()
        pl.show()
        
    #need options for displaying Pxx, freqs and bins
    
    
if __name__ == '__main__':
    spectrogram()   
