def gaussianNoise(signal=None, txtout=False):
    '''Gaussian-noise-generator
   
    Program applies random gaussian noise to input signal. From numpy, the 
    random.normal function allows for random samples to be drawn from a Normal
    distribution.
    
    
    Parameters
    ----------
    signal: np.array 
        The clean signal that noise is to be added to, formatted as a two
        column array with time in RH column and strain in LH column. If not 
        present program asks for a text file.
    
    txtout: bool
        Toggles text file output, defaults to false.
    
    Returns
    -------
    noisyData: np.array
        Gaussian Noise applied to input data / signal.
        
    '''
    
    import numpy as np
    
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
    
    N = len(h)
    m = 0
    sd = 10
    
    noise = np.random.normal(m,sd,N) #Gaussian Noise 
    ny = h + noise #noisy signal 
    
    if txtout:
        np.savetxt('noisyData-col.txt', np.array([t, ny]).T)
    else:
        return np.array([t, ny]).T
        
    
    
    
    
if __name__ == '__main__':
    gaussianNoise(txtout=1)   
