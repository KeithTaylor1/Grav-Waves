def WhitenFunc(data, Sn_interp, Ts):
    '''
    Created on Tue Feb 20 15:53:05 2018
    @author: njg573, MNS543
    
    Whitening  - transform to freq domain, divide by asd, transform back
    also plotting the whitened data in the time domain
    
    
    Parameters
    ----------
    
    data: float64 numpy array
        LIGO strain data in time domain
        
    Sn_interp: float64 numpy array
        interpolated PSD of data, should have same length
    
    Ts: float64
        sample spacing of data
        
        
    Returns
    -------
    
    White_d: float64 numpy array
        whitened strain data
        
    '''
    
    import numpy as np
    df = np.fft.rfft(data)
    white_df = df / np.sqrt(Sn_interp) * np.sqrt(Ts*2)
    white_d = np.fft.irfft(white_df)
    return white_d
