def templateGen(mA, mB, Ts, zeroPad=False, makePlots=False, savePlots=False, txtout=False):
    '''TemplateGen
    Author: Emma, Michael
    
    Produces model of GW from binary black hole merger at a distance of 3e24m
    signal includes inspiral and ringdown.
    
    Parameters
    ----------
    
    mA: double
        Mass of BH A
        
    mB: double
        Mass of BH B
        
    Ts: double
        Sample spacing
        
    zeroPad: bool
        pads with zeros up to 32 sec, randomly locates signal in this range, to
        be used before colouredNoiseGenerator, defaults to false
    
    makePlots: bool
        toggles output plot of signal, defaults to false
        
    txtout: bool
        toggle output text file, defaults to false
        
    Return
    ------
    
    signal: float64 numpy array
        two column array, LH column is time, RH column is strain
    
    '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    from numpy import pi, sin
    
    
    #Defining constants
    c = 3.0e8 #ms-1
    G = 6.667e-11 #m3 kg-1 s-2
    SM = 2.e30 #kg
    f_GW = 20 #Hz
    r = 8.e24 #m
    t = 0.0 #s
    
    mA *= SM
    mB *= SM
    # complex ringdown frequency for a BH spin 0.7
    Fx = 0.5326
    Fy = 0.0808
 
    T=[]
    S=[]
    
    #Minimum separation
    a0 = (G * (mA + mB) / (pi*f_GW)**2) ** (1./3.)
    
    #Cut-off time
    Tcutoff = ((5 * c**5) / (256 * G**3 * mA * mB * (mA + mB))) * (a0**4 - (((6 * G * (mA + mB)) / (c**2)) ** 4))
    
    #Separation as a function of time
    def a(mA, mB, G, c, t):
        return ( (a0)**4 - (256 * G**3 * mA * mB * (mA+mB) * (t))/(5 * c**5) )**(1./4.)
    
    #Calculates amplitude with increasing time
    
    while t <= Tcutoff:
        f = ( (G * (mA+mB) * a(mA, mB, G, c, t)**(-3.)) / (pi**2) ) ** (1./2.)
        h = ( (G**2 * mA * mB) * (a(mA, mB, G, c, t)**(-1)) ) / (c**4 * r)
        s = h * (sin(2*pi*t*f))
        T.append(t)
        S.append(s)
        t += Ts
    #print t
    
    N = len(T) - 1
    fstitch = S[N]
    t=0.0
    
    #Ringdown
    f_GW_RD = (c**3 * Fx * 1j) / (G * (mA + mB))
    damp = (G * (mA + mB)) / (c**3 * Fy)
    
    while t <= 3*damp:
        e1 = f_GW_RD * t
        e2 = -t / damp
        z = np.real(fstitch * np.exp(e1) * np.exp(e2))
        ts = t + Tcutoff
        t += Ts
        T.append(ts)
        S.append(z)
    
    
    
    #Zero pad data if noise going to be added
    if zeroPad:
        while t <= 32:
            T.append(t)
            t += Ts
        
        S = np.pad(S, (int((32-T[-1])*np.random.rand()/Ts),0), 'constant', constant_values=(0,0))
        S = np.pad(S, (0,len(T)-len(S)), 'constant', constant_values=(0,0))
    
    if txtout:
        np.savetxt(f'Template{int(round(mA/SM))}_{int(round(mB/SM))}.txt', np.array([T, S]).T)

    if makePlots:
        plt.figure()
        plt.plot(T, S, label='Template')
        #plt.xscale('log')
        plt.xlabel('Time/s')
        plt.ylabel('Strain')
        if savePlots: plt.savefig(f'Template{int(round(mA/SM))}_{int(round(mB/SM))}.jpg')
        plt.legend(loc='best')
        plt.grid()
        plt.show()
        
    return np.array([T,S]).T







if __name__ == '__main__':
    templateGen(56., 30., 1/4096, zeroPad=0, makePlots=1, txtout=0)
