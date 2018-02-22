def templateGen(mA, mB, Ts, zeroPad=False, makePlots=False, savePlots=False, txtout=False):
    '''TemplateGen
    Author: Emma
    
    Produces model of GW from binary black hole merger at a distance of 3e24m
    
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
        be used before colouredNoiseGenerator
    
    makePlots: bool
        output plot of signal
        
    txtout: bool
        output text file if true, else return signal
        
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
    r = 3.e24 #m
    t = 0.0 #s
    
    mA = 50. * SM
    mB = 50. * SM
    
    T=[]
    f=[]
    
    #Minimum separation
    a0 = (G * (mA + mB) / (pi*f_GW)**2) ** (1./3.)
    #print a0
    
    #Merge Time
    Tmerge = (5 * c**5 * a0**4)/(256 * G**3 * mA * mB * (mA +mB))
    #print Tmerge
    
    #Separation as a function of time
    def a(mA, mB, G, c, t):
        return ( (a0)**4 - (256 * G**3 * mA * mB * (mA+mB) * (t))/(5 * c**5) )**(1./4.)
    
    #Calculates amplitude with increasing time
    
    while t <= Tmerge:
        y = ( (G * (mA+mB) * a(mA, mB, G, c, t)**(-3.)) / (pi**2) ) ** (1./2.)
        z = ( (G**2 * mA * mB) * (a(mA, mB, G, c, t)**(-1)) ) / (c**4 * r)
        s = z * (sin(2*pi*t*y))
        T.append(t)
        f.append(s)
        t += Ts
    #print t
    
    #print len(T)
    
    #Zero pad data if noise going to be added
    if zeroPad:
        while t <= 32:
            T.append(t)
            t += Ts
        
        f = np.pad(f, (int((32-Tmerge)*np.random.rand()/Ts),0), 'constant', constant_values=(0,0))
        f = np.pad(f, (0,len(T)-len(f)), 'constant', constant_values=(0,0))
    
    if txtout:
        file = open(f'Template{int(mA/SM)}_{int(mB/SM)}.txt',"w")
        #file.write("Time/s" + " " + "Strain" + "\n")
        for i in range(len(T)):
            file.write(str(T[i]) + " " + str(f[i]) + "\n")
        file.close
    
    if makePlots:    
        plt.plot(T, f)
        #plt.xscale('log')
        plt.xlabel('Time/s')
        plt.ylabel('Strain')
        if savePlots: plt.savefig(f'Template{int(mA/SM)}_{int(mB/SM)}.jpg')
        plt.grid
        plt.show()
        
    return np.array([T,f]).T







if __name__ == '__main__':
    templateGen(50, 50, 1/4096, zeroPad=1, makePlots=1, txtout=1)