import numpy as np
import scipy as sp
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as pl
import matplotlib.mlab as ml
from scipy import fftpack as fftp

"""
variable suffix meaning: f_1 frequency domain 1 sided (only +ve values)
                         f_2 frequency domain 2 sided (+ve & -ve values)
                         f_2s frequency domain 2 sided but with f=0 centred
"""


'''Q1a'''

data1 = np.loadtxt('sinusoid.txt')

t = data1[:,0]
h = data1[:,1]


"""plot data"""
if(0):
    fig1 = pl.figure()
    
    ax1 = fig1.add_subplot(2,1,1)
    ax1.plot(t,h)
    ax1.set_xlabel('time, s')
    ax1.set_ylabel('amplitude')
    ax1.grid()
    
    ax2 = fig1.add_subplot(2,1,2)
    ax2.plot(t,h)
    ax2.set_xlim((0,0.2))
    ax2.set_xlabel('time, s')
    ax2.set_ylabel('amplitude')
    ax2.grid()
    
    pl.tight_layout()
    pl.show()



N = len(t) #number sample points
Ts = (t[len(t)-1]-t[0])/N #sample spacing
fs = 1./Ts #sampling frequence

'''Default frequency rannge between minus half the sampling frequency and half 
the sampling frequency. Number of bins for fft equal to number of sampling 
points, only looking at +ve frequencies so only half the number of points'''

hf_2 = fftp.fft(h) #2 sided fft
hf_1 = hf_2[0:N//2]#1 sided fft

xf_2 = fftp.fftfreq(N,Ts)#gets 2 sided frequencies
xf_1 = xf_2[0:N//2]

'''############## maybe plot full fft? ############'''

'''plot data in frequency domain'''
if (0):
    pl.figure()
    pl.plot(xf_1, 2./N*np.abs(hf_1))
    pl.xlabel('frequency, Hz')
    pl.ylabel('amplitude')
    pl.grid()
    pl.show()


f_spike = xf_1[np.argmax(2./N*np.abs(hf_1))]
Amp = max(2./N*np.abs(hf_1))
print(f'\nPlain sinusoidal signal:\nThere is a spike with amplitude {Amp:.2f} at frequency {f_spike:.2f}Hz.')

'''Parseval's theorem'''
print(f'Power in the time domain is {sum(abs(h)**2)} and power in the  frequency domain is {1/N*sum(abs(hf_2)**2)}.')


#%%

'''Q1b'''

data2 = np.loadtxt('noisy_sinusoid.txt')

t = data2[:,0]
d = data2[:,1]


"""plot data"""
if(0):
    fig1 = pl.figure()
    
    ax1 = fig1.add_subplot(2,1,1)
    ax1.plot(t,d)
    ax1.set_xlabel('time, s')
    ax1.set_ylabel('amplitude')
    ax1.grid()
    
    ax2 = fig1.add_subplot(2,1,2)
    ax2.plot(t,d)
    ax2.set_xlim((0,0.2))
    ax2.set_xlabel('time, s')
    ax2.set_ylabel('amplitude')
    ax2.grid()
    
    pl.tight_layout()
    pl.show()



N = len(t) #number sample points
Ts = (t[len(t)-1]-t[0])/N #sample spacing
fs = 1./Ts #sampling frequence


df_2 = fftp.fft(d) #2 sided fft
df_1 = df_2[0:N//2] #1 sided fft

xf_2 = fftp.fftfreq(N,Ts) #gets 2 sided frequency bins
xf_1 = xf_2[0:N//2]


'''plot data in frequency domain'''
if(0):
    pl.figure()
    pl.plot(xf_1, 2./N*abs(df_1))
    pl.xlabel('frequency, Hz')
    pl.ylabel('amplitude')
    pl.grid()
    pl.show()


f_spike = xf_1[np.argmax(2./N*np.abs(df_1))]
Amp = max(2./N*np.abs(df_1))
print(f'\nNoisy sinusoidal signal:\nThere is a spike with amplitude {Amp:.2f} at frequency {f_spike:.2f}Hz.')

#%%

'''Q1c'''

data2 = np.loadtxt('noisy_sinusoid.txt')

t = data2[:,0]
d = data2[:,1]

window = np.blackman(len(d)) #blackman window, reduces 

''' put data and template in frequency domain'''

N = len(t) #number sample points
Ts = (t[len(t)-1]-t[0])/N #sample spacing
fs = 1./Ts #sampling frequency

NFFT = int(4*fs)

'''data'''
df_2 = fftp.fft(d)/fs #2 sided fft
df_1 = df_2[0:N//2+1] #1 sided fft
#df_2s = fftp.fftshift(df_2) #centre f=0

#xf_2 = fftp.fftfreq(N,Ts) #gets 2 sided frequency bins
#xf_1 = xf_2[0:N//2+1]
#xf_2s = fftp.fftshift(xf_2) #centre f=0

"""
'''noise'''
nf_1 = df_1.copy()
nf_1[np.argmax(nf_1)] = 400*np.random.rand()
n = np.fft.irfft(nf_1) #data minus the signal to give just noise
"""

'''Noise power spectral density'''
#psd, freq = ml.psd(n, Fs=fs, NFFT=NFFT, window=np.blackman(NFFT), noverlap=NFFT/2) #PSD, args taken from LOSC to minimise spectral leakage
Sn_1 = 1 ####only ok because noise is white#####
#np.interp(np.abs(xf_1), freq, psd) #PSD, 1 sided with correct no of points
'''index = abs(np.arange(-N//2,N//2)) #calculate 2 sided PSD with f=0 at centre
Sn_2s = list()                       #and |f|
for i in index:
    
    Sn_2s.append(Sn_1[i])
    
Sn_2s = np.array(Sn_2s)
'''

'''limits:'''
f_min = 40
f_max = 50
N_points = 1000

f = np.linspace(f_min,f_max,N_points)

SNR = []
for i in f:

    '''template'''
    h = sin(2.*pi*i*t)
    hf_2 = fftp.fft(h)/fs
    hf_1 = hf_2[0:N//2+1]
    #hf_2s = fftp.fftshift(hf_2) #centre f=0
    
    '''inner products'''
    #<d|h>    
    inner_dh = 4*np.real(sum(Ts*df_1*np.conj(hf_1)/Sn_1))
    #<h|h>
    inner_hh = 4*np.real(sum(Ts*hf_1*np.conj(hf_1)/Sn_1))
    
    '''SNR'''
    SNR.append(inner_dh/sqrt(inner_hh))

SNR = abs(np.array(SNR))

f_spike = f[np.argmax(SNR)]
print(f'\nMatched filtering of noisy data:\nThere is a spike at frequency {f_spike:.2f}Hz.')


'''plot SNR'''
if(0):
    pl.figure()
    pl.plot(f, SNR)
    pl.xlabel('frequency, Hz')
    pl.ylabel('SNR')
    pl.grid()
    pl.show()

#%%
    
'''Q1c simple: in time domain'''

data2 = np.loadtxt('noisy_sinusoid.txt')

t = data2[:,0]
d = data2[:,1]

'''limits:'''
f_min = 30
f_max = 60
N_points = 4000

f = np.linspace(f_min,f_max,N_points)

SNR = []
for i in f:

    '''template'''
    h = sin(2.*pi*i*t)
    
    '''inner products, time domain'''
    #<d|h>    
    inner_dh = np.inner(d, h)
    #<h|h>
    inner_hh = np.inner(h, h)
    
    '''SNR'''
    SNR.append(inner_dh/sqrt(inner_hh))

SNR = abs(np.array(SNR))



'''plot SNR'''
if(0):
    pl.figure()
    pl.plot(f, SNR)
    pl.xlabel('frequency, Hz')
    pl.ylabel('SNR')
    pl.grid()
    pl.show()
    
    
#%%
    
'''Q1d'''


t = np.linspace(0, 4*pi, 2000)
h = sin(2*pi*40*t)
N = len(t)
Ts = (t[N-1]-t[0])/N

window1 = np.piecewise(t, [t>=t[N//4], t<t[N//4], t>t[3*N//4]], [ 1, 0, 0])
window2 = np.piecewise(t, [t>=t[3*N//8], t<t[3*N//8], t>t[5*N//8]], [ 1, 0, 0])
window3 = np.piecewise(t, [t>=t[7*N//16], t<t[7*N//16], t>t[9*N//16]], [ 1, 0, 0])

xf = np.fft.rfftfreq(N, d=Ts)
hf = np.fft.rfft(h)*2./N
hf1 = np.fft.rfft(h*window1)*2./N
hf2 = np.fft.rfft(h*window2)*2./N
hf3 = np.fft.rfft(h*window3)*2./N



if(0):
    pl.figure()
    pl.plot(xf, np.abs(hf), label='No Window')
    pl.plot(xf, np.abs(hf1), label='Half Duration')
    pl.plot(xf, np.abs(hf2), label='Quater Duration')
    pl.plot(xf, np.abs(hf3), label='Eighth Duration')
    pl.xlabel('frequency, Hz')
    pl.ylabel('Amplitude')
    pl.legend()
    pl.grid()
    pl.show()


