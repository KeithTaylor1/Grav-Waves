"""coloured-noise-generator
Created on Tue Feb  6 16:23:46 2018

@author: michael stephens

Program to take signal in time domain from *.txt file in same directory and
apply noise. Output is a .txt file in time domain. 

Input signal should contain time and signal amplitudes in same file and be over
16500 data points long.

Program assumes that the file H1-PSD.txt is in same directory as 
coloured-noise-generator.py, this is generated by noise-psd-ligoH1-GW150914
"""

import numpy as np
import numpy.fft as fft


# load signal from *.txt file, checking user input
filename = input('Enter file name (including \".txt\"): ')
while 1:
    try:
        signal = np.loadtxt(filename)
        break
    except FileNotFoundError:
        filename = input('File not found, please try again... ')

t = signal[0, :]
h = signal[1, :]

N = len(h)
Ts = t[1]-t[0]

xf_h = fft.rfftfreq(N, Ts)
hf = fft.rfft(h)  # convert signal to freq domain


#load noise and and scale to signal, also randomise phase
ligo_psd = np.loadtxt('H1-PSD.txt')
xf_n = ligo_psd[1, :]
nf = np.sqrt(ligo_psd[0, :])*np.exp(2j*np.pi*np.random.rand())


# ensure both have same length in the frequency domain
if len(hf) < len(nf):
    hf = np.interp(xf_n, xf_h, hf)
    xf = xf_n
elif len(hf) > len(nf):
    nf = np.interp(xf_h, xf_n, nf)
    xf = xf_h


# combine noise and signal and convert back to time domain
df = nf#+nf  
d = fft.irfft(df)






