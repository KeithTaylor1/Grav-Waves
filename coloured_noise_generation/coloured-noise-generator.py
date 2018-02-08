"""
Created on Tue Feb  6 16:23:46 2018

@author: michael stephens

Program to take signal in time domain from *.txt file in same directory and 
apply noise. Output is a .txt file in time domain. 

Input signal should contain time and signal amplitudes in same file.

Program assumes that the file H1-PSD.txt is in same directory as 
coloured-noise-generator.py
"""

import numpy as np
import numpy.fft as fft

filename = input('Enter file name (including \".txt\"): ')
while 1:
    try:
        signal = np.loadtxt(filename)
        break
    except FileNotFoundError:
        print('File not found, please try again... ')




ligo_psd = np.loadtxt('H1-PSD.txt')
nf = ligo_psd[0, :]
xf = ligo_psd[1, :]

