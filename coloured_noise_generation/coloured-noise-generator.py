"""
Created on Tue Feb  6 16:23:46 2018

@author: michael

Program to take signal in time domain from *.txt file or numpy array and apply 
noise. Output is a .txt file in time domain. 
"""

import numpy as np

ligo_psd = np.loadtxt('H1-PSD.txt')

Pxx = ligo_psd[0,:]
freq = ligo_psd[1,:]
