# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 14:19:08 2018

@author: njg573
"""

import numpy as np
import scipy as sp
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as pl
import matplotlib.mlab as ml
from scipy import fftpack as fftp
from scipy import signal

t = np.linspace(0,100,1000)
h = sin(2.*pi*t)
window1 = np.blackman(1000)
'''np.blackman cannot be applied to any 
data with a different # of points than itself'''
data1 = h * window1
window1fft = fftp.fft(window1) / (len(window1)/2.0)
freq = np.linspace(-0.5, 0.5, len(window1fft))
response = 20 * np.log10(np.abs(fftp.fftshift(window1fft / abs(window1fft).max())))

fig1 = pl.figure()

ax1 = fig1.add_subplot(3,1,1)
ax1.plot(window1)
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.grid()
'''shows spectral leakage of the window in the F domain'''
ax2 = fig1.add_subplot(3,1,2)
ax2.plot(freq, response)
ax2.set_xlabel('Frequency')
ax2.set_ylabel('response')
ax2.grid()

ax3 = fig1.add_subplot(3,1,3)
ax3.plot(data1)
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.grid()

pl.show()

#%%
window2 = sp.signal.tukey(1000)

data2 = h * window2
'''sp.signal.tukey cannot be applied to any 
data with a different # of points than itself'''
data2fft = fftp.fft(window2)
A = data2fft / (len(window2)/2.0)
freq = np.linspace(-0.5, 0.5, len(A))
response = 20 * np.log10(np.abs(fftp.fftshift(A / abs(A).max())))

fig2 = pl.figure()

ax1 = fig2.add_subplot(3,1,1)
ax1.plot(window2)
ax1.set_title('Window plot')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.grid()
'''shows spectral leakage of the window in the F domain'''
ax2 = fig2 .add_subplot(3,1,2)
ax2.plot(freq, response)
ax2.set_xlabel('Frequency')
ax2.set_ylabel('response')
ax2.grid()

ax3 = fig2.add_subplot(3,1,3)
ax3.plot(data2)
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.grid()

pl.show()

#%%

h = ['noisySinusoid.npy']
t = ['time.npy']
DataWindow = sp.signal.tukey(len(h))

data3 = h * DataWindow
'''sp.signal.tukey cannot be applied to any 
data with a different # of points than itself'''
data3fft = fftp.fft(DataWindow)
A = data3fft / (len(DataWindow)/2.0)
freq = np.linspace(-0.5, 0.5, len(A))
response = 20 * np.log10(np.abs(fftp.fftshift(A / abs(A).max())))

fig3 = pl.figure()

ax1 = fig3.add_subplot(3,1,1)
ax1.plot(h,t)
ax1.set_title('Window plot')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.grid()
'''shows spectral leakage of the window in the F domain'''
ax2 = fig3.add_subplot(3,1,2)
ax2.plot(freq, response)
ax2.set_xlabel('Frequency')
ax2.set_ylabel('response')
ax2.grid()

ax3 = fig3.add_subplot(3,1,3)
ax3.plot(data3)
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.grid()

pl.show()



pl.show()
