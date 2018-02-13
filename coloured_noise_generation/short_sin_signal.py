
"""
Created on Thu Feb  8 15:33:37 2018

@author: michael stephens

Creates 3 sec 60Hz sine signal and outputs as *.txt

"""

import numpy as np


t = np.arange(0, 3, 0.0001)
h = np.sin(2*np.pi*60*t)*1e-19

signal = np.array([t,h]).T

np.savetxt('signal.txt', signal)