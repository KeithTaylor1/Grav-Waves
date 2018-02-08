'''
Code taken from LOSC Turorial, edited by Michael Stephens on 08/02/2018

Takes raw LIGO data from LOSC totorial *.zip file and outputs text file 
containing frequencies and PSD for 32s around H1 detection of GW150914

'''
# Standard python numerical analysis imports:
import numpy as np
import json

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# LIGO-specific readligo.py 
import readligo as rl



eventname = 'GW150914'

# Read the event properties from a local json file
fnjson = "BBH_events_v3.json"
try:
    events = json.load(open(fnjson,"r"))
except IOError:
    print("Cannot find resource file "+fnjson)
    print("You can download it from https://losc.ligo.org/s/events/"+fnjson)
    print("Quitting.")
    quit()

# did the user select the eventname ?
try: 
    events[eventname]
except:
    print('You must select an eventname that is in '+fnjson+'! Quitting.')
    quit()

# Extract the parameters for the desired event:
event = events[eventname]
fn_H1 = event['fn_H1']              # File name for H1 data
fn_template = event['fn_template']  # File name for template waveform
fs = event['fs']                    # Set sampling rate
tevent = event['tevent']            # Set approximate event GPS time
fband = event['fband']              # frequency band for bandpassing signal
#print("Reading in parameters for event " + event["name"])
#print(event)


#----------------------------------------------------------------
# Load LIGO data from a single file.
# FIRST, define the filenames fn_H1 and fn_L1, above.
#----------------------------------------------------------------
try:
    # read in data from H1 and L1, if available:
    strain_H1, time_H1, chan_dict_H1 = rl.loaddata(fn_H1, 'H1')
except:
    print("Cannot find data files!")
    print("You can download them from https://losc.ligo.org/s/events/"+eventname)
    print("Quitting.")
    quit()

time = time_H1
# the time sample interval (uniformly sampled!)
dt = time[1] - time[0]

    
make_psds = 1
if make_psds:
    # number of sample for the fast fourier transform:
    NFFT = 4*fs
    Pxx_H1 = mlab.psd(strain_H1, Fs = fs, NFFT = NFFT)
    
np.savetxt('H1-PSD.txt', Pxx_H1)

if 0:
    f_min = 20.
    f_max = 2000. 
    plt.figure(figsize=(10,8))
    plt.loglog(freqs, Pxx_H1,'r',label='H1 strain')
    plt.xlim([f_min, f_max])#, 1e-24, 1e-19])
    plt.grid('on')
    plt.ylabel('PSD ($strain^2/Hz$)')
    plt.xlabel('Freq (Hz)')
