import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import matplotlib.animation as animation



################################################
# Actuation waveform generation:

def genWaveform1(t,omega):
    # Standard waveform used in NatureComm paper
    m = np.sin(omega/2.0*t)
    m = (np.abs(m) - 0.3)/0.7
    m[np.where(m<0)[0]] = 0
    return m

def genWaveform2(t,omega):
    # Pefect sinusoidal, 2-sided driving function
    m = np.sin(omega*t)
    #m = (np.abs(m) - 0.3)/0.7
    #m[np.where(m<0)[0]] = 0
    return m

def genWaveform3(t,omega):
    # Nominally tweaked waveform form conference proceedings etc...
    m = np.sin(omega/2.0*t)
    m = (np.abs(m) - 0.25)/0.75
    m[np.where(m<0)[0]] = 0
    return m


####################################################

class name(object):
    def __init__(self,params):
        print('Initializing system')
        self.LT1 = params['LT1']



    def actuator(self):
        print('Calculating normal driving forces w(x,t)')
