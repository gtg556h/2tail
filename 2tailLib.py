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
        self.LT1 = params['LT1']; self.LT2 = params['LT2']
        self.theta = params['theta']
        self.dx = params['dx']
        self.omega1 = params['']
        self.omega2 = params['']
        self.dOmegaInit = params['']
        self.nT = params['']
        self.dt = params['']
        self.drivingFunction = params['']
        self.E = params['']; self.A1 = params['']; self.A2 = params['']
        self.zetaN1 = params['zetaN1']; self.zetaT1 = params['zetaT1']
        self.zetaN2 = params['zetaN2']; self.zetaT2 = params['zetaT2']
        self.moment1 = params['moment1']; self.moment2 = params['moment2']
        self.mStart1 = params['mStart1']; self.mEnd1 = params['mEnd1']
        self.mStart2 = params['mStart2']; self.mEnd2 = params['mEnd2']

        self.painter()


    def painter(self):
        self.s = np.concatenate((np.arange(-LT1,0,dx), np.arange(0,LT2+dx,dx)), axis=0)
        x = np.cos(self.theta/2) * np.concatenate(np.arange(LT1,0,-dx), np.arange(0,LT2+dx,dx), axis=0)
        y = np.sin(self.theta/2) * np.concatenate(np.arange(-LT1,0,dx), np.arange(0,LT2+dx,dx), axis=0)
        self.x = np.array([x[:],y[:]])


        
    def actuator(self):
        print('Calculating normal driving forces w(x,t)')

       
