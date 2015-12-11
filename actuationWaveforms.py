import numpy as np
import numpy.linalg

def genWaveform(drivingFunction, t, omega):
    if drivingFunction == 2:
        return genWaveform2(t, omega)
    elif drivingFunction == 3:
        return genWaveform3(t, omega)
    elif drivingFunction == 4:
        return genWaveform4(t, omega)
    else:
        return genWaveform1(t, omega)

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

def genWaveform4(t,omega):
    # Nominally tweaked waveform form conference proceedings etc...

    m = np.sin(omega/2.0*t)
    m = (np.abs(m) - 0.25)/0.75
    m[np.where(m<0)[0]] = 0

    m = np.ones_like(t)
    
    return m
