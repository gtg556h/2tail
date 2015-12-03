import numpy as np
import matplotlib.pyplot as plt
import 2tailLib as 2t

LT1 = 1200       # Length of tail 1
LT2 = LT1        # Length of tail 2
theta = 50/180*np.pi    # Angle between tails
dx = 3.0        

omega1 = 3.25*2*np.pi     # Frequency of actuator 1
omega2 = omega1           # Frequency of actuator 2
dOmegaInit = 0   # Difference in *initial* phase of actuators
nT = 1.0                  # Number of actuation periods to compute
dt = 0.001       
 
drivingFunction = 3       # Shape of actuation (see 2tailLib)

E = 3.86                  # Modulus of elasticity
IT1 = 7.9**3 * 22.0 / 12.0
IT2 = IT1
A1 = E*IT1                # Bending stiffness of tail 1
A2 = E*IT2                # Bending stiffness of tail 2

zetaN1 = 9.25E-9          # Normal drap coeff, tail 1
zetaT1 = zetaN1 / 2.0     # Tangential drag coeff, tail 1
zetaN2 = zetaN1           # Tail 2
zetaT2 = zetaT1           # Tail 2

moment1 = 37.37            # Peak bending moment, actuator 1
mStart1 = 200.0            # Distance of start of actuator 1 from tail1 root
mEnd1 = mStart1 + 60.0     # Distance of end of actuator 1 from tail root
moment2 = moment1          # Peak bending moment, actuator 2
mStart2 = 200.0            # Distance of start of actuator 2 from tail2 root
mEnd2 = mStart2 + 60.0     # Distance of end of actuator 2 from tail2 root

params = {'LT1':LT1, 'LT2':LT2, 'theta':theta, 'dx':dx, 'omega1':omega1, 'omega2':omega2, 'dOmegaInit':dOmegaInit, 'nT':nT, 'dt':dt, 'drivingFunction':drivingFunction, 'E':E, 'A1':A1, 'A2':A2, 'zetaN1':zetaN1, 'zetaT1':zetaT1, 'zetaN2':zetaN2, 'zetaT2':zetaT2, 'moment1':moment1, 'mStart1':mStart1, 'mEnd1':mEnd1, 'moment2':moment2, 'mStart2':mStart2, 'mEnd2':mEnd2}


###########################################

# Initialize and process swimmer:

s = 2t.2tail(params)
# s.actuator()
# s.numSolve()
# s.propulsionCalc()
# print(np.mean(s.Ux))
# s.plotDisp(16,1)
