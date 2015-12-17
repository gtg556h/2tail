import numpy as np
import matplotlib.pyplot as plt
import twoTailLib as tt


ID1 = 'flagella1'
ID2 = 'flagella2'

LT1 = 1503.0       # Length of tail 1
LT2 = 424.0#LT1        # Length of tail 2
alpha1 = 0#np.pi/4
alpha2 = np.pi#-alpha1
dx = 3.0
origin1 = np.array([0,0])
origin2 = origin1

omega1 = 3.25*2*np.pi     # Frequency of actuator 1
omega2 = omega1          # Frequency of actuator 2
phi1 = 0   # initial phase of actuator 1
phi2 = 0       
nT = 8.0                  # Number of actuation periods to compute
dt = 0.001
dt = 0.004
 
drivingFunction = 1       # Shape of actuation (see 2tailLib)

E = 3.86                  # Modulus of elasticity
IT1 = 7.9**3 * 22.0 / 12.0
IT2 = IT1
A1 = E*IT1                # Bending stiffness of tail 1
A2 = E*IT2                # Bending stiffness of tail 2

zetaN1 = 9.25E-9          # Normal drap coeff, tail 1
zetaT1 = zetaN1 / 2.0     # Tangential drag coeff, tail 1
zetaN2 = 11.19E-9           # Tail 2
zetaT2 = zetaN2 / 2.0           # Tail 2

moment1 = -37.37            # Peak bending moment, actuator 1
mStart1 = 45.0            # Distance of start of actuator 1 from tail1 root
mEnd1 = mStart1 + 60.0     # Distance of end of actuator 1 from tail root
moment2 = 0          # Peak bending moment, actuator 2
mStart2 = mStart1            # Distance of start of actuator 2 from tail2 root
mEnd2 = mEnd1     # Distance of end of actuator 2 from tail2 root

tMax = 2*np.pi/omega1 * nT

params1 = {'ID':ID1, 'LT':LT1, 'alpha':alpha1, 'origin':origin1, 'dx':dx, 'omega':omega1, 'phi':phi1, 'drivingFunction':drivingFunction, 'E':E, 'A':A1, 'zetaN':zetaN1, 'zetaT':zetaT1, 'moment':moment1, 'mStart':mStart1, 'mEnd':mEnd1}

params2 = {'ID':ID2, 'LT':LT2, 'alpha':alpha2, 'origin':origin2, 'dx':dx, 'omega':omega2, 'phi':phi2, 'drivingFunction':drivingFunction, 'E':E, 'A':A2, 'zetaN':zetaN2, 'zetaT':zetaT2, 'moment':moment2, 'mStart':mStart2, 'mEnd':mEnd2}



f1 = tt.flagella(params1)
f2 = tt.flagella(params2)
structures = [f1, f2]
#structures = [flagella1]

s = tt.swimmer(structures, tMax, dt)
s.numSolve(hydrodynamicCoupling=1, propulsion=1)
s.assemble()
plt.plot(s.t,s.X)
plt.show()
s.plotDisp(DF=4)




###########################################

def plotFrame(i):
    plt.plot(s.x[:,i], s.y[:,i])
    plt.show()









############################################

# Initialize and process swimmer:

# s = tt.swimmer(params)
# s.actuator()
# s.numSolve()
# s.propulsionCalc()
# print(np.mean(s.Ux))
# s.plotDisp(16,1)

if 0:
    plt.plot(s.x[0,:],s.x[1,:])
    plt.axis('equal')
    plt.show()
    

if 0:
    for i in s.structures:
        plt.plot(i.x[:,0], i.y[:,0])

    plt.show()
