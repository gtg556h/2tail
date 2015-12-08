import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import matplotlib.animation as animation
import actuationWaveforms

##########################################################

# Some basic notes:

# {x y z} lab coordinate system, {xi eta zeta} element (flagella) coordinate system.

# Class swimmer has attributes x,y, which are lists of x, y arrays for each constituent system component

# flagella.x is 2D array with form (x,t)

##########################################################

class flagella(object):
    def __init__(self, params):
        self.params = params
        self.LT = params['LT']
        self.alpha = params['alpha']
        self.origin = params['origin']
        self.dx = params['dx']
        self.omega = params['omega']
        self.phi = params['phi']        # Initial phase
        self.drivingFunction = params['drivingFunction']
        self.E = params['E']
        self.A = params['A']
        self.zetaN = params['zetaN']
        self.zetaT = params['zetaT']
        self.moment = params['moment']
        self.mStart = params['mStart']  # distance from root to start of actuator
        self.mEnd = params['mEnd']     # distance from root to end of actuator
        self.xi0 = np.arange(0, self.LT+self.dx, self.dx)
        self.eta0 = np.zeros_like(self.xi0)
        self.x0 = self.origin[0] + self.xi0 * np.cos(self.alpha) - self.eta0 * np.sin(self.alpha)
        self.y0 = self.origin[1] + self.xi0 * np.sin(self.alpha) + self.eta0 * np.cos(self.alpha)


    def processActuator(self, t):
        self.t = t
        self.m = np.zeros_like(self.xi0)
        mStartIndex = np.round(self.params['mStart'] / self.params['dx'])
        mEndIndex = np.round(self.params['mEnd'] / self.params['dx'])
        self.m[mStartIndex:mEndIndex + 1] = self.params['moment']
        self.mFunc = actuationWaveforms.genWaveform(self.drivingFunction, t, self.omega)

        self.momentCalc()


    #######

    def momentCalc(self):

        self.m[0:2] = np.zeros([2]); self.m[-2:] = np.zeros([2])
        a = np.zeros([self.m.shape[0],self.m.shape[0]])
        self.w = np.zeros([self.x.shape[0],self.t.shape[0]])

        for i in range(0, self.x.shape[0]-1):
            a[i,i+1:self.m.shape[0]] = self.dx*np.arange(1,self.m.shape[0]-i,1)

        a[self.m.shape[0]-1,:] = np.ones([1,self.m.shape[0]])

        wBase = np.linalg.solve(a,self.m)

        for i in range(0,self.t.shape[0]):
            self.w[:,i] = np.multiply(wBase,self.mFunc[i])


            
#########################################################


        
class swimmer(object):
    def __init__(self, structures, tMax, dt):
        self.structures = structures
        self.t = np.arange(0,tMax+dt,dt)

        for i in structures:
            i.x = np.zeros([i.x0.shape[0], self.t.shape[0]])
            i.y = np.zeros([i.y0.shape[0], self.t.shape[0]])
            i.xi = np.zeros([i.xi0.shape[0], self.t.shape[0]])
            i.eta = np.zeros([i.eta0.shape[0], self.t.shape[0]])

            i.x[:,0] = i.x0
            i.y[:,0] = i.y0
            i.xi[:,0] = i.xi0
            i.eta[:,0] = i.eta0

        print('Calculating normal driving forces w(x,t)')
        for i in structures:
            i.processActuator(self.t)


        self.theta = np.zeros_like(self.t)
        self.X = np.zeros_like(self.t)
        self.Y = np.zeros_like(self.t)


    def numSolve(self):
        print('Solving for y(x,t)')

        # for each element, generate D4, Dt, a, c

        # Solve for displacement for each element, *including* current swimmer translation,
        # rotation, and contribution from integrating green's functions from all other elements.

        # Calculate required origin translation and rotation to enforce zero force and torque

        # Recompute element displacements and iterate until correction translations and
        # rotations are less than some convergence threshold

            
####################################################











