import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import matplotlib.animation as animation




####################################################

class swimmer(object):
    def __init__(self,params):
        print('Initializing system')
        self.LT1 = params['LT1']; self.LT2 = params['LT2']
        self.theta = params['theta']
        self.dx = params['dx']
        self.omega1 = params['omega1']
        self.omega2 = params['omega2']
        self.dOmegaInit = params['dOmegaInit']
        self.nT = params['nT']
        self.dt = params['dt']
        self.drivingFunction = params['drivingFunction']
        self.E = params['E']; self.A1 = params['A1']; self.A2 = params['A2']
        self.zetaN1 = params['zetaN1']; self.zetaT1 = params['zetaT1']
        self.zetaN2 = params['zetaN2']; self.zetaT2 = params['zetaT2']
        self.moment1 = params['moment1']; self.moment2 = params['moment2']
        self.mStart1 = params['mStart1']; self.mEnd1 = params['mEnd1']
        self.mStart2 = params['mStart2']; self.mEnd2 = params['mEnd2']

        #############################
        # Draw relaxed swimmer geometry:
        self.painter()

        #############################

        
    def painter(self):
        self.s = np.concatenate((np.arange(-self.LT1,0,self.dx), np.arange(0,self.LT2+self.dx,self.dx)), axis=0)
        x = np.cos(self.theta/2) * np.concatenate((np.arange(self.LT1,0,-self.dx), np.arange(0,self.LT2+self.dx,self.dx)), axis=0)
        y = np.sin(self.theta/2) * np.concatenate((np.arange(-self.LT1,0,self.dx), np.arange(0,self.LT2+self.dx,self.dx)), axis=0)
        self.x = np.array([x[:],y[:]])


    ###############################
    
    def actuator(self):
        print('Calculating normal driving forces w(x,t)')


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
        
        
        if self.drivingFunction == 2:
            mFunc1 = genWaveform2(self.t, self.omega1)
            mFunc2 = genWaveform2(self.t, self.omega2)
        elif self.drivingFunction == 3:
            mFunc1 = genWaveform3(self.t, self.omega1)
            mFunc2 = genWaveform3(self.t, self.omega2)
        else:
            mFunc1 = genWaveform1(self.t, self.omega1)
            mFunc2 = genWaveform1(self.t, self.omega2)

        self.m = np.zeros(self.x.shape[1])
        mStartIndex1 = np.round(self.mStart1/self.dx)
        mEndIndex1 = np.round(self.mEnd1/self.dx)
        self.m[mStartIndex1:mEndIndex1 + 1] = self.moment1

        mStartIndex2 = np.round(self.mStart2 / self.dx)
        mEndIndex2 = np.round(self.mEnd2 / self.dx)
        self.m[mStartIndex2:mEndIndex2 + 1] = self.moment2

        self.mFunc1 = mFunc1
        self.mFunc2 = mFunc2

        self.w = np.zeros([self.x.shape[1], self.t.shape[0]])
        self.momentCalc(mFunc1, 1)
        self.momontCalc(mFunc2, 1


    #######

    def momentCalc(self, mFunc, tailIndex):
        self.m[0:2] = np.zeros([2])
        self.m[self.m.shape[0]-2:self.m.shape[0]] = np.zeros([2])

        a = np.zeros([self.m.shape[0],self.m.shape[0]])
        self.w = np.zeros([self.x.shape[0],self.t.shape[0]])

        for i in range(0, self.x.shape[0]-1):
            a[i,i+1:self.m.shape[0]] = self.dx*np.arange(1,self.m.shape[0]-i,1)

        a[self.m.shape[0]-1,:] = np.ones([1,self.m.shape[0]])

        wBase = np.linalg.solve(a,self.m)

        for i in range(0,self.t.shape[0]):
            if tailIndex == 1:
                self.w1[:,i] = np.multiply(wBase,mFunc[i])
            elif tailIndex == 2:
                self.w2[:,i] = np.multiply(wBase,mFunc[i])
                

        #return self.w

