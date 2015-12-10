import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import matplotlib.animation as animation
import actuationWaveforms
import pdb

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
        self.A = np.ones_like(self.x0)*self.A
        self.zetaN = np.ones_like(self.x0)*self.zetaN
        self.zetaT = np.ones_like(self.x0)*self.zetaT


    #####################################################

        
    def processActuator(self, t, dt):
        self.t = t
        self.dt = dt
        self.m = np.zeros_like(self.xi0)
        mStartIndex = np.round(self.params['mStart'] / self.params['dx'])
        mEndIndex = np.round(self.params['mEnd'] / self.params['dx'])
        self.m[mStartIndex:mEndIndex + 1] = self.params['moment']
        self.mFunc = actuationWaveforms.genWaveform(self.drivingFunction, t, self.omega)

        self.momentCalc()


    #####################################################
    

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


    ######################################################
            

    def differentialOperator(self):
        nx = self.x.shape[0]
        nt = self.t.shape[0]

        self.D4 = np.zeros([nx, nx])
        self.Dt = np.zeros_like(self.D4)
        self.a = np.zeros_like(self.D4)
        self.c = np.zeros([nx, nt])
        self.c2 = np.zeros([nx, nt])

        for i in range(2, nx-2):
            self.D4[i, i-2:i+3] = self.A[i] * np.array([1, -4, 6, -4, 1]) / self.dx**3

        self.Dt[2:nx-2, 2:nx-2] = np.diag(self.zetaN[2:nx-2] * self.dx / self.dt)

        # Apply *fixed* angular and position LHS BC:
        self.a[0, 0] = 1
        self.a[1, 1] = 1

        self.c[0,:] = np.zeros(nt)
        self.c[1,:] = self.c[0,:]

        # Apply *free* BC to RHS:
        self.a[nx-2, nx-4:nx] = np.array([-1, 3, -3, 1]) / self.dx**3
        self.a[nx-1, nx-4:nx] = np.array([-1, 4, -5, 2]) / self.dx**2
        self.c[nx-2:nx, :] = np.zeros([2, nt])

        # Fill out c:
        self.c[2:nx-2, :] = self.c[2:nx-2, :] + self.w[2:nx-2, :]

        # Build differential operator:
        self.a = self.a + self.D4 + self.Dt


    ########################################################
    

    def stepElement(self, i, v, theta, X, Y):
        nx = self.x.shape[0]
        nt = self.t.shape[0]

        self.c2[2:nx-2, i] = self.c[2:nx-2, i] + np.multiply(self.zetaN[2:nx-2], self.eta[2:nx-2, i-1]) * self.dx/self.dt
        
        # pdb.set_trace()
        
        # Add contribution of local fluid flow in eta direction:
        self.c2[2:nx-2, i] = self.c2[2:nx-2, i] + np.multiply(self.zetaN[2:nx-2], v[2:nx-2]) * self.dx/self.dt

        self.eta[:,i] = np.linalg.solve(self.a, self.c2[:,i])
        self.xi[:,i] = self.xi0

        self.x[:,i] = X + self.xi0 * np.cos(self.alpha + theta) - self.eta[:,i] * np.sin(self.alpha + theta)
        self.y[:,i] = Y + self.xi0 * np.sin(self.alpha + theta) + self.eta[:,i] * np.cos(self.alpha + theta)


    #########################################################

        
    def computeDragConstants(self, theta):
        # Following *assumes* that flagella is radial WRT to origin.  Modify!!
        # Rotational moment = self.cxy * thetaDot
        self.cxy = np.sum(sef.xi0 * self.zetaN * self.dx)    # Rotational drag coeff

        # x-drag = self.cx * Xdot
        self.cx = self.dx * np.sum(self.zetaN * np.sin(self.alpha + theta))
        self.cx = self.cx + self.dx * np.sum(self.zetaT * np.cos(self.alpha + theta))

        # y-drag = self.cy * Ydot
        self.cy = self.dx * np.sum(self.zetaN * np.cos(self.alpha + theta))
        self.cy = self.cy + self.dx * np.sum(self.zetaT * np.sin(self.alpha + theta))


    #############################################################

        
    def calcDrag(self, theta, thetaDot, XDot, YDot, timestep):
        
        
            
#########################################################
#########################################################
#########################################################

        
class swimmer(object):
    def __init__(self, structures, tMax, dt):
        self.structures = structures
        self.t = np.arange(0,tMax+dt,dt)
        self.dt = dt

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
        for i in self.structures:
            i.processActuator(self.t, dt)

        self.theta = np.zeros_like(self.t); self.thetaDot = np.zeros_like(self.t)
        self.X = np.zeros_like(self.t); self.XDot = np.zeros_like(self.t)
        self.Y = np.zeros_like(self.t); self.YDot = np.zeros_like(self.t)


    ####################################################    
        

    def numSolve(self):
        print('Solving for y(x,t)')

        # for each element, generate D4, Dt, a, c

        # Solve for displacement for each element, *including* current swimmer translation,
        # rotation, and contribution from integrating green's functions from all other elements.

        # Calculate required origin translation and rotation to enforce zero force and torque

        # Recompute element displacements and iterate until correction translations and
        # rotations are less than some convergence threshold

        self.computeIZZ()

        for j in self.structures:
            j.differentialOperator()

        for i in range(1,self.t.shape[0]):
            if np.mod(i, 50) == 0:
                print(i/self.t.shape[0])

            self.stepStructure(i)

            
    ###################################################

    
    def computeDragConstants(self, theta):
        #####
        self.cxy = 0
        
        for j in self.structures:
            j.IZZ = np.sum(self.xi0 * self.zetaN * self.dx)
            self.IZZ = self.IZZ + j.IZZ

            
    ####################################################

            
    def stepStructure(self, i):

        thetaThresh = 100
        XThresh = 100
        YThresh = 100
        
        theta = self.theta[i-1] + self.thetaDot[i-1]*self.dt
        X = self.X[i-1] + self.XDot[i-1]*self.dt
        Y = self.Y[i-1] + self.YDot[i-1]*self.dt

        thetaOld = theta; XOld = X; YOld = Y
        
        iteration = 0

        while iteration < 2 or (theta - thetaOld > thetaThresh) or (X - XOld > XThresh) or (Y - YOld > YThresh):

            iteration +=1
            print(iteration)

            XDot = (X - self.X[i-1])/self.dt
            YDot = (Y - self.Y[i-1])/self.dt
            thetaDot = (theta - self.theta[i-1])/self.dt

            M = 0
            Fx = 0
            Fy = 0
            
            for j in self.structures:

                
                # v = self.computeVelocity(theta, XDot, YDot, j)
                
                v = -2*np.ones_like(j.x0)

                theta = 0
                X = 0
                Y = 0
                j.stepElement(i, v, theta, X, Y)

                # The following moment calculation is a bit simplistic....
                # Should be updated to accomodate truly generic assemblies
                # (With roots not necessarily fixed at origin)
                M = M - np.sum((j.eta[:,i] - j.eta[:,i-1])/dt * j.xi[:,i] * j.zetaN)*j.dx
                M = M - thetaDot * j.IZZ

                
                Fx = Fx - np.sum(XDot * j.zetaN * j.dx * sin(j.alpha + theta))
                Fx = Fx - np.sum(XDot * j.zetaT * j.dx * cos(j.alpha + theta))
                Fx = Fx + np.sum(j.eta[:,i] * j.zetaN * j.dx * sin(j.alpha + theta))

                Fy = Fy - np.sum(YDot * j.zetaN * j.dx * cos(j.alpha + self.theta))
                Fy = Fy - np.sum(YDot * j.zetaT * j.dx * sin(j.alpha + self.theta))
                Fy = Fy - np.sum(j.eta[:,i] * j.zetaN * j.dx * cos(j.alpha + self.theta))

            XOld = X; YOld = Y; thetaOld = theta
            X, Y, theta = self.balanceForces(XOld, YOld, thetaOld, XDot, YDot, thetaDot)


    ###############################################################
            

    def balanceForces(XOld, YOld, thetaOld, XDot, YDot, thetaDot):
        cx = 0; yx = 0; IZZ = self.IZZ

        for j in self.structures:
            cx = cx + np.sum(j.zetaN * j.dx * sin(j.alpha + thetaOld))
            cx = cx + np.sum(j.zetaT * j.dx * cos(j.alpha + thetaOld))

            cy = cy + np.sum(j.zetaN * j.dx * cos(j.alpha + thetaOld))
            cy = cy + np.sum(j.zetaT * j.dx * sin(j.alpha + thetaOld))

        Fx = (f(thetaDot, XDot))
        Fy = f(thetaDot, YDot)
        M = f()
    

        

                
    ###########################################################


    def assemble(self):
        # Build global displacement from set of individual element displacement functions:

        self.x = None
        self.y = None

        for j in self.structures:
            if self.x == None:
                self.x = j.x
                self.y = j.y
            else:
                self.x = np.concatenate((self.x, j.x), axis=0)
                self.y = np.concatenate((self.y, j.y), axis=0)

        ## !!!!! Fix me later
        self.X = self.x
            
                
    ###############################################            


    def plotDisp(self,DF=1,plotFrac=1):

        print('Displaying solution shapes')
        nx = self.x.shape[0]
        nt = self.t.shape[0]

        nFrames = np.int(np.floor(nt/DF*plotFrac))

        yRange = np.max(self.y)-np.min(self.y)
        xRange = np.max(self.X)-np.min(self.X)


        figScale = 8.0
        fig = plt.figure(figsize=[figScale,yRange/np.max(self.x)*figScale])
        ax = plt.axes()

        ax.set_xlim([np.min(self.X)-.1*xRange,np.max(self.X)+.1*xRange])
        ax.set_ylim([np.min(self.y)-.25*yRange,np.max(self.y)+.25*yRange])
        ax.set_yticklabels([])
        
        line, = ax.plot([], [],'.',lw=2, markersize=4)

        # Initialization function: plot the background of each frame
        def init():
            line.set_data([], [])
            return line,

        # animation function, called sequentially:
        def animate(i):
            y2 = self.y[:,DF*i]
            X2 = self.X[:,DF*i]
            #line.set_data(self.x,y2)
            line.set_data(X2,y2)
            return line,

        # Call the animator:
        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=nFrames, interval=50, blit=False, repeat=True)

        plt.show()
                
            
####################################################











