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
    

    def stepElement(self, i, theta, X, Y, XDot, YDot, thetaDot):
        nx = self.x.shape[0]
        nt = self.t.shape[0]

        v = self.calcLocalVelocity(theta, XDot, YDot, thetaDot)

        self.c2[2:nx-2, i] = self.c[2:nx-2, i] + np.multiply(self.zetaN[2:nx-2], self.eta[2:nx-2, i-1]) * self.dx/self.dt
        
        # pdb.set_trace()
        
        # Add contribution of local fluid flow in eta direction:
        self.c2[2:nx-2, i] = self.c2[2:nx-2, i] - np.multiply(self.zetaN[2:nx-2], v[2:nx-2]) * self.dx

        self.eta[:,i] = np.linalg.solve(self.a, self.c2[:,i])
        self.xi[:,i] = self.xi0

        self.x[:,i] = X + self.xi0 * np.cos(self.alpha + theta) - self.eta[:,i] * np.sin(self.alpha + theta)
        self.y[:,i] = Y + self.xi0 * np.sin(self.alpha + theta) + self.eta[:,i] * np.cos(self.alpha + theta)



    def calcLocalVelocity(self, theta, XDot, YDot, thetaDot):
        # For given element, calculate the local velocity
        # generated by translation and rotation of swimmer
        # As well as element-element coupling:
                
        v = -self.xi0 * thetaDot
        v = v - YDot * np.cos(self.alpha + theta)
        v = v + XDot * np.sin(self.alpha + theta)

        return v

    #########################################################

        
    def computeDragConstants(self, theta):
        # Formulas compute drag on body provided associated velocities:

        # Returs 3x3 matrix c, such that c * {XDot, YDot, thetaDot} returns FX, FY, M

        c = np.zeros([3,3])
        

        # x-drag = self.cx * Xdot
        c[0,0] = -self.dx * np.sum(self.zetaN) * np.abs(np.sin(self.alpha + theta))
        c[0,0] = c[0,0] - self.dx * np.sum(self.zetaT) * np.abs(np.cos(self.alpha + theta))
        # Assuming flagella radial to origin:
        c[2,0] = self.dx * np.sum(self.xi0 * self.zetaN) * np.sin(self.alpha + theta)
                                  
        # y-drag = self.cy * Ydot
        c[1,1] = -self.dx * np.sum(self.zetaN) * np.abs(np.cos(self.alpha + theta))
        c[1,1] = c[1,1] - self.dx * np.sum(self.zetaT) * np.abs(np.sin(self.alpha + theta))
        # Assuming flagella radial to origin:
        c[2,1] = -self.dx * np.sum(self.xi0 * self.zetaN) * np.cos(self.alpha + theta)

        # Following *assumes* that flagella is radial WRT to origin.  Modify!!
        # Rotational moment = self.cxy * thetaDot
        c[2,2] = -self.dx * np.sum(self.xi0 * self.xi0 * self.zetaN)
        c[0,2] = self.dx * np.sum(self.xi0 * self.zetaN) * np.sin(self.alpha + theta)
        c[1,2] = -self.dx * np.sum(self.xi0 * self.zetaN) * np.cos(self.alpha + theta)
                

        return c

    #############################################################
        
    def calcDrag(self, theta, thetaDot, XDot, YDot, timestep):
        # Compute drag forces and moments ON swimmer elements:

        c = self.computeDragConstants(theta)
        #pdb.set_trace()
        etaDot = (self.eta[:,timestep] - self.eta[:, timestep-1]) / self.dt

        XDotVec = np.array([[XDot],[YDot],[thetaDot]])
        
        fx = self.dx * np.sum(etaDot * self.zetaN) * np.sin(self.alpha + theta)
        fx = fx + np.dot(c, XDotVec)[0,0]

        fy = self.dx * np.sum(-etaDot * self.zetaN) * np.cos(self.alpha + theta)
        fy = fy + np.dot(c, XDotVec)[1,0]

        
        # Following calculation assumes that flagella is radial to origin.
        # Create scaling array that maps flagella location with radial component
        # in future!!!
        
        moment = -self.dx * np.sum(etaDot * self.xi0 * self.zetaN)
        moment = moment + np.dot(c, XDotVec)[2,0]

        return fx, fy, moment

    #def calcCorrectionDrag(self, theta, thetaDot, XDot, YDot):
        # Compute drag forces for
        
        
        
            
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


        for j in self.structures:
            j.differentialOperator()

        for i in range(1,self.t.shape[0]):
        #for i in range(1,3):
            if np.mod(i, 50) == 0:
                print(i/self.t.shape[0])

            self.stepStructure(i)

            
    ####################################################

            
    def stepStructure(self, i):

        thetaThresh = .001
        XThresh = .001
        YThresh = .001
        
        theta = self.theta[i-1] + self.thetaDot[i-1]*self.dt
        X = self.X[i-1] + self.XDot[i-1]*self.dt
        Y = self.Y[i-1] + self.YDot[i-1]*self.dt

        #############
        # DEBUGGING:
        #X = self.X[i-1]
        #Y = self.Y[i-1]
        #############

        thetaOld = theta; XOld = X; YOld = Y
        
        iteration = 0
 
        while iteration < 10 and (iteration < 2 or (theta - thetaOld > thetaThresh) or (X - XOld > XThresh) or (Y - YOld > YThresh)):
                                  

            iteration +=1
            #print(iteration)


            XDot = (X - self.X[i-1])/self.dt
            YDot = (Y - self.Y[i-1])/self.dt
            thetaDot = (theta - self.theta[i-1])/self.dt

            
            for j in self.structures:
                j.stepElement(i, theta, X, Y, XDot, YDot, thetaDot)

                
            XOld = X; YOld = Y; thetaOld = theta

            # Calculate *corrections* to velocities in order to balance forces/moments
            # with current flagella deformations:
            cXDot, cYDot, cthetaDot = self.balanceForces(X, Y, theta, XDot, YDot, thetaDot, i)
            #cXDot = cXDot / 10
            #cYDot = cYDot / 10
            #cthetaDot = cthetaDot / 10
            #print(cXDot,cYDot,cthetaDot)
            X = X + cXDot * self.dt;  XDot = XDot + cXDot
            Y = Y + cYDot * self.dt;  YDot = YDot + cYDot
            theta = theta + cthetaDot * self.dt; thetaDot = thetaDot + cthetaDot

        
          
            #pdb.set_trace()

        # print(iteration)
        self.X[i] = X; self.XDot[i] = XDot
        self.Y[i] = Y; self.YDot[i] = YDot
        self.theta[i] = theta; self.thetaDot[i] = thetaDot

        # Log net drag and moment:
        # self.FX[i] = self.computeDrag
        
        

    ###############################################################


    def calcForces(self, theta, XDot, YDot, thetaDot, timestep):
        # Do me next.

        FX = 0; FY = 0; M = 0;
        
        for j in self.structures:
            fx, fy, moment = j.calcDrag(theta, thetaDot, XDot, YDot, timestep)
            FX = FX + fx
            FY = FY + fy
            M = M + moment

        return FX, FY, M


    ###############################################################

    def balanceForces(self, X, Y, theta, XDot, YDot, thetaDot, timestep):
        # Returns cXDot, cYDot, cthetaDot, the correction terms
        # to XDot, etc..., req'd to yield zero net force and moment
        # With current flagella deformations:


        C = self.computeDragConstants(theta)

        # Calculate forces applied ON the swimmer:
        FX, FY, M = self.calcForces(theta, XDot, YDot, thetaDot, timestep)
        #print(FX,FY,M)
        # pdb.set_trace()

        cXDotVec = np.linalg.solve(C, np.array([[-FX],[-FY],[-M]]))
        cXDot = cXDotVec[0,0]
        cYDot = cXDotVec[1,0]
        cthetaDot = cXDotVec[2,0]
        #pdb.set_trace()
        return cXDot, cYDot, cthetaDot
    

        
    ###################################################

    
    def computeDragConstants(self, theta):
        # CX, CY, etc... are drag coeff for structure
        # cx, cy, etc... are drag coeff for individual elements


        C = np.zeros([3,3])
        
        for j in self.structures:

            C = C + j.computeDragConstants(theta)

        return C

                
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
        #self.X = self.x
            
                
    ###############################################            


    def plotDisp(self,DF=1,plotFrac=1):

        x = self.x
        y = self.y
        
        print('Displaying solution shapes')
        #nx = self.x.shape[0]
        nx = x.shape[0]
        
        nt = self.t.shape[0]

        nFrames = np.int(np.floor(nt/DF*plotFrac))

        yRange = np.max(y)-np.min(y)
        xRange = np.max(x)-np.min(x)


        figScale = 8.0
        fig = plt.figure(figsize=[figScale,yRange/np.max(x)*figScale])
        ax = plt.axes()

        ax.set_xlim([np.min(x)-.1*xRange,np.max(x)+.1*xRange])
        ax.set_ylim([np.min(y)-.25*yRange,np.max(y)+.25*yRange])
        ax.set_yticklabels([])
        
        line, = ax.plot([], [],'.',lw=2, markersize=4)

        # Initialization function: plot the background of each frame
        def init():
            line.set_data([], [])
            return line,

        # animation function, called sequentially:
        def animate(i):
            y2 = y[:,DF*i]
            x2 = x[:,DF*i]
            #line.set_data(self.x,y2)
            line.set_data(x2,y2)
            return line,

        # Call the animator:
        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=nFrames, interval=50, blit=False, repeat=True)

        plt.show()
                
            
####################################################











