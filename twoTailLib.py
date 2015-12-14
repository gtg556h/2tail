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
        self.ID = params['ID']
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
        self.A = np.ones_like(self.xi0)*self.A
        self.zetaN = np.ones_like(self.xi0)*self.zetaN
        self.zetaT = np.ones_like(self.xi0)*self.zetaT


    #####################################################

    def local2global(self, xi, eta, theta):

        x = xi * np.cos(theta + self.alpha) - eta * np.sin(theta + self.alpha)
        y = xi * np.sin(theta + self.alpha) + eta * np.cos(theta + self.alpha)
        
        return x, y

    ##########################################################

    # If operator c computes drag *on swimmer* in local basis, then
    # operator np.dot(l2g, np.dot(c, g2l)) computes drag *on swimmer* in global basis
    
    def l2g(self, theta):
        return np.array([[np.cos(theta + self.alpha), -np.sin(theta + self.alpha)], [np.sin(theta + self.alpha), np.cos(theta + self.alpha)]])

    def g2l(self, theta):
        return np.array([[np.cos(theta + self.alpha), np.sin(theta + self.alpha)], [-np.sin(theta + self.alpha), np.cos(theta + self.alpha)]])


    ##########################################################


    def global2local(self, x, y, theta):

        xi = x * np.cos(theta + self.alpha) + y * np.sin(theta + self.alpha)
        eta = -x * np.sin(theta + self.alpha) + y * np.cos(theta + self.alpha)

        return xi, eta
    

    ##########################################################

        
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
        self.w = np.zeros([self.xi.shape[0],self.t.shape[0]])

        for i in range(0, self.xi.shape[0]-1):
            a[i,i+1:self.m.shape[0]] = self.dx*np.arange(1,self.m.shape[0]-i,1)

        a[self.m.shape[0]-1,:] = np.ones([1,self.m.shape[0]])

        wBase = np.linalg.solve(a,self.m)

        for i in range(0,self.t.shape[0]):
            self.w[:,i] = np.multiply(wBase,self.mFunc[i])


    ######################################################
            

    def differentialOperator(self):
        nx = self.xi.shape[0]
        nt = self.t.shape[0]

        self.D4 = np.zeros([nx, nx])
        self.Dt = np.zeros_like(self.D4)
        self.a = np.zeros_like(self.D4)
        self.c = np.zeros([nx, nt])
        self.c2 = np.zeros([nx, nt])
        self.Dx = np.zeros([nx, nx])

        for i in range(2, nx-2):
            self.D4[i, i-2:i+3] = self.A[i] * np.array([1, -4, 6, -4, 1]) / self.dx**3

        for i in range(1, nx-1):
            self.Dx[i, i-1:i+2] = np.array([-1, 0, 1]) / (2*self.dx)
        self.Dx[0,0:2] = np.array([-1, 1]) / (self.dx)
        self.Dx[-1, nx-2:nx] = np.array([-1,1]) / (self.dx)

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

        hydrodynamicCouplingVelocity = np.zeros_like(self.xi0)
        v = self.calcLocalVelocity(theta, XDot, YDot, thetaDot, hydrodynamicCouplingVelocity)

        self.c2[2:nx-2, i] = self.c[2:nx-2, i] + np.multiply(self.zetaN[2:nx-2], self.eta[2:nx-2, i-1]) * self.dx/self.dt
        
        # Add contribution of local fluid flow in eta direction, caution with sign!!
        self.c2[2:nx-2, i] = self.c2[2:nx-2, i] + np.multiply(self.zetaN[2:nx-2], v[2:nx-2]) * self.dx

        self.eta[:,i] = np.linalg.solve(self.a, self.c2[:,i])
        self.xi[:,i] = self.xi0

        x, y = self.local2global(self.xi0, self.eta[:,i], theta)
        self.x[:,i] = X + x
        self.y[:,i] = Y + y

        self.f_x[i], self.f_y[i] = self.propulsionCalc(theta, v, i)


    ########################################################


    def calcLocalVelocity(self, theta, XDot, YDot, thetaDot, hydrodynamicCoupling):
        # For given element, calculate the local velocity
        # generated by translation and rotation of swimmer
        # As well as element-element coupling:

        v =  -self.xi0 * thetaDot - self.global2local(XDot,YDot,(self.alpha+theta))[1]
        # Add hydrodynamic coupling:
        v = v + hydrodynamicCoupling

        return v

    
    #########################################################

    
    def propulsionCalc(self, theta, vFlow, i):

        eta_x = np.dot(self.Dx, self.eta[:,i])
        eta_t = (self.eta[:, i] - self.eta[:, i-1]) / self.dt

        # Compute propulsive force *on swimmer* in xi-direction:
        f_xi = self.dx * np.sum((eta_t - vFlow) * eta_x * (self.zetaN - self.zetaT))

        f_x, f_y = self.local2global(f_xi, 0, theta)

        return f_x, f_y

    
    ##########################################################

    
    def computeDragConstants(self, theta):
        # Formulas compute drag on body provided associated velocities:

        # Returs 3x3 matrix c, such that c * {XDot, YDot, thetaDot} returns FX, FY, M

        # SUSPECTED BUG IN THIS FUNCTION!!!
        
        c = np.zeros([3,3])

        c_xi = self.dx * np.sum(self.zetaT)
        c_eta = self.dx * np.sum(self.zetaN)

        

        c[0,0] = self.dx * self.local2global(c_xi, 0, theta)
        
        c[0,0] = -self.dx * np.sum(self.zetaN) * np.abs(np.sin(self.alpha + theta))
        c[0,0] = c[0,0] - self.dx * np.sum(self.zetaT) * np.abs(np.cos(self.alpha + theta))
        # Assuming flagella radial to origin:
        c[2,0] = self.dx * np.sum(self.xi0 * self.zetaN) * np.sin(self.alpha + theta)
                                  
        c[1,1] = -self.dx * np.sum(self.zetaN) * np.abs(np.cos(self.alpha + theta))
        c[1,1] = c[1,1] - self.dx * np.sum(self.zetaT) * np.abs(np.sin(self.alpha + theta))
        # Assuming flagella radial to origin:
        c[2,1] = -self.dx * np.sum(self.xi0 * self.zetaN) * np.cos(self.alpha + theta)


        # Following *assumes* that flagella is radial WRT to origin.  Modify!!
        c[2,2] = -self.dx * np.sum(self.xi0 * self.xi0 * self.zetaN)
        c[0,2] = self.dx * np.sum(self.xi0 * self.zetaN) * np.sin(self.alpha + theta)
        c[1,2] = -self.dx * np.sum(self.xi0 * self.zetaN) * np.cos(self.alpha + theta)
                

        
        return c

    #############################################################
        
    def calcDrag(self, theta, thetaDot, XDot, YDot, timestep):
        # Compute drag forces and moments ON swimmer elements:

        # Generate array of drag coefficients:
        c = self.computeDragConstants(theta)

        # Calculate eta velocity, map to global coordinate system:
        etaDot = (self.eta[:,timestep] - self.eta[:, timestep-1]) / self.dt
        xDot, yDot = self.local2global(0,etaDot,theta)

        # Compose vector of system motion:
        XDotVec = np.array([XDot,YDot,thetaDot])

        pdb.set_trace()
        # Compute forces and moments on element:
        fx = -self.dx * np.sum(xDot * self.zetaN) 
        fx = fx + np.dot(c, XDotVec)[0,0]

        fy = -self.dx * np.sum(yDot * self.zetaN) 
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

        for j in structures:
            j.x0 = j.origin[0] + j.local2global(j.xi0, j.eta0, 0)[0]
            j.y0 = j.origin[1] + j.local2global(j.xi0, j.eta0, 0)[1]
            j.x = np.zeros([j.x0.shape[0], self.t.shape[0]])
            j.y = np.zeros([j.y0.shape[0], self.t.shape[0]])
            j.xi = np.zeros([j.xi0.shape[0], self.t.shape[0]])
            j.eta = np.zeros([j.eta0.shape[0], self.t.shape[0]])

            j.x[:,0] = j.x0
            j.y[:,0] = j.y0
            j.xi[:,0] = j.xi0
            j.eta[:,0] = j.eta0

            j.f_x = np.zeros(self.t.shape[0])
            j.f_y = np.zeros(self.t.shape[0])

        print('Calculating normal driving forces w(x,t)')
        for j in self.structures:
            j.processActuator(self.t, dt)

        self.theta = np.zeros_like(self.t); self.thetaDot = np.zeros_like(self.t)
        self.X = np.zeros_like(self.t); self.XDot = np.zeros_like(self.t)
        self.Y = np.zeros_like(self.t); self.YDot = np.zeros_like(self.t)


    ####################################################    
        

    def numSolve(self, hydrodynamicCoupling=1, propulsion=1):
        print('Solving for y(x,t)')

        self.hydrodynamicCoupling = hydrodynamicCoupling
        self.propulsion = propulsion

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
 
        while iteration < 10:# and (iteration < 2 or (theta - thetaOld > thetaThresh) or (X - XOld > XThresh) or (Y - YOld > YThresh)):
                                  

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
            # Add propulsive forces j.f_x, j.f_y:
            FX = FX + fx + j.f_x[timestep]
            FY = FY + fy + j.f_y[timestep]

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











