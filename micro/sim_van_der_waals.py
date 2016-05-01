import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from math import sqrt

import time

def rand_unif( n,m ):
    return (2*np.random.rand(n,m) - 1)

class CParticleSystem:
    
    
    def __init__( self , dim ):
        self.dim = dim
        self.X = np.array([])
        self.DX = np.array([])
        
        self.plt_X = None
        self.plt_DX = None
        
        # gradient of the potential
        self.DV = None
        
        self.epsilon = 1.4e-23
        self.sigma = 2.56e-10
        #self.m = 1.66e-27
        self.m = 1
        
        self.particleCount = 0
        
        self.bPlotVelocities = True
        self.bBrownMotion = True
        
        self.boundaryConditions = 'open'
        
        self.domain = None        
        
    def V( self , r ):
        if r == 0.:
            return 0
        else:
            r6 = (self.sigma/r)**6
            return 4*self.epsilon* r6 * ( r6 - 1)
            
            
    def gradV( self , x ):
        #r = np.linalg.norm(x)
        r2 = 0
        for i in xrange(0,self.dim):
            r2 += x[i]*x[i]
        if r2 == 0.:
            return 0
        else:
            #return 4*self.epsilon* ( -12. * self.sigma**12 / r**14 + 6. * self.sigma**6 / (r**8) ) * x
            r2 = 1/r2
            return 4*self.epsilon* r2**(4) * ( -12. * self.sigma**12 * r2**3 + 6. * self.sigma**6 ) * x
    
    def rand( self , particleCount ):
        self.X = rand_unif( self.dim , particleCount ) * self.sigma
        #self.DX = rand_unif( self.dim , particleCount )
        self.DX = np.zeros_like( self.X )
        self.DV = np.zeros_like( self.DX )
        self.particleCount = particleCount
        self.domain = np.array( [[-self.sigma,self.sigma],[-self.sigma,self.sigma]] )
        
    def createGrid( self , particleCount , distance ):
        d = sqrt( particleCount ) / 2 * distance
        #  only valid for dim == 2
        self.X = np.reshape( np.mgrid[ -d:d:distance , -d:d:distance ] , ( self.dim , -1 ) )
        self.DX = np.zeros_like( self.X )
        self.DV = np.zeros_like( self.DX)
        self.particleCount = self.X.shape[1]
        self.domain = np.array( [[-d,d],[-d,d]] )  # (xmin,xmax ; ymin,ymax)

    def updateBoundaryConditions(self):
    
        for i in xrange(0,self.particleCount):
            for j in xrange( 0 , self.dim ):
                if( self.X[j,i] < self.domain[j,0] ):
                    self.X[j,i] = 2*self.domain[j,0] - self.X[j,i]
                    self.DX[j,i] *= -1
                if( self.X[j,i] > self.domain[j,1] ):
                    self.X[j,i] = 2*self.domain[j,1] - self.X[j,i]
                    self.DX[j,i] *= -1
                    
    def updateForce(self):
        for i in xrange(0,self.particleCount):
            self.DV[:,i] = 0            
            for j in xrange(0,i):
                gradV = self.gradV( self.X[:,i] - self.X[:,j]  )
                self.DV[:,i] += gradV
                self.DV[:,j] -= gradV
    
    def verletStep( self , dt ):
        
        self.DX += 0.5 * dt * (-self.DV)
        if self.bBrownMotion :
            self.DX +=  sqrt( 0.5*dt) * rand_unif( self.dim , self.particleCount ) * self.sigma
        self.X += dt* self.DX / self.m
        
                                        
        self.updateBoundaryConditions()
        self.updateForce()

        
        self.DX += 0.5 * dt * (-self.DV)
        if self.bBrownMotion :
            self.DX +=  sqrt( 0.5*dt) * rand_unif( self.dim , self.particleCount ) * self.sigma
    
            
    
    def eulerStep( self , dt ):

        self.updateForce()
        
        self.X  += dt*self.DX / self.m
        
        self.DX += -dt*self.DV
        if self.bBrownMotion :
            self.DX += sqrt(dt) * rand_unif( self.dim , self.particleCount ) * self.sigma
        #self.DX = self.DX - dt*self.DV
        #print( self.DX 
        
    #def totalEnergy( self ):
        
    
    def plot( self ):
        if( self.plt_X == None ):
            self.plt_X, = plt.plot( self.X[0,:] , self.X[1,:] , linewidth = 0 , marker = '.' )
            if self.bPlotVelocities:
                self.plt_DX = plt.quiver( self.X[0,:] , self.X[1,:] , self.DX[0,:] , self.DX[1,:] )
            
        else:
            self.plt_X.set_xdata( self.X[0,:] )
            self.plt_X.set_ydata( self.X[1,:] )
            
            if self.bPlotVelocities:
                self.plt_DX.remove()
                self.plt_DX = plt.quiver( self.X[0,:] , self.X[1,:] , self.DX[0,:] , self.DX[1,:] , linewidth = 0.01)
        
        
fig = plt.figure()

P = CParticleSystem(2)

#P.rand(30)
#P.X = P.X * 5

P.bPlotVelocities = False
P.bBrownMotion = False
P.epsilon = 1.4e-23
P.sigma = 2.56e-10
#self.m = 1.66e-27
P.m = 1
P.boundaryConditions = 'reflecting'

P.createGrid( 200 , 8 * P.sigma )
P.DX = rand_unif( 2 , P.X.shape[1] ) *10*P.sigma
P.plot()

dt = 0.01

def update_plot(dt):
    #P.eulerStep( dt )
    P.verletStep( dt ) # performance: O(particleCount**2)
    P.plot()

anim = animation.FuncAnimation( fig , update_plot , [dt] , interval = 10 , blit=False )
