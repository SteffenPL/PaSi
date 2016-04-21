import numpy as np
import matplotlib.pyplot as plt


#
#%x = np.random.rand(20)
#%y = 1e7*np.random.rand(20)
#
#fig, ax = plt.subplots()
#ax.fmt_ydata = millions
#plt.plot(x, y, 'o')
#
#plt.show()


class CParticleSystem:
    
    
    def __init__( self , dim ):
        self.dim = dim
        self.Z = np.array([])
        self.plt_X = None
        self.plt_Y = None
        
    
    def rand( self , particleCount ):
        self.Z = np.random.rand( 2* self.dim , particleCount )
    
    def plot( self ):
        if( self.plt_X == None ):
            self.plt_X = plt.plot( self.Z[0,:] , self.Z[1,:] , linewidth = 0 , marker = '*' )
            self.plt_Y = plt.quiver( self.Z[0,:] , self.Z[1,:] , self.Z[2,:] , self.Z[3,:] )
            
            plt.show();
        
        else:
            self.plt_X.set_xdata = self.Z[0,:]
            self.plt_X.set_ydata = self.Z[1,:]
            
            self.plt_V.set_xdata = self.Z[0,:]
            self.plt_V.set_ydata = self.Z[1,:]
            self.plt_V.set_udata = self.Z[2,:]
            self.plt_V.set_vdata = self.Z[3,:]
        
        
plt.clf()

P = CParticleSystem(2)

P.rand(10)

P.plot()
    