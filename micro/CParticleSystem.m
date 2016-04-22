classdef CParticleSystem < handle
   
    % I will try to use the notations given in the lecture 'interaction
    % particle systems' 2016
    properties
        
        % dimension
        dim = 2;
        
        % Z denotes a matrix containing both: [positions ; velocities]
        % the dimension is 2*dim x N
        Z = [];

        
        % variable holding the current plot handle
        pltPos = [];
        pltVelocities = [];
    end
    
    
    
    methods
        
        
        function self = CParticleSystem( dim )
            self.dim = dim;
            self.Z = [];
        
        end
        
        
        % generate normal distributed particles (pos & normals)
        function self = rand( self, N )
           self.Z = rand( 2*self.dim , N );
        end
        
        
        function self = euler_step( self , dt )
           
            self.Z 
            
        end
        
        % plot the current state
        function self = plot(self)
            
            if( isempty(self.pltPos ) )
               
                % create plot for the first time
                
                if( self.dim == 2 )
                   self.pltPos = plot( self.Z(1,:) , self.Z(2,:) );
                   self.pltVelocities = quiver( self.Z(1,:) , self.Z(2,:) , ...
                                                self.Z(3,:) , self.Z(4,:) );
                end
                
                
                if( self.dim == 3 )
                    
                end
                
            else
                if( self.dim == 2 )
                    self.pltPos.XData = self.Z(1,:);
                    self.pltPos.YData = self.Z(2,:);
                    self.pltVelocities.XData = self.Z(1,:);
                    self.pltVelocities.YData = self.Z(2,:);
                    self.pltVelocities.UData = self.Z(3,:);
                    self.pltVelocities.VData = self.Z(4,:);
                end
                
                if( self.dim == 3 )
                    
                end
            end
            
        end
    end
    
end