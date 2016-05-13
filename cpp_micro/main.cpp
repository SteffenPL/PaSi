#include <vtkVector.h>
#include <vtkMath.h>

// std
#include <random>
#include <vector>
#include <functional>
#include <sstream>
#include <algorithm>

// own includes
#include "ParticleSystem.hpp"
#include "ParticleRenderer.hpp"

#include "ConfigManager.hpp"


template< typename Real >
void Central_Potential(const Real* x , Real* y , const CParticleSystem& )
{
    y[0] = -x[0]*x[0] - x[1]*x[1] - x[2]*x[2];
}

template< typename Real >
void Quader_Potential(const Real* x , Real* y , const CParticleSystem& p )
{
    y[0] = -x[0]*x[0] - x[1]*x[1] - x[2]*x[2];
    if( y[0] > -0.1 * p.getSigma() )
        y[0] = 0.;
}

template< typename Real >
void Gravity_Potential(const Real* x , Real* y , const CParticleSystem& )
{
    y[0] = -x[1]*1e-15;
}

int main(int argc , char ** argv)
{
    CConfigManager config;
    config.parse(argc,argv);

    if( config.hasKey("configFile"))
    {
        ifstream file( config.getValue("configFile") );
        if( file.is_open() )
            config.parse(file);
    }

    std::cout << config;

    CParticleSystem particles( 0 );

    particles.loadConfig( config );

    particles.setPotential1( Gravity_Potential<CParticleSystem::TCodiVec3d> );
    particles.setPotential2( VanDerWaals_Potential<CParticleSystem::TCodiVec3d> );


    particles.setRadius( 1000. * particles.getSigma() );

    //particles.setNumberOfNeighbourCells( 6,6,6 );


    //particles.generateRandomVelocities(1);

    CParticleRenderer renderer;

    renderer.init( &particles );
    renderer.start();


    return EXIT_SUCCESS;
}
