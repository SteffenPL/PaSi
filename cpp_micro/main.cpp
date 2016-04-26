#include <vtkVector.h>
#include <vtkMath.h>

// std
#include <random>
#include <vector>
#include <functional>
#include <sstream>
#include <algorithm>

// own includes
#include "CParticleSystem.hpp"
#include "CParticleRenderer.hpp"

template< typename Real >
void Central_Potential(const Real* x , Real* y , double eps , double sigma )
{
    y[0] = -x[0]*x[0] - x[1]*x[1] - x[2]*x[2];
}

template< typename Real >
void Quader_Potential(const Real* x , Real* y , double eps , double sigma )
{
    y[0] = -x[0]*x[0] - x[1]*x[1] - x[2]*x[2];
    if( y[0] > -0.1 * sigma )
        y[0] = 0.;
}

template< typename Real >
void Gravity_Potential(const Real* x , Real* y , double eps , double sigma )
{
    y[0] = -x[2]*1e-14;
}

int main(int argc , char ** argv)
{
    int N = 25;
    double dt = 0.01;


    CParticleSystem particles( N );
    particles.setDimension( 3 );

    particles.parseParameters(argc,argv);

    particles.setPotential1( Zero_Potential<CParticleSystem::TCodiVec3d> );
    particles.setPotential2( VanDerWaals_Potential<CParticleSystem::TCodiVec3d> );


    //particles.generateRandomPositions( particles.getSigma()*particles.getSigma() );

    particles.generateGriddedPositions( 2*particles.getSigma()*sqrt(2) );
    particles.setDomain( 10. * particles.getSigma() );
    particles.setRadius( 1000. * particles.getSigma() );



    //particles.generateRandomVelocities(1);

    CParticleRenderer renderer;

    renderer.init( &particles );
    renderer.start();


    return EXIT_SUCCESS;
}
