#include <vtkVector.h>
#include <vtkMath.h>

// std
#include <random>
#include <vector>
#include <functional>
#include <sstream>

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

int main(int argc , char ** argv)
{
    int cParticles = 25;
    double dt = 0.01;

    switch( argc )
    {

        case 3:
        // input parameters
        {
            std::istringstream ss(argv[2]);
            if( !(ss >> dt) || dt < 0. )
            {
                std::cout << "Invalid second argument, should be the time step size ("
                          << dt << ")" << std::endl;
                dt = 0.01;
            }
        }

        case 2:
            // input parameters
            {
                std::istringstream ss(argv[1]);
                if( !(ss >> cParticles) )
                {
                    std::cout << "Invalid first argument, should be the number of particles! ("
                              << cParticles << ")" << std::endl;
                    cParticles = 25;
                }
            }
        break;

        default:
        break;
    }

    CParticleSystem particles( cParticles );
    particles.setDimension( 2 );
    particles.setTimeStepSize( dt );

    particles.setPotential1( Zero_Potential<CParticleSystem::TCodiVec3d> );
    particles.setPotential2( VanDerWaals_Potential<CParticleSystem::TCodiVec3d> );

    //particles.generateRandomPositions( 30.* particles.getSigma() );
    particles.generateGriddedPositions( particles.getSigma() * 2. );
    particles.setRadius( 20. * particles.getSigma());
    particles.generateRandomVelocities(particles.getSigma()*particles.getSigma());

    CParticleRenderer renderer;

    renderer.init( &particles );
    renderer.start();


    return EXIT_SUCCESS;
}
