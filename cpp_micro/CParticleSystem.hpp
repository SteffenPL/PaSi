#ifndef CPARTICLESYSTEM_HPP
#define CPARTICLESYSTEM_HPP

// vtk includes
#include <vtkVector.h>
#include <vtkVectorOperators.h>
#include <vtkMath.h>

// std includes
#include <random>
#include <vector>
#include <functional>

// CoDiPack
#include "codi.hpp"

template< typename Real >
void Zero_Potential(const Real* x , Real* y , double eps , double sigma )
{
    y[0] = 0.;
}

template< typename Real >
void VanDerWaals_Potential(const Real* x , Real* y , double eps , double sigma )
{
    Real d2 = sigma*sigma / (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    y[0] = 4 * eps * codi::pow(d2,3) * (codi::pow(d2,3) - 1.);
}


class CParticleSystem
{
public:
    CParticleSystem(size_t cParticles , int dim = 3);

    void generateRandomPositions( double variance = 1.);
    void generateRandomVelocities(double variance = 1.);

    void generateGriddedPositions(double distance );
    void generateGriddedPositions(double distance , int nx , int ny , int nz );

    void updateForces();

    vtkVector3d getPotential1Gradient( vtkVector3d x ) const;
    vtkVector3d getPotential2Gradient( vtkVector3d r ) const;

    void eulerStep( double dt );
    void verletStep(double dt);


    inline void update( double dt )
    {
        verletStep(dt);
    }

    inline void update() { update(m_dt); }

    void updateBoundaryConditions();


    // access functions

    inline const vtkVector3d& getPosition( size_t i ) const;
    inline void               setPosition( size_t i , const vtkVector3d& p );

    inline double getEpsilon() const;
    inline void   setEpsilon( double epsilon );

    inline double getSigma() const;
    inline void   setSigma( double sigma );

    inline int    getDimension() const;
    inline void   setDimension( int dim );

    inline size_t getNumberOfParticles() const;
    inline void   setNumberOfParticles( size_t new_size );

    inline double getTimeStepSize() const;
    inline void   setTimeStepSize( double dt );

    inline double getRadius() const;
    inline void   setRadius( double radius );

    using TCodiVec3d = codi::RealForwardVec<3>;
    using TPotentialFunctional  = std::function<void(const TCodiVec3d* x , TCodiVec3d* y , double , double)>;


    // Potential depending on the space position of a particle
    void setPotential1( TPotentialFunctional potential );

    // Potential between two particles, depending on the
    // distance of two particles
    void setPotential2( TPotentialFunctional potential );

private:
    size_t m_cParticles;
    int    m_dim;

    // Positions
    std::vector<vtkVector3d> m_pos;

    // Velocities
    std::vector<vtkVector3d> m_vel;

    // Velocities
    std::vector<vtkVector3d> m_forces;

    // Mass
    double m_mass;

    // atom constants
    double m_epsilon;
    double m_sigma;

    // Potential depending on the space position of a particle
    TPotentialFunctional m_potential1;

    // Potential between two particles, depending on the
    // distance of two particles
    TPotentialFunctional m_potential2;

    // space ares
    double m_radius;

    // simulations time step
    double m_dt;


public:
    // Helper functions for random number/vector generation
    std::mt19937 generator;
    inline double rand_norm( double variance = 1. );
    inline vtkVector3d rand_norm3d( double variance = 1. );
};

//
// inline implementations:
//

inline const vtkVector3d& CParticleSystem::getPosition( size_t i ) const
{
    return m_pos[i];
}

inline void CParticleSystem::setPosition( size_t i , const vtkVector3d& p )
{
    m_pos[i] = p;
}

inline double CParticleSystem::getEpsilon() const
{
    return m_epsilon;
}

inline void CParticleSystem::setEpsilon( double epsilon )
{
    m_epsilon = epsilon;
}

inline double CParticleSystem::getSigma() const
{
    return m_sigma;
}

inline void CParticleSystem::setSigma( double sigma )
{
    m_sigma = sigma;
}

inline size_t CParticleSystem::getNumberOfParticles() const
{
    return m_cParticles;
}

inline void CParticleSystem::setNumberOfParticles( size_t new_size )
{
    m_cParticles = new_size;
    m_pos.resize( new_size );
    m_vel.resize( new_size );
    m_forces.resize( new_size );
}


inline int CParticleSystem::getDimension() const
{
    return m_dim;
}

inline void CParticleSystem::setDimension( int dim )
{
    m_dim = dim;
}

inline double CParticleSystem::getTimeStepSize() const
{
    return m_dt;
}

inline void   CParticleSystem::setTimeStepSize( double dt )
{
    m_dt = dt;
}


inline double CParticleSystem::getRadius() const
{
    return m_radius;
}

inline void   CParticleSystem::setRadius( double radius )
{
    m_radius = radius;
}


#endif // CPARTICLESYSTEM_HPP
