#include "CParticleSystem.hpp"

// CoDiPack includes
#include "codi.hpp"



CParticleSystem::CParticleSystem( size_t cParticles , int dim ):
    m_cParticles( cParticles ),
    m_dim(dim),
    m_pos( cParticles , vtkVector3d(0.) ),
    m_vel( cParticles , vtkVector3d(0.) ),
    m_forces( cParticles , vtkVector3d(0.) ),
    m_mass( 1. ),
    m_epsilon( 1.4e-23 ),
    m_sigma( 2.56e-10 ),
    m_potential1(Zero_Potential<TCodiVec3d>),
    m_potential2(VanDerWaals_Potential<TCodiVec3d>),
    m_radius(20. * m_sigma ),
    m_dt( 0.1 ),
    generator()
{
}

void CParticleSystem::generateRandomPositions(double variance)
{
    for( auto& p : m_pos )
        p = rand_norm3d(variance);
}


void CParticleSystem::generateRandomVelocities(double variance)
{
    for( auto& v : m_vel )
        v = rand_norm3d(variance);
}

void CParticleSystem::generateGriddedPositions(double distance , int nx , int ny , int nz )
{
    setNumberOfParticles( nx*ny*nz );
    size_t k = 0;

    for( int x = 0 ; x < nx ; ++x )
        for( int y = 0 ; y < ny ; ++y )
            for( int z = 0 ; z < nz ; ++z )
            {
                setPosition( k++ , vtkVector3d(x*distance,y*distance,z*distance) );
            }

}

void CParticleSystem::generateGriddedPositions(double distance )
{
    int gridWidth = floor( pow( static_cast<float>(m_cParticles) , 1/ (double)m_dim ) );
    generateGriddedPositions(distance,gridWidth,gridWidth, ((m_dim==3)? gridWidth : 1) );
}

vtkVector3d CParticleSystem::getPotential1Gradient( vtkVector3d x ) const
{
    codi::RealForwardVec<3> in_x[3];
    codi::RealForwardVec<3> out_y[1];
    in_x[0] = x.GetX();
    in_x[1] = x.GetY();
    in_x[2] = x.GetZ();

    for( int i = 0 ; i < 3 ; ++i)
        in_x[i].gradient()[i] = 1.;

    m_potential1( in_x , out_y , m_epsilon , m_sigma );

    vtkVector3d grad;
    for( int i = 0 ; i < 3 ; ++i)
        grad[i] = out_y[0].getGradient()[i];
    return grad;
}


vtkVector3d CParticleSystem::getPotential2Gradient(vtkVector3d r ) const
{
    codi::RealForwardVec<3> in_x[3];
    codi::RealForwardVec<3> out_y[1];
    in_x[0] = r.GetX();
    in_x[1] = r.GetY();
    in_x[2] = r.GetZ();

    for( int i = 0 ; i < 3 ; ++i)
        in_x[i].gradient()[i] = 1.;

    m_potential2( in_x , out_y , m_epsilon , m_sigma );

    vtkVector3d grad;
    for( int i = 0 ; i < 3 ; ++i)
        grad[i] = out_y[0].getGradient()[i];
    return grad;
}


void CParticleSystem::updateForces()
{
    static vtkVector3d grad;
    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {
        m_forces[i] = getPotential1Gradient( m_pos[i] );
        for( size_t j = 0 ; j < i ; ++j )
        {
            grad = getPotential2Gradient( m_pos[j] - m_pos[i] );
            m_forces[i] = m_forces[i] - grad;
            m_forces[j] = m_forces[j] + grad;
            //std::cout << grad[0] << std::endl;
        }
    }
}

void CParticleSystem::verletStep( double dt )
{
    double sqrtdt = sqrt( dt/2. );

    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {
        m_vel[i] = m_vel[i] + 0.5 * dt / m_mass * m_forces[i];
        //m_vel[i] = m_vel[i] + sqrtdt * rand_norm3d( m_sigma )  ;

        m_pos[i] = m_pos[i] + dt * m_vel[i];
    }

    updateForces();

    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {

        m_vel[i] = m_vel[i] + 0.5 * dt / m_mass * m_forces[i];
        //m_vel[i] = m_vel[i] + sqrtdt * rand_norm3d( m_sigma )  ;

    }

    updateBoundaryConditions();
}

void CParticleSystem::eulerStep(double dt)
{
    updateForces();
    double sqrtdt = sqrt( dt );
    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {
        m_pos[i] = m_pos[i] + dt * m_vel[i];
        m_vel[i] = m_vel[i] + dt / m_mass * m_forces[i];
        //m_vel[i] = m_vel[i] + sqrtdt * rand_norm3d()  ;

        updateBoundaryConditions();
    }


}

void CParticleSystem::updateBoundaryConditions()
{
    // boundary conditions
    for( size_t i = 0 ; i < m_cParticles ; ++i )
        for( int d = 0 ; d < m_dim ; ++d )
        {
            if( fabs( m_pos[i][d] ) >= m_radius  )
            {
                double sign = ((m_pos[i][d] >= 0.) ? 1. : -1.);
                m_pos[i][d] = 2*sign*m_radius - m_pos[i][d];
                m_vel[i][d] = -m_vel[i][d];
            }
        }
}

void CParticleSystem::setPotential1(CParticleSystem::TPotentialFunctional potential)
{
    m_potential1 = potential;
}

void CParticleSystem::setPotential2(CParticleSystem::TPotentialFunctional potential)
{

    m_potential2 = potential;
}


double CParticleSystem::rand_norm(double variance)
{
    std::normal_distribution<double> dist( 0. , variance );
    return dist( generator );
}

vtkVector3d CParticleSystem::rand_norm3d(double variance)
{
    return vtkVector3d( rand_norm() , rand_norm() , (( m_dim >= 3 ) ? rand_norm() : 0.) );
}
