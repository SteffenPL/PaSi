#include "CParticleSystem.hpp"

// CoDiPack includes
#include "codi.hpp"

// std includes
#include <sstream>


// Bolzmann constant
namespace constants
{

    constexpr double k_B = 1.3806504e-23;
    constexpr double dalton = 1.660539e-27; //kg
    constexpr double helium_mass = 4.002602*dalton;

}

CParticleSystem::CParticleSystem( size_t cParticles , int dim ):
    m_cParticles( cParticles ),
    m_dim(dim),
    m_pos( cParticles , vtkVector3d(0.) ),
    m_vel( cParticles , vtkVector3d(0.) ),
    m_forces( cParticles , vtkVector3d(0.) ),
    m_mass( constants::helium_mass ),
    m_epsilon( 1.4e-23 ),
    m_sigma( 2.56e-10 ),
    m_temperatur( 0.0001 ),
    m_potential1(Zero_Potential<TCodiVec3d>),
    m_potential2(VanDerWaals_Potential<TCodiVec3d>),
    m_radius(20. * m_sigma ),
    m_quader( vtkVector3d(0.,0.,0.) , vtkVector3d(1.,1.,1.) ),
    m_dt( 0.1 ),
    m_time( 0. ),
    m_boundaryConditions( EBoundaryConditions::FreeBoundary ),
    generator()
{
    double quaderSize = 2*m_sigma*pow( static_cast<double>(cParticles)
                                       , 1./ static_cast<double>(m_dim ));
    m_quader.second = m_quader.second * quaderSize;
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

void CParticleSystem::generateGriddedPositions(double distance , size_t nx , size_t ny , size_t nz )
{
    setNumberOfParticles( nx*ny*nz );
    size_t k = 0;

    for( size_t x = 0 ; x < nx ; ++x )
        for( size_t y = 0 ; y < ny ; ++y )
            for( size_t z = 0 ; z < nz ; ++z )
            {
                setPosition( k++ , vtkVector3d(x*distance,y*distance,z*distance) );
            }

}

void CParticleSystem::generateGriddedPositions(size_t nx, size_t ny, size_t nz)
{
    vtkVector3d diagonal = m_quader.second - m_quader.first;
    const double dx = diagonal[0] / static_cast<double>(nx+1);
    const double dy = diagonal[1] / static_cast<double>(ny+1);
          double dz = diagonal[2] / static_cast<double>(nz+1);

    if( m_dim < 3 )
    {
        nz = 1;
        dz = 0;
    }


    setNumberOfParticles( nx*ny*nz );
    size_t k = 0;

    for( size_t x = 0; x < nx ; ++x )
        for( size_t y = 0; y < ny ; ++y )
            for( size_t z = 0; z < nz ; ++z )
            {
                m_pos[k++] = m_quader.first + vtkVector3d(   0.5*dx + x*dx
                                                           , 0.5*dy + y*dy
                                                           , 0.5*dz + z*dz );
            }

}

void CParticleSystem::generateGriddedPositions(double distance )
{
    int gridWidth = floor( pow( static_cast<float>(m_cParticles) , 1/ static_cast<double>(m_dim) ) );
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
    static codi::RealForwardVec<3> in_x[3];
    static codi::RealForwardVec<3> out_y[1];
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


void CParticleSystem::updateForcesNaiv()
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

void CParticleSystem::updateForcesNeighbourlist()
{

}

void CParticleSystem::verletStep( double dt )
{
    double sqrtdt = sqrt( dt/2. );

    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {
        m_vel[i] = m_vel[i] + 0.5 * dt / m_mass * m_forces[i];
        m_vel[i] = m_vel[i] + 2*constants::k_B*m_temperatur / m_mass * sqrtdt * rand_norm3d( 1. )  ;

        m_pos[i] = m_pos[i] + dt * m_vel[i];
    }

    updateForces();

    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {

        m_vel[i] = m_vel[i] + 0.5 * dt / m_mass * m_forces[i];
        m_vel[i] = m_vel[i] + 2*constants::k_B*m_temperatur / m_mass * sqrtdt * rand_norm3d( 1. )  ;

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
        m_vel[i] = m_vel[i] + 2*constants::k_B*m_temperatur / m_mass * sqrtdt * rand_norm3d( 1. )  ;

    }

    updateBoundaryConditions();
}

void CParticleSystem::updateBoundaryReflecting()
{

    // note: for the stochasic part this isn't a practical abroach
    // better set the brownian motion to zero if in conflict with the boundaries
    for( size_t i = 0 ; i < m_cParticles ; ++i )
        for( int d = 0 ; d < m_dim ; ++d )
        {
            if( m_pos[i][d] < m_quader.first[d]  )
            {
                m_pos[i][d] = 2*m_quader.first[d] - m_pos[i][d];
                m_vel[i][d] = -m_vel[i][d];
            }

            if( m_pos[i][d] >= m_quader.second[d]  )
            {
                m_pos[i][d] = 2*m_quader.second[d] - m_pos[i][d];
                m_vel[i][d] = -m_vel[i][d];
            }
        }
}



void CParticleSystem::updateBoundaryPeriodic()
{
    vtkVector3d diagonal = m_quader.second - m_quader.first;

    for( size_t i = 0 ; i < m_cParticles ; ++i )
        for( int d = 0 ; d < m_dim ; ++d )
        {
            if( m_pos[i][d] < m_quader.first[d] )
                m_pos[i][d] += diagonal[d];

            if( m_pos[i][d] > m_quader.second[d] )
                m_pos[i][d] -= diagonal[d];
        }
}

void CParticleSystem::updateBoundarySticky()
{
    for( size_t i = 0 ; i < m_cParticles ; ++i )
        for( int d = 0 ; d < m_dim ; ++d )
        {
            if( m_pos[i][d] < m_quader.first[d] )
            {
                m_pos[i][d] = m_quader.first[d];
                m_vel[i][d] = 0.;
            }

            if( m_pos[i][d] > m_quader.second[d] )
            {
                m_pos[i][d] = m_quader.second[d];
                m_vel[i][d] = 0.;
            }
        }
}

void CParticleSystem::updateBoundaryConditions()
{
    switch( m_boundaryConditions )
    {
        case EBoundaryConditions::FreeBoundary:
            updateBoundaryFree();
        break;

        case EBoundaryConditions::PeriodicBoundary:
            updateBoundaryPeriodic();
        break;

        case EBoundaryConditions::ReflectingBoundary:
            updateBoundaryReflecting();
        break;

        case EBoundaryConditions::StickyBoundary:
            updateBoundarySticky();
        break;
    }
}

namespace
{

template<typename T>
void parseParameter( const std::string& key , char* argv , const std::string& name
                  , void(CParticleSystem::*fPtr)(T) , CParticleSystem* obj)
{
    if( key == name )
    {
        std::stringstream ss(argv);
        T value;

        if( !(ss>>value) )
        {
            std::cout << "Couldn't parse the value for key " << key << "." << std::endl;
        }
        else
        {
            std::cout << key << " = " << value << std::endl;
            (obj->*fPtr)(value);
        }
    }
}


void parseParameter( const std::string& key , const std::string& name
                     , void(CParticleSystem::*fPtr)() , CParticleSystem* obj)
{
    if( key == name )
    {
        (obj->*fPtr)();
        std::cout << key << std::endl;
    }
}

}

void CParticleSystem::parseParameters(int argc, char **argv)
{
    // parse arguments
    for( int i = 1 ; i < argc ; i += 1 )
    {
        std::istringstream ss( argv[i] );

        std::string key;

        ss.ignore(2);
        if( !(ss >> key) )
            std::cout << "Couldn't parse parameter argument " << i << ".\n";
        else
        {
            parseParameter( key, "freeBoundary" , &CParticleSystem::setFreeBoundary , this );
            parseParameter( key, "periodicBoundary" , &CParticleSystem::setPeriodicBoundary , this );
            parseParameter( key, "reflectingBoundary" , &CParticleSystem::setReflectingBoundary , this );
            parseParameter( key, "stickyBoundary" , &CParticleSystem::setStickyBoundary , this );


            if( i+1 < argc )
            {
                parseParameter(key,argv[i+1],"stepSize"
                              , &CParticleSystem::setTimeStepSize , this );

                parseParameter(key,argv[i+1],"numberOfParticles"
                              , &CParticleSystem::setNumberOfParticles , this );

                parseParameter(key,argv[i+1],"temperature"
                              , &CParticleSystem::setTemparature , this );

                parseParameter(key,argv[i+1],"dimension"
                              , &CParticleSystem::setDimension , this );
                ++i;
            }

        }


    }
}

void CParticleSystem::setDomain(double outerDistance)
{
    vtkVector3d min( std::numeric_limits<double>::infinity() );
    vtkVector3d max = -1. * min;

    for( auto& p : m_pos )
        for( int d = 0 ; d < m_dim ; ++d )
        {
            if( p[d] < min[d] ) min[d] = p[d];
            if( p[d] > max[d] ) max[d] = p[d];
        }

    setDomain(min - vtkVector3d(outerDistance) , max + vtkVector3d(outerDistance));
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
