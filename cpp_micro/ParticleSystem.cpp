#include "ParticleSystem.hpp"

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
    constexpr double helium_epsilon = 1.4e-23;
    constexpr double helium_sigma = 2.56e-10;
}

CParticleSystem::CParticleSystem( size_t cParticles , int dim ):
    m_cParticles( cParticles ),
    m_dim(dim),
    m_pos( cParticles , vtkVector3d(0.) ),
    m_vel( cParticles , vtkVector3d(0.) ),
    m_forces( cParticles , vtkVector3d(0.) ),
    m_mass( constants::helium_mass ),
    m_epsilon( constants::helium_epsilon ),
    m_sigma( constants::helium_sigma ),
    m_temperatur( 0.0 ),
    m_friction( 0. ),
    m_potential1(Zero_Potential<TCodiVec3d>),
    m_potential2(VanDerWaals_Potential<TCodiVec3d>),
    m_radius(20. * m_sigma ),
    m_quader( vtkVector3d(0.,0.,0.) , vtkVector3d(1.,1.,1.) ),
    m_dt( 0.1 ),
    m_time( 0. ),
    m_boundaryConditions( EBoundaryConditions::FreeBoundary ),
    m_solver(ESolver::Verlet),
    m_forceUpdateMethod(EForceUpdateMethod::Naiv),
    m_neighbourQuaderSize( std::numeric_limits<double>::infinity() ),
    m_expDecayOfTemperatur(0.),
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

void CParticleSystem::insertPosition(size_t i)
{
    size_t nx = static_cast<size_t>((m_pos[i][0] - m_quader.first[0]) / m_neighbourQuaderSize[0]);
    size_t ny = static_cast<size_t>((m_pos[i][1] - m_quader.first[1]) / m_neighbourQuaderSize[1]);
    size_t nz = 0;

    if( m_dim == 3 )
        nz = static_cast<size_t>((m_pos[i][2] - m_quader.first[2]) / m_neighbourQuaderSize[2]);

    auto& vec = m_neighbours;
    if( nx < vec.size() && ny < vec[nx].size() && nz < vec[nx][ny].size() )
        m_neighbours[nx][ny][nz].push_back( i );
}

void CParticleSystem::initNeighbourlist()
{
    using std::vector;
    using std::list;

    vtkVector3d diag = getDomainDiagonal();
    size_t nx = diag[0] / m_neighbourQuaderSize[0];
    size_t ny = diag[1] / m_neighbourQuaderSize[1];
    size_t nz = 1;

    if( m_dim == 3 )
        nz = diag[2] / m_neighbourQuaderSize[2];

    // wtf :X sry for this lines of code...
    m_neighbours = vector<vector<vector<list<size_t>>>>( nx ,
                          vector<vector<list<size_t>>>( ny ,
                                 vector<list<size_t>>( nz , list<size_t>() ) ) );

    for( size_t i = 0 ; i < m_pos.size() ; ++i )
        insertPosition(i);
}

void CParticleSystem::updateNeighbourlist()
{
    using std::vector;
    using std::list;

    size_t counter = 0;

    vector<vector<vector<list<size_t>>>>& vecXYZ = m_neighbours;
    for( size_t ix = 0 ; ix < vecXYZ.size() ; ++ix )
    {
        vector<vector<list<size_t>>>& vecYZ = vecXYZ[ix];
        for( size_t iy = 0 ; iy < vecYZ.size() ; ++iy )
        {


            vector<list<size_t>>& vecZ = vecYZ[iy];
            for( size_t iz = 0 ; iz < vecZ.size() ; ++iz )
            {
                list<size_t>& l = vecZ[iz];
                auto it = l.begin();
                while( it != l.end() )
                {
                    // check if the particle has moved to an other domain
                    ++counter;
                    size_t nx = static_cast<size_t>((m_pos[*it][0] - m_quader.first[0]) / m_neighbourQuaderSize[0]);
                    size_t ny = static_cast<size_t>((m_pos[*it][1] - m_quader.first[1]) / m_neighbourQuaderSize[1]);

                    size_t nz = 0;
                    if( m_dim == 3)
                        nz = static_cast<size_t>((m_pos[*it][2] - m_quader.first[2]) / m_neighbourQuaderSize[2]);

                    bool hasRightPos = (nx == ix && ny == iy && nz == iz);

                    if( !hasRightPos && nx < vecXYZ.size() && ny < vecYZ.size() && nz < vecZ.size() )
                    {
                        m_neighbours[nx][ny][nz].push_back( *it );
                        l.erase( it++ );
                    }
                    else
                        ++it;
                }
            }
        }
    }

    if( counter < m_cParticles )
        initNeighbourlist();
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

    m_potential1( in_x , out_y , *this );

    vtkVector3d grad(0.);
    for( int i = 0 ; i < m_dim ; ++i)
        grad[i] = out_y[0].getGradient()[i];
    return grad;
}


vtkVector3d CParticleSystem::getPotential2Gradient(const vtkVector3d& x1 ,const vtkVector3d x2) const
{
    vtkVector3d diagonal = m_quader.second - m_quader.first;


    codi::RealForwardVec<3> in_x[3];
    codi::RealForwardVec<3> out_y[1];

    in_x[0] = x2.GetX() - x1.GetX();
    in_x[1] = x2.GetY() - x1.GetY();
    in_x[2] = x2.GetZ() - x1.GetZ();

    if( m_boundaryConditions == EBoundaryConditions::PeriodicBoundary )
    {
        for( int i = 0 ; i < m_dim ; ++i )
        {
            if( in_x[i] >  diagonal[i] / 2.)
                in_x[i] = in_x[i] - diagonal[i];
            if( in_x[i] < -diagonal[i] / 2. )
                in_x[i] = in_x[i] + diagonal[i];
        }

    }

    for( int i = 0 ; i < 3 ; ++i)
        in_x[i].gradient()[i] = 1.;

    m_potential2( in_x , out_y , *this );

    vtkVector3d grad(0.);
    for( int i = 0 ; i < m_dim ; ++i)
        grad[i] = out_y[0].getGradient()[i];
    return grad;
}


void CParticleSystem::updateForcesNaiv()
{
    static vtkVector3d grad;
    for( size_t i = 0 ; i < m_pos.size() ; ++i )
    {
        m_forces[i] = getPotential1Gradient( m_pos[i] );
        for( size_t j = 0 ; j < i ; ++j )
        {
            grad = getPotential2Gradient( m_pos[i] , m_pos[j] );
            m_forces[i] = m_forces[i] - grad;
            m_forces[j] = m_forces[j] + grad;
        }
    }
}

void CParticleSystem::updateForcesNeighbourlist()
{
    static vtkVector3d grad;

    updateNeighbourlist();

    // first of all we need to calculate the global potential
    // for each particle
    for( size_t i = 0 ; i < m_pos.size() ; ++i )
        m_forces[i] = getPotential1Gradient( m_pos[i] );

    // create the cells:

    using std::vector;
    using std::list;

    // we will not only search in
    vector<const list<size_t>*> nearCells( (m_dim==3)? 7 : 3 , nullptr );

    // reminder: m_neighbours is a 3d vector containting list<size_t>
    // where these lists hold the indices of these particles which are currently
    // near to each other.
    const vector<vector<vector<list<size_t>>>>& vecXYZ = m_neighbours;
    for( size_t ix = 0 ; ix < vecXYZ.size() ; ++ix )
    {
        const vector<vector<list<size_t>>>& vecYZ = vecXYZ[ix];
        for( size_t iy = 0 ; iy < vecYZ.size() ; ++iy )
        {
            const vector<list<size_t>>& vecZ = vecYZ[iy];
            for( size_t iz = 0 ; iz < vecZ.size() ; ++iz )
            {
                const list<size_t>& currentCell = vecZ[iz];
                // now we have a list of the points,
                // between which we to calculate the forces between them
                for( auto i = currentCell.begin() ; i != currentCell.end() ; ++i )
                {
                    for( auto j = currentCell.begin() ; j != i ; ++j )
                    {
                        // *i,*j are indices wrt m_pos,m_vel,m_forces
                        grad = getPotential2Gradient( m_pos[*i] , m_pos[*j] );
                        m_forces[*i] = m_forces[*i] - grad;
                        m_forces[*j] = m_forces[*j] + grad;
                    }

                    // now move to the neighbour cells

                    int ixNext = (int)(ix+1) % (int)vecXYZ.size() ;
                    int iyNext = (int)(iy+1) % (int)vecYZ.size() ;
                    int izNext = (int)(iz+1) % (int)vecZ.size() ;

                    nearCells[0] = &(vecXYZ[ixNext]   [iy]        [iz]);
                    nearCells[1] = &(vecXYZ[ix]       [iyNext]    [iz]);
                    nearCells[2] = &(vecXYZ[ixNext]   [iyNext]    [iz]);

                    if( ixNext == (int)ix )
                    {
                        nearCells[0] = nullptr;
                        nearCells[2] = nullptr;
                    }

                    if( iyNext == (int)iy )
                    {
                        nearCells[1] = nullptr;
                        nearCells[2] = nullptr;
                    }

                    if( m_dim == 3 )
                    {
                        nearCells[3] = &(vecXYZ[ix]       [iy]        [izNext]);
                        nearCells[4] = &(vecXYZ[ixNext]   [iy]        [izNext]);
                        nearCells[5] = &(vecXYZ[ix]       [iyNext]    [izNext]);
                        nearCells[6] = &(vecXYZ[ixNext]   [iyNext]    [izNext]);

                        if( ixNext == (int)ix )
                        {
                            nearCells[4] = nullptr;
                            nearCells[6] = nullptr;
                        }

                        if( iyNext == (int)iy )
                        {
                            nearCells[5] = nullptr;
                            nearCells[6] = nullptr;
                        }

                        if( izNext == (int)iz )
                        {
                            nearCells[3] = nullptr;
                            nearCells[4] = nullptr;
                            nearCells[5] = nullptr;
                            nearCells[6] = nullptr;
                        }

                    }


                    for( size_t iCell = 0 ; iCell < nearCells.size() ; ++iCell )
                    {
                        if( nearCells[iCell] != nullptr )
                            for( auto j = nearCells[iCell]->begin() ; j != nearCells[iCell]->end() ; ++j )
                            {
                                // *i,*j are indices wrt m_pos,m_vel,m_forces
                                grad = getPotential2Gradient( m_pos[*i] , m_pos[*j] );
                                m_forces[*i] = m_forces[*i] - grad;
                                m_forces[*j] = m_forces[*j] + grad;
                            }
                    }


                }
            }
        }
    }
}

void CParticleSystem::updateForces()
{
    switch( m_forceUpdateMethod)
    {
        case EForceUpdateMethod::Naiv:
            updateForcesNaiv();
        break;

        case EForceUpdateMethod::NeighbourList:
            updateForcesNeighbourlist();
        break;
    }
}

void CParticleSystem::velocityVerletStep( double dt )
{
    double sqrtdt = sqrt( dt/2. );

    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {
        m_vel[i] = m_vel[i] + 0.5 * dt / m_mass * m_forces[i];
        m_vel[i] = m_vel[i] + 2*constants::k_B*m_temperatur / m_mass * sqrtdt * rand_norm3d( 1. )  ;

        m_pos[i] = m_pos[i] + dt * m_vel[i] + dt*dt/(2. * m_mass) * m_forces[i];
    }

    updateBoundaryConditions();
    updateForces();

    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {

        m_vel[i] = m_vel[i] + 0.5 * dt / m_mass * m_forces[i];
        m_vel[i] = m_vel[i] + 2*constants::k_B*m_temperatur / m_mass * sqrtdt * rand_norm3d( 1. )  ;

    }

    updateBoundaryConditions();
}

void CParticleSystem::update(double dt)
{
    switch (m_solver) {
    case ESolver::ExplicitEuler:
        eulerStep(dt);
        break;

    case ESolver::Verlet:
        verletStep(dt);
        break;


    case ESolver::VelocityVerlet:
        velocityVerletStep(dt);
        break;


    default:
        verletStep(dt);
        break;
    }

    if( m_expDecayOfTemperatur != 0. )
        m_temperatur += m_dt * m_expDecayOfTemperatur * m_temperatur;
}

void CParticleSystem::update()
{
    update(m_dt);
}

void CParticleSystem::verletStep( double dt )
{
    double sqrtdt = sqrt( dt/2. );

    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {
        m_vel[i] = m_vel[i] + 0.5 * dt / m_mass * m_forces[i]
                            - m_friction * m_vel[i]
                            + 2*constants::k_B*m_temperatur / m_mass * sqrtdt * rand_norm3d( 1. )  ;

        m_pos[i] = m_pos[i] + dt * m_vel[i];
    }

    updateBoundaryConditions();
    updateForces();

    for( size_t i = 0 ; i < m_cParticles ; ++i )
    {

        m_vel[i] = m_vel[i] + 0.5 * dt / m_mass * m_forces[i]
                            - m_friction * m_vel[i]
                            + 2*constants::k_B*m_temperatur / m_mass * sqrtdt * rand_norm3d( 1. )  ;

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
        m_vel[i] = m_vel[i] + dt / m_mass * m_forces[i]
                            - m_friction * m_vel[i]
                            + 2*constants::k_B*m_temperatur / m_mass * sqrtdt * rand_norm3d( 1. )  ;

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

vtkVector3d CParticleSystem::getCamCenter() const
{
    vtkVector3d camCenter(0.,0.,0.);

    //for( auto& p : m_pos )
    //{
    //    camCenter = camCenter + 1/m_mass * p;
    //}
    return m_mass / m_cParticles * camCenter ;
}

void CParticleSystem::loadConfig(const CConfigManager &config)
{

    if( config.hasKey("boundary") )
    {
        const std::string cond = config.getValue<std::string>("boundary");

        if( cond == "free" )
            setBoundaryConditions( EBoundaryConditions::FreeBoundary );

        if( cond == "periodic")
            setBoundaryConditions( EBoundaryConditions::PeriodicBoundary );

        if( cond == "reflecting")
            setBoundaryConditions( EBoundaryConditions::ReflectingBoundary );

        if( cond == "sticky" )
            setBoundaryConditions( EBoundaryConditions::StickyBoundary );
    }

    if( config.hasKey("solver") )
    {
        std::string solver = config.getValue<std::string>("solver");

        if( solver == "explicitEuler" )
            setSolver( ESolver::ExplicitEuler );

        if( solver == "verlet")
            setSolver( ESolver::Verlet );

        if( solver == "velocityVerlet")
            setSolver( ESolver::VelocityVerlet );
    }

    if( config.hasKey("forceUpdate") )
    {
        std::string method = config.getValue<std::string>("forceUpdate");

        if( method == "naiv" )
            setForceUpdateMethod( EForceUpdateMethod::Naiv );

        if( method == "neighbourList" )
            setForceUpdateMethod( EForceUpdateMethod::NeighbourList );
    }

    if( config.hasKey<double>("dt") )
        setTimeStepSize( config.getValue<double>("dt") );

    if( config.hasKey<size_t>("N") )
        setNumberOfParticles( config.getValue<size_t>("N") );

    if( config.hasKey<double>("temperature") )
        setTemparature( config.getValue<double>("temperature") );

    if( config.hasKey<int>("dim") )
        setDimension( config.getValue<int>("dim") );

    if( config.hasKey<double>("friction") )
        setFriction( config.getValue<double>("friction") );

    if( config.hasKey<double>("exponetialDecayOfTemperatur") )
        m_expDecayOfTemperatur = config.getValue<double>("exponetialDecayOfTemperatur");

    if( config.hasKey<double>("genGrid") )
    {
        const double dist  = config.getValue<double>("genGrid");
        generateGriddedPositions( dist*pow(2,1/2.)*m_sigma );
    }

    if( config.hasKey<double>("gridBoundary") )
    {
        const double dist  = config.getValue<double>("gridBoundary");
        setDomain( dist*pow(2,1/2.)*m_sigma );
    }

    if( config.hasKey<int>("subcells") )
    {
        const int n  = config.getValue<int>("subcells");
        setNumberOfNeighbourCells( n, n, n );
    }
}

void CParticleSystem::setDomain(double outerDistance)
{
    vtkVector3d min( std::numeric_limits<double>::infinity() );
    vtkVector3d max( -std::numeric_limits<double>::infinity() );

    for( auto& p : m_pos )
        for( int d = 0 ; d < 3 ; ++d )
        {
            if( p[d] < min[d] ) min[d] = p[d];
            if( p[d] > max[d] ) max[d] = p[d];
        }

    setDomain(min - vtkVector3d(outerDistance) , max + vtkVector3d(outerDistance));
}

vtkVector3d CParticleSystem::getDomainDiagonal() const
{
    return m_quader.second - m_quader.first;
}

void CParticleSystem::setNumberOfNeighbourCells(int nx, int ny, int nz)
{
    vtkVector3d diagonal = m_quader.second - m_quader.first;
    m_neighbourQuaderSize.Set(  diagonal[0] / static_cast<double>(nx) ,
                                diagonal[1] / static_cast<double>(ny) ,
                                diagonal[2] / static_cast<double>(nz) );

    initNeighbourlist();
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

