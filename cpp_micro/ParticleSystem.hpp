#ifndef CPARTICLESYSTEM_HPP
#define CPARTICLESYSTEM_HPP

// vtk includes
#include <vtkVector.h>
#include <vtkVectorOperators.h>
#include <vtkMath.h>

// eigen includes
#include <eigen3/Eigen/Sparse>

// std includes
#include <random>
#include <vector>
#include <functional>
#include <list>
#include <numeric>

// CoDiPack
#include "codi.hpp"

#include "ConfigManager.hpp"

class CParticleSystem;



template< typename Real >
void Zero_Potential(const Real* x , Real* y , const CParticleSystem& )
{
    y[0] = 0.;
}


/*template< typename Real , double R >
void Bond2_Potential(const Real* x , Real* y , const CParticleSystem& p )
{
    x[0] = x[3] - x[0];
    x[1] = x[4] - x[1];
    x[2] = x[5] - x[2];

    Real d2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    y[0] = (d2-R)*(d2-R);
}*/


struct SAtom
{
    double mass;
};

struct SMolecule
{

    using TCodiVec3d = codi::RealForwardVec<3>;

    std::vector<size_t>         atoms;
    Eigen::SparseMatrix<bool>   bonds;

    using TPotentialFunctional = std::function<void( const TCodiVec3d* x , TCodiVec3d* y , const CParticleSystem& )>;

    TPotentialFunctional pot2;  // potential between two atoms
    TPotentialFunctional pot3;  // potential between three atoms
    TPotentialFunctional pot4;  // potential between four atoms


public:

    void loadFromFile( const std::string& filename )
    {

    }

    void initPotentials()
    {
        //pot2 = Bond2_Potential<TCodiVec3d,1.>;
        //pot3 = Bond3_Potential;
        //pot4 = Bond4_Potential;

    }
};

class CParticleSystem
{
public:

    // used typedefs
    using TQuader = std::pair<vtkVector3d,vtkVector3d>;
    using TNeighbourList = std::vector< std::vector< std::vector< std::list< size_t >>>>;
    using TCodiVec3d = SMolecule::TCodiVec3d;
    using TPotentialFunctional  = std::function<void(const TCodiVec3d* x , TCodiVec3d* y , const CParticleSystem&)>;

    enum EBoundaryConditions
    {
        FreeBoundary,
        PeriodicBoundary,
        ReflectingBoundary,
        StickyBoundary
    };

    enum ESolver
    {
        ExplicitEuler,
        Verlet,
        VelocityVerlet,
    };

    enum EForceUpdateMethod
    {
        Naiv,
        NeighbourList,
    };

    using TBond = std::pair<size_t,size_t>;

public:
    CParticleSystem(size_t cParticles , int dim = 3);

    void generateRandomPositions( double variance = 1.);
    void generateRandomVelocities(double variance = 1.);

    void generateGriddedPositions(double distance );
    void generateGriddedPositions(double distance , size_t nx , size_t ny , size_t nz );
    void generateGriddedPositions(size_t nx , size_t ny , size_t nz = 1);


    void initNeighbourlist();
    void updateNeighbourlist();
    void insertPosition( size_t i );

    void updateForcesNaiv();
    void updateForcesNeighbourlist();


    void updateForces();

    vtkVector3d getPotential1Gradient( vtkVector3d x ) const;
    vtkVector3d getPotential2Gradient(const vtkVector3d& x1 ,const vtkVector3d x2) const;

    void eulerStep( double dt );
    void verletStep(double dt);
    void velocityVerletStep(double dt);


    void update( double dt );

    void update();


    void updateBoundaryReflecting();
    void updateBoundaryPeriodic();
    void updateBoundarySticky();
    void updateBoundaryFree(){}

    void updateBoundaryConditions();

    vtkVector3d getCamCenter() const;

    // access functions

    inline const vtkVector3d& getPosition( size_t i ) const;
    inline void               setPosition( size_t i , const vtkVector3d& p );

    void addMolecule( std::vector<vtkVector3d> pos , Eigen::SparseMatrix<bool> bonds )
    {
        auto firstAtom = m_pos.size();

        m_pos.insert( m_pos.end() , pos.begin() , pos.end() );

        SMolecule molecule;
        molecule.atoms.resize(pos.size());
        std::iota( molecule.atoms.begin() , molecule.atoms.end() , firstAtom );

        molecule.bonds = bonds;

        m_molecules.push_back( molecule );
    }

    inline double getEpsilon() const;
    inline void   setEpsilon( double epsilon );

    inline double getSigma() const;
    inline void   setSigma( double sigma );

    inline double getTemperature() const;
    inline void   setTemparature(double temperatur);

    inline double getFriction() const;
    inline void   setFriction( double friction );

    inline int    getDimension() const;
    inline void   setDimension( int dim );

    inline size_t getNumberOfParticles() const;
    inline void   setNumberOfParticles( size_t new_size );

    inline double getTimeStepSize() const;
    inline void   setTimeStepSize( double dt );

    inline double getTime() const;
    inline void   setTime( double time = 0. );

    inline double getRadius() const;
    inline void   setRadius( double radius );

    void loadConfig(const CConfigManager &config );

    // the domain is given by two vectors containing the minimal (x,y,z) coords [minEdge]
    // and the maximal (x,y,z) coords [maxEdge]
    inline const TQuader& getDomain() const;
    inline void          setDomain( vtkVector3d minEdge , vtkVector3d maxEdge );

    // set a domain with the given outerDistance
    void setDomain( double outerDistance );
    vtkVector3d getDomainDiagonal() const;


    // set number of subcells
    void setNumberOfNeighbourCells( int nx , int ny , int nz );

    // Potientials are of the form void(const TCodiVec3d* x , TCodiVec3d* y , double eps , double sigma )

    // Potential depending on the space position of a particle
    void setPotential1( TPotentialFunctional potential );

    // Potential between two particles, depending on the
    // distance of two particles
    void setPotential2( TPotentialFunctional potential );

    inline void setBoundaryConditions( EBoundaryConditions cond );
    inline void setSolver( ESolver solver );
    inline void setForceUpdateMethod( EForceUpdateMethod method );


private:
    size_t m_cParticles;
    int    m_dim;

    // Positions
    std::vector<vtkVector3d>            m_pos;
    std::vector<SMolecule>    m_molecules;

    // Velocities
    std::vector<vtkVector3d> m_vel;

    // Velocities
    std::vector<vtkVector3d> m_forces;

    // Mass
    double m_mass;

    // atom constants
    double m_epsilon;
    double m_sigma;

    // temperatur
    double m_temperatur;

    // friction
    double m_friction;

    // Potential depending on the space position of a particle
    TPotentialFunctional m_potential1;

    // Potential between two particles, depending on the
    // distance of two particles
    TPotentialFunctional m_potential2;

    // space area
    // radial
    double m_radius;

    // quader
    std::pair<vtkVector3d,vtkVector3d> m_quader;

    // simulations time step
    double m_dt;

    // elapsed time
    double m_time;

    // boundary conditions
    EBoundaryConditions m_boundaryConditions;

    // solver
    ESolver m_solver;

    // solver
    EForceUpdateMethod m_forceUpdateMethod;


    // neighbour list
    TNeighbourList m_neighbours;
    vtkVector3d    m_neighbourQuaderSize;

    // exponential decay
    float m_expDecayOfTemperatur;

public:
    // Helper functions for random number/vector generation
    std::mt19937 generator;
    inline double rand_norm( double variance = 1. );
    inline vtkVector3d rand_norm3d( double variance = 1. );
};



template< typename Real >
void VanDerWaals_Potential(const Real* x , Real* y , const CParticleSystem& p )
{
    Real d2 = p.getSigma()*p.getSigma() / (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    y[0] = 4 * p.getEpsilon() * codi::pow(d2,3) * (codi::pow(d2,3) - 1.);
}


//
// inline implementations:
//

inline const vtkVector3d& CParticleSystem::getPosition( size_t i ) const
{    return m_pos[i];}

inline void CParticleSystem::setPosition( size_t i , const vtkVector3d& p )
{    m_pos[i] = p;}

inline double CParticleSystem::getEpsilon() const
{    return m_epsilon;}

inline void CParticleSystem::setEpsilon( double epsilon )
{    m_epsilon = epsilon;}

inline double CParticleSystem::getSigma() const
{    return m_sigma;}

inline void CParticleSystem::setSigma( double sigma )
{    m_sigma = sigma;}


inline double CParticleSystem::getTemperature() const
{   return m_temperatur;}

inline void   CParticleSystem::setTemparature( double temperatur )
{   m_temperatur = temperatur; }

inline size_t CParticleSystem::getNumberOfParticles() const
{    return m_cParticles;}


inline double CParticleSystem::getFriction() const
{    return m_friction;}

inline void   CParticleSystem::setFriction( double friction )
{   m_friction = friction; }

inline void CParticleSystem::setNumberOfParticles( size_t new_size )
{
    m_cParticles = new_size;
    m_pos.resize( new_size );
    m_vel.resize( new_size );
    m_forces.resize( new_size );
}


inline int CParticleSystem::getDimension() const
{    return m_dim;}

inline void CParticleSystem::setDimension( int dim )
{    m_dim = dim;}

inline double CParticleSystem::getTimeStepSize() const
{    return m_dt;}

inline void   CParticleSystem::setTimeStepSize( double dt )
{    m_dt = dt;}


inline double CParticleSystem::getRadius() const
{    return m_radius;}

inline void   CParticleSystem::setRadius( double radius )
{    m_radius = radius;}


inline const CParticleSystem::TQuader&   CParticleSystem::getDomain() const
{
    return m_quader;
}

inline void CParticleSystem::setDomain( vtkVector3d minEdge , vtkVector3d maxEdge )
{
    m_quader.first = minEdge;
    m_quader.second = maxEdge;

    if( !m_neighbours.empty() )
    {
        auto& vec = m_neighbours;
        setNumberOfNeighbourCells( vec.size() , vec[0].size() , vec[0][0].size() );
    }
}

inline double CParticleSystem::getTime() const
{    return m_time;}

inline void   CParticleSystem::setTime( double time )
{   m_time = time;}


inline void CParticleSystem::setBoundaryConditions( EBoundaryConditions cond )
{
    m_boundaryConditions = cond;
}

inline void CParticleSystem::setSolver( ESolver solver )
{
    m_solver = solver;
}


inline void CParticleSystem::setForceUpdateMethod( EForceUpdateMethod method )
{
    m_forceUpdateMethod = method;
}

#endif // CPARTICLESYSTEM_HPP
