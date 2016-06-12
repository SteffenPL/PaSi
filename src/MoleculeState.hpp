#ifndef CMOLECULE_STATE_HPP
#define CMOLECULE_STATE_HPP

// eigen includes
#include <Eigen/Dense>

//vtk includes
#include <vtkNew.h>

// own includes
#include "ForceField.hpp"

// std includes
#include <map>
#include <vector>
#include <set>
#include <string>

// foreward delcaration
class vtkMolecule;

struct SAtom
{
    Eigen::Vector3d  p;
    Eigen::Vector3d  v;
    int              atomicNumber;
    std::set< int >  neighbours_fwd;
    std::set< int >  neighbours;
};

class CMoleculeState
{

    using TVec3 = Eigen::Vector3d;

public:

    CMoleculeState();
    ~CMoleculeState();


    void loadMoleculeFromFile( const std::string& filename );
    void loadForceFieldFromFile( const std::string& filename );

    void deepCopyMolecule( vtkMolecule* pMolecule );
    vtkMolecule* getMolecule();

    void printInfo();

    void calculateForces();

    void update()
    {
        double dt = 0.001;
        double sqrtdt = sqrt( dt/2. );

        for( size_t i = 0 ; i < numberOfAtoms() ; ++i )
        {
            vel(i) = vel(i) + 0.5 * dt * m_forces[i];

            pos(i) = pos(i) + dt * vel(i) + dt*dt/(2.) * m_forces[i];
        }

        calculateForces();

        for( size_t i = 0 ; i < numberOfAtoms() ; ++i )
        {
            vel(i) = vel(i) + 0.5 * dt  * m_forces[i];
        }
    }

    size_t numberOfAtoms() const
    {
        return m_atoms.size();
    }

    inline Eigen::Vector3d& pos( size_t i )
    {
        return m_atoms[i].p;
    }

    inline Eigen::Vector3d& vel( size_t i )
    {
        return m_atoms[i].v;
    }

    inline int atomicNumber( int i )
    {
        return m_atoms[i].atomicNumber;
    }

protected:


private:

    CForceField     m_forceField;

    vtkNew<vtkMolecule>   m_molecule;


    std::vector< SAtom >            m_atoms;
    std::vector< Eigen::Vector2i >    m_bonds;

    std::vector< Eigen::Vector3d > m_forces;



};



#endif // CMOLECULE_STATE_HPP
