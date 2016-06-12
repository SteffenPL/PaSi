#include "MoleculeState.hpp"


// vtk includes
#include <vtkPeriodicTable.h>
#include <vtkStringArray.h>
#include <vtkFloatArray.h>
#include <vtkPDBReader.h>
#include <vtkBlueObeliskData.h>

#include <vtkMolecule.h>


CMoleculeState::CMoleculeState()
{
    // load atomic names
    vtkNew<vtkPeriodicTable> table;
    auto atomicNames   = table->GetBlueObeliskData()->GetNames();
    auto atomicSymbols = table->GetBlueObeliskData()->GetSymbols();

    for( size_t i = 0 ; i < atomicNames->GetSize() ; ++i )
        m_forceField.setAtomicNumber( atomicSymbols->GetValue(i) ,
                                      atomicNames->GetValue(i) ,
                                      i );

}

CMoleculeState::~CMoleculeState()
{

}

void CMoleculeState::loadMoleculeFromFile(const std::string &filename)
{
    // since vtk should not be a strong dependency, we here import all
    // data given by vtk to our own data structure,
    // this read function should be the only dependency on vtk
    // (and the constructor of cause).

    vtkNew<vtkPDBReader> pdb;
    pdb->SetFileName( filename.c_str()   );
    pdb->SetHBScale(1.0);
    pdb->SetBScale(1.0);
    pdb->Update();

    auto mol = vtkMolecule::SafeDownCast( pdb->GetOutputDataObject(1) );

    vtkNew<vtkPeriodicTable> table;
    auto atomicSymbols = table->GetBlueObeliskData()->GetSymbols();

    // offsets
    const size_t oAtoms = m_atoms.size();
    const size_t oBonds = m_bonds.size();

    // new elements
    const size_t cAtoms = (size_t)mol->GetNumberOfAtoms();
    const size_t cBonds = (size_t)mol->GetNumberOfBonds();

    m_atoms.resize(  oAtoms + cAtoms );
    m_forces.resize( oAtoms + cAtoms );
    m_bonds.resize( oBonds + cBonds );


    // init position, velocities , bonds etc
    {
        auto it = m_atoms.begin() + oAtoms;
        auto itForce = m_forces.begin() + oAtoms;
        for( vtkIdType i = 0 ; i < (vtkIdType)cAtoms ; ++i )
        {
            vtkVector3f pos = mol->GetAtomPosition(i);
            int id          = mol->GetAtomAtomicNumber(i);

            it->p << pos[0] , pos[1] , pos[2];
            it->v << 0.     , 0.     , 0.    ;
            it->atomicNumber = m_forceField.symbolToNumber( atomicSymbols->GetValue(id) );

            *itForce << 0. , 0. , 0.;
            ++it;
            ++itForce;
        }
    }

    m_bonds.resize( mol->GetNumberOfBonds() );
    {
        auto it = m_bonds.begin() + oBonds;
        for( vtkIdType i = 0 ; i < (vtkIdType)cBonds ; ++i )
        {
            vtkBond bond = mol->GetBond(i);
            int a = (int)bond.GetBeginAtomId();
            int b = (int)bond.GetEndAtomId();
            *it++ << a,b;

            m_atoms[a].neighbours.insert(b);
            m_atoms[b].neighbours.insert(a);

            m_atoms[ std::min(a,b) ].neighbours_fwd.insert( std::max(a,b) );
        }
    }

    // TODO, find better way
    m_molecule->DeepCopy( mol );
}

void CMoleculeState::loadForceFieldFromFile(const std::string &filename)
{
    m_forceField.loadFromFile(filename);
}

void CMoleculeState::deepCopyMolecule(vtkMolecule *pMolecule)
{
    for( size_t i = 0 ; i < numberOfAtoms(); ++i )
    {
        Eigen::Vector3d& p = pos(i);
        m_molecule->SetAtomPosition( (vtkIdType)i , vtkVector3f( p[0] , p[1] , p[2] ) );
        ++i;
    }

    pMolecule->DeepCopy( m_molecule.GetPointer() );
}



vtkMolecule *CMoleculeState::getMolecule()
{
    return m_molecule.GetPointer();
}

void CMoleculeState::printInfo()
{
    using std::cout;
    using std::endl;

    for( size_t i = 0; i < numberOfAtoms() ; ++i )
    {
    }
}

void CMoleculeState::calculateForces()
{
    const size_t cAtoms = numberOfAtoms();

    // Monopotential
    for( int i = 0; i < cAtoms ; ++i )
        m_forces[i] << 0. , 0. , 0.;

    // Bipotential
    Eigen::Vector2i atomicNumbers;
    Eigen::Matrix<size_t,2,1> i;
    for( i[0] = 0 ; i[0] < cAtoms ; ++i[0] )
    {
        for( i[1] = i[0]+1 ; i[1] < cAtoms ; ++i[1] )
        {
            atomicNumbers << atomicNumber(i[0]) , atomicNumber(i[1]);
            TVec3 F = m_forceField.force_lenard_jones( pos(i[0]) , pos(i[1]) , atomicNumbers );

            m_forces[i[0]] += F;
            m_forces[i[1]] -= F;
        }
    }

    // Bounded potentials
    // iterate through all bonds

    // Bipotentials

    for( auto& b : m_bonds )
    {
        atomicNumbers << atomicNumber(b[0]) , atomicNumber(b[1]);

        TVec3 F = m_forceField.force_bond( pos(b[0]) , pos(b[1]) , atomicNumbers );
        m_forces[b[0]] += F;
        m_forces[b[1]] -= F;
    }
}
