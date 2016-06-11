#include <vtkVector.h>
#include <vtkMath.h>

// std
#include <random>
#include <vector>
#include <functional>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <sstream>

// boost
#include "boost/functional/hash.hpp"

// own includes
#include "ParticleSystem.hpp"
#include "ParticleRenderer.hpp"

#include "ConfigManager.hpp"

#include "vtkMolecule.h"
#include "vtkMoleculeMapper.h"
#include "vtkProperty.h"
#include "vtkLight.h"
#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkPolyDataMapper.h"
#include "vtkProteinRibbonFilter.h"
#include "vtkPeriodicTable.h"
#include "vtkBlueObeliskData.h"
#include "vtkFloatArray.h"
#include "vtkStringArray.h"

#include "vtkPDBReader.h"
#include "vtkGaussianCubeReader2.h"


#include "eigen3/Eigen/Dense"

class CSolver
{

};


void appendMolecule( vtkMolecule& mol , vtkMolecule& newMol )
{
    for( vtkIdType i = 0 ; newMol.GetNumberOfAtoms() ; ++i )
    {
        vtkAtom atom = newMol.GetAtom(i);
        mol.AppendAtom( atom.GetAtomicNumber() , atom.GetPosition() );
    }

    vtkIdType offset = mol.GetNumberOfAtoms();

    for( vtkIdType i = 0 ; i < newMol.GetNumberOfBonds() ; ++i )
    {
        vtkBond bond = newMol.GetBond( i );
        mol.AppendBond( offset + bond.GetBeginAtomId() ,
                        offset + bond.GetEndAtomId() ,
                        bond.GetOrder() );
    }
}



namespace std {
  template < typename T, int N >
  struct hash<Eigen::Matrix<T,N,1>>
  {
    std::size_t operator()(const Eigen::Matrix<T,N,1>& k) const
    {
        size_t seed = 0;
        for( int i = 0 ; i < N ; ++i )
        boost::hash_combine(seed , k[i]);
        return seed;
    }
  };
}

class CForceField
{

    using TVec3  = Eigen::Vector3d;

public:

    void setAtomicNumber( const std::string& name , int id )
    {
        m_atomicNumbers[name] = id;
    }

    bool loadFromFile( const std::string& filename )
    {
        setAtomicNumber("X", -1);

        std::ifstream file( filename );

        enum EParameterType{
            EUnknown,
            EMass,
            EBonds,
            EAngles,
            EDihedrals,
            EImproper,
            ENonBonded
        } state;

        state = EUnknown;

        std::string line;
        std::string value;

        while( file.good() )
        {
            std::getline( file , line );

            auto it = line.find('!');
            if( it != std::string::npos )
                line = line.substr( 0 , --it );

            std::stringstream ss;
            ss << line;
            ss >> value;


            bool newCategory = true;

            // there are different cases

            // is the current value a new category?

            if( value == "BONDS" )
                state = EBonds;
            else if( value == "ANGLES" )
                state = EAngles;
            else if( value == "DIHEDRALS" )
                state == EDihedrals;
            else if( value == "IMPROPER" )
                state = EImproper;
            else if( value == "NONBONDED")
                state = ENonBonded;
            else
                newCategory = false;

            // ok, if the current line is not the title of a new category,
            // then we exect it to contain parameters or a comment

            if( !newCategory )
            {
                switch(state)
                {
                    case EBonds:
                    {
                        std::string a1,a2;
                        double kb,b0;

                        a1 = value;
                        ss >> a2 >> kb >> b0;

                        Eigen::Vector2i v;
                        v[0] = atomID(a1);
                        v[1] = atomID(a2);

                        m_kb[v] = kb;
                        m_b0[v] = b0;
                    }
                    break;


                    case EAngles:
                    {
                        std::string a1,a2,a3;
                        double k_theta,theta0;

                        a1 = value;
                        ss >> a2 >> a3 >> k_theta >> theta0;

                        Eigen::Vector3i v;
                        v[0] = atomID(a1);
                        v[1] = atomID(a2);
                        v[2] = atomID(a3);

                        m_k_theta[v] = k_theta;
                        m_theta0[v] = theta0;
                    }
                    break;


                    case EDihedrals:
                    {
                        std::string a1,a2,a3,a4;
                        double k_chi,delta;
                        int n;

                        a1 = value;
                        ss >> a2 >> a3 >> a4 >> k_chi >> n >> delta;

                        Eigen::Vector4i v;
                        v[0] = atomID(a1);
                        v[1] = atomID(a2);
                        v[2] = atomID(a3);
                        v[3] = atomID(a4);

                        m_k_chi[v] = k_chi;
                        m_n[v] = n;
                        m_delta[v] = delta;
                    }
                    break;


                    case EImproper:
                    {
                        std::string a1,a2,a3,a4;
                        double k_psi,psi0;

                        a1 = value;
                        ss >> a2 >> a3 >> a4 >> k_psi >> psi0;

                        Eigen::Vector4i v;
                        v[0] = atomID(a1);
                        v[1] = atomID(a2);
                        v[2] = atomID(a3);
                        v[3] = atomID(a4);

                        m_k_psi[v] = k_psi;
                        m_psi0[v] = psi0;
                    }
                    break;

                    case ENonBonded:
                    {
                        std::string a1;
                        double epsilon,r_min_half;

                        a1 = value;
                        ss >> epsilon >> r_min_half;

                        int v = atomID(a1);

                        m_epsilon[v] = epsilon;
                        m_r_min_half[v] = r_min_half;
                    }
                    break;

                    case EMass:
                    {
                        std::string a1;
                        double mass;

                        a1 = value;
                        ss >> mass;

                        int v = atomID(a1);

                        m_mass[v] = mass;
                    }
                    break;


                    default: // do nothing
                    break;
                }
            }
        }


        return true;
    }


    Eigen::Vector3d force_bond( const Eigen::Vector3d& a1 , const Eigen::Vector3d& a2 , const Eigen::Vector2i& atomicNumbers )
    {
        // V_bond = kb(r-b0)**2
        // grad V = 2*kb*(r-b0) * vec r / norm r;

        TVec3 r = a2 - a1;

        return 2. * m_kb[atomicNumbers] * ( r.norm() - m_b0[atomicNumbers] ) * r.normalized();
    }

    Eigen::Vector3d force_lenard_jones( const Eigen::Vector3d& a1 , const Eigen::Vector3d& a2 , const Eigen::Vector2i& atomicNumbers )
    {
        //V = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
        //
        //epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
        //Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
        //
        // grad V = epsij * ...

        double eps  = sqrt( m_epsilon[atomicNumbers[0]] + m_epsilon[atomicNumbers[1]] );
        double Rmin = m_r_min_half[atomicNumbers[0]] + m_r_min_half[atomicNumbers[1]];

        TVec3 r = a2 - a1;

        double d = Rmin / r.squaredNorm();

        return 2. * m_kb[atomicNumbers] * ( r.norm() - m_b0[atomicNumbers] ) * r.normalized();
    }


private:

    inline int atomID( const std::string& ref )
    {
        return m_atomicNumbers[ref];
    }


    std::map<std::string,int>   m_atomicNumbers;

    // parameters...

    std::unordered_map<int,double> m_mass;

    std::unordered_map<int,double> m_epsilon;
    std::unordered_map<int,double> m_r_min_half;

    std::unordered_map<Eigen::Vector2i,double> m_kb;
    std::unordered_map<Eigen::Vector2i,double> m_b0;

    std::unordered_map<Eigen::Vector3i,double> m_k_theta;
    std::unordered_map<Eigen::Vector3i,double> m_theta0;

    std::unordered_map<Eigen::Vector4i,double> m_k_chi;
    std::unordered_map<Eigen::Vector4i,int>    m_n;
    std::unordered_map<Eigen::Vector4i,double> m_delta;

    std::unordered_map<Eigen::Vector4i,double> m_k_psi;
    std::unordered_map<Eigen::Vector4i,double> m_psi0;

};

class CState
{

    using TVec3 = Eigen::Vector3d;

public:

    CState()
    {

        vtkNew<vtkPeriodicTable> table;
        m_blueObeliskData = table->GetBlueObeliskData();

        auto atomicNames = m_blueObeliskData->GetNames();
        for( int i = 0 ; i < atomicNames->GetSize() ; ++i )
        {
            m_forceField.setAtomicNumber( atomicNames->GetValue(i) , i );
        }

    }

    void loadMoleculeFromFile( const std::string& filename )
    {
        auto pdb = vtkSmartPointer<vtkPDBReader>::New();
        pdb->SetFileName( filename.c_str()   );
        pdb->SetHBScale(1.0);
        pdb->SetBScale(1.0);
        pdb->Update();
        m_molecule = vtkMolecule::SafeDownCast(
                    pdb->GetOutputDataObject(1) );

        m_vel.resize( m_molecule->GetNumberOfAtoms() );
        m_pos.resize( m_molecule->GetNumberOfAtoms() );


        // init position

        auto it = m_pos.begin();
        for( vtkIdType i = 0 ; i < numberOfAtoms() ; ++i )
        {
            vtkVector3f pos = m_molecule->GetAtomPosition(i);
            *it = Eigen::Vector3d( pos[0] , pos[1] , pos[2] );
            ++it;
        }


    }

    void loadForceFieldFromFile( const std::string& filename )
    {
        m_forceField.loadFromFile(filename);
    }

    vtkMolecule* getMolecule()
    {
        // copy positions to molecule

        vtkIdType i = 0;

        for( auto& pos : m_pos )
        {
            m_molecule->SetAtomPosition( i , vtkVector3f( pos[0] , pos[1] , pos[2] ) );
            ++i;
        }

        return m_molecule;
    }

    void printInfo()
    {
        using std::cout;
        using std::endl;

        for( int i = 0; i < numberOfAtoms() ; ++i )
        {
            cout << i  << "   id: " << atomicNumber(i)
                 << "   name: " << name(i) << "   mass: " << mass(i)
                 << "   vdw radius: " << vdwRadius(i) << endl;
        }
    }

    void calculateForces()
    {
        const int cAtoms = numberOfAtoms();

        // Monopotential
        for( int i = 0; i < cAtoms ; ++i )
            m_forces[i] = Force1( i );

        // Bipotential
        for( int i = 0 ; i < cAtoms ; ++i )
        {
            for( int j = 0 ; j < i ; ++j )
            {
                TVec3 F = Force2( i , j );

                m_forces[i] += F;
                m_forces[j] -= F;
            }
        }

        // Bounded potentials
        // iterate through all bonds

        const int cBonds = m_molecule->GetNumberOfBonds();

        // Bipotentials

        for( int k = 0 ; k < cBonds ; ++k )
        {
            vtkBond bond = m_molecule->GetBond(k);

            auto i = bond.GetBeginAtomId();
            auto j = bond.GetEndAtomId();

            TVec3 F = Force2B( i , j );
            m_forces[i] += F;
            m_forces[j] -= F;
        }
    }

    void update()
    {

    }

    int numberOfAtoms() const
    {
        return m_molecule->GetNumberOfAtoms();
    }

    inline Eigen::Vector3d& pos( int i )
    {
        return m_pos[i];
    }

    inline Eigen::Vector3d& vel( int i )
    {
        return m_vel[i];
    }

    inline unsigned int atomicNumber( int i )
    {
        return m_molecule->GetAtomAtomicNumber(i);
    }

    // in dalton
    inline double mass( int i )
    {
        return blueObeliskData()->GetMasses()->GetTuple1(atomicNumber(i));
    }

    inline std::string name( int i )
    {
        return blueObeliskData()->GetNames()->GetValue(atomicNumber(i));
    }

    // in pico meter 10^-12
    inline double vdwRadius( int i )
    {
        return blueObeliskData()->GetVDWRadii()->GetValue( atomicNumber(i) );
    }

protected:

    inline vtkBlueObeliskData* blueObeliskData()
    {
        return m_blueObeliskData;
    }

private:

    CForceField     m_forceField;

    vtkSmartPointer<vtkBlueObeliskData> m_blueObeliskData;
    vtkSmartPointer<vtkMolecule>   m_molecule;

    std::vector< Eigen::Vector3d > m_pos;
    std::vector< Eigen::Vector3d > m_vel;

    std::vector< Eigen::Vector3d > m_forces;

    std::list< Eigen::Vector2i > m_bonds;
    std::list< Eigen::Vector3i > m_tribond;
    std::list< Eigen::Vector4i > m_quadbond;


};


struct SParams
{

    CState* state;
    vtkSmartPointer<vtkProgrammableFilter> filter;
};

// far to complex construction to get a periodic update
// of the particles
void TimerCallbackFunction ( vtkObject* caller
                             , long unsigned int vtkNotUsed(eventId)
                             , void* clientData
                             , void* vtkNotUsed(callData) )
{

    auto para = static_cast<SParams*>(clientData);

    // update molecule
    para->state->update();

    // update plot
    para->filter->Modified();
    static_cast<vtkRenderWindowInteractor*>(caller)->Render();
}



void Update( void* args )
{
    auto para = static_cast<SParams*>(args);
    auto newmol = vtkSmartPointer<vtkMolecule>::New();
    newmol->DeepCopy( para->state->getMolecule() );
    para->filter->GetOutput()->ShallowCopy( newmol );
}


int main( int argc , char ** argv )
{
    CState state;
    state.loadMoleculeFromFile( "1skl.pdb");

    state.printInfo();

    /*
    auto pdb = vtkSmartPointer<vtkPDBReader>::New();
    pdb->SetFileName( "5a2m.pdb"   );
    pdb->SetHBScale(1.0);
    pdb->SetBScale(1.0);
    pdb->Update();*/

    auto programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
    programmableFilter->SetInputData( state.getMolecule() );


    SParams params = { &state , programmableFilter };

    programmableFilter->SetExecuteMethod( Update , &params);

    // Add Callback function
    auto callback = vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback( TimerCallbackFunction );
    callback->SetClientData( &params );

    /*vtkNew<vtkElectronMapper> eleMapper;
    eleMapper->SetInputConnection( programmableFilter->GetOutputPort() );
    eleMapper->UseMolecularOrbita(4);

    vtkNew<vtkElectronActor> eleActor;
    eleActor->SetMapper( eleMapper.GetPointer() );*/


    /*vtkNew<vtkProteinRibbonFilter> ribbonFilter;
    ribbonFilter->SetInputConnection( pdb->GetOutputPort() );
    ribbonFilter->Update();

    vtkNew<vtkPolyDataMapper> polyDataMapper;
    polyDataMapper->SetInputData( ribbonFilter->GetOutput() );
    polyDataMapper->Update();

    vtkNew<vtkActor> ribbonActor;
    ribbonActor->SetMapper( polyDataMapper.GetPointer() );*/


    // Visualization
    auto mapper = vtkSmartPointer<vtkMoleculeMapper>::New();
    mapper->SetInputConnection( programmableFilter->GetOutputPort());
    //mapper->UseBallAndStickSettings();
    //mapper->UseVDWSpheresSettings();
    mapper->UseLiquoriceStickSettings();
    //mapper->UseFastSettings();

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper( mapper );

    // Create a renderer, render window, and interactor
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetMultiSamples(2);

    // Add the actor to the scene
    renderer->AddActor(actor);
    //renderer->AddActor(ribbonActor.GetPointer());
    //renderer->AddActor(eleActor);
    renderer->SetBackground(.1,.1,.1); // Background color green

    auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindowInteractor->Initialize();
    renderWindowInteractor->CreateRepeatingTimer(1);
    renderWindowInteractor->AddObserver( vtkCommand::TimerEvent , callback );


    // Render and interact
    renderWindow->Render();
    renderWindowInteractor->Start();

    return 0;
}

/*
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

    //particles.setPotential1( Gravity_Potential<CParticleSystem::TCodiVec3d> );
    //particles.setPotential2( VanDerWaals_Potential<CParticleSystem::TCodiVec3d> );


    particles.setRadius( 1000. * particles.getSigma() );

    //particles.setNumberOfNeighbourCells( 6,6,6 );


    //particles.generateRandomVelocities(1);

    CParticleRenderer renderer;

    renderer.init( &particles );
    renderer.start();


    return EXIT_SUCCESS;
}
*/
