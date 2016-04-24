// vtk
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataWriter.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCommand.h>
#include <vtkProgrammableFilter.h>
#include <vtkCallbackCommand.h>

#include <vtkVector.h>
#include <vtkMath.h>

// codi pack
#include <codi.hpp>

// std
#include <random>
#include <vector>
#include <functional>

vtkVector3d operator*( double f , const vtkVector3d& vec )
{
    return vtkVector3d( f * vec.GetX() , f*vec.GetY() , f*vec.GetZ()  );
}

vtkVector3d operator+( const vtkVector3d& a ,const vtkVector3d& b)
{
    return vtkVector3d(   a.GetX() + b.GetX()
                        , a.GetY() + b.GetY()
                        , a.GetZ() + b.GetZ() );
}
vtkVector3d operator-( const vtkVector3d& a ,const vtkVector3d& b)
{
    return vtkVector3d(   a.GetX() - b.GetX()
                        , a.GetY() - b.GetY()
                        , a.GetZ() - b.GetZ() );
}


template< typename Real >
void Potential(const Real* x , Real* y , double eps , double sigma )
{

    Real d2 = sigma*sigma / (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    y[0] = 4 * eps * (d2*d2*d2*d2*d2*d2 - d2*d2*d2);

}

/*template< typename T , int N , int M >
class CMathFunction
{
public:
    CMathFunction
};*/



class CParticleSystem
{

public:
    CParticleSystem( size_t cParticles );

    void generateRandomPositions();
    void generateRandomVelocities();

    void updateForces()
    {
        static vtkVector3d grad;
        for( size_t i = 0 ; i < m_cParticles ; ++i )
        {
            m_forces[i] = vtkVector3d(0.);
            for( size_t j = 0 ; j < i ; ++j )
            {
                grad = PotentialGradient( m_pos[j] - m_pos[i] );
                m_forces[i] = m_forces[i] - grad;
                m_forces[j] = m_forces[j] + grad;
                //std::cout << grad[0] << std::endl;
            }
        }
    }

    vtkVector3d PotentialGradient( vtkVector3d x )
    {
        codi::RealForwardVec<3> in_x[3];
        codi::RealForwardVec<3> out_y[1];
        in_x[0] = x.GetX();
        in_x[1] = x.GetY();
        in_x[2] = x.GetZ();

        for( int i = 0 ; i < 3 ; ++i)
            in_x[i].gradient()[i] = 1.;

        Potential( in_x , out_y , m_epsilon , m_sigma );

        vtkVector3d grad;
        for( int i = 0 ; i < 3 ; ++i)
            grad[i] = out_y[0].getGradient()[i];
        return grad;
    }

    void eulerStep( double dt )
    {
        updateForces();
        double sqrtdt = sqrt( dt );
        for( size_t i = 0 ; i < m_cParticles ; ++i )
        {
            m_pos[i] = m_pos[i] + dt * m_vel[i];
            m_vel[i] = m_vel[i] + dt / m_mass * m_forces[i];
            m_vel[i] = m_vel[i] + sqrtdt * rand_norm3d()  ;

            // bc
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


    }



    // access functions

    inline double* GetPosition( size_t i )
    {
        return m_pos[i].GetData();
    }

    inline double GetEpsilon()
    {
        return m_epsilon;
    }

    inline double GetSigma()
    {
        return m_sigma;
    }

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

    // space ares
    double m_radius;


public:
    std::mt19937 generator;
    std::normal_distribution<double> normal_distirbution;

    inline double rand_norm()
    {
        return normal_distirbution( generator );
    }

    inline vtkVector3d rand_norm3d()
    {
        return vtkVector3d( rand_norm() , rand_norm() , (( m_dim >= 3 ) ? rand_norm() : 0.) );
    }
};

//std::normal_distribution<double> CParticleSystem::normal_distirbution(0.f,1.f);
//std::mt19937 CParticleSystem::generator;

CParticleSystem::CParticleSystem( size_t cParticles ):
    m_cParticles( cParticles ),
    m_dim(2),
    m_pos( cParticles , vtkVector3d(0.) ),
    m_vel( cParticles , vtkVector3d(0.) ),
    m_forces( cParticles , vtkVector3d(0.) ),
    m_mass( 1. ),
    m_epsilon( 1.4e-23 ),
    m_sigma( 2.56e-10 ),
    generator(),
    normal_distirbution( 0. , m_sigma ),
    m_radius(20. * m_sigma )
{
}

void CParticleSystem::generateRandomPositions()
{
    for( auto& p : m_pos )
        p = 4. * rand_norm3d();
}


void CParticleSystem::generateRandomVelocities()
{
    for( auto& v : m_vel )
        v = rand_norm3d();
}


// far to complex construction to get a periodic update
// of the particles
void TimerCallbackFunction ( vtkObject* caller
                             , long unsigned int vtkNotUsed(eventId)
                             , void* clientData
                             , void* vtkNotUsed(callData) )
{
    static_cast<vtkProgrammableFilter*>(clientData)->Modified();
    static_cast<vtkRenderWindowInteractor*>(caller)->Render();
}

struct SUpdateParams
{
    CParticleSystem*        particles;
    vtkProgrammableFilter*  progammableFilter;
    double dt;
};

void Update( void* args )
{
    SUpdateParams* para = static_cast<SUpdateParams*>(args);

    CParticleSystem*        particles = para->particles;
    vtkProgrammableFilter*  filter    = para->progammableFilter;

    // simulate

    particles->eulerStep( para->dt );

    // assume the the number of points is consant!
    vtkPoints* inPts = filter->GetPolyDataInput()->GetPoints();
    vtkIdType numPts = inPts->GetNumberOfPoints();


    double f = 1. / particles->GetSigma();
    double* v;
    double  w[3];
    for( vtkIdType i = 0 ; i < numPts ; ++i )
    {
        v = particles->GetPosition(i);
        w[0] = f*v[0];
        w[1] = f*v[1];
        w[2] = f*v[2];
        if( (w[0]*w[0] + w[1]*w[1] + w[2]*w[2]) <= 10000. )
            inPts->SetPoint( i , w[0] , w[1] , w[2] );
    }

}

int main(int, char *[])
{
    int cParticles = 10;


    CParticleSystem particles( cParticles );

    particles.generateRandomPositions();
    //particles.generateRandomVelocities();



    // Create points
    auto points = vtkSmartPointer<vtkPoints>::New();

    points->SetNumberOfPoints( cParticles );

    for( int i = 0 ; i < cParticles ; ++i )
        points->SetPoint( i , particles.GetPosition(i) );


    auto pointsPolyData = vtkSmartPointer<vtkPolyData>::New();
    pointsPolyData->SetPoints( points );

    auto vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->SetInputData( pointsPolyData );
    vertexFilter->Update();

    auto polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->ShallowCopy(vertexFilter->GetOutput());



    auto programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
    programmableFilter->SetInputData( polydata );

    SUpdateParams params{ &particles , programmableFilter , 0.0001 };
    programmableFilter->SetExecuteMethod( Update , &params);


    // Visualization
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(programmableFilter->GetOutputPort());

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(3);

    // Create a renderer, render window, and interactor
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Add the actor to the scene
    renderer->AddActor(actor);
    renderer->SetBackground(.1,.1,.1); // Background color green


    auto timerCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    timerCallback->SetCallback( TimerCallbackFunction );
    timerCallback->SetClientData( programmableFilter );

    renderWindowInteractor->Initialize();
    renderWindowInteractor->CreateRepeatingTimer(1);
    renderWindowInteractor->AddObserver( vtkCommand::TimerEvent , timerCallback );

    // Render and interact
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
