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


class CParticleSystem
{

public:
    CParticleSystem( size_t cParticles );

    void generateRandomPositions();
    void generateRandomVelocities();

    void updateForces()
    {
        double d;
        for( size_t i = 0 ; i < m_cParticles ; ++i )
        {
            m_forces[i] = vtkVector3d(0.);
            for( size_t j = 0 ; j < m_cParticles ; ++j )
            {
                m_forces[i] = m_forces[i];
            }
        }
    }

    void eulerStep( double dt )
    {
        for( size_t i = 0 ; i < m_cParticles ; ++i )
        {
            m_pos[i] = m_pos[i] + dt * m_vel[i];
            m_vel[i] = m_vel[i] + dt / m_mass * m_forces[i];
        }
    }



    // access functions

    inline double* GetPosition( size_t i )
    {
        return m_pos[i].GetData();
    }

private:
    size_t m_cParticles;

    // Positions
    std::vector<vtkVector3d> m_pos;

    // Velocities
    std::vector<vtkVector3d> m_vel;

    // Velocities
    std::vector<vtkVector3d> m_forces;

    // Mass
    double m_mass;


public:
    static std::mt19937 generator;
    static std::normal_distribution<double> normal_distirbution;

    static inline double rand_norm()
    {
        return normal_distirbution( generator );
    }

    static inline vtkVector3d rand_norm3d()
    {
        return vtkVector3d( rand_norm() , rand_norm() , rand_norm() );
    }
};

std::normal_distribution<double> CParticleSystem::normal_distirbution(0.f,1.f);
std::mt19937 CParticleSystem::generator;

CParticleSystem::CParticleSystem( size_t cParticles ):
    m_cParticles( cParticles ),
    m_pos( cParticles , vtkVector3d(0.) ),
    m_vel( cParticles , vtkVector3d(0.) ),
    m_forces( cParticles , vtkVector3d(0.) ),
    m_mass( 1. )
{
}

void CParticleSystem::generateRandomPositions()
{
    for( auto& p : m_pos )
        p = rand_norm3d();
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


    for( vtkIdType i = 0 ; i < numPts ; ++i )
    {
        inPts->SetPoint( i , particles->GetPosition(i) );
    }

}

int main(int, char *[])
{
    int cParticles = 100000;


    CParticleSystem particles( cParticles );

    particles.generateRandomPositions();
    particles.generateRandomVelocities();



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

    SUpdateParams params{ &particles , programmableFilter , 0.01 };
    programmableFilter->SetExecuteMethod( Update , &params);


    // Visualization
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(programmableFilter->GetOutputPort());

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(1);

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
    renderWindowInteractor->CreateRepeatingTimer(10);
    renderWindowInteractor->AddObserver( vtkCommand::TimerEvent , timerCallback );

    // Render and interact
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
