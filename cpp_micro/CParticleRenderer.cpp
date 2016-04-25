#include "CParticleRenderer.hpp"

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

CParticleRenderer::CParticleRenderer()
{

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



void Update( void* args )
{
    auto para = static_cast<CParticleRenderer::SUpdateParams*>(args);

    CParticleSystem*        particles = para->particleSystem;
    vtkProgrammableFilter*  filter    = para->progammableFilter;

    // simulate

    particles->update();

    // assume the the number of points is consant!
    vtkPoints* inPts = filter->GetPolyDataInput()->GetPoints();
    vtkIdType numPts = inPts->GetNumberOfPoints();


    double factor = 1. / particles->getSigma();
    double r = 2. * factor * particles->getRadius() ;
    double rr = r*r;

    vtkVector3d v;
    for( vtkIdType i = 0 ; i < numPts ; ++i )
    {
        v = factor * particles->getPosition(i);

        if( v.SquaredNorm() > rr )
            v = r/v.Norm()* v;

        inPts->SetPoint( i , v.GetData() );

    }

}

 // end of anonymous namespace

void CParticleRenderer::init( CParticleSystem* particleSystem )
{
    // Create points
    m_points = vtkSmartPointer<vtkPoints>::New();

    size_t cParticles = particleSystem->getNumberOfParticles();
    m_points->SetNumberOfPoints( cParticles );

    for( size_t i = 0 ; i < cParticles ; ++i )
        m_points->SetPoint( i , 0. , 0. , 0. );


    m_pointsPolyData = vtkSmartPointer<vtkPolyData>::New();
    m_pointsPolyData->SetPoints( m_points );

    // Create vertices for these points
    m_vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    m_vertexFilter->SetInputData( m_pointsPolyData );
    m_vertexFilter->Update();

    m_polydata = vtkSmartPointer<vtkPolyData>::New();
    m_polydata->ShallowCopy(m_vertexFilter->GetOutput());


    // Set up a programmableFilter to manipulate the points while rendering
    m_programmableFilter = vtkSmartPointer<vtkProgrammableFilter>::New();
    m_programmableFilter->SetInputData( m_polydata );

    m_params.particleSystem = particleSystem;
    m_params.progammableFilter = m_programmableFilter;

    m_programmableFilter->SetExecuteMethod( Update , &m_params);


    // Visualization
    m_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_mapper->SetInputConnection(m_programmableFilter->GetOutputPort());

    m_actor = vtkSmartPointer<vtkActor>::New();
    m_actor->SetMapper(m_mapper);
    m_actor->GetProperty()->SetPointSize(3);

    // Create a renderer, render window, and interactor
    m_renderer = vtkSmartPointer<vtkRenderer>::New();
    m_renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    m_renderWindow->AddRenderer(m_renderer);
    m_renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    m_renderWindowInteractor->SetRenderWindow(m_renderWindow);

    // Add the actor to the scene
    m_renderer->AddActor(m_actor);
    m_renderer->SetBackground(.1,.1,.1); // Background color green


    m_timerCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    m_timerCallback->SetCallback( TimerCallbackFunction );
    m_timerCallback->SetClientData( m_programmableFilter );

    m_renderWindowInteractor->Initialize();
    m_renderWindowInteractor->CreateRepeatingTimer(20);
    m_renderWindowInteractor->AddObserver( vtkCommand::TimerEvent , m_timerCallback );
}

void CParticleRenderer::start()
{
    // Render and interact
    m_renderWindow->Render();
    m_renderWindowInteractor->Start();
}
