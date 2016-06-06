#include "ParticleRenderer.hpp"

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
#include <vtkTextProperty.h>
#include <vtkProperty2D.h>

#include <vtkVector.h>
#include <vtkMath.h>

CParticleRenderer::CParticleRenderer()
{

}

namespace
{

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

    CParticleSystem*           particles = para->particleSystem;
    vtkProgrammableFilter*     filter    = para->progammableFilter;
    vtkSliderRepresentation2D* slider    = para->sliderTemperature;

    // simulate

    particles->update();

    // TODO: share memory!!!
    // assume the the number of points is consant!
    vtkPoints* inPts = filter->GetPolyDataInput()->GetPoints();
    vtkIdType numPts = inPts->GetNumberOfPoints();


    vtkVector3d camCenter = particles->getCamCenter();
    double factor = 1. / particles->getSigma();
    double r = 2. * factor * particles->getRadius() ;
    double rr = r*r;

    vtkVector3d v;
    for( vtkIdType i = 0 ; i < numPts ; ++i )
    {
        v = factor * (particles->getPosition(i) - camCenter);

        if( v.SquaredNorm() > rr )
            v = r/v.Norm()* v;

        inPts->SetPoint( i , v.GetData() );

    }

    // update view / sliders

    slider->SetValue( pow(particles->getTemperature(),1/5.) );
}

// Slider Callback


} // end of anonymous namespace


// see: http://www.vtk.org/Wiki/VTK/Examples/Cxx/Widgets/Slider2D
class CSliderCallback : public vtkCommand
{
public:
    static CSliderCallback *New()
    {
        return new CSliderCallback;
    }

    virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
        vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
        double temperature = static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue();

        this->particleSystem->setTemparature( pow(temperature,5) );
    }
    CSliderCallback():particleSystem(0) {}

    CParticleSystem* particleSystem;
};

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

    // create sphere sources for each points
    m_sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    m_glyph3D = vtkSmartPointer<vtkGlyph3D>::New();

    m_sphereSource->SetRadius(1);
    m_glyph3D->SetSourceConnection(m_sphereSource->GetOutputPort());
    m_glyph3D->SetInputConnection(m_programmableFilter->GetOutputPort());
    m_glyph3D->Update();

    // Visualization
    m_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_mapper->SetInputConnection(m_glyph3D->GetOutputPort());

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
    // Add slider

    m_sliderRep =
      vtkSmartPointer<vtkSliderRepresentation2D>::New();

    m_sliderRep->SetMinimumValue(0.0);
    m_sliderRep->SetMaximumValue(10.0);
    m_sliderRep->SetValue( pow( particleSystem->getTemperature() , 1./5. ) );
    m_sliderRep->SetTitleText("Temperature");


    m_sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
    m_sliderRep->GetPoint1Coordinate()->SetValue( 10 , 10);
    m_sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
    m_sliderRep->GetPoint2Coordinate()->SetValue( 110, 10);

    m_params.particleSystem = particleSystem;
    m_params.progammableFilter = m_programmableFilter;
    m_params.sliderTemperature = m_sliderRep;

    m_programmableFilter->SetExecuteMethod( Update , &m_params);


    m_sliderWidget =
      vtkSmartPointer<vtkSliderWidget>::New();
    m_sliderWidget->SetInteractor(m_renderWindowInteractor);
    m_sliderWidget->SetRepresentation(m_sliderRep);
    m_sliderWidget->SetAnimationModeToAnimate();
    m_sliderWidget->EnabledOn();

    m_sliderCallback = vtkSmartPointer<CSliderCallback>::New();
    m_sliderCallback->particleSystem = particleSystem;

    m_sliderWidget->AddObserver(vtkCommand::InteractionEvent,m_sliderCallback);

    m_renderWindowInteractor->Initialize();
    m_renderWindowInteractor->CreateRepeatingTimer(1);
    m_renderWindowInteractor->AddObserver( vtkCommand::TimerEvent , m_timerCallback );
}

void CParticleRenderer::start()
{
    // Render and interact
    m_renderWindow->Render();
    m_renderWindowInteractor->Start();
}
