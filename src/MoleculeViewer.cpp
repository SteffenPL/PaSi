#include "MoleculeViewer.hpp"

// includes
#include "MoleculeState.hpp"

// vtk includes
#include <vtkNew.h>
#include <vtkProgrammableFilter.h>

#include "vtkMolecule.h"
#include "vtkMoleculeMapper.h"
#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"


CMoleculeViewer::CMoleculeViewer():
    m_state( nullptr )
{

}

CMoleculeViewer::~CMoleculeViewer()
{

}

void CMoleculeViewer::init(CMoleculeState *state)
{
    m_state = state;

    m_programmableFilter->SetInputData( state->getMolecule() );
    m_programmableFilter->SetExecuteMethod( Update , this);

    m_callback->SetCallback( TimerCallbackFunction );
    m_callback->SetClientData( this );

    m_moleculeMapper->SetInputConnection( m_programmableFilter->GetOutputPort() );
    //m_moleculeMapper->UseBallAndStickSettings();
    //m_moleculeMapper->UseVDWSpheresSettings();
    m_moleculeMapper->UseLiquoriceStickSettings();
    //m_moleculeMapper->UseFastSettings();

    m_moleculeActor->SetMapper( m_moleculeMapper.GetPointer() );

    m_renWin->AddRenderer( m_ren.GetPointer() );

    m_ren->AddActor( m_moleculeActor.GetPointer() );
    m_ren->SetBackground( 0.1 , 0.1 , 0.1 );

    //m_renWin->SetMultiSamples(2);

    m_renInteractor->SetRenderWindow( m_renWin.GetPointer() );
    m_renInteractor->Initialize();

    m_renInteractor->CreateRepeatingTimer( 1 );
    m_renInteractor->AddObserver( vtkCommand::TimerEvent , m_callback.GetPointer() );

}

void CMoleculeViewer::start()
{
    m_renWin->Render();
    m_renInteractor->Start();
}


// implementations

// far to complex construction to get a periodic update
// of the particles
void CMoleculeViewer::TimerCallbackFunction ( vtkObject* caller
                             , long unsigned int vtkNotUsed(eventId)
                             , void* clientData
                             , void* vtkNotUsed(callData) )
{

    auto viewer = static_cast<CMoleculeViewer*>(clientData);

    // update molecule
    for( int i = 0 ; i < 10 ; ++i)
        viewer->m_state->update();

    // update plot
    viewer->m_programmableFilter->Modified();
    static_cast<vtkRenderWindowInteractor*>(caller)->Render();
}


void CMoleculeViewer::Update( void* args )
{
    auto viewer = static_cast<CMoleculeViewer*>(args);
    auto newmol = vtkSmartPointer<vtkMolecule>::New();
    viewer->m_state->deepCopyMolecule( newmol );
    viewer->m_programmableFilter->GetOutput()->ShallowCopy( newmol );
}

