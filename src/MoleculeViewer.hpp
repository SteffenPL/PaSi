#ifndef CMOLECULEVIEWER_HPP
#define CMOLECULEVIEWER_HPP

// vtk includes
#include <vtkNew.h>
#include <vtkSmartPointer.h>

class CMoleculeState;
class vtkProgrammableFilter;
class vtkCallbackCommand;
class vtkMoleculeMapper;
class vtkOpenGLMoleculeMapper;
class vtkMolecule;
class vtkActor;
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkObject;

class CMoleculeViewer
{
public:

    CMoleculeViewer();
    ~CMoleculeViewer();

    /// set up the vtk pipeline
    void init( CMoleculeState* state );
    void start();

public:

    enum ERenderStyle
    {
        EBallAndSticks,
        EVDW,
        ESticks,
        EFast
    } renderStyle;

protected:

    virtual void initRenderWindow();
    vtkSmartPointer<vtkRenderWindow>    m_renWin;

private:

    static void Update( void* args );
    static void TimerCallbackFunction ( vtkObject* caller
                                 , long unsigned int eventId
                                 , void* clientData
                                 , void* callData );
    friend void Update( void* args );
    friend void TimerCallbackFunction ( vtkObject* caller
                                 , long unsigned int eventId
                                 , void* clientData
                                 , void* callData );

    CMoleculeState*                     m_state;
    vtkNew<vtkProgrammableFilter>       m_programmableFilter;

    vtkNew<vtkCallbackCommand>          m_callback;
    vtkNew<vtkMoleculeMapper>     m_moleculeMapper;
    vtkNew<vtkActor>                    m_moleculeActor;
    vtkNew<vtkRenderer>                 m_ren;
    vtkNew<vtkRenderWindowInteractor>   m_renInteractor;

};

#endif // CMOLECULEVIEWER_HPP
