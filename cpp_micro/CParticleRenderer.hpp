#ifndef CPARTICLERENDERER_HPP
#define CPARTICLERENDERER_HPP

// vtk includes
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProgrammableFilter.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCallbackCommand.h>

// own headers
#include "CParticleSystem.hpp"




class CParticleRenderer
{
public:
    struct SUpdateParams
    {
        CParticleSystem*        particleSystem;
        vtkProgrammableFilter*  progammableFilter;
    };

public:
    CParticleRenderer();

    void init(CParticleSystem *particleSystem);

    void start();

    void registerUpdateCall();


private:

    vtkSmartPointer<vtkPoints>                  m_points;
    vtkSmartPointer<vtkPolyData>                m_pointsPolyData;
    vtkSmartPointer<vtkPolyData>                m_polydata;
    vtkSmartPointer<vtkVertexGlyphFilter>       m_vertexFilter;
    vtkSmartPointer<vtkPolyDataMapper>          m_mapper;
    vtkSmartPointer<vtkActor>                   m_actor;
    vtkSmartPointer<vtkRenderer>                m_renderer;
    vtkSmartPointer<vtkRenderWindow>            m_renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor>  m_renderWindowInteractor;
    vtkSmartPointer<vtkProgrammableFilter>      m_programmableFilter;
    vtkSmartPointer<vtkCallbackCommand>         m_timerCallback;

    SUpdateParams m_params;
};

#endif // CPARTICLERENDERER_HPP
