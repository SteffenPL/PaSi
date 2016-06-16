#include "MoleculeEditor.hpp"

#include "ui_MoleculeEditor.h"

// vtk includes
#include <vtkNew.h>
#include <vtkProgrammableFilter.h>

#include "vtkMolecule.h"
#include "vtkMoleculeMapper.h"
#include "vtkOpenGLMoleculeMapper.h"
#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"


CMoleculeEditor::CMoleculeEditor(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::CMoleculeEditor)
{
    ui->setupUi(this);

    this->renderStyle = EBallAndSticks;
}

CMoleculeEditor::~CMoleculeEditor()
{
    delete ui;
}

void CMoleculeEditor::initRenderWindow()
{
    m_renWin.TakeReference( ui->renderWidget->GetRenderWindow() );
}
