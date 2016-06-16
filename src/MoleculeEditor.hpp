#ifndef MOLECULEEDITOR_HPP
#define MOLECULEEDITOR_HPP

#include <QMainWindow>

#include "MoleculeViewer.hpp"

namespace Ui {
class CMoleculeEditor;
}

class CMoleculeEditor : public QMainWindow , public CMoleculeViewer
{
    Q_OBJECT

public:
    explicit CMoleculeEditor(QWidget *parent = 0);
    ~CMoleculeEditor();


protected:

    virtual void initRenderWindow();

private:
    Ui::CMoleculeEditor *ui;
};

#endif // MOLECULEEDITOR_HPP
