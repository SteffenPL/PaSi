#include <QApplication>
#include <QCommandLineParser>
#include <QCommandLineOption>


// includes
#include "ConfigManager.hpp"

#include "MoleculeViewer.hpp"
#include "MoleculeState.hpp"

// user interface
#include "MoleculeEditor.hpp"


using namespace std;


int main(int argc, char *argv[])
{


    CConfigManager config;
    config.parse(argc , argv );

    std::ifstream configFile("../data/config.txt");
    config.parse( configFile );

    cout << config << "\n" << std::flush;

    // init molecule

    CMoleculeState state;
    if( config.hasKey("molecule") )
        state.loadMoleculeFromFile( "../data/" + config.getValue<std::string>("molecule") );
    else
        cout << "No molecule file given...\n";

    if( config.hasKey("forceField") )
        state.loadForceFieldFromFile( "../data/" + config.getValue("forceField") );
    else
        cout << "No force field field given...\n";


    QApplication app(argc, argv);


    CMoleculeEditor mainWin;

    mainWin.init( &state );

    mainWin.show();
    return app.exec();
}
