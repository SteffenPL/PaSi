#include <iostream>


// includes
#include "ConfigManager.hpp"

#include "MoleculeViewer.hpp"
#include "MoleculeState.hpp"

using namespace std;


int main( int argc , char** argv )
{

    CConfigManager config;
    config.parse(argc , argv );

    // init molecule

    CMoleculeState state;
    if( config.hasKey("molecule") )
        state.loadMoleculeFromFile( config.getValue<std::string>("molecule") );
    else
        cout << "No molecule file given...\n";

    if( config.hasKey("forceField") )
        state.loadForceFieldFromFile( config.getValue("forceField") );
    else
        cout << "No force field field given...\n";


    CMoleculeViewer viewer;

    viewer.init( &state );
    viewer.start();


}
