#!/bin/bash

# Todo: this should be using CMake to compile the programs
# These compilation scripts have been moved here from the tests, and will do for now

RESET_STYLE='\033[0m'
INFO_STYLE='\033[0;33m\033[1m'

GMML_ROOT_DIR=$(pwd)
BIN_PATH="${GMML_ROOT_DIR}/bin/"

echo -e "\n${INFO_STYLE}###### Compiling internal programs ######${RESET_STYLE}"
echo -e "Glycoprotein builder"
g++ -std=c++17 -g -I ./ -L"${BIN_PATH}" -Wl,-rpath,"${BIN_PATH}" internalPrograms/GlycoproteinBuilder/gpBuilder_main.cpp -lgmml2 -pthread -o bin/gpBuilder
echo -e "Glycoprotein builder table"
g++ -std=c++17 -g -I ./ -L"${BIN_PATH}" -Wl,-rpath,"${BIN_PATH}" internalPrograms/createGlycosylationTables.cpp -lgmml2 -pthread -o bin/gpBuilderTable
echo -e "Carbohydrate builder"
g++ -std=c++17 -g -I ./ -L"${BIN_PATH}" -Wl,-rpath,"${BIN_PATH}" internalPrograms/CarbohydrateBuilder/main.cpp -lgmml2 -pthread -o bin/carbohydrateBuilder
echo -e "Wiggle to site"
g++ -std=c++17 -g -I ./ -L"${BIN_PATH}" -Wl,-rpath,"${BIN_PATH}" internalPrograms/WiggleToSite/wiggleToSiteDriver.cpp -lgmml2 -pthread -o bin/wiggleToSite
echo -e "Draw glycan"
g++ -std=c++17 -g -I ./ -L"${BIN_PATH}" -Wl,-rpath,"${BIN_PATH}" internalPrograms/DrawGlycan/drawGlycanMain.cpp -lgmml2 -pthread -o bin/drawGlycan
