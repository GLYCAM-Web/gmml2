#!/bin/bash

# Todo: this should be using CMake to compile the programs

RESET_STYLE='\033[0m'
INFO_STYLE='\033[0;33m\033[1m'

GMML_ROOT_DIR=$(pwd)
LIB_PATH="${GMML_ROOT_DIR}/lib/"

mkdir -p "${GMML_ROOT_DIR}/bin/"

echo -e "\n${INFO_STYLE}###### Compiling internal programs ######${RESET_STYLE}"
echo -e "Glycoprotein builder"
g++ -std=c++17 -g -I ./ -L"${LIB_PATH}" -Wl,-rpath,"${LIB_PATH}" programs/GlycoproteinBuilder/gpBuilder_main.cpp -lgmml2 -pthread -o bin/glycoproteinBuilder
echo -e "Glycosylation site finder"
g++ -std=c++17 -g -I ./ -L"${LIB_PATH}" -Wl,-rpath,"${LIB_PATH}" programs/createGlycosylationTables.cpp -lgmml2 -pthread -o bin/glycosylationSiteFinder
echo -e "Carbohydrate builder"
g++ -std=c++17 -g -I ./ -L"${LIB_PATH}" -Wl,-rpath,"${LIB_PATH}" programs/CarbohydrateBuilder/carbBuilderMain.cpp -lgmml2 -pthread -o bin/carbohydrateBuilder
echo -e "Wiggle to site"
g++ -std=c++17 -g -I ./ -L"${LIB_PATH}" -Wl,-rpath,"${LIB_PATH}" programs/WiggleToSite/wiggleToSiteDriver.cpp -lgmml2 -pthread -o bin/wiggleToSite
echo -e "Draw glycan"
g++ -std=c++17 -g -I ./ -L"${LIB_PATH}" -Wl,-rpath,"${LIB_PATH}" programs/DrawGlycan/drawGlycanMain.cpp -lgmml2 -pthread -o bin/drawGlycan
echo -e "PDB Analyzer"
g++ -std=c++17 -g -I ./ -L"${LIB_PATH}" -Wl,-rpath,"${LIB_PATH}" programs/pdbAnalyzer.cpp -lgmml2 -pthread -o bin/pdbAnalyzer
