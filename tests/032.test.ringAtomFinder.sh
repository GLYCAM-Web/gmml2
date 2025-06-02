#!/bin/bash
printf "Testing 032.graph..."

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 032 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/lib/ -Wl,-rpath,"${GMML_ROOT_DIR}"/lib/ tests/032.test.ringAtomFinder.cpp -lgmml2 -pthread -o 032.ringAtomFinder.exe

./032.ringAtomFinder.exe tests/inputs/032.4q1p.pdb > 032.output_graph.txt 2>&1

if ! cmp 032.output_graph.txt tests/correct_outputs/032.output_graph.txt >/dev/null 2>&1; then
    printf "Test FAILED! Output file different. Try\ndiff %s %s\n" 032.output_graph.txt tests/correct_outputs/032.output_graph.txt
    echo "Exit Code: 1"
    return 1
    exit 1
fi
#if ! [ -f finished.pdb ]; then
#    echo "Test FAILED! Did not create finished.pdb"
#    echo "Exit Code: 1"
#    return 1
#    exit 1
#fi
printf "Test passed.\n"
rm -r 032.output_graph.txt ./032.ringAtomFinder.exe 
echo "Exit Code: 0"
return 0
