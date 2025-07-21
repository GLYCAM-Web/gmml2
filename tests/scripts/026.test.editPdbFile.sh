#!/bin/bash

printf "Testing 026.test.editPDB.cpp... ~2 seconds. "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 026 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/lib/ -Wl,-rpath,"${GMML_ROOT_DIR}"/lib/ src/026.test.editPDB.cpp -lgmml2 -pthread -o 026.editPdb
shopt -s nullglob
for filepath in inputs/026.*.pdb; do
    file=$(basename "${filepath}")
    ./026.editPdb inputs/"${file}" >026.output.txt
    if ! cmp 026.outputPdbFile.pdb correct_outputs/"${file}"-output.pdb >/dev/null 2>&1; then
        echo -e "Test FAILED!. 026.outputPdbFile.pdb different from correct_outputs/${file}-output.pdb\n Compare using diff or VMD\n"
        echo "Exit Code: 1"
        return 1
    elif ! cmp 026.outputSelection.pdb correct_outputs/"${file}"-outputSelection.pdb >/dev/null 2>&1; then
        echo -e "Test FAILED!. 026.outputSelection.pdb different from correct_outputs/${file}-outputSelection.pdb\n Compare using diff or VMD\n"
        echo "Exit Code: 1"
        return 1
    fi
done
printf "Test passed.\n"
rm ./026.editPdb 026.outputPdbFile.pdb 026.output.txt 026.outputSelection.pdb
echo "Exit Code: 0"
return 0
