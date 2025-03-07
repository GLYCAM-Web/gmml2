#!/bin/bash

printf "Testing 030.test.gmPreProcessor.cpp... ~3 seconds. "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 030 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/lib/ -Wl,-rpath,"${GMML_ROOT_DIR}"/lib/ ../internalPrograms/glycomimeticPreprocessor/glycomimeticPreprocessor.cpp -lgmml2 -pthread -o gmPreProcessor.exe
shopt -s nullglob
for filepath in tests/inputs/030.*.pdb; do
    file=$(basename "${filepath}")
    ./gmPreProcessor.exe tests/inputs/"${file}" 030.outputPdbFile.pdb >030.output.txt
    if ! cmp 030.output.txt tests/correct_outputs/"${file}"-output.txt >/dev/null 2>&1; then
        echo -e "Test FAILED!. 030.output.txt different from tests/correct_outputs/${file}-output.txt\n Compare using diff\n"
        echo "Exit Code: 1"
        return 1
        #exit 1
    elif ! cmp 030.outputPdbFile.pdb tests/correct_outputs/"${file}"-output.pdb >/dev/null 2>&1; then
        echo -e "Test FAILED!. 030.outputPdbFile.pdb different from tests/correct_outputs/${file}-output.pdb\n Compare using diff or VMD\n"
        echo "Exit Code: 1"
        return 1
        #exit 1
    fi
done
printf "Test passed.\n"
rm ./gmPreProcessor.exe 030.outputPdbFile.pdb 030.output.txt
echo "Exit Code: 0"
return 0
