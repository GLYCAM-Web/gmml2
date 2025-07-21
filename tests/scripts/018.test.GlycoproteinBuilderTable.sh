#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH="${GMML_ROOT_DIR}/bin/gpBuilderTable"

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 018 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

directory="018"
output="output/${directory}"
expected="correct_outputs/${directory}"
rm -r "${output}" >/dev/null 2>&1
mkdir -p "${output}"

printf "Testing 018.test.createGlycosylationTables.cpp... "
eval "${BIN_PATH} inputs/018.4mbzEdit.pdb --format list > ${output}/GlycoproteinBuilderTableList.txt 2>&1"
eval "${BIN_PATH} inputs/018.4mbzEdit.pdb --format csv > ${output}/GlycoproteinBuilderTable.csv 2>&1"
eval "${BIN_PATH} inputs/018.4mbzEdit.pdb --format txt > ${output}/GlycoproteinBuilderTable.txt 2>&1"
eval "${BIN_PATH} inputs/018.4mbzEdit.pdb --format html > ${output}/GlycoproteinBuilderTable.html 2>&1"

if [ ! -d "${expected}" ]; then
    echo "Test FAILED"
    echo "directory ${expected} does not exist"
    echo "Exit Code: 1"
    return 1
fi
DIFF=$(diff -qr "${output}" "${expected}")
if [ "$DIFF" ]
then
    echo "Test FAILED"
    echo "${DIFF}"
    echo "Exit Code: 1"
    return 1
fi

rm -r "${output}" >/dev/null 2>&1
printf "Test passed.\n"
echo "Exit Code: 0"
return 0
