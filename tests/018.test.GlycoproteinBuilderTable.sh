#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH="${GMML_ROOT_DIR}/bin/gpBuilderTable"

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 018 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

printf "Testing 018.test.createGlycosylationTables.cpp... "
eval "${BIN_PATH} tests/inputs/018.4mbzEdit.pdb >018.output_GlycoproteinBuilderTable.txt 2>&1"
if ! cmp 018.output_GlycoproteinBuilderTable.txt tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt >/dev/null 2>&1; then
    printf "Test FAILED!. tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt different from 018.output_GlycoproteinBuilderTable.txt\ndiff tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt 018.output_GlycoproteinBuilderTable.txt\n"
    diff tests/correct_outputs/018.output_GlycoproteinBuilderTable.txt 018.output_GlycoproteinBuilderTable.txt
    echo "Exit Code: 1"
    return 1
else
    printf "Test passed.\n"
    rm 018.output_GlycoproteinBuilderTable.txt
    echo "Exit Code: 0"
    return 0
fi
