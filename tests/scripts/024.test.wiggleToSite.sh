#!/bin/bash
printf "Testing 024.wiggleToSite..."

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH=${GMML_ROOT_DIR}/bin/wiggleToSite

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 024 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

eval "${BIN_PATH} inputs/024.input.txt >024.output_wiggleToSite.txt 2>&1"

if ! cmp 024.output_wiggleToSite.txt correct_outputs/024.output_wiggleToSite.txt >/dev/null 2>&1; then
    printf "Test FAILED! Output file %s different from %s \n" 024.output_wiggleToSite.txt correct_outputs/024.output_wiggleToSite.txt
    echo "Exit Code: 1"
    return 1
    exit 1
fi
if ! [ -f finished.pdb ]; then
    echo "Test FAILED! Did not create finished.pdb"
    echo "Exit Code: 1"
    return 1
    exit 1
fi
printf "Test passed.\n"
rm -r 024.output_wiggleToSite.txt initial.pdb initial.off superimposed.off superimposed.pdb finished.pdb finished.off
echo "Exit Code: 0"
return 0
