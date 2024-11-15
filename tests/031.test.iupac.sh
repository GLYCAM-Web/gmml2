#!/bin/bash

printf "Testing 031.test.iupac.cpp... "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 031 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/031.test.iupac.cpp -lgmml2 -pthread -o iupac
./iupac > 031.output.txt

if ! cmp -s 031.output.txt tests/correct_outputs/031.output.txt; then
    echo -e "Test FAILED!\n 031.output.txt is different from tests/correct_outputs/031.output.txt\n"
    echo "Exit Code: 1"
    return 1
fi
printf "Test passed.\n"
rm iupac 031.output.txt
echo "Exit Code: 0"

return 0
