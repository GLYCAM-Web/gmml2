#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH="${GMML_ROOT_DIR}/bin/gpBuilder"

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 017b failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

printf "Testing 017b.test.GlycoproteinBuilderFailure.cpp... "
output=017.failure.txt
correct=tests/correct_outputs/$output
eval "${BIN_PATH} tests/inputs/017.GlycoproteinBuilderInputDosAndError.txt > $output 2>&1"
if ! cmp "${output}" "${correct}" >/dev/null 2>&1; then
    echo -e "Failure test failed to fail as expected!\n ${output} is different from ${correct}\n"
    echo "Exit Code: 1"
    return 1
fi
rm ${output}

printf "Test passed.\n"
echo "Exit Code: 0"
return 0
