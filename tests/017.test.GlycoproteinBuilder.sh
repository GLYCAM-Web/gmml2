#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH="${GMML_ROOT_DIR}/bin/gpBuilder"

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 017 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

prefix=$1
input=$2
directory=$3
output="output/${directory}"

echo -n "Testing 017.test.GlycoproteinBuilder.cpp ${prefix}..."
rm -r "${output}" >/dev/null 2>&1
mkdir -p "${output}"
"${BIN_PATH}" "${input}" "${output}" "--test-mode" > "${output}"/GlycoproteinBuilder.txt 2>&1
expected="tests/correct_outputs/${directory}"

if [ ! -d "${expected}" ]; then
    echo "Test FAILED"
    echo "directory ${expected} does not exist"
    echo "Exit Code: 1"
    return 1
fi
# note: if diff starts being slow, consider comparing checksums of directory contents instead
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
