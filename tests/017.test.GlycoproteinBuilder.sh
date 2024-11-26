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

echo -n "Testing 017.test.GlycoproteinBuilder.cpp ${prefix}..."
rm -r "${directory}" >/dev/null 2>&1
mkdir -p "${directory}"
"${BIN_PATH}" "${input}" "${directory}" > "${directory}"/GlycoproteinBuilder.txt 2>&1
expected="tests/correct_outputs/${directory}"

if [ ! -d "${expected}" ]; then
    echo "Test FAILED"
    echo "directory ${expected} does not exist"
    echo "Exit Code: 1"
    return 1
fi
# note: if diff starts being slow, consider comparing checksums of directory contents instead
DIFF=$(diff -qr "${directory}" "${expected}")
if [ "$DIFF" ]
then
    echo "Test FAILED"
    echo "Difference between actual and expected output"
    echo "Run the following command for details:"
    echo "    diff -qr ${directory} ${expected}"
    echo "Exit Code: 1"
    return 1
fi

rm -r "${directory}" >/dev/null 2>&1
printf "Test passed.\n"
echo "Exit Code: 0"
return 0
