#!/bin/bash

printf "Testing 030.test.gmPreProcessor.cpp... ~3 seconds. "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH="${GMML_ROOT_DIR}/bin/glycomimeticPreprocessor"

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 030 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

directory="030"
output="output/${directory}"

rm -r "${output}" >/dev/null 2>&1
mkdir -p "${output}"

shopt -s nullglob
for filepath in inputs/"${directory}"/*.pdb; do
    file=$(basename "${filepath}")
    filename="${file%.*}"
    "${BIN_PATH}" "inputs/${directory}/${file}" "${output}/${file}" > "${output}/${filename}.txt"
done

expected="correct_outputs/${directory}"

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
