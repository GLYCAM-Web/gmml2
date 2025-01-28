#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH="${GMML_ROOT_DIR}/bin/drawGlycan"

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 016 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

filename=$1
sequence=$2
directory="016"
output="output/${directory}/${filename}"

echo -n "Testing 016 drawGlycan ${sequence}... "

rm "${output}" >/dev/null 2>&1
mkdir -p "output/${directory}"
"${BIN_PATH}" "${sequence}" "${output}" --base-dir "${GMML_ROOT_DIR}/" --relative-paths 2>&1

if ! cmp "${output}" "tests/correct_outputs/${directory}/${filename}" >/dev/null 2>&1; then
    echo "Test FAILED! Output file different"
    echo "Exit Code: 1"
    return 1
else
    echo "Test passed."
    rm "${output}"
    echo "Exit Code: 0"
    return 0
fi
