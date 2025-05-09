#!/bin/bash

printf "Testing 019.test.newPDBClass.cpp... ~30 seconds. "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 019 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

directory="019"
output="output/${directory}"

rm -r "${output}" >/dev/null 2>&1
mkdir -p "${output}"

g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/lib/ -Wl,-rpath,"${GMML_ROOT_DIR}"/lib/ tests/019.test.newPDBClass.cpp -lgmml2 -pthread -o newPdbClass
shopt -s nullglob
for filepath in tests/inputs/"${directory}"/*.pdb; do
    file=$(basename "${filepath}")
    filename="${file%.*}"
    ./newPdbClass "tests/inputs/${directory}/${file}" "${output}/${file}" > "${output}/${filename}.txt"
done

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
rm ./newPdbClass
printf "Test passed.\n"
echo "Exit Code: 0"
return 0
