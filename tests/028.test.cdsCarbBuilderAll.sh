#!/bin/bash
echo -n "Testing 028 carbBuilderAll..."

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 028 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

input=$1
directory=$2
output="output/${directory}"

g++ -std=c++17 -I "${GMML_ROOT_DIR}"/ -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ tests/028.test.cdsCarbBuilderAll.cpp -lgmml2 -pthread -o 028.carbBuilder.exe
rm -r "${output}" >/dev/null 2>&1
mkdir -p ${output}
./028.carbBuilder.exe "${input}" "_" "${output}" > "${output}"/carbohydrateBuilder.txt 2>&1

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

rm -r "${output}" 028.carbBuilder.exe >/dev/null 2>&1
printf "Test passed.\n"
echo "Exit Code: 0"
return 0
