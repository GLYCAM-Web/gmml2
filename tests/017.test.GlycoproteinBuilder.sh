#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 017 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

test_case()
{
    prefix=$1
    input=$2
    directory=$3
    mkdir -p ${directory}
    ./gpBuilder ${input} ${directory} >${directory}/GlycoproteinBuilder.txt 2>&1
    fileList=("glycoprotein_initial.pdb" "glycoprotein.pdb" "0_glycoprotein.pdb" "1_glycoprotein.pdb" "glycoprotein.off" "glycoprotein_serialized.pdb" "GlycoproteinBuilder.txt")
    for file in "${fileList[@]}"; do
        output="${directory}/${file}"
        if [ ! -f "${output}" ]; then
            echo -e "${prefix} Test FAILED!\n ${output} does not exist\n"
            echo "Exit Code: 1"
            failed=1
            return 1
        fi
        correct="tests/correct_outputs/${output}"
        if ! cmp "${output}" "${correct}" >/dev/null 2>&1; then
            echo -e "${prefix} Test FAILED!\n ${output} is different from ${correct}\n"
            echo "Exit Code: 1"
            failed=1
            return 1
        fi
    done
    failed=0
}

printf "Testing 017.test.GlycoproteinBuilder.cpp... "
g++ -std=c++17 -g -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ "${GMML_ROOT_DIR}"/internalPrograms/GlycoproteinBuilder/gpBuilder_main.cpp -lgmml2 -pthread -o gpBuilder
rm -r 017/ >/dev/null 2>&1

test_case "" "tests/inputs/017.GlycoproteinBuilderInput.txt" "017/standard"
if [ "$failed" == 1 ]; then
    return 1
fi
test_case "SkipMdPrep " "tests/inputs/017.GlycoproteinBuilderInputNoMDPrep.txt" "017/skipMdPrep"
if [ "$failed" == 1 ]; then
    return 1
fi
test_case "FreezeGSConformation " "tests/inputs/017.GlycoproteinBuilderInputFreezeGSConformation.txt" "017/freezeGSConformation"
if [ "$failed" == 1 ]; then
    return 1
fi

rm -r 017/ >/dev/null 2>&1
printf "Test passed.\n"
rm gpBuilder
echo "Exit Code: 0"
return 0
