#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH="${GMML_ROOT_DIR}/bin/gpBuilder"

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 017 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

test_case()
{
    prefix=$1
    input=$2
    directory=$3
    rm -r ${directory} >/dev/null 2>&1
    mkdir -p ${directory}
    eval "${BIN_PATH} ${input} ${directory} >${directory}/GlycoproteinBuilder.txt 2>&1"
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
    rm -r ${directory} >/dev/null 2>&1
    failed=0
}

printf "Testing 017.test.GlycoproteinBuilder.cpp... "

test_case "Standard " "tests/inputs/017.GlycoproteinBuilderInput.txt" "017/standard"
if [ "$failed" == 1 ]; then
    return 1
fi
test_case "Skip MD prep " "tests/inputs/017.GlycoproteinBuilderInputNoMDPrep.txt" "017/skipMdPrep"
if [ "$failed" == 1 ]; then
    return 1
fi
test_case "Freeze glycosite conformation " "tests/inputs/017.GlycoproteinBuilderInputFreezeGSConformation.txt" "017/freezeGSConformation"
if [ "$failed" == 1 ]; then
    return 1
fi
test_case "Single glycan " "tests/inputs/017.GlycoproteinBuilderInputSingleGlycan.txt" "017/singleGlycan"
if [ "$failed" == 1 ]; then
    return 1
fi
test_case "Delete incompatible sites " "tests/inputs/017.GlycoproteinBuilderInputDeleteSites.txt" "017/deleteSites"
if [ "$failed" == 1 ]; then
    return 1
fi

printf "Test passed.\n"
echo "Exit Code: 0"
return 0
