#!/bin/bash

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 017 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

printf "Testing 017.test.GlycoproteinBuilder.cpp... "
g++ -std=c++17 -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ "${GMML_ROOT_DIR}"/internalPrograms/GlycoproteinBuilder/gpBuilder_main.cpp -lgmml2 -pthread -o gpBuilder
rm -r 017/ >/dev/null 2>&1
directory="017/standard"
mkdir -p ${directory}
./gpBuilder tests/inputs/017.GlycoproteinBuilderInput.txt ${directory} >${directory}/GlycoproteinBuilder.txt 2>&1
fileList=("glycoprotein_initial.pdb" "glycoprotein.pdb" "0_glycoprotein.pdb" "1_glycoprotein.pdb" "glycoprotein.off" "glycoprotein_serialized.pdb" "GlycoproteinBuilder.txt")
for file in "${fileList[@]}"; do
    output="${directory}/${file}"
    if [ ! -f "${output}" ]; then
        echo -e "Test FAILED!\n ${output} does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi
    correct="tests/correct_outputs/${output}"
    if ! cmp "${output}" "${correct}" >/dev/null 2>&1; then
        echo -e "Test FAILED!\n ${output} is different from ${correct}\n"
        echo "Exit Code: 1"
        return 1
    fi
done
# Do it again but skip preprocessing
echo "Dingo"
directory="017/skipMdPrep"
mkdir -p ${directory}
./gpBuilder tests/inputs/017.GlycoproteinBuilderInputNoMDPrep.txt ${directory} >${directory}/GlycoproteinBuilder.txt 2>&1
fileList=("glycoprotein_initial.pdb" "glycoprotein.pdb" "0_glycoprotein.pdb" "1_glycoprotein.pdb" "glycoprotein.off" "glycoprotein_serialized.pdb" "GlycoproteinBuilder.txt")
for file in "${fileList[@]}"; do
    output="${directory}/${file}"
    if [ ! -f "${output}" ]; then
        echo -e "SkipMdPrep Test FAILED!\n ${output} does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi
    correct="tests/correct_outputs/${output}"
    if ! cmp "${output}" "${correct}" >/dev/null 2>&1; then
        echo -e "SkipMdPrep Test FAILED!\n ${output} is different from ${correct}\n"
        echo "Exit Code: 1"
        return 1
    fi
done
directory="017/freezeGSConformation"
mkdir -p ${directory}
./gpBuilder tests/inputs/017.GlycoproteinBuilderInputFreezeGSConformation.txt ${directory} >${directory}/GlycoproteinBuilder.txt 2>&1
fileList=("glycoprotein_initial.pdb" "glycoprotein.pdb" "0_glycoprotein.pdb" "1_glycoprotein.pdb" "glycoprotein.off" "glycoprotein_serialized.pdb" "GlycoproteinBuilder.txt")
for file in "${fileList[@]}"; do
    output="${directory}/${file}"
    if [ ! -f "${output}" ]; then
        echo -e "FreezeGSConformation Test FAILED!\n ${output} does not exist\n"
        echo "Exit Code: 1"
        return 1
    fi
    correct="tests/correct_outputs/${output}"
    if ! cmp "${output}" "${correct}" >/dev/null 2>&1; then
        echo -e "FreezeGSConformation Test FAILED!\n ${output} is different from ${correct}\n"
        echo "Exit Code: 1"
        return 1
    fi
done
rm -r 017/ >/dev/null 2>&1
printf "Test passed.\n"
rm gpBuilder
echo "Exit Code: 0"
return 0
