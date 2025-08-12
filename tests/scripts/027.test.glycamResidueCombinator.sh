#!/bin/bash

printf "Testing 027 glycam residue combinator... "

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
BIN_PATH="${GMML_ROOT_DIR}/bin/glycamResidueCombinator"
PREP="${GMML_ROOT_DIR}/dat/prep"
BACT="${PREP}/bacterial"

if [[ "${GMML_ROOT_DIR}" != *"gmml2" ]]; then
    echo -e "Test 027 failed, we think our GMML root directory is:\t${GMML_ROOT_DIR}\n"
    exit 1
fi

outputFile=027.output.txt
"${BIN_PATH}" "${PREP}"/GLYCAM_06j-1_GAGS_KDN.prep "${BACT}"/0aN-alfa_final.prep "${BACT}"/0an-beta_final.prep "${BACT}"/0DH-alfa_final.prep "${BACT}"/0Dh-beta_final.prep "${BACT}"/0eC-alfa_final.prep "${BACT}"/0ec-beta_final.prep "${BACT}"/0FC-alfa_final.prep "${BACT}"/0Fc-beta_final.prep "${BACT}"/0gF-alfa_final.prep "${BACT}"/0gf-beta_final.prep "${BACT}"/0KX-alfa_final.prep "${BACT}"/0Kx-beta_final.prep "${BACT}"/0LD_final.prep "${BACT}"/0LG-alfa_final.prep "${BACT}"/0Lg-beta_final.prep "${BACT}"/0LH-alfa_final.prep "${BACT}"/0Lh-beta_final.prep "${BACT}"/0LU_final.prep "${BACT}"/0mP-alfa_final.prep "${BACT}"/0mp-beta_final.prep "${BACT}"/0MR-alfa_final.prep "${BACT}"/0Mr-beta_final.prep "${BACT}"/0QF-alfa_final.prep "${BACT}"/0Qf-beta_final.prep "${BACT}"/0ZF-alfa_final.prep "${BACT}"/0Zf-beta_final.prep "${BACT}"/0KO-alfa_final.prep "${BACT}"/0Ko-beta_final.prep > $outputFile

file="GLYCAM_06k.lib"
if [ ! -f "${file}" ]; then
    echo -e "Test FAILED!\n ${file} does not exist\n"
    echo "Exit Code: 1"
     return 1
fi
if ! cmp -s "${file}" ../dat/parameters/"${file}"; then
    echo -e "Test FAILED!\n ${file} is different from ../dat/parameters/${file}\n"
    echo "Exit Code: 1"
     return 1
fi
  
sed -n '/\!entry/q;p' GLYCAM_06k.lib | sed 's/ "//g' | sed 's/"//g' | grep -v "index" >libEntries.txt
grep INT ../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep | cut -d \  -f1 >prepEntriesList.txt
for id in `cat libEntries.txt`
do
    if ! grep -q "^$id" ../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep;then
        echo "$id is new in lib file" >>$outputFile; 
    fi;
done

uniq -D libEntries.txt >>$outputFile

for id in `cat prepEntriesList.txt`; 
do 
    if ! grep -q " \"$id\"$" GLYCAM_06k.lib; 
    then 
        echo "$id not found in lib file" >>$outputFile  
    fi;     
done

if ! cmp -s $outputFile correct_outputs/$outputFile; then
    echo -e "Test FAILED!\nOutput files different. Try:\ndiff $outputFile correct_outputs/$outputFile"
    echo "Exit Code: 1"
    return 1    
fi
    
printf "Test passed.\n"
rm $outputFile libEntries.txt prepEntriesList.txt $file
echo "Exit Code: 0"

return 0
