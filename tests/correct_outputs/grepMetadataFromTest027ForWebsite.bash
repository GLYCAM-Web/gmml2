#!/bin/bash
for name in $(grep "Found name" 027.output.txt | cut -d _ -f2 | sort --unique )
do
    pPositions=$(grep -A2 "_${name}_p" 027.output.txt | tail -n1)
    fPositions=$(grep -A2 "_${name}_f" 027.output.txt | tail -n1)
    ano=$(grep -A1 "_${name}_" 027.output.txt | tail -n1 | cut -d \  -f2)
    if [ ! -z "${pPositions}" ]; then p="p"; else p=""; fi;
    if [ ! -z "${fPositions}" ]; then f="f"; else f=""; fi;
    grep -B2 "_${name}_" 027.output.txt > info.txt
    if grep -q "is D_._." info.txt; then d="D"; else d=""; fi;
    if grep -q "is L_._." info.txt; then l="L"; else l=""; fi;
    if grep -q "is ._._a" info.txt; then a="a"; else a=""; fi;
    if grep -q "is ._._b" info.txt; then b="b"; else b=""; fi;
    echo "$name|$d|$l|$p|$f|$a|$b|$ano|$pPositions|$fPositions|"
done
rm info.txt

# Apr 2025
# Oliver made this one and done script just to grab out the metadata Dan needs in the website. This does not set defaults for any of the values. Oliver takes this and pastes it into a google drive sheet called Residue Metadata.
# There are some spurious "pa" at the end of some of the names and "p" in the middle of some due to the whacky state of the gmml metadata. I manually remove those in the google sheet.
# I want a record of this, and would need to run this again if new residues get added to the website. 
