# Oliver wrote this 2024-10-14 to update the metadata with the residues from the new bacterial sugars.
# This used to be updated manually, but now we have a record:
# !/bin/bash


localFile=All_Residues.txt
libFile=../../../CurrentParams/GLYCAM_06k.lib

echo "Checking for entries in $localFile that are not present in $libFile. If there are any this is bad:"
for id in $(cat $localFile); 
do 
    if ! grep -q "\"$id\"" $libFile; then 
        echo "$libFile does not contain $id found in $localFile"; 
    fi 
done
echo "Finished. If there were outputs you need to figure that out."
echo "Now checking for entries in $libFile that are not in $localFile, meaning that $localFile needs updating"
>addMeToBottomOfLocalFile.txt
for id in $(sed '/\!entry/q' $libFile | grep "^ " | sed 's/"//g' | sed 's/ //g' )
do
    if ! grep -q "^$id$" $localFile; then
	echo "$localFile does not contain $id found in $libFile"
	echo "$id" >> addMeToBottomOfLocalFile.txt
    fi
done
echo "Done. If there were outputs add them to the bottom of $localFile and run the scripts here."
