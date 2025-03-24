Usage instructions:
The programs like gpBuilder are compiled and available in bin/.

Setup:

# Create a directory to mount into the container
mkdir ioData

Connect and run the internal test:

# Mount the folder and connect to the container
docker run -it -v ./ioData:/app/data gmml2 bash

# Run the test: For now we are unfortunately very picky about where you run this from.
cd tests/
./../bin/gpBuilder tests/inputs/017.GlycoproteinBuilderInput.txt ../../data/
# It should finish in less than 1 second. Here is what a success looks like:
root@c30afa15eafa:/app/gmml2/tests# ./../bin/gpBuilder tests/inputs/017.GlycoproteinBuilderInput.txt ../../data/
Input file is tests/inputs/017.GlycoproteinBuilderInput.txt
Reading input file complete, on to construction
Resolving overlaps
Program got to end ok
root@c30afa15eafa:/app/gmml2/tests# ls ../../data/
0_glycoprotein.pdb  1_glycoprotein.pdb	glycoprotein.off  glycoprotein.pdb  glycoprotein_initial.pdb  glycoprotein_serialized.pdb  tmp.txt
root@c30afa15eafa:/app/gmml2/tests#

More realistic user example run from the mounted folder:

Step 1. Place a file called exampleInput.txt with the following content into the shared folder inputs/

Protein:/app/data/1rbj.pdb

persistCycles:10

seed:0

deleteIncompatibleSites:true

ProteinResidue, GlycanName:
A_34|DManpa1-2DManpa1-6[DManpa1-2DManpa1-3]DManpa1-6[DManpa1-2DManpa1-2DManpa1-3]DManpb1-4DGlcpNAcb1-4DGlcpNAcb1-OH
END

Step 2. Download the pdb file for 1rbj from the protein data bank (https://www.rcsb.org/structure/1RBJ⁠

) and save it into the mounted folder.

Step 3. Connect to the container and run the gp builder:

docker run -it -v ./inputs:/app/data gmml2 bash
cd /app/data/
/app/gmml2/bin/gpBuilder exampleInput.txt

Step 4. Your results should be in a file called glycoprotein.pdb.
