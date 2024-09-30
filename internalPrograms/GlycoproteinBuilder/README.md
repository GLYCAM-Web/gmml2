# GlycoProteinBuilder
Uses GMML2 to add and adapt 3D structures of N-glycans and O-glycans onto glycoproteins. It can do this for Asn, Ser, Thr and Tyr.

## General Concept
![schematic](schematic/schematic.png)

### Notes
Project is under development, contact olivercgrant "at" gmail.com with queries. 
This code is intended to replace the glycoprotein builder currently available on glycam.org/gp. You can compile and run it locally.
Has been tested on linux, but might install on both Mac and Windows.

### Prerequisites

You'll need GMML2: [Click here for installation instructions](https://github.com/GLYCAM-Web/gmml2#readme)

### Installation of GlycoProteinBuilder
    cd gmml2/
    GMML_ROOT_DIR=$(git rev-parse --show-toplevel) # assumes you have cloned the gmml repo using git.
    g++ -std=c++17 -g -I "${GMML_ROOT_DIR}" -L"${GMML_ROOT_DIR}"/bin/ -Wl,-rpath,"${GMML_ROOT_DIR}"/bin/ "${GMML_ROOT_DIR}"/internalPrograms/GlycoproteinBuilder/gpBuilder_main.cpp -lgmml2 -pthread -o gpBuilder

### Testing
Once compiled, you can run:

./gpBuilder tests/inputs/017.GlycoproteinBuilderInput.txt 

### Setup
Edit or create an input.txt file. See tests/inputs/017.GlycoproteinBuilderInput.txt for an example.

You must provide:

    A protein 3D structure
    A Glycan 3D structure(s) or sequences in GLYCAM condensed nomenclature (just like the carb builder here: glycam.org/cb)
    An input.txt, which contains:
        protein file name
        the protein residue numbers you want to attach to (no automatic detection of sequons)
        the glycan you want to attach in Glycam condensed sequence format.

NumberOfOutputStructures -> produces different shapes for the glycan at each site.
maxThreads -> currently does nothing.
persistCycles -> how long to keep trying to find an improvement in overlap. A higher value takes more time, but may better resolve overlaps.
seed:0 -> means you will get the same result every time you run the program with the same inputs. Useful for testing and reproducibility. You can delete this line if you want it to be random.
