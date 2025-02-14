# Glycoprotein Builder
Uses GMML2 to add and adapt 3D structures of N-glycans and O-glycans onto glycoproteins. It can do this for Asn, Ser, Thr and Tyr.

## General Concept
![schematic](schematic/schematic.png)

### Notes
Project is under development, contact olivercgrant "at" gmail.com with queries. 
This code is intended to replace the glycoprotein builder currently available on glycam.org/gp. You can compile and run it locally.
Only tested on Linux.

### Prerequisites

You'll need GMML2: [Click here for installation instructions](https://github.com/GLYCAM-Web/gmml2#readme)

### Installation
The Glycoprotein Builder will be compiled to gmml2/bin/gpBuilder after running the gmml2 make.sh script

### Testing
Once compiled, a call to the program can look as follows
```
./bin/gpBuilder inputFile.txt outputDir
```
Or using a test input file
```
cd gmml2/tests
../bin/gpBuilder tests/inputs/017.GlycoproteinBuilderInput.txt outputDir
```

You can run `./bin/gpBuilder --help` for a list of options

```
usage: ./bin/gpBuilder [-h | --help]
                       [-v | --version]
                       [-t <value> | --num-threads <value>]
                       [--overwrite-existing-files]
                       input-file
                       [output-directory]
```

### Setup
Edit or create an input.txt file. See gmml2/tests/tests/inputs/017.GlycoproteinBuilderInput.txt for an example.

Required input:
```
A protein 3D structure in .pdb format
An input.txt, containing:
  Protein: file name of the protein structure

  ProteinResidue, GlycanName:
  protein residue number | glycan sequence
  ...
  protein residue number | glycan sequence
  END

Any number of protein residues | glycan pairs can be entered
Glycans are defined in Glycam condensed sequence format (see the carb builder here: glycam.org/cb)
```

Additional settings:

| Setting  | Default | Explanation |
| --- | --- | --- |
| NumberOfOutputStructures | 1 | Each output structure is an independent sample of possible glycan shapes. Rotamers are randomized according to their statistical likelihood |
| persistCycles | 5 | How long the algorithm persists in looking for a solution with a lower number of overlaps. Higher values take more time, but may yield better results |
| seed | none | Using a specific seed will provide reproducible output across different computers and different runs. Delete this line for a random seed to be used |
| residueOverlapRejectionThreshold | none | Structures with one or more glycan residues having overlaps with this many atoms in other residues will be placed in a /rejected directory. Delete this line to not reject any structures |
| atomOverlapTolerance | 0.6 | How much any pair of atoms are allowed to overlap into their vdW radii before the algorithm tries to separate them. Higher values are more permissive |
| freezeGlycositeResidueConformation | false | Enable to preserve chi1 and chi2 angles of protein-glycan linkages according to their initial shape |
| allowSidechainAdjustment | false | Enable to adjust protein sidechains from their initial shape if and only if it allows for a greater range of glycan sampling |
| deleteIncompatibleSites | false | Enable to delete glycans where overlaps fail to resolve in order to produce a structure with no overlaps |
