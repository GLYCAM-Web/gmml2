## File overview
### Lib files
Lib files contain templates for constructing the atoms of a residue given its name.
The data used is atom names, types, charges, numbers, coordinates, and connectivity.
There is also data for protein residue atom types and charges, used when preparing for molecular dynamics simulations.

### Prep files
[Format description](https://ambermd.org/doc/prep.html)

Prep files are not used directly by gmml2, but are used by the residue combinator in order to generate lib file templates of the residues in the prep file with slight variations in linkages. The lib file templates are then used in carbohydrate construction.

### Pdb files
[Format description](http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)

The protein data bank format is used both as output for writing carbohydrates and glycoproteins, as well as being used as input for the protein component of glycoproteins.

[VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD), [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/download.html) and [Mol*](https://molstar.org/viewer/) are all excellent tools for visualizing pdb files.

#### Pdb preprocessing
After being read, pdb files are processed in order to be easier to work with.
While amino acids with their correct residue and atom names can be bonded unproblematically, there sometimes occur residues that are unknown, misnamed, or malformed. Since these are impossible to account for in every way, we bond atoms by distance as a backup.
When gaps are present in the chain, capping groups are inserted to neutralize charges at the respective terminals. This is important for MD simulations.

### Off files
[Format description](https://ambermd.org/doc/OFF_file_format.txt)

Off files are used as input to molecular dynamics simulations. They contain atom charges, which pdb files do not.
### Dot files
[Format description](https://graphviz.org/doc/info/lang.html)

Dot files are an input format to graphviz 
Paired with the svg images in the SNFG folder, they are used to generate glycan images.
