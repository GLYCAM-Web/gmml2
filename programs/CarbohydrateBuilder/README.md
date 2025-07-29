# CarbohydrateBuilder
Uses GMML to create 3D structures of carbohydrates based on the GLYCAM forcefield.

### Notes
Project is under development, contact olivercgrant "at" gmail.com with queries. 
The underlying code is used by the builder currently available on glycam.org/cb.

### Prerequisites
You'll need GMML2: [Click here for installation instructions](https://github.com/GLYCAM-Web/gmml2#readme)

### Installation
The Carbohydrate Builder will be compiled to gmml2/bin/carbohydrateBuilder after running the gmml2 make.sh script

### Testing
Once compiled, a call to the program can look as follows
```
./bin/carbohydrateBuilder inputFile.txt separator outputDir
```
Or using a test input file
```
cd gmml2/tests
../bin/carbohydrateBuilder inputs/023.smallLibrary.txt _ outputDir
```

You can run `./bin/carbohydrateBuilder --help` for a list of options

```
usage: ./bin/carbohydrateBuilder [-h | --help]
                                 [-v | --version]
                                 [--overwrite-existing-files]
                                 input-file
                                 list-delimiter
                                 output-directory
```