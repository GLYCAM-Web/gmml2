# GMML2
The GLYCAM Molecular Modeling Library (GMML) is a C++ library created by the Woods Group. This is the scientific code underlying much of GLYCAM-Web (glycam.org). 
GMML2 is a new version of GMML that has replaced much of the old codebase. Some tools like Glyfinder (glycam.org/gf) and the glycomimetic builder use the original GMML.
Note that the GMML2 code was originally developed within GMML1. We then forked to this repo and deleted the old code from here.

[Overview](#overview)

[Prerequisites](#prerequisites)

[Obtaining the software](#obtaining-the-software)

[Compiling the Library](#compiling-the-library)

[Using the Internal Programs](#using-the-internal-programs)

[Testing the Library](#testing-the-library)

[For Developers](#for-developers)

[Coding Standards](#coding-standards)

---
## Overview

GMML2 provides a library for common molecular modeling tasks.  It is 
particularly well-tuned for models of carbohydrates and systems that
contain carbohydrates.

For more detailed descriptions of how it works, have a look at our [developer guide](https://glycam-web.github.io/gmml2/). 

### Used by [GLYCAM-Web](https/glycam.org)

This code also serves as the main molecular modeling engine for GLYCAM-Web.  

### Funding Sources

We are very grateful to our funders.  
[Please check them out!](https://github.com/GLYCAM-Web/website/blob/master/funding.md)

### Acknowledgements

GMML2 includes the following third-party data and source code:

A curated data set of the [Smooth Backbone-Dependent Rotamer Library 2010](http://dunbrack.fccc.edu/lab/bbdep2010) licensed under [Creative Commons CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode)

The [PCG Library](https://www.pcg-random.org/index.html), licensed under the [Apache License](https://www.apache.org/licenses/LICENSE-2.0)

## Prerequisites

### Building GMML2

In order to build GMML2, you are required to have the following software available on your system:

* `cmake` (Version >= `3.13.4`)
* `g++` (Version >= `7.0`)
* `openmp`
* `make`
* `git`

Installation instructions will vary according to what package manager your distro uses. If you are using apt as a package manager on a linux system, you should be able to use a command like this:

```bash
sudo apt-get update &&\
sudo apt-get install git git-all cmake g++ libgomp1
```
For other linux distros, please follow the instructions for the package managment software included with your system.

---
## Obtaining the software
The following does not require root access, but it does require one has `git` installed.

1. Navigate to the directory that you would like to have `gmml2` live.

2. Clone `gmml2` from the git repo
```bash
git clone https://github.com/GLYCAM-Web/gmml2.git
```

**NOTE:** There are non-git ways to obtain gmml2. Don't do this as our tests won't compile outside of a git repo.

---
## Compiling the Library

Make sure you are in the gmml2 folder. To control the number of processors used during the *`make`* process, use the `-j` flag for our `make.sh`, so to run with 8 cores we would run `./make.sh -j 8`.

This will create the needed `cmake` files and will add the following directories within the `gmml2` directory:

* `lib` (Contains the `gmml2` shared object libary, `libgmml2.so`)
* `bin` (Contains the compiled gmml2-dependent programs: carbohydrateBuilder, gpBuilder (glycoprotein builder), gpBuilderTable (glycosylation site finder), wiggleToSite)
* `cmakeBuild` (Contains all files produced by the `cmake` command, a `compile_commands.json` file to be used with tools of your choice, and all files contained within the directories listed above)

You can either use the `libgmml2.so` file within the `lib` directory or the `libgmml2.so` file within the `cmakeBuild` directory. They are the exact same.

Please enter `./make.sh -h` for help regarding the make script.

---
## Using the Internal Programs

The GMML2 internal programs will be available in gmml2/bin/ after running the make.sh script

It is recommended that you test them (see next section) before use

[Glycoprotein Builder Instructions](programs/GlycoproteinBuilder/README.md)

[Carbohydrate Builder Instructions](programs/CarbohydrateBuilder/README.md)

---
## Testing the Library

From within the `gmml2` directory, you must change your current working directory to the `gmml2/tests` directory. Note that `<NUM_JOBS>` is however many tests you want to run at once.

```bash
gmml2$ cd tests/
gmml2/tests$ ./compile_run_tests.bash -j <NUM_JOBS>
```
The output will tell you whether or not the library is behaving appropriately and if all tests are passed the output will look similar to the following:

```bash
$ bash compile_run_tests.bash -j4

#### Beginning GMML tests ####
Number of tests found:	14
Number of testing jobs:	4

mkdir: created directory './tempTestOutputs'

Beginning test: ./016.test.DrawGlycan.sh
Beginning test: ./017.test.GlycoproteinBuilder.sh
Beginning test: ./017b.test.GlycoproteinBuilderFailure.sh
Beginning test: ./018.test.GlycoproteinBuilderTable.sh

Testing 017b.test.GlycoproteinBuilderFailure.cpp... Test passed.
Exit Code: 0

Beginning test: ./019.test.newPDBClass.sh

Testing 018.test.createGlycosylationTables.cpp... Test passed.
Exit Code: 0

Beginning test: ./020.test.parameterFiles.sh

Testing 016.test.DrawGlycan.cc...0.svg correct_outputs/016.output_SVGs/0.svg differ: byte 15325, line 70
Test FAILED! Output file 0.svg different to correct_outputs/016.output_SVGs/0.svg
Exit Code: 1

Beginning test: ./022.test.libraryFileReader.sh

Testing 020.test.parameterFiles.cpp... Test passed.
Exit Code: 0

Beginning test: ./023.test.carbohydrateBuilder.sh

Testing 022.test.libraryFileReader.cpp... Test passed.
Exit Code: 0

Beginning test: ./024.test.wiggleToSite.sh

Testing 024.wiggleToSite...Test passed.
Exit Code: 0

Beginning test: ./026.test.editPdbFile.sh

Testing 017.test.GlycoproteinBuilder.cpp... Test passed.
Exit Code: 0

Beginning test: ./027.test.glycamResidueCombinator.sh

Testing 026.test.editPDB.cpp... ~2 seconds. Test passed.
Exit Code: 0

Beginning test: ./028.test.cdsCarbBuilderAll.sh

Testing 023.carbohydrateBuilder... Test passed.
Exit Code: 0

Beginning test: ./030.test.gmPreProcessor.sh

Testing 030.test.gmPreProcessor.cpp... ~3 seconds. Test passed.
Exit Code: 0

Beginning test: ./031.test.iupac.sh

Testing 031.test.iupac.cpp... Test passed.
Exit Code: 0

Testing 028.test.cdsCarbBuilderAll.cpp...Test passed.
Exit Code: 0

Testing 027.test.glycamResidueCombinator.cpp... Test passed.
Exit Code: 0

Testing 019.test.newPDBClass.cpp... ~30 seconds. Test passed.
Exit Code: 0

######## GMML TESTS COMPLETED ########
Required tests:	14
Passed tests:	13
Failed tests:	1
Time taken:	11 seconds
######################################

!!! OUTPUT OF THE 1 GMML TEST(S) THAT FAILED !!!

Testing 016.test.DrawGlycan.cc...0.svg correct_outputs/016.output_SVGs/0.svg differ: byte 15325, line 70
Test FAILED! Output file 0.svg different to correct_outputs/016.output_SVGs/0.svg
Exit Code: 1

!!! FINISHED PRINTING FAILED TESTS !!!
```

If any tests fail then something is wrong.
