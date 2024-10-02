# GMML2
The GLYCAM Molecular Modeling Library (GMML) is a C++ library created by the Woods Group. This is the scientific code underlying much of GLYCAM-Web (glycam.org). 
GMML2 is a new version of GMML that has replaced much of the old codebase. Some tools like Glyfinder (glycam.org/gf) and the glycomimetic builder use the original GMML.
Note that the GMMl2 code was originally developed within GMML1. We then forked to this repo and deleted the old code from here.

[Overview](#overview)

[Prerequisites](#prerequisites)

[Obtaining the software](#obtaining-the-software)

[Compiling the Library](#compiling-the-library)

[Testing the Library](#testing-the-library)

[Coding Standards](#coding-standards)

[Documentation](#documentation)

---
## Overview

GMML2 provides a library for common molecular modeling tasks.  It is 
particularly well-tuned for models of carbohydrates and systems that
contain carbohydrates.

### Used by [GLYCAM-Web](https/glycam.org)

This code also serves as the main molecular modeling engine for GLYCAM-Web.  

### Funding Sources

We are very grateful to our funders.  
[Please check them out!](https://github.com/GLYCAM-Web/website/blob/master/funding.md)


## Prerequisites

### Building GMML

In order to build GMML, you are required to have the following software available on your system:

* `libssl1.1`
* `libssl-dev`
* `cmake` (Version >= `3.13.4`)
* `g++` (Version >= `7.0`)
* `make`
* `python3.9` (Version `3.9.12`)
* `python3.9-dev` (Version `3.9.12`)
* `git`
* `swig` (Version `4.0.2`)
* `libeigen3-dev` (Version >= `3.3.7`)

Installation instructions will vary according to what package manager your distro uses. If you are using apt as a package manager on a linux system, you should be able to use a command like this:

```bash
sudo apt-get update &&\
sudo apt-get install libssl1.1 libssl-dev git python3.9 python3.9-dev cmake g++ git-all libeigen3-dev
```
For other linux distros, please follow the instructions for the package managment software included with your system.

Most linux distros have swig 4.0. Otherwise swig 4.0.2 must be installed from [their website](https://www.swig.org/download.html). 

### Contributing to GMML

If you want to contribute to `gmml2` you will also need to install the following packages:

* `clang-tidy-15`
* `shellcheck`
* `shfmt`

---
## Obtaining the software
The following does not require root access, but it does require one has `git` installed.

1. Navigate to the directory that you would like to have `gmml2` live. Please note that in order to use the produced library with [`gems`](https://github.com/glycam-web/gems/) the `gmml2` directory must be placed within the `gems` directory.

2. Clone `gmml2` from the git repo
```bash
git clone https://github.com/GLYCAM-Web/gmml2.git -b feature_ChangesForFork
```

**NOTE:** There are non-git ways to obtain gmml2. Don't do this as our tests won't compile outside of a git repo.
**NOTE:** The -b feature_ChangesForFork should be temporary. If it's not working we likely forgot to update this doc. Try remove "-b feature_ChangesForFork" from the command.

---
## Compiling the Library
```
Make sure you are in the gmml2 folder. To control the number of processors used during the *`make`* process, use the `-j` flag for our `make.sh`, so to run with 8 cores we would run `./make.sh -j 8`.

Also, we have the option to wrap our code into python using `swig`. If you are not using GEMS you can skip this step.
There are two methods to do this.

1. Once the makefile is generated using `cmake`, you can go into the `cmakeBuild` directory (or wherever you threw the makefile) and use the `gmml_wrapped` make target.

2. You can just call, from the base `GMML` directory, `./make.sh -w` 

Now all one must do is run the make script.

```bash
$./make.sh
```

This will create the needed `cmake` files and will add the following directories within the `gmml2` directory:

* `lib` (Contains: the `gmml2` shared object libary, `libgmml2.so`)
* `cmakeBuild` (Contains: all files produced by the `cmake` command, a `compile_commands.json` file to be used with tools of your choice, and all files contained within the directories listed above)

You can either use the `libgmml2.so` file within the `lib` directory or the `libgmml2.so` file within the `cmakeBuild` directory. They are the exact same.

Both the `build` and `lib` directories must remain untouched because `gems` utilizes them both and expects both to be in the state that `./make.sh` leaves them.

Please enter `./make.sh -h` for help regarding the make script.

---
## Testing the Library

From within the `gmml2` directory, you must change your current working directory to the `gmml2/tests` directory. Note that `<NUM_JOBS>` is however many tests you want to run at once.

```bash
gmml2$ cd tests/
gmml2/tests$ ./compile_run_tests.bash -j <NUM_JOBS>
```

Please note that running GMML bare metal will cause test 016 (svg drawing) to fail, this is due to not setting the created svgs correctly and will eventually be fixed but for now don't worry if `016.test.DrawGlycan.sh` fails while running on bare metal; if you are utilizing the dev enviroment all tests are expected to pass. This is of no concern because these tests need some extra things running to check, but those are internal for now.

The output will tell you whether or not the library is behaving appropriately and if all tests are passed the output will look similar to the following:

```bash
$ bash compile_run_tests.bash -j4

#### Beginning GMML tests ####
Number of tests found:	13
Number of testing jobs:	4

mkdir: created directory './tempTestOutputs'

Beginning test: ./016.test.DrawGlycan.sh
Beginning test: ./017.test.GlycoproteinBuilder.sh
Beginning test: ./018.test.GlycoproteinBuilderTable.sh
Beginning test: ./019.test.newPDBClass.sh

Testing 016.test.DrawGlycan.cc...0.svg tests/correct_outputs/016.output_SVGs/0.svg differ: byte 15325, line 70
Test FAILED! Output file 0.svg different to tests/correct_outputs/016.output_SVGs/0.svg
Exit Code: 1

Beginning test: ./020.test.parameterFiles.sh

Testing 018.test.createGlycosylationTables.cpp... Test passed.
Exit Code: 0

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

Beginning test: ./029.test.graph.sh

Testing 029.graph...Test FAILED! Output file different. Try
diff 029.output_graph.txt tests/correct_outputs/029.output_graph.txt
Exit Code: 1

Beginning test: ./030.test.gmPreProcessor.sh

Testing 030.test.gmPreProcessor.cpp... ~3 seconds. Test passed.
Exit Code: 0

Testing 028.test.cdsCarbBuilderAll.cpp...Test passed.
Exit Code: 0

Testing 027.test.glycamResidueCombinator.cpp... Test passed.

Exit Code: 0

Testing 019.test.newPDBClass.cpp... ~30 seconds. Test passed.
Exit Code: 0

######## GMML TESTS COMPLETED ########
Required tests:	13
Passed tests:	11
Failed tests:	2
Time taken:	11 seconds
######################################

!!! OUTPUT OF THE 2 GMML TEST(S) THAT FAILED !!!

Testing 016.test.DrawGlycan.cc...0.svg tests/correct_outputs/016.output_SVGs/0.svg differ: byte 15325, line 70
Test FAILED! Output file 0.svg different to tests/correct_outputs/016.output_SVGs/0.svg
Exit Code: 1

Testing 029.graph...Test FAILED! Output file different. Try
diff 029.output_graph.txt tests/correct_outputs/029.output_graph.txt
Exit Code: 1

!!! FINISHED PRINTING FAILED TESTS !!!

```
Note that both test 016 and 029 fail outside of the developer environment and that's ok. If any other tests fail the something is wrong.


## Using the Glycoprotein Builder:
[Glycoprotein Builder Instructions](internalPrograms/GlycoproteinBuilder/README.md)


## Developers only (other users can ignore below here): 

###Updating file lists and whatnot

**DO NOT JUST FIRE THE `updateCmakeFileList.sh` SCRIPT AND NOT KNOW WHAT IS GOING ON. The method implemented is done in order to avoid a taxing typical cmake pattern; if the script is just fired off too many times we will have to remove it in order to avoid possible undefined behavior**. Please note that not only for cmake, but for all compilers, one should not just grab every file present and compile; these type of things must have some thought to them. The reason why one should never just glob files that one *thinks* are what one needs to compile is due to the huge increase in chances of introducing unknown behavior.

Basically treat this the same way as one treats using `git add --all` as bad practice due to priming the code base to have a bunch of random files (that should not be pushed) added to the repo; instead of being able to directly avoid `git add --all` and using `git add <YOUR_SPECIFIC_FILES>` instead, **YOU** must be the difference between that logic if you call the script check the git.

The `cmakeFileLists` directory contains the ouput from our `./updateCmakeFileList.sh` script. This script goes through and grabs all our files that we want to compile. There are 3 types:

* `cFileList.txt` - this contains all of our cpp files and where they be

* `hDirectoryList.txt` - this contains all of the directories that OUR source headers are. In the compiler this gets passed `-I` flag

---
## Coding Standards

In order to make deving on the library consistent, we must enforce coding standards. They will be added piecewise, including the appropriate tests (be them pre-commit/push hooks, ci/cd hooks, etc.) and will be outlined below.

### Git Branches that define our workflow:
release: Updated from stable where commits get tagged with the semver e.g. 2.1.0 for use in the test, dev and actual websites.
stable: Updated from dev. Hotfix from here. Before pushing: (spin up dev env, do tests. precommit, no skipping).
dev: Updated from feature branches created from here. Before pushing: (pre-commit & pre-push so gmml2 level tests need to pass).

### Branch Naming

All branch names must take the form of `<branchType>_<descriptiveName>`. Be sure that you have a good descriptive name. The branch types are as follows:

- `feature` - any new feature that is being created. This will be the most commonly used branch type used.
- `bugfix` - branches that are fixing bugs
- `hotfix` - This is to fix critical bugs, etc. that is blocking people and that need to be handled immediately
- `playground` - For when you are not creating an actual feature or something that will be integrated into the main branches. Branches with this name would be for playing around with ideas.
- `juggle` - this is for branches that are intermediates between others, think if you have some horrible merge conflicts or something and want to create a to help make your life a bunch easier, that is what this branch is for. For instance, you could make the `juggle_zingusMerge` branch that could be made. This one will be used rarely, but it is nice to have in our back pocket. For instance lets say two devs are working on branches that need to be merged together before everything is placed into `gmml-test`, that is when this naming will be used.

Some examples of good branch names are:

- `feature_basicGithubActions`
- `bugfix_addKDNOparams`
- `playground_llvmTooling`

### Pre-Commit Hooks

We run various pre-commit hooks to help ensure `gmml`'s commit history remains as clean and readable as possible.

#### Hooks for C-Files

All code must follow the format described in the `.clang-format` file, and the pre-commit hook will ensure the commited format is correct. The precommit hook will ensure all files you want to commit are correctly formatted. Any files that are not correctly formatted will be listed in the terminal you tried to commit from, if you are using something like `gitflow` or `gitkraken` check the logs. Many code editors, IDEs or text editors, have the ability to apply a specific format on save of the file, so save yourself headaches and set that up.

Now, how do you format a specific file?

```bash
user@host:.../gmml2$ clang-tidy-15 -i path/to/bad/file.cpp 
```

What if you did a bunch of files and want to be lazy? This can miss a couple bits that need to be changed so run it a couple times, it also will use all your cores but hey it is pretty quick.

```bash
user@host:.../gmml2$ find . -not -path "./cmakeBuild/*" -type f -iname "*.cpp" -o -iname "*.hpp" -o -iname "*.h" -o -iname "*.cc" | xargs -P $(nproc --all --ignore=2)  -I % sh -c 'clang-format-15 -i %'
```

#### Hooks for Shell Scripts

In order to commit any shell scripts, the files must adhear to both our formatting and linting commands. Do not worry, linting shell scripts is very quick. Our checks are defined directly below:

* Formatting a shell script (NOTE: the formatting will be immediately applied to the script in question):
```bash
user@host:.../gmml2$ shfmt -i 4 -ci -fn -w path/to/bad/script.sh
```

* Linting a shell script (NOTE: you will have to edit your script so it no longer has the issues the linter displays):
```bash
user@host:.../gmml2$ shellcheck --enable=require-variable-braces,quote-safe-variables,add-default-case path/to/bad/script.sh
```

---
## Depreciated Instructions

The GLYCAM Molecular Modeling Library (GMML) was designed to be used as a library accessed by GEMS (GLYCAM Extensible Modeling Script), but can be used as a standalone library.

More information about GEMS can be found here:

Website:  [http://glycam.org/gems](http://glycam.org/gems "GLYCAM GEMS")  
Github:  [https://github.com/GLYCAM-Web/gems](https://github.com/GLYCAM-Web/gems "GLYCAM GEMS - Github")

To get started, follow the [Download and Install](http://glycam.org/docs/gems/download-and-install/ "Download and Install") instructions. These instructions will walk you through the steps to obtain and configure the software, and also test the installation.
  
To compile and use the programs that are based on gmml2 (e.g. the carbohydrate or glycoprotein builders) go to their subfolders (e.g. internalPrograms/GlycoproteinBuilder/) and follow the compilation instructions in the readme there.
