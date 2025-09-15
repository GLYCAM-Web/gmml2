## Scripts and coding standards

If you want to contribute to `gmml2` you will also need to install the following packages:

* `clang-tidy-15`
* `shellcheck`
* `shfmt`

### Updating file lists and whatnot

**DO NOT JUST FIRE THE `updateCmakeFileList.sh` SCRIPT AND NOT KNOW WHAT IS GOING ON. The method implemented is done in order to avoid a taxing typical cmake pattern; if the script is just fired off too many times we will have to remove it in order to avoid possible undefined behavior**. Please note that not only for cmake, but for all compilers, one should not just grab every file present and compile; these type of things must have some thought to them. The reason why one should never just glob files that one *thinks* are what one needs to compile is due to the huge increase in chances of introducing unknown behavior.

Basically treat this the same way as one treats using `git add --all` as bad practice due to priming the code base to have a bunch of random files (that should not be pushed) added to the repo; instead of being able to directly avoid `git add --all` and using `git add <YOUR_SPECIFIC_FILES>` instead, **YOU** must be the difference between that logic if you call the script check the git.

The `cmakeFileLists` directory contains the ouput from our `./updateCmakeFileList.sh` script. This script goes through and grabs all our files that we want to compile. There are 3 types:

* `cFileList.txt` - this contains all of our cpp files and where they be

* `hDirectoryList.txt` - this contains all of the directories that OUR source headers are. In the compiler this gets passed `-I` flag

### Clang

GMML2 uses GCC per default, but will also compile with Clang if you wish to use associated tools. To do so, openmp will have to be installed for Clang aswell. Using apt, this can be done as follows
```bash
sudo apt-get update &&\
sudo apt-get install libomp-dev
```

---
## Coding Standards

In order to make deving on the library consistent, we must enforce coding standards. They will be added piecewise, including the appropriate tests (be them pre-commit/push hooks, ci/cd hooks, etc.) and will be outlined below.

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

We run various pre-commit hooks to help ensure `gmml2`'s commit history remains as clean and readable as possible.

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
