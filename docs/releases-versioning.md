## Releases and versioning

The two important branches in gmml2 are `dev` and `main`. Most work is done either on `dev` or on feature branches. `main` is expected to be stable.

#### Versioning

`version.txt` contains the current library version, which follows a semantic versioning scheme of `major.minor.patch`. Upon compilation, CMake generates a `version.h` file which can be included in the code. This file is not committed.

Each release has a tag associated with it, pinned to a specific commit. This is useful for debugging or reproducing older results.


#### Releases

The general pattern when releasing new features is to

1. Create a release branch from `dev`
2. Increment the major or minor number in `version.txt` in a commit with no other changes
3. Merge the release branch into `main`. A tag named after the current version will be created by a git hook
4. Merge `main` back into `dev`
5. Push both branches. Include the flag `--follow-tags`, or the tag will remain local

#### Hotfixes

When fixing a bug on `main`, steps 1-2 should be substituted by creating a hotfix branch from `main`, where fixes are first committed, followed by incrementing the patch number in `version.txt`. Then resume normally with steps 3-5.
