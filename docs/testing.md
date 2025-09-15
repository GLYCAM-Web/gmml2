## Testing
Most of the tests are functional tests, verifying that entire programs or features work as expected. It’s a good idea to add a new test when an optional feature is added.
Some tests code is compiled by the test scripts, others are run from compiled programs.

### Running tests

From the test directory, tests can be run through a script
```
gmml2/tests$ ./compile_run_tests.bash -j <NUM_JOBS>
```

Each test is run as a separate job by a multithreaded job scheduler, and the tests themselves are split among the input files `testsBase.txt`, `testsDrawGlycan.txt`, `testsGpBuilder.txt`. The script runs them all by default, but if we only want to run a subset of them we can do so, e.g
```
gmml2/tests$ ./compile_run_tests.bash -i testsBase.txt -i testsDrawGlycan.txt
```
How far the tests should be split up is a tradeoff between multithreading speedups and readability. Note that some programs, like the carbohydrate builder, spend a non-trivial amount of time loading parameter files compared to the time spent constructing small and medium-sized carbohydrates. Splitting this from one job into too many would risk slowing the tests down overall.

### test mode
The output of certain programs contain the gmml2 version number and / or absolute file paths.
In order to keep test output consistent across multiple versions and installations, the (hidden) flag `--test-mode` is used in these cases. When provided, version numbers are omitted, and file paths are written as relative.
Note that file paths do need to be absolute to be useful for other programs, but that is not the concern of the tests.
