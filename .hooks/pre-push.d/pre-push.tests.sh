#!/bin/bash

##### PRETTY PRINT STUFF #####
#lazy and dont want to have to type all these color variables a bunch, also improves readability
RESET_STYLE='\033[0m'
GREEN_BOLD='\033[0;32m\033[1m'
YELLOW_BOLD='\033[0;33m\033[1m'
RED_BOLD='\033[0;31m\033[1m'

#Just a global var to check if GEMS or GMML is out of date. Here we are using string existence
#instead of a wannabe bool type flag
BRANCH_IS_BEHIND=""

#this is to enforce our branch naming scheme
branch_regex="^(feature|bugfix|hotfix|playground|juggle)_[a-zA-Z0-9]{2,36}$"

#How many commits a feature branch can be missing from the gmml-test branch
MAX_FEATURE_BEHIND_TEST=15

if [[ $(pwd) == */gems/gmml2 ]]; then
    GMML_DIR=$(pwd)
    #checkin to be sure
    if [ "${GMML_DIR}" != "$(git rev-parse --show-toplevel)" ]; then
        echo -e "Prepush hook failed, we think our GMML root directory is:\t${GMML_DIR}\n"
        exit 1
    fi
    GEMS_DIR=$(cd .. && pwd)
    if [ "${GEMS_DIR}" != "$(cd .. && git rev-parse --show-toplevel)" ]; then
        echo -e "Prepush hook failed, we think our GEMS root directory is:\t${GEMS_DIR}\n"
        exit 1
    fi
else
    echo -e "${RED_BOLD}ERROR: Trying to run from incorrect directory being: $(pwd)\nRun from the base GMML directory.\nABORTING${RESET_STYLE}"
    exit 1
fi

if [ -z "${GEMSHOME}" ]; then
    echo -e "${YELLOW_BOLD}WARNING: Your GEMSHOME environment variable is not set! It should be set to the GEMS directory\nthat is the parent of the GMML directory. This can cause some issues as some of the codebase\nstill relies upon the GEMSHOME variable. Continuing but you have been warned.${RESET_STYLE}"
fi

## OG Oct 2021 have the hooks update themselves.
#TODO: Do this more auto like, if this script is updated the first run the next time will not
#reflect the made changes due to the old script calling the copy then continuing.
cp -r "${GMML_DIR}"/.hooks/* "${GMML_DIR}"/.git/hooks/


TEST_SKIP=0
#### Allow skipping tests ####
branch=$(git rev-parse --abbrev-ref HEAD)
if [[ "${branch}" != "dev" ]] && [[ "${branch}" != "main" ]]; then

    echo -e "Branch is ${branch}\nSkipping tests is allowed.\nDo you want to skip them?\ns=skip\na=abort\nEnter anything to run tests.\n"
    read -p "Enter response: " response </dev/tty
    if [[ "${response}" == [sS] ]]; then
        echo -e "Skipping tests!\n"
        TEST_SKIP=1
    elif [[ "${response}" == [aA] ]]; then
        printf "Abort!\n"
        exit 1
    else
        printf "Running tests.\n"
    fi
fi

if [ "${TEST_SKIP}" == 1 ]; then
    echo "Skipping tests, you are good to push."
    exit 0
else
    echo "Beginning tests"
fi

cd "${GEMS_DIR}" || {
    echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}"
    echo "Exiting..."
    exit 1
}

#Add these removes so the tests don't pass on an old version of the library
#This below is commented out. I need to ensure that GEMs doesnt do anything funky with the generated files like move
#them within the gmml directory. End goal for gmml2 is that all generated files should remain where they are
#generated and ref to correct location. Maybe just delete the cmakeBuild dir unsure...
#rm -f "${GMML_DIR}/gmml.py" "${GMML_DIR}/_gmml.so"
rm -rf "${GMML_DIR:?}/lib"
if [ -d "${GMML_DIR}/cmakeBuild" ]; then
    echo "Removing the libgmml.so from our cmakeBuild directory"
    rm "${GMML_DIR}/cmakeBuild/libgmml2.so"
    rm "${GMML_DIR}/cmakeBuild/_gmml2.so"
    rm "${GMML_DIR}/cmakeBuild/gmml2.py"
fi

echo "Compiling gmml2 using GEMS ./make.sh, no wrap flag cause it auto wraps"

./make.sh -j "$(nproc --all --ignore=2)"

echo "Running mandatory GMML tests..."
cd "${GMML_DIR}"/tests/ || {
    echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GMML_DIR}/tests/"
    echo "Exiting..."
    exit 1
}

#assuming > 2 cores
nice -10 ./compile_run_tests.bash -j "$(nproc --all --ignore=2)"
result=$? # record the exit status from compile_run_tests.bash

cd "${GMML_DIR}" || {
    echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GMML_DIR}"
    echo "Exiting..."
    exit 1
}

if [ "${result}" -eq 0 ]; then
    echo "GMML level tests have passed. Doing gems level tests."
    cd "${GEMS_DIR}"/tests/ || {
        echo -e "${RED_BOLD}failed...${RESET_STYLE} We could not change directory to the following:\n\t ${GEMS_DIR}/tests"
        echo "Exiting..."
        exit 1
    }
    bash run_tests.sh
    gems_tests_result=$? # record the exit status of previous command
    if [ "${gems_tests_result}" -ne 0 ]; then
        echo "GEMS level tests have failed. Make sure you have pulled the latest version and are on the appropriate branch. "
        echo "If you are up-to-date, this failure indicates that you have caused the outputs of ${GEMS_DIR}/tests to change. You can open the ${GEMS_DIR}/tests/run_tests.sh file and run the test line by line to get an output file. Compare it to the saved \"correct\" version in ${GEMS_DIR}/tests/correct_outputs."
        echo "Sometimes the changes you make are fine, and you just need to update what the correct output is by overwriting the old output. Make sure it is ok though, or you will be mur-didely-urdered."
        exit 1
    else
        echo -e "${GREEN_BOLD}All tests have passed. Pushing allowed.${RESET_STYLE}"
        exit 0
    fi
else
    echo -e "${RED_BOLD}
         *****************************************************************
         The GMML level tests have failed! 
         Push cancelled.
         *****************************************************************
         ${RESET_STYLE}"
    exit 1
fi
