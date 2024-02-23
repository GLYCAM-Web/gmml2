#!/bin/bash

#flags for behav control
FORCE_RUN=0
KEEP_BUILD=0
KEEP_LOG=0
NMP=1

#just init em
VALGRIND_COMM=""
CACHEGRIND_COMM=""

#gets full path of all directories matching <NUM><NUM><NUM>.test.<testName>/
TEST_DIRS=$(find . -maxdepth 1 -type d -regex "\./[0-9][0-9][0-9].test.*" -exec realpath {} \; 2>/dev/null)
printHelp()
{
    echo -e "\t========= GMML TEST RUNNER SCRIPT =========
This script compiles and runs each test inside each directory
that follows the pattern of NNN.test.<dirName> where N is some
number.
*************************************************************
Options are as follows:
\t-f\t\t\tForce all tests to be rebuilt, deleting all previously
\t\t\t\t\tbuilt test code inside its respective directory
\t-j <NUM_JOBS>\t\tRun <NUM_JOBS> test at a time
\t-p\t\t\tPreserve specific files within each test
\t\tlogs\t\t\tPreserve build log, test output, etc.
\t\tbuild\t\t\tPreserve all build files, including the binary
\t-t\t\t\tRun a tool on each test, and save the output to a file
\t\tcache\t\t\tRun \"cool command\" on each test
\t\tval\t\t\tRun valgrind with sane-ish options on each test
\t-h\t\t\tPrint this help message and exit
*************************************************************
Exiting."
    exit 1
}

while getopts "j:fp:t:h" option; do
    case "${option}" in
        j)
            if [[ "${OPTARG}" =~ ^[1-9][0-9]*$ ]]; then
                NMP="${OPTARG}"
            else
                printHelp
            fi
            ;;
        f)
            FORCE_RUN=1
            ;;
        p)
            if [ "${OPTARG}" == "build" ]; then
                KEEP_BUILD=1
            elif [ "${OPTARG}" == "logs" ]; then
                KEEP_LOG=1
            else
                printHelp
            fi
            ;;
        t)
            if [ "${OPTARG}" == "cache" ]; then
                CACHEGRIND_COMM="cool command"
            elif [ "${OPTARG}" == "val" ]; then
                #may want to remove track origins, its slow
                VALGRIND_COMM="valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes"
            else
                printHelp
            fi
            #This assumes that cacherun or val run is 1. Was gonna xor but w/e
            if [ -n "${VALGRIND_COMM}" ] && [ -n "${CACHEGRIND_COMM}" ]; then
                echo "ERROR: Cannot run both valgrind and cachegrind at the same time. Exiting."
                exit 1
            fi
            ;;
        h)
            printHelp
            ;;
        *)
            echo -e "Incorrect flag passed"
            printHelp
            ;;
    esac
done

for DIR in ${TEST_DIRS}; do

    cd "${DIR}" || {
        echo "ERROR: Cannot cd into ${DIR}"
        exit 1
    }

    mkdir build || {
        if [[ ${FORCE_RUN} == 0 ]]; then
            #possibly switch to -rf
            echo "ERROR: build directory for ${DIR} already exists"
            exit 1

        fi
        rm -r "${DIR}/build"
        mkdir build
    }

    #If we can make the build directory and continue, lets actually care about making a log
    #file. We do this so we can have all logs in one file
    LOG_FILE="${DIR}/$(basename "$(pwd)").log"
    LOG_CMD=" &>> ${LOG_FILE}"

    if [ -f "${LOG_FILE}" ]; then
        if [[ ${FORCE_RUN} == 0 ]]; then
            echo "ERROR: ${LOG_FILE} Already exists"
            exit 1
        fi
        rm "${LOG_FILE}"
    fi

    #gotta do this so the redirect works :/
    COOL_COMMAND="cmake -S . -B ./build ${LOG_CMD}"

    eval "${COOL_COMMAND}" || {
        echo "ERROR: Cmake command for ${DIR} test failed"
        exit 1
    }
    cd build || {
        echo "ERROR: Cannot cd into ${DIR}/build"
        exit 1
    }

    COOL_COMMAND="make -j ${NMP} ${LOG_CMD}"

    eval "${COOL_COMMAND}" || {
        echo "ERROR: Make command for ${DIR} test failed"
        exit 1
    }
    #TODO: Figure out the executable naming stuff. Also what about different
    #       args....
    #Kinda funky looking, but this assumes either cachegrind or valgrind are empty. Check done in options parsing
    COOL_COMMAND="${CACHEGRIND_COMM}${VALGRIND_COMM} ./test ${LOG_CMD}"
    #TODO: Check it doesnt bork
    eval "${COOL_COMMAND}" ||
        {
            echo "ERROR: The test in ${DIR} failed!"
            #shorts to next loop iteration
            continue
        }
    echo "PASSED: $(basename ${DIR})"

    if [[ ${KEEP_LOG} == 0 ]]; then
        rm "${LOG_FILE}"
    fi

    if [[ ${KEEP_BUILD} == 0 ]]; then
        rm -r "${DIR}"/build
    fi
done
