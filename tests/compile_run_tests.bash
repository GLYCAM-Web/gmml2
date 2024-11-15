#!/bin/bash

################################################################
##########               Cool variables             ############
################################################################
#Tests that dont work on baremetal. Need to make this a list but lol
TEST_SKIP_LIST=("./016.test.DrawGlycan.sh")
SKIP_TIME=0

#Get base list of all files we want so we dont need to deal with
#figuring out file list more than once, we want to take advantage of output splitting
#here so we can hit the list ez pz
# shellcheck disable=SC2207=
readarray -t GMML_TEST_FILE_LIST < ./testList.txt
#This is mostly to make sure that we actually keep track of our failed tests.
GMML_FAILED_TESTS=0
#Just to keep track for our cool output at the end.
GMML_PASSED_TESTS=0
#amount of tests we are actually running, thus amount of tests we need to pass
#GMML_REQUIRED_PASSING=$(find . -type f -name "*.test.*.sh" | wc -l)

##### PRETTY PRINT STUFF #####
#lazy and dont want to have to type all these color variables a bunch
RESET_STYLE='\033[0m'
PASSED_STYLE='\033[0;32m\033[1m'
INFO_STYLE='\033[0;33m\033[1m'
ERROR_STYLE='\033[0;31m\033[1m'

################################################################
##########               Print help message         ############
################################################################

printHelp()
{
    echo -e "
===== GMML TEST RUNNING SCRIPT =====
$0 is used to allow us to run multiple of our
test scripts at once. This is not to be confused with compiling
the test code, we will be running multiple of our scripts at once.
There exists some GNU utils that could be used to accomplish this
but we want to keep these bash scripts as bare bones as possible.

We find tests to run by globbing \"*.test.*.sh\", to disable a
test just change the file name of test you would like to disable
so it does not match that glob.
*************************************************************
Options are as follows:
\t-d bare_metal\t\tSkip tests that dont work on bare metal and need
\t\t\t\tbe run in the docker container
\t-h \t\t\tPrint this msg
*************************************************************
Exiting"
    exit 1
}

################################################################
#########               COMMAND LINE INPUTS            #########
################################################################
#get opts allows for users to pass flags to the script, i.e. calling
#./compile_run_tests.bash -j 11 will allow us to run the script with 11 jobs.
#The colon after j (in the while decleration) means that said flag expects an "input"
#The lack of a colon after h means that said flag does not accept any "input"
while getopts "j:hd:" option; do
    case "${option}" in
        h)
            printHelp
            ;;
        #we need a wildcard case so that our script will print out the help msg
        #if a user passes an unexpected flag
        d)
            dIn="${OPTARG}"
            if [ "${dIn}" == "bare_metal" ]; then
                for TEST_SKIP in "${TEST_SKIP_LIST[@]}"; do
                    for ((SKIP_INDEX = ${#GMML_TEST_FILE_LIST[@]}; SKIP_INDEX >= 0; SKIP_INDEX--)); do
                        if [ "${GMML_TEST_FILE_LIST["${SKIP_INDEX}"]}" == "${TEST_SKIP}" ]; then
                            unset "GMML_TEST_FILE_LIST[${SKIP_INDEX}]"
                            GMML_TEST_FILE_LIST=("${GMML_TEST_FILE_LIST[@]}")
                            SKIP_TIME=1
                        fi
                    done
                done
            else
                printHelp
            fi
            ;;
        *)
            printHelp
            ;;
    esac
done
#$(....) runs the command in the subshell, strips trailing whitespaces, newlines, etc.
#and once all completes in the subshell it outputs the final info that would be sent to
#the stdout buffer and sets the variable to said value
START_TIME=$(date +%s)

echo -e "
${INFO_STYLE}#### Beginning GMML tests ####${RESET_STYLE}
Number of tests found: ${#GMML_TEST_FILE_LIST[@]}
"

mkdir -v ./tempTestOutputs
echo ""
#when this script wants to exit, trap will run the code before the exit code is actually returned
#this allows us to automatically delete the temp files no matter what happens. They are deleted if
#you use ctrl-c to exit, an error happens and the code exits, the script completes and the code
#exits, etc. etc.
trap "rm -rf ./tempTestOutputs" EXIT

#Gonna be used in the main loop, here so its close to first usage.
#These are set up as empty arrays
JOB_OUTPUT_FILES=()
FAILED_JOB_OUTPUTS=()

for CURRENT_TEST in "${GMML_TEST_FILE_LIST[@]}"; do
    #let user know whats going on
    echo -e "${INFO_STYLE}Beginning test:${RESET_STYLE} ${CURRENT_TEST}"

    #Wanted to do below but couldnt figure out force check that the only .sh
    #we would be replacing would be at the very end like I can with sed....
    #CURRENT_TEST_OUTPUT_FILENAME="${CURRENT_TEST//.sh/.output.txt}"

    #do some sed stuff to make our filename, this is a POSIX compliant method
    #Basically we first take the dir path and append the current_test name to the end of it after
    #replacing .sh with .output.txt, the sed -e means that we treat the next thing as an expression, "s" means we
    #are dealing with a substitute with sed, then / is a delimeter to the next part for the pattern we are matching,
    #we then have to escape the period with \. so it actually knows we are trying to match a period. Thus, what we are
    #matching is .sh, the $ as the end means that we are only matching .sh when it is at the END of the string, we then
    #have our delimeter / again in order to change to the bit where we define what the .sh is going to be replaced with,
    #finally the final / just means we are done with the pattern dealio
    CURRENT_TEST_OUTPUT_FILENAME="./tempTestOutputs/"$(echo "${CURRENT_TEST}" | sed -e "s/\.sh$/.output.txt/")

    #actually run the script we currently want to run, the &> redirects the stderr output and the stdout
    #output from our script to the file we define. Finally the & at the end runs the command in a subshell
    #that will not block this script. This means that instead of waiting for the command to complete, we just
    #spawn it in a subshell and continue running this script while the subshell does its thing.
    source "${CURRENT_TEST}" &>"${CURRENT_TEST_OUTPUT_FILENAME}" &
    
    #Add corresponding output file to our output file array
    JOB_OUTPUT_FILES+=("${CURRENT_TEST_OUTPUT_FILENAME}")
done

#wait for all jobs to complete....
wait

for JOB_OUTPUT in "${JOB_OUTPUT_FILES[@]}"; do
    #now we need to extract the exit code from the file itself, this will
    #be the last line of the file. I store as a variable cause we know we
    #gonna have to run the tail command at least once, and possibly more,
    #so just take the hit once and we are good. Gotta use 2 cause include newline
    DUMB_EXIT_CODE=$(tail -c -2 "${JOB_OUTPUT}")
    #Since we are checking arithmatic values, we need to use [[ .... ]] so bash doesnt
    #treat whats in the brackets as a string
    if [[ ${DUMB_EXIT_CODE} == 0 ]]; then
        #When we actually DO arithmatic we need to surround our code with ((....)) in order
        #to let bash know we are dealing with numbers and adding them
        ((GMML_PASSED_TESTS++))
        echo -e "${PASSED_STYLE}"
    elif [[ ${DUMB_EXIT_CODE} == 1 ]]; then
        ((GMML_FAILED_TESTS++))
        FAILED_JOB_OUTPUTS+=("${JOB_OUTPUT}")
        echo -e "${ERROR_STYLE}"
    fi
    cat "${JOB_OUTPUT}"
    echo -ne "${RESET_STYLE}"
done

case "${GMML_PASSED_TESTS}" in
    #if we passed as many tests as we have, we know we passed all tests thus we
    #go ahead and prepare the output to be green
    "${#GMML_TEST_FILE_LIST[@]}")
        RESULT_COLOR=${PASSED_STYLE}
        ;;
        #if our passed tests does not equal our number of tests, we know something is
        #wrong so we go ahead and do the spooky red output
    *)
        RESULT_COLOR=${ERROR_STYLE}
        ;;
esac

if [ "${SKIP_TIME}" == 1 ]; then
    echo -e "${INFO_STYLE}
######## WARNING YOU SKIPPED TESTS ########
Tests Skipped:\t${TEST_SKIP_LIST[*]}
######################################"
fi

#dont want to make another variable
#START_TIME=$(( $(date +%s) - START_TIME ))
echo -e "${RESULT_COLOR}
######## GMML TESTS COMPLETED ########
Required tests:\t${#GMML_TEST_FILE_LIST[@]}
Passed tests:\t${GMML_PASSED_TESTS}
Failed tests:\t${GMML_FAILED_TESTS}
Time taken:\t$(($(date +%s) - START_TIME)) seconds
######################################
"

#Ps paranoid programming
if [ "${GMML_PASSED_TESTS}" == "${#GMML_TEST_FILE_LIST[@]}" ]; then
    echo -e "ALL GMML TESTS PASSED${RESET_STYLE}\n"
    exit 0
elif [[ $((GMML_PASSED_TESTS + GMML_FAILED_TESTS)) != "${#GMML_TEST_FILE_LIST[@]}" ]]; then
    echo -e "!!! ERROR WE DIDNT GET EXIT CODES FROM ALL NEEDED SCRIPTS !!!${RESET_STYLE}\n"
    exit 1
elif [ "${GMML_FAILED_TESTS}" != 0 ]; then
    echo -e "!!! OUTPUT OF THE ${GMML_FAILED_TESTS} GMML TEST(S) THAT FAILED !!!\n"
    for BORKED_FILE in "${FAILED_JOB_OUTPUTS[@]}"; do
        cat "${BORKED_FILE}"
        echo ""
    done
    echo -e "!!! FINISHED PRINTING FAILED TESTS !!!${RESET_STYLE}\n"
    exit 1
else
    echo -e "!!! SOMETHING BORKED OH NO NOT GOOD !!!${RESET_STYLE}\n"
    exit 1
fi
