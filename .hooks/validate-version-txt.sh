#!/bin/bash

RESET_STYLE='\033[0m'
GREEN_BOLD='\033[0;32m\033[1m'
YELLOW_BOLD='\033[0;33m\033[1m'
RED_BOLD='\033[0;31m\033[1m'

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
read -r version < ${GMML_ROOT_DIR}/version.txt

# version.txt has unstaged changes when git diff returns 1
$(git diff --diff-filter=MD --quiet -- ${GMML_ROOT_DIR}/version.txt)
# bash doesn't let me capture the exit code of git diff directly for some reason
# capturing it indirectly as below is the hack solution
unstaged_version_changes=$?

$(git diff --diff-filter=MD --quiet --cached HEAD -- ${GMML_ROOT_DIR}/version.txt)
staged_version_changes=$?
number_of_staged_files=$(git diff --cached --numstat | wc -l)

if [[ $unstaged_version_changes != 0 ]]; then
    echo -e "${RED_BOLD}ERROR: ${RESET_STYLE}changes to version.txt should not be left uncommitted, or the scripts will break"
    echo "Exiting"
    exit 1
elif [[ $staged_version_changes != 0 && $number_of_staged_files != 1 ]]; then
    echo -e "${RED_BOLD}ERROR: ${RESET_STYLE}please commit nothing but version.txt on version updates"
    echo "Exiting"
    exit 1
elif [[ ! $version =~ ^([0-9]+)\.([0-9]+)\.([0-9]+)(\-[A-Z]+\.[0-9]+)?$ ]]; then
    echo -e "${RED_BOLD}ERROR: ${RESET_STYLE}version.txt '$version' does not conform to semantic versioning format MAJOR.MINOR.PATCH"
    echo "Exiting"
    exit 1
fi
