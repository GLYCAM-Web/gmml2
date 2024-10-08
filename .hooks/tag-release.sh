#!/bin/bash

RESET_STYLE='\033[0m'
GREEN_BOLD='\033[0;32m\033[1m'
YELLOW_BOLD='\033[0;33m\033[1m'

GMML_ROOT_DIR=$(git rev-parse --show-toplevel)
read -r version < ${GMML_ROOT_DIR}/version.txt
branch=$(git rev-parse --abbrev-ref HEAD)

if [ "main" = "$branch" ]; then
    # validate-version.sh should have been run before this
    # assume that version.txt is up to date and we're good
    if git tag -a "v$version" -m "`date +%Y-%m-%d`"; then
        echo -e "Successfully created tag ${GREEN_BOLD}v$version${RESET_STYLE}"
    # tag couldn't be created
    else
        echo -e "${YELLOW_BOLD}Couldn't create tag v$version. Might already exist?${RESET_STYLE}"
    fi
fi
