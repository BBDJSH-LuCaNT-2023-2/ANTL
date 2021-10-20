#!/bin/bash
# this script changes the current shell. To run it use the command 'source' followed by this file's name
# The purpose of this script is to customize the environment (currently a Docker container running a linux base image) where a CI job is being run

# use 256 color codes. See https://robotmoon.com/256-colors/ 
colors="error=01;38;5;200" # magenta
colors+=":warning=01;38;5;11" # yellow
export GCC_COLORS=$colors
