#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

# list output files to detect new or missing files
ls -1

#check md5sum for metrics json file
cat *data.json | md5sum
