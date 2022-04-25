#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
#omitting outputData.txt since it contains time-based info
#fastqs are sorted before hashed

files=$(ls -I "*.txt")
for file in $files
do
    if [[ $file = *.fastq.gz ]]
    then
        md5sums+=$(zcat $file | sort | md5sum | awk '{print $1" "}')
    else
        md5sums+=$(md5sum $file | awk '{print $1" "}')
    fi
done

printf '%s\n' ${md5sums[*]} | sort