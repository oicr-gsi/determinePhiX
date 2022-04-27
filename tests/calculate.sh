#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
#find . -xtype f ! -name "*.txt" -exec md5sum {} + | sort

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

#do
#    if [[ $file = *.fastq.gz ]]
#    then
#        test+=$(zcat $file | sort | md5sum | awk '{print $1" "}')
#    fi
#done

printf '%s\n' ${md5sums[*]} | sort




#fastqs=$(ls *.fastq.gz)
#or file in $fastqs
#do
#    md5sums+=$(zcat "$file" | sort | md5sum | awk '{print $1" "}')
#done
#md5sums+=$(md5sum ./*.json | awk '{print $1" "}')

#printf '%s\n' ${md5sums[*]} | sort
