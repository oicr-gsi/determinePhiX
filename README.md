# determinePhiX

Workflow to determine PhiX contamination of undetermined reads

## Overview

![](./determinePhiX.svg?raw=true "Workflow diagram")

## Dependencies

* [bcl2fastq 2.20.0.422](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
* [htslib 1.9](https://github.com/samtools/htslib)


## Usage

### Cromwell
```
java -jar cromwell.jar run determinePhiX.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`runDirectory`|String|Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE).
`lanes`|Array[Int]|A single lane or a list of lanes for no lane splitting (merging lanes).
`basesMask`|String|The bases mask to produce the index reads (e.g. single 8bp index = "Y1N*,I8,N*", dual 8bp index = "Y1N*,I8,I8,N*").


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String?|None|Output prefix to prefix output file names with.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`generateFastqs.bcl2fastq`|String|"bcl2fastq"|bcl2fastq binary name or path to bcl2fastq.
`generateFastqs.modules`|String|"bcl2fastq/2.20.0.422"|Environment module name and version to load (space separated) before command execution.
`generateFastqs.mem`|Int|32|Memory (in GB) to allocate to the job.
`generateFastqs.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`getPhiXData.modules`|String|"bbmap/38.75"|Environment module name and version to load (space separated) before command execution.
`getPhiXData.mem`|Int|32|Memory (in GB) to allocate to the job.
`getPhiXData.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`formatData.modules`|String|""|Environment module name and version to load (space separated) before command execution.
`formatData.mem`|Int|32|Memory (in GB) to allocate to the job.
`formatData.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.


### Outputs

Output | Type | Description
---|---|---
`matchedRead1`|File|Read 1 with contaminated bases removed
`matchedRead2`|File|Read 2 with contaminated bases removed
`dataFile`|File|All output data contained in TXT file
`metricsJson`|File|Contamination data contained in a JSON file


## Commands
 This section lists command(s) run by determinePhiX workflow
 
 * Running WORKFLOW
 
 Workflow to determine PhiX contamination of undetermined reads.
 
 ### Gathering undetermined reads from BCL data
 ```
     module load ~{modules}
     ~{bcl2fastq} \
     --runfolder-dir "~{runDirectory}" \
     --intensities-dir "~{runDirectory}/Data/Intensities/" \
     --processing-threads 8 \
     --output-dir "~{outputDirectory}" \
     --create-fastq-for-index-reads \
     --sample-sheet "/dev/null" \
     --tiles "s_[~{sep='' lanes}]" \
     --use-bases-mask "~{basesMask}" \
     --no-lane-splitting \
     --interop-dir "~{outputDirectory}/Interop"
 ```
 ### Obtaining decontaminated reads and metrics
 ```
 $BBMAP_ROOT/share/bbmap/bbduk.sh \
 threads=6 \
 in1=~{read1} \
 in2=~{read2} \
 outm1=matched1.fastq.gz \
 outm2=matched2.fastq.gz \
 ref=$BBMAP_ROOT/share/bbmap/resources/phix174_ill.ref.fa.gz \
 k=31 \
 hdist=1
 
 ```
 ### Formatting output data into TXT and JSON format
 ```
 
 cp ~{data} ./outputData.txt
 
 python3<<CODE
 import re
 import json
 
 with open(r'~{data}') as f:
     lines = f.readlines()
 for i in lines:
     if "Contaminants" in i:
         data = i
 
 # Create dictionary with contamination data
 metricsJson = {}
 
 # Adding read data
 reads = re.compile("(\d+)(?=\s*reads)")
 metricsJson["reads"] = reads.findall(data)[0]
 
 # Adding bases data
 bases = re.compile("(\d+)(?=\s*bases)")
 metricsJson["bases"] = bases.findall(data)[0]
 
 # Adding contamination data (formatted to 2 decimal places)
 contamination = re.compile("(\d\.\d+)(?=\s*%)")
 metricsJson["phiX_contamination"] = float("{:.4}".format(contamination.findall(data)[0]))
 
 # Write dictionary out to JSON file
 out = open("data.json","w")
 json.dump(metricsJson, out)
 out.close()
 
 CODE
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
