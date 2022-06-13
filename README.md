# determinePhiX

Workflow to determine PhiX content of a sequencing run lane, by alignment of PhiX reference genome to undetermined reads and using the PhiX indices to pull PhiX reads, also runs fastQC on undetermined reads.

## Overview

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
`runDirectory`|String|{'description': 'Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE).', 'vidarr_type': 'directory'}
`lane`|String|A single lane to get metrics from.
`basesMask`|String|The bases mask to produce the index reads (e.g. single 8bp index = "Y1N*,I8,N*", dual 8bp index = "Y1N*,I8,I8,N*").
`outputFileNamePrefix`|String|Output prefix to prefix output file names with.


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`phiXindices`|Array[String]?|None|List of PhiX index or indices to generate PhiX fastqs as additional method to estimate PhiX content, eg. ["GGGGGGGG", "AGATCTCG"] (Optional).


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`generateFastqs.modules`|String|"bcl2fastq/2.20.0.422"|Environment module name and version to load (space separated) before command execution.
`generateFastqs.mem`|Int|32|Memory (in GB) to allocate to the job.
`generateFastqs.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`generateFastqs.threads`|Int|8|Requested CPU threads.
`fastQC.secondMateZip_timeout`|Int|1|Timeout, in hours, needed to override imposed limits.
`fastQC.secondMateZip_jobMemory`|Int|2|Memory allocated to this task.
`fastQC.secondMateHtml_timeout`|Int|1|Timeout, in hours, needed to override imposed limits.
`fastQC.secondMateHtml_jobMemory`|Int|2|Memory allocated to this task.
`fastQC.secondMateFastQC_modules`|String|"perl/5.28 java/11 fastqc/0.11.9"|Names and versions of required modules.
`fastQC.secondMateFastQC_threads`|Int?|None|Threads param for fastqc
`fastQC.secondMateFastQC_javaHeap`|Int|4|Memory allocated to java heap, in G.
`fastQC.secondMateFastQC_timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`fastQC.secondMateFastQC_jobMemory`|Int|6|Memory allocated to fastqc.
`fastQC.firstMateZip_timeout`|Int|1|Timeout, in hours, needed to override imposed limits.
`fastQC.firstMateZip_jobMemory`|Int|2|Memory allocated to this task.
`fastQC.firstMateHtml_timeout`|Int|1|Timeout, in hours, needed to override imposed limits.
`fastQC.firstMateHtml_jobMemory`|Int|2|Memory allocated to this task.
`fastQC.firstMateFastQC_modules`|String|"perl/5.28 java/11 fastqc/0.11.9"|Names and versions of required modules.
`fastQC.firstMateFastQC_threads`|Int?|None|Threads param for fastqc
`fastQC.firstMateFastQC_javaHeap`|Int|4|Memory allocated to java heap, in G.
`fastQC.firstMateFastQC_timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`fastQC.firstMateFastQC_jobMemory`|Int|6|Memory allocated to fastqc.
`fastQC.r1Suffix`|String|"_R1"|Suffix for R1 file.
`fastQC.r2Suffix`|String|"_R2"|Suffix for R2 file.
`getPhiXData.modules`|String|"bbmap/38.75"|Environment module name and version to load (space separated) before command execution.
`getPhiXData.mem`|Int|32|Memory (in GB) to allocate to the job.
`getPhiXData.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`getPhiXData.threads`|Int|6|Requested CPU threads.
`formatData.mem`|Int|12|Memory (in GB) to allocate to the job.
`formatData.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`generatePhixFastqs.modules`|String|"bcl2fastq/2.20.0.422"|Environment module name and version to load (space separated) before command execution.
`generatePhixFastqs.mem`|Int|32|Memory (in GB) to allocate to the job.
`generatePhixFastqs.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`generatePhixFastqs.threads`|Int|8|Requested CPU threads.
`formatPhiXdata.mem`|Int|12|Memory (in GB) to allocate to the job.
`formatPhiXdata.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.


### Outputs

Output | Type | Description
---|---|---
`metricsJson`|File|Collection of metrics for the sequencing run and PhiX content in JSON file.
`fastqcResultsR1`|File|FastQC results for undetermined reads 1, zipped.
`fastqcResultsR2`|File?|FastQC results for undetermined reads 2, zipped.


## Commands
 This section lists command(s) run by WORKFLOW workflow
 
 * Running WORKFLOW
 
 Workflow to estimate the PhiX content of a sequencing run, by performing an alignment of all the sequencing reads agains the PhiX genome. The workflow also runs fastQC on the all reads for QC and has an optional task to extract the PhiX reads using the PhiX indices.
 
 ### Obtain undetermined reads from sequencing run BCL data
 ```
     bcl2fastq \
     --runfolder-dir "~{runDirectory}" \
     --intensities-dir "~{runDirectory}/Data/Intensities/" \
     --processing-threads 8 \
     --output-dir "~{outputDirectory}" \
     --create-fastq-for-index-reads \
     --sample-sheet "/dev/null" \
     --tiles "s_[~{lane}]" \
     --use-bases-mask "~{basesMask}" \
     --no-lane-splitting \
     --interop-dir "~{outputDirectory}/Interop"
 
     #rename files to include run name
     mv ~{outputDirectory}/Undetermined_S0_R1_001.fastq.gz ~{outputDirectory}/~{outputFileNamePrefix}_Undetermined_S0_R1_001.fastq.gz
     mv ~{outputDirectory}/Undetermined_S0_R2_001.fastq.gz ~{outputDirectory}/~{outputFileNamePrefix}_Undetermined_S0_R2_001.fastq.gz
 ```
 
 ### Align undetermined reads to PhiX genome
 ```
   $BBMAP_ROOT/share/bbmap/bbmap.sh \
   threads=6 \
   in=~{fastqSingleRead} \
   out=~{outputFileNamePrefix}.sam \
   ref=$BBMAP_ROOT/share/bbmap/resources/phix174_ill.ref.fa.gz \
   nodisk -Xmx32g
 
   echo "read=~{read}" > ~{outputFileNamePrefix}_data.txt
   cat stderr >>  ~{outputFileNamePrefix}_data.txt
 ```
 
 ### Format output metrics into JSON
 ```
     python3<<CODE
     import re
     import json
 
     with open(r'~{data[0]}') as f1:
         lines_f1 = f1.readlines()
 
     with open(r'~{data[1]}') as f2:
         lines_f2 = f2.readlines()
 
     stats_file = open("~{stats}")
     data = json.load(stats_file)
 
     # Create dictionary for json file
     metricsJson = {}
     metricsJson["run_id"] = data["RunId"]
     metricsJson["lane_number"] = data["ConversionResults"][0]["LaneNumber"]
     metricsJson["total_clusters"] = data["ConversionResults"][0]["TotalClustersPF"]
     #print("TotalClustersRaw ", data["ConversionResults"][0]["TotalClustersRaw"])
 
     #get metrics for each read
     metricsJson["read_1"] = {}
     metricsJson["read_2"] = {}
 
     for file in [lines_f1, lines_f2]:
         d = {}
         for line in file:
 
             line_split = line.split("\t")
             if "Reads Used" in line:
                 d["total_reads"] = line_split[1]
                 d["total_bases"] = line_split[2].split(" ")[0][1:]
 
             if "mapped" in line:
                 d["pct_reads_align_phiX"] = round(float(line_split[1].strip()[:-1]),2)
                 d["number_reads_align_phiX"] = line_split[2].strip()
 
         if file[0].strip() == "read=1":
             metricsJson["read_1"] = d
         elif file[0].strip() == "read=2":
             metricsJson["read_2"] = d
 
     #Write dictionary out to JSON file
     out = open("~{outputFileNamePrefix}_data.json","w")
     json.dump(metricsJson, out)
     out.close()
 
     CODE
 ```
 
 ### Obtain PhiX reads from sequencing run BCL data
 - (OPTIONAL, only runs if phiXindices argument is provided)
 
 ```
     #create sample sheet for bcl2fastq with phiX indices
     echo [Header],,,,, > SampleSheet.csv
     echo Date,2022-05-26-04:00,,,, >> SampleSheet.csv
     echo ,,,,, >> SampleSheet.csv
     echo [Reads],,,,, >> SampleSheet.csv
     echo 75,,,,, >> SampleSheet.csv
     echo 75,,,, >> SampleSheet.csv
     echo ,,,,, >> SampleSheet.csv
     echo [Data],,,,, >> SampleSheet.csv
 
     if [ "~{indexTotal}" -eq 1 ]; then
       echo Sample_ID,Sample_Name,I7_Index_ID >> SampleSheet.csv
       echo PHIX,PHIX_0001,PhiX Adapter,"~{phiXindices[0]}" >> SampleSheet.csv
     fi
 
     if [ "~{indexTotal}" -eq 2 ]; then
       echo Sample_ID,Sample_Name,I7_Index_ID,index,I5_Index_ID,index2 >> SampleSheet.csv
       echo PHIX,PHIX_0001,PhiX Adapter,"~{phiXindices[0]}",PhiX Adapter,"~{phiXindices[1]}" >> SampleSheet.csv
     fi
 
     #run bcl2fastq using sample sheet created above
     bcl2fastq \
     --runfolder-dir "~{runDirectory}" \
     --processing-threads 8 \
     --output-dir "~{outputDirectory}" \
     --create-fastq-for-index-reads \
     --sample-sheet SampleSheet.csv \
     --tiles "s_[~{lane}]" \
     --use-bases-mask "~{basesMask}" \
     --no-lane-splitting \
     --interop-dir "~{outputDirectory}/Interop"
 
     #rename files to include run name
     mv ~{outputDirectory}/PHIX_0001_S1_R1_001.fastq.gz ~{outputDirectory}/~{outputFileNamePrefix}_PHIX_0001_S1_R1_001.fastq.gz
     mv ~{outputDirectory}/PHIX_0001_S1_R2_001.fastq.gz ~{outputDirectory}/~{outputFileNamePrefix}_PHIX_0001_S1_R2_001.fastq.gz
 ```
 
 ### Format PhiX reads metrics into JSON
 - (OPTIONAL, only runs if phiXindices argument is provided)
 ```
 
     python3<<CODE
     import re
     import json
 
     #stats from generatePhixFastqs
     stats_file = open("~{phixReadsStats}")
     phiXdata = json.load(stats_file)
 
     #metrics from formatData
     data_file = open("~{data}")
     metricsJson = json.load(data_file)
 
     #phix reads
     metricsJson["phixReads"] = {}
     metricsJson["phixReads"]["MismatchCounts"] = phiXdata["ConversionResults"][0]["DemuxResults"][0]['IndexMetrics'][0]["MismatchCounts"]
     metricsJson["phixReads"]["IndexSequence"] = phiXdata["ConversionResults"][0]["DemuxResults"][0]['IndexMetrics'][0]["IndexSequence"]
     metricsJson["phixReads"]["NumberofReads"] = phiXdata["ConversionResults"][0]["DemuxResults"][0]["NumberReads"]
 
     #Write dictionary out to JSON file
     out = open("~{outputFileNamePrefix}_data.json","w")
     json.dump(metricsJson, out)
     out.close()
 
     CODE
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
