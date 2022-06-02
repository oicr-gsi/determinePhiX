version 1.0

import "imports/fastqc.wdl" as fastqc

workflow determinePhiX {
  input {
    String runDirectory
    String lane
    String basesMask
    String outputFileNamePrefix
  }

  call generateFastqs {
    input:
      runDirectory = runDirectory,
      lane = lane,
      basesMask = basesMask,
      outputFileNamePrefix = outputFileNamePrefix
  }

  call fastqc.fastQC {
    input:
      fastqR1 = generateFastqs.read1,
      fastqR2 = generateFastqs.read2,
      basesMask = basesMask,
      outputFileNamePrefix = outputFileNamePrefix
  }

  call generatePhixFastqs {
    input:
      runDirectory = runDirectory,
      lane = lane,
      outputFileNamePrefix = outputFileNamePrefix
  }

  scatter (read in generateFastqs.reads) {
    call getPhiXData {
      input:
        fastqSingleRead = read.right,
        read = read.left,
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

  call formatData {
    input:
      data = getPhiXData.data,
      phixReadsStats = generatePhixFastqs.phixStats,
      outputFileNamePrefix = outputFileNamePrefix
  }

  output {
    File metricsJson = formatData.metrics
    File fastqcResultsR1 = fastQC.zip_bundle_R1
    File? fastqcResultsR2 = fastQC.zip_bundle_R2
  }

  parameter_meta {
    runDirectory: {
      description: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE).",
      vidarr_type: "directory"
    }
    lane: "A single lane to get metrics from."
    basesMask: "The bases mask to produce the index reads (e.g. single 8bp index = \"Y1N*,I8,N*\", dual 8bp index = \"Y1N*,I8,I8,N*\")."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
  }

  meta {
    author: "Murto Hilali & Michael Laszloffy, updated by Beatriz Lujan"
    email: "michael.laszloffy@oicr.on.ca, beatriz.lujantoro@oicr.on.ca"
    description: "Workflow to determine PhiX content of a sequencing run lane, by alignment of PhiX reference genome to undetermined reads and using the PhiX indices to pull PhiX reads, also runs fastQC on undetermined reads."
    dependencies: [
      {
        name: "bcl2fastq/2.20.0.422",
        url: "https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html"
      },
      {
        name: "htslib/1.9",
        url: "https://github.com/samtools/htslib"
      }
    ]
    output_meta: {
      metricsJson: "Collection of metrics for the sequencing run and PhiX content in JSON file.",
      fastqcResultsR1: "FastQC results for undetermined reads 1, zipped.",
      fastqcResultsR2: "FastQC results for undetermined reads 2, zipped."
    }
  }
}

task generateFastqs {
  input {
    String runDirectory
    String lane
    String basesMask
    String outputFileNamePrefix
    String modules = "bcl2fastq/2.20.0.422"
    Int mem = 32
    Int timeout = 6
    Int threads = 8
  }

  String outputDirectory = "out"

  command <<<
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
  >>>

  output {
    File read1 = "~{outputDirectory}/~{outputFileNamePrefix}_Undetermined_S0_R1_001.fastq.gz"
    File read2 = "~{outputDirectory}/~{outputFileNamePrefix}_Undetermined_S0_R2_001.fastq.gz"
    Array[Pair[String, File]] reads = [("1","~{outputDirectory}/~{outputFileNamePrefix}_Undetermined_S0_R1_001.fastq.gz"),("2","~{outputDirectory}/~{outputFileNamePrefix}_Undetermined_S0_R2_001.fastq.gz")]
  }

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
    threads: "{threads}"
  }

  parameter_meta {
    runDirectory: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE)."
    lane: "A single lane to produce fastqs from."
    basesMask: "The bases mask to produce the index reads (e.g. single 8bp index = \"Y1N*,I8,N*\", dual 8bp index = \"Y1N*,I8,I8,N*\")."
    outputFileNamePrefix: "Prefix to name output files."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      reads: "An array of fastq files with the read number associated with it.",
      read1: "Undetermined fastq reads 1.",
      read2: "Undetermined fastq reads 2."
    }
  }
}

task generatePhixFastqs {
  input {
    String runDirectory
    String lane
    String basesMask
    String outputFileNamePrefix
    String modules = "bcl2fastq/2.20.0.422"
    Int mem = 32
    Int timeout = 6
    Int threads = 8
  }

  String outputDirectory = "out"

  command <<<
    #create sample sheet for bcl2fastq with phiX indices
    echo [Header],,,,, > SampleSheet.csv
    echo Date,2022-05-26-04:00,,,, >> SampleSheet.csv
    echo ,,,,, >> SampleSheet.csv
    echo [Reads],,,,, >> SampleSheet.csv
    echo 75,,,,, >> SampleSheet.csv
    echo ,,,,, >> SampleSheet.csv
    echo [Data],,,,, >> SampleSheet.csv
    echo Sample_ID,Sample_Name,I7_Index_ID,index,I5_Index_ID,index2 >> SampleSheet.csv
    echo PHIX,PHIX_0001,PhiX Adapter,GGGGGGGG,PhiX Adapter,AGATCTCG >> SampleSheet.csv

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
  >>>

  output {
    File phixStats = "~{outputDirectory}/Stats/Stats.json"
  }

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
    threads: "{threads}"
  }

  parameter_meta {
    runDirectory: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE)."
    lane: "A single lane to get fastqs from."
    outputFileNamePrefix: "Prefix to name output files."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      phixStats: "bcl2fastq stats for PhiX reads."
    }
  }
}

task getPhiXData {
  input {
    File fastqSingleRead
    String read
    String outputFileNamePrefix
    String modules = "bbmap/38.75"
    Int mem = 32
    Int timeout = 6
    Int threads = 6
  }

  command <<<

  $BBMAP_ROOT/share/bbmap/bbmap.sh \
  threads=6 \
  in=~{fastqSingleRead} \
  out=~{outputFileNamePrefix}.sam \
  ref=$BBMAP_ROOT/share/bbmap/resources/phix174_ill.ref.fa.gz \
  nodisk -Xmx32g

  echo "read=~{read}" > ~{outputFileNamePrefix}_data.txt
  cat stderr >>  ~{outputFileNamePrefix}_data.txt

  >>>

  output {
    File data = "~{outputFileNamePrefix}_data.txt"
  }

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
    cpu: "~{threads}"
  }

  parameter_meta {
    fastqSingleRead: "Undetermined read from bcl2fastq"
    read: "Indicates the read number passed to bbmap, either 1 or 2."
    outputFileNamePrefix: "Prefix to name output files."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    threads: "Requested CPU threads"

  }

  meta {
    output_meta: {
      data: "File containing ouput from alignment."
    }
  }

}

task formatData {
  input {
    Array[File] data
    String outputFileNamePrefix
    File phixReadsStats
    Int mem = 32
    Int timeout = 6
  }

  command <<<

    python3<<CODE
    import re
    import json

    with open(r'~{data[0]}') as f1:
        lines_f1 = f1.readlines()

    with open(r'~{data[1]}') as f2:
        lines_f2 = f2.readlines()

    stats_file = open("~{phixReadsStats}")
    data = json.load(stats_file)

    # Create dictionary for json file
    metricsJson = {}
    metricsJson["run_id"] = data["RunId"]
    metricsJson["lane_number"] = data["ConversionResults"][0]["LaneNumber"]
    metricsJson["total_clusters"] = data["ConversionResults"][0]["TotalClustersPF"]
    #print("TotalClustersRaw ", data["ConversionResults"][0]["TotalClustersRaw"])

    #phix reads
    metricsJson["phix_reads"] = {}
    metricsJson["phix_reads"]["index_sequence"] = data["ConversionResults"][0]["DemuxResults"][0]['IndexMetrics'][0]["IndexSequence"]
    metricsJson["phix_reads"]["number_of_reads"] = data["ConversionResults"][0]["DemuxResults"][0]["NumberReads"]
    metricsJson["phix_reads"]["mismatch_counts"] = data["ConversionResults"][0]["DemuxResults"][0]['IndexMetrics'][0]["MismatchCounts"]

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
  >>>

  output {
    File metrics = "~{outputFileNamePrefix}_data.json"
  }

  runtime {
    memory: "~{mem} GB"
    timeout: "~{timeout}"
  }

  parameter_meta {
    data: "File containing output metrics from PhiX alignment."
    phixReadsStats: "bcl2fastq stats for PhiX reads.."
    outputFileNamePrefix: "Prefix to name output files."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      metrics: "Contamination data contained in a JSON file."
    }
  }
}
