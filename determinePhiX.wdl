version 1.0

import "imports/fastqc.wdl" as fastQC

workflow determinePhiX {
  input {
    String runDirectory
    Array[Int] lanes
    String basesMask
    String outputFileNamePrefix
    Array[String]? phiXindices
  }

  call generateFastqs {
    input:
      runDirectory = runDirectory,
      lanes = lanes,
      basesMask = basesMask,
      outputFileNamePrefix = outputFileNamePrefix
  }

  call fastQC.fastQC {
    input:
      fastqR1 = generateFastqs.read1,
      fastqR2 = generateFastqs.read2,
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
      stats = generateFastqs.stats,
      outputFileNamePrefix = outputFileNamePrefix
  }

  if(defined(phiXindices)){
    Array[String] indices = select_first([phiXindices,[]])
    call generatePhixFastqs {
      input:
        runDirectory = runDirectory,
        lanes = lanes,
        basesMask = basesMask,
        phiXindices = indices,
        outputFileNamePrefix = outputFileNamePrefix
    }

    call formatPhiXdata {
      input:
        phixReadsStats = generatePhixFastqs.phixStats,
        data = formatData.metrics,
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

  output {
    File metricsJson = select_first([formatPhiXdata.metrics,formatData.metrics])
    File fastqcResultsR1 = fastQC.zip_bundle_R1
    File? fastqcResultsR2 = fastQC.zip_bundle_R2
  }

  parameter_meta {
    runDirectory: {
      description: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE).",
      vidarr_type: "directory"
    }
    lanes: "A single lane or a list of lanes for no lane splitting (merging lanes)."
    basesMask: "The bases mask to produce the index reads (e.g. single 8bp index = \"Y1N*,I8,N*\", dual 8bp index = \"Y1N*,I8,I8,N*\")."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    phiXindices: "List of PhiX index or indices to generate PhiX fastqs as additional method to estimate PhiX content, eg. [\"GGGGGGGG\", \"AGATCTCG\"] (Optional)."
  }

  meta {
    author: "Murto Hilali & Michael Laszloffy, updated by Beatriz Lujan"
    email: "michael.laszloffy@oicr.on.ca, beatriz.lujantoro@oicr.on.ca"
    description: "Workflow to determine PhiX content of a sequencing run lane, by alignment of PhiX reference genome to undetermined reads and using the PhiX indices to pull PhiX reads, also runs fastQC on undetermined reads. The PhiX Control Library (commonly referred to as PhiX, FC-110-3001) is derived from the small, well-characterized bacteriophage genome, PhiX. It is a concentrated Illumina library (10 nM in 10 µl) that has an average size of 500 bp and consists of balanced base composition at ~45% GC and ~55% AT. The PhiX library provides a quality control for cluster generation, sequencing, and alignment, and a calibration control for cross-talk matrix generation, phasing, and prephasing. It can be rapidly aligned to estimate relevant sequencing by synthesis (SBS) metrics such as phasing and error rate."
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
    metricsJson: {
        description: "Collection of metrics for the sequencing run and PhiX content in JSON file.",
        vidarr_label: "metricsJson"
    },
    fastqcResultsR1: {
        description: "FastQC results for undetermined reads 1, zipped.",
        vidarr_label: "fastqcResultsR1"
    },
    fastqcResultsR2: {
        description: "FastQC results for undetermined reads 2, zipped.",
        vidarr_label: "fastqcResultsR2"
    }
}
  }
}

task generateFastqs {
  input {
    String runDirectory
    Array[Int] lanes
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
    --tiles "s_[~{sep='' lanes}]" \
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
    File stats = "~{outputDirectory}/Stats/Stats.json"
  }

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
    threads: "{threads}"
  }

  parameter_meta {
    runDirectory: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE)."
    lanes: "A single lane or a list of lanes for no lane splitting (merging lanes)."
    basesMask: "The bases mask to produce the index reads (e.g. single 8bp index = \"Y1N*,I8,N*\", dual 8bp index = \"Y1N*,I8,I8,N*\")."
    outputFileNamePrefix: "Prefix to name output files."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    threads: "Requested CPU threads."
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
    Array[Int] lanes
    String basesMask
    String outputFileNamePrefix
    Array[String] phiXindices
    String modules = "bcl2fastq/2.20.0.422"
    Int mem = 32
    Int timeout = 6
    Int threads = 8
  }

  String outputDirectory = "out"
  Int indexTotal = length(phiXindices)

  command <<<
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
    --tiles "s_[~{sep='' lanes}]" \
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
    lanes: "A single lane or a list of lanes for no lane splitting (merging lanes)."
    outputFileNamePrefix: "Prefix to name output files."
    phiXindices: "List of PhiX index or indices to generate PhiX fastqs as additional method to estimate PhiX content."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    threads: "Requested CPU threads."
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
    threads: "Requested CPU threads."

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
    File stats
    String outputFileNamePrefix
    Int mem = 12
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
    outputFileNamePrefix: "Prefix to name output files."
    stats: "Sequencing stats from the generateFastqs task."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      metrics: "Contamination data contained in a JSON file."
    }
  }
}

task formatPhiXdata {
  input {
    File data
    String outputFileNamePrefix
    File phixReadsStats
    Int mem = 12
    Int timeout = 6
  }

  command <<<

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
  >>>

  output {
    File metrics = "~{outputFileNamePrefix}_data.json"
  }

  runtime {
    memory: "~{mem} GB"
    timeout: "~{timeout}"
  }

  parameter_meta {
    data: "File containing sequencing and alginment metrics from previous tasks."
    phixReadsStats: "bcl2fastq stats for PhiX reads."
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
