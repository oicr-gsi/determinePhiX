version 1.0

workflow determinePhiX {
  input {
    String runDirectory
    Array[Int] lanes
    String basesMask
    String? outputFileNamePrefix
  }

  call generateFastqs {
    input:
      runDirectory = runDirectory,
      lanes = lanes,
      basesMask = basesMask
  }

  call getPhiXData {
    input:
      read1 = generateFastqs.read1,
      read2 = generateFastqs.read2
  }

  call formatData {
    input:
      data = getPhiXData.data
  }

  output {
    File matchedRead1 = getPhiXData.matchedRead1
    File matchedRead2 = getPhiXData.matchedRead2
    File dataFile = formatData.outputData
    File metricsJson = formatData.metrics
  }

  parameter_meta {
    runDirectory: {
      description: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE).",
      vidarr_type: "directory"
    }
    lanes: "A single lane or a list of lanes for no lane splitting (merging lanes)."
    basesMask: "The bases mask to produce the index reads (e.g. single 8bp index = \"Y1N*,I8,N*\", dual 8bp index = \"Y1N*,I8,I8,N*\")."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
  }

  meta {
    author: "Murto Hilali & Michael Laszloffy"
    email: "murto.hilali@oicr.on.ca OR michael.laszloffy@oicr.on.ca"
    description: "Workflow to determine PhiX contamination of undetermined reads"
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
      matchedRead1: "Read 1 with contaminated bases removed",
      matchedRead2: "Read 2 with contaminated bases removed",
      dataFile: "All output data contained in TXT file",
      metricsJson: "Contamination data contained in a JSON file"
    }
  }
}

task generateFastqs {
  input {
    String runDirectory
    Array[Int] lanes
    String basesMask
    String bcl2fastq = "bcl2fastq"
    String modules = "bcl2fastq/2.20.0.422"
    Int mem = 32
    Int timeout = 6
  }

  String outputDirectory = "out"

  command <<<
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
  >>>

  output {
    File read1 = "~{outputDirectory}/Undetermined_S0_R1_001.fastq.gz"
    File read2 = "~{outputDirectory}/Undetermined_S0_R2_001.fastq.gz"
  }

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  parameter_meta {
    runDirectory: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE)."
    lanes: "A single lane or a list of lanes for no lane splitting (merging lanes)."
    basesMask: "The bases mask to produce the index reads (e.g. single 8bp index = \"Y1N*,I8,N*\", dual 8bp index = \"Y1N*,I8,I8,N*\")."
    bcl2fastq: "bcl2fastq binary name or path to bcl2fastq."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      read1: "Read 1 fastq.gz.",
      read2: "Read 2 fastq.gz."
    }
  }
}

task getPhiXData {
  input {
    File read1
    File read2
    String modules = "bbmap/38.75"
    Int mem = 32
    Int timeout = 6
  }

  command <<<

$BBMAP_ROOT/share/bbmap/bbduk.sh \
threads=6 \
in1=~{read1} \
in2=~{read2} \
outm1=matched1.fastq.gz \
outm2=matched2.fastq.gz \
ref=$BBMAP_ROOT/share/bbmap/resources/phix174_ill.ref.fa.gz \
k=31 \
hdist=1

  >>>

  output {
    File data = stderr()
    File matchedRead1 = "matched1.fastq.gz"
    File matchedRead2 = "matched2.fastq.gz"
  }

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  parameter_meta {
    read1: "Undetermined read 1 from bcl data"
    read2: "Undetermined read 2 from bcl data"
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."

  }

  meta {
    output_meta: {
      data: "File containing stdout from PhiX determination step",
      matchedRead1: "Read 1 with contaminated bases removed",
      matchedRead2: "Read 2 with contaminated bases removed"
    }
  }

}

task formatData {
  input {
    File data
    String modules = ""
    Int mem = 32
    Int timeout = 6
  }

  command <<<

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

  >>>

  output {
    File metrics = "data.json"
    File outputData = "outputData.txt"
  }

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  parameter_meta {
    data: "File containing stdout from PhiX determination step"
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      metrics: "Contamination data contained in a JSON file",
      outputData: "All output data contained in TXT file"
    }
  }
}