version 1.0

# ======================================================
# Workflow accepts two fastq files, with R1 and R2 reads
# ======================================================
workflow fastQC {
input {
        Int secondMateZip_timeout = 1
        Int secondMateZip_jobMemory = 2
        Int secondMateHtml_timeout = 1
        Int secondMateHtml_jobMemory = 2
        String secondMateFastQC_modules = "perl/5.28 java/11 fastqc/0.11.9"
        Int? secondMateFastQC_threads
        Int secondMateFastQC_javaHeap = 4
        Int secondMateFastQC_timeout = 20
        Int secondMateFastQC_jobMemory = 6
        Int firstMateZip_timeout = 1
        Int firstMateZip_jobMemory = 2
        Int firstMateHtml_timeout = 1
        Int firstMateHtml_jobMemory = 2
        String firstMateFastQC_modules = "perl/5.28 java/11 fastqc/0.11.9"
        Int? firstMateFastQC_threads
        Int firstMateFastQC_javaHeap = 4
        Int firstMateFastQC_timeout = 20
        Int firstMateFastQC_jobMemory = 6
        File fastqR1 
        File? fastqR2
        String outputFileNamePrefix = ""
        String r1Suffix = "_R1"
        String r2Suffix = "_R2"
}
Array[File] inputFastqs = select_all([fastqR1,fastqR2])
String outputPrefixOne = if outputFileNamePrefix == "" then basename(inputFastqs[0], '.fastq.gz') + "_fastqc"
                                                       else outputFileNamePrefix + r1Suffix

call runFastQC as firstMateFastQC { input: inputFastq = inputFastqs[0], jobMemory = firstMateFastQC_jobMemory, timeout = firstMateFastQC_timeout, javaHeap = firstMateFastQC_javaHeap, threads = firstMateFastQC_threads, modules = firstMateFastQC_modules }
call renameOutput as firstMateHtml { input: inputFile = firstMateFastQC.html_report_file, extension = "html", customPrefix = outputPrefixOne, jobMemory = firstMateHtml_jobMemory, timeout = firstMateHtml_timeout }
call renameOutput as firstMateZip { input: inputFile = firstMateFastQC.zip_bundle_file, extension = "zip", customPrefix = outputPrefixOne, jobMemory = firstMateZip_jobMemory, timeout = firstMateZip_timeout }

if (length(inputFastqs) > 1) {
 String outputPrefixTwo = if outputFileNamePrefix=="" then basename(inputFastqs[1], '.fastq.gz') + "_fastqc"
                                                      else outputFileNamePrefix + r2Suffix
 call runFastQC as secondMateFastQC { input: inputFastq = inputFastqs[1], jobMemory = secondMateFastQC_jobMemory, timeout = secondMateFastQC_timeout, javaHeap = secondMateFastQC_javaHeap, threads = secondMateFastQC_threads, modules = secondMateFastQC_modules }
 call renameOutput as secondMateHtml { input: inputFile = secondMateFastQC.html_report_file, extension = "html", customPrefix = outputPrefixTwo, jobMemory = secondMateHtml_jobMemory, timeout = secondMateHtml_timeout }
 call renameOutput as secondMateZip { input: inputFile = secondMateFastQC.zip_bundle_file, extension = "zip", customPrefix = outputPrefixTwo, jobMemory = secondMateZip_jobMemory, timeout = secondMateZip_timeout }
}

parameter_meta {
    secondMateZip_timeout: "Timeout, in hours, needed to override imposed limits."
    secondMateZip_jobMemory: "Memory allocated to this task."
    secondMateHtml_timeout: "Timeout, in hours, needed to override imposed limits."
    secondMateHtml_jobMemory: "Memory allocated to this task."
    secondMateFastQC_modules: "Names and versions of required modules."
    secondMateFastQC_threads: "Threads param for fastqc"
    secondMateFastQC_javaHeap: "Memory allocated to java heap, in G."
    secondMateFastQC_timeout: "Timeout in hours, needed to override imposed limits."
    secondMateFastQC_jobMemory: "Memory allocated to fastqc."
    firstMateZip_timeout: "Timeout, in hours, needed to override imposed limits."
    firstMateZip_jobMemory: "Memory allocated to this task."
    firstMateHtml_timeout: "Timeout, in hours, needed to override imposed limits."
    firstMateHtml_jobMemory: "Memory allocated to this task."
    firstMateFastQC_modules: "Names and versions of required modules."
    firstMateFastQC_threads: "Threads param for fastqc"
    firstMateFastQC_javaHeap: "Memory allocated to java heap, in G."
    firstMateFastQC_timeout: "Timeout in hours, needed to override imposed limits."
    firstMateFastQC_jobMemory: "Memory allocated to fastqc."
  fastqR1: "Input file with the first mate reads."
  fastqR2: " Input file with the second mate reads (if not set the experiments will be regarded as single-end)."
  outputFileNamePrefix: "Output prefix, customizable. Default is the first file's basename."
  r1Suffix: "Suffix for R1 file."
  r2Suffix: "Suffix for R2 file."
}

meta {
    author: "Peter Ruzanov"
    email: "peter.ruzanov@oicr.on.ca"
    description: "Niassa-wrapped Cromwell (widdle) workflow for running FastQC tools on paired or unpaired reads.\n![fastqc flowchart](docs/fastqc-wf.png)\n"
    dependencies: [
      {
        name: "fastqc/0.11.9",
        url: "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
      }
    ]
    output_meta: {
      html_report_R1: "HTML report for the first mate fastq file.",
      zip_bundle_R1: "zipped report from FastQC for the first mate reads.",
      html_report_R2: "HTML report for read second mate fastq file.",
      zip_bundle_R2: "zipped report from FastQC for the second mate reads."
    }
}

output {
 File html_report_R1  = firstMateHtml.renamedOutput
 File zip_bundle_R1   = firstMateZip.renamedOutput
 File? html_report_R2 = secondMateHtml.renamedOutput
 File? zip_bundle_R2  = secondMateZip.renamedOutput
}

}

# ===================================
#            MAIN STEP
# ===================================
task runFastQC {
input {
        Int    jobMemory = 6
        Int    timeout   = 20
        Int    javaHeap  = 4
        Int?   threads
        File   inputFastq
        String modules = "perl/5.28 java/11 fastqc/0.11.9"
}

command <<<
 set -euo pipefail
 FASTQC=$(which fastqc)
 JAVA=$(which java)
 export _JAVA_OPTIONS=-Xmx~{javaHeap}G
 perl $FASTQC ~{inputFastq} --java=$JAVA ~{"-t" + threads}--noextract --outdir "."
>>>

parameter_meta {
 jobMemory: "Memory allocated to fastqc."
 javaHeap: "Memory allocated to java heap, in G."
 threads: "Threads param for fastqc"
 inputFastq: "Input fastq file, gzipped."
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File html_report_file = "~{basename(inputFastq, '.fastq.gz')}_fastqc.html"
  File zip_bundle_file  = "~{basename(inputFastq, '.fastq.gz')}_fastqc.zip"
}
}

# =================================================
#      RENAMING STEP - IF WE HAVE CUSTOM PREFIX
# =================================================
task renameOutput {
input {
  Int  jobMemory = 2
  File inputFile
  String extension
  String customPrefix
  Int timeout    = 1
}

parameter_meta {
 inputFile: "Input file, html or zip."
 extension: "Extension for a file (without leading dot)."
 customPrefix: "Prefix for making a file."
 jobMemory: "Memory allocated to this task."
 timeout: "Timeout, in hours, needed to override imposed limits."
}

command <<<
 set -euo pipefail
 if [[ ~{basename(inputFile)} != "~{customPrefix}.~{extension}" ]];then 
   ln -s ~{inputFile} "~{customPrefix}.~{extension}"
 else
   ln -s ~{inputFile}
 fi
>>>

runtime {
  memory:  "~{jobMemory} GB"
  timeout: "~{timeout}"
}


output {
  File renamedOutput = "~{customPrefix}.~{extension}"
}
}
