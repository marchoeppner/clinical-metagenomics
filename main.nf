// Pipeline variables

OUTDIR = params.outdir 

PATHOSCOPE_INDEX_DIR=file(params.pathoscope_index_dir)
PATHOSCOPE=file(params.pathoscope)

METAPHLAN=file(params.metaphlan)
METAPHLAN_PKL=file(params.metaphlan_pkl)
METAPHLAN_DB=params.metaphlan_db

FASTQC=file(params.fastqc)

TRIMMOMATIC = params.trimmomatic
leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen
adapters = params.adapters

FOLDER=file(params.folder)

// Logging and reporting

logParams(params, "nextflow_parameters.txt")

VERSION = "1.0" 
// Header log info 

log.info "=========================================" 
log.info "IKMB pipeline version v${VERSION}" 
log.info "Nextflow Version: $workflow.nextflow.version" 
log.info "Command Line: $workflow.commandLine" 
log.info "=========================================" 


// Starting the workflow

Channel
        .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
        .map { prefix, file1, file2 ->  tuple(prefix.substring(0,6), file1, file2) }
        .groupTuple()
        .into { inputMerge }

process Merge {

	tag "${id}"
        publishDir("${OUTDIR}/Data/${id}")

        input:
        set id,forward_reads,reverse_reads from inputMerge

        output:
        set id,file(left_merged),file(right_merged) into inputTrimmomatic

        script:
        left_merged = id + "_R1.fastq.gz"
        right_merged = id + "_R2.fastq.gz"

        """
                zcat ${forward_reads.join(" ")} | gzip > $left_merged
		zcat ${reverse_reads.join(" ")} | gzip > $right_merged
        """
}

process Trimmomatic {

   tag "${id}"
   publishDir "${OUTDIR}/trimmomatic", mode: 'copy'

   input:
   set id,file(left_reads),file(right_reads) from inputTrimmomatic

   output:
   set id,file("${id}_R1_paired.fastq.gz"), file("${id}_R2_paired.fastq.gz") into inputFastqc,inputPathoscopeMap,inputMetaphlan

   script:

    """
        java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads 8 $left_reads $right_reads \
        ${id}_R1_paired.fastq.gz ${id}_1U.fastq.gz ${id}_R2_paired.fastq.gz ${id}_2U.fastq.gz \
        ILLUMINACLIP:${TRIMMOMATIC}/adapters/${adapters}:2:30:10:3:TRUE\
        LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen} && sleep 5
   """

}

process Fastqc {

   tag "${id}"
   publishDir "${OUTDIR}/Data/${id}/FastQC/", mode: 'copy'

    input:
    set id, file(left_reads), file(right_reads) from inputFastqc

    output:
    set file("*.zip"), file("*.html") into outputFastqc

    script:
    """
    $FASTQC -t 1 -o . ${left_reads} ${right_reads}
    """

}

process runPathoscopeMap {

   tag "${id}"
   publishDir "${OUTDIR}/Data/${id}/Pathoscope"

   input:
   set id,file(left_reads),file(right_reads) from inputPathoscopeMap

   output:
   set id,file(pathoscope_sam) into inputPathoscopeId

   script:
   pathoscope_sam = id + ".sam"

   """
	$PATHOSCOPE MAP -1 $left_reads -2 $right_reads -indexDir $PATHOSCOPE_INDEX_DIR -filterIndexPrefixes hg19_rRNA \
	-targetIndexPrefix A-Lbacteria.fa,M-Zbacteria.fa,virus.fa -outAlign $pathoscope_sam -expTag $id -numThreads 8
   """

}

process runPathoscopeId {

   tag "${id}"
   publishDir "${OUTDIR}/Data/${id}/Pathoscope"

   input:
   set id,file(samfile) from inputPathoscopeId

   output:
   set id,file(pathoscope_tsv),file(pathoscope_sam) into outputPathoscopeId

   script:

   pathoscope_sam = "updated_" + samfile
   pathoscope_tsv = id + "-sam-report.tsv"

   """
	$PATHOSCOPE ID -alignFile $samfile -fileType sam -expTag $id
   """

}

process runMetaphlan {

   tag "${id}"
   publishDir "${OUTDIR}/Data/${id}/Metaphlan2"

   input:
   set id,file(left_reads),file(right_reads) from inputMetaphlan

   output:
   set id,file(metaphlan_out) into outputMetaphlan

   script:

   metaphlan_out = id + "_metaphlan.out"

   """
     $METAPHLAN --mpa_pkl $METAPHLAN_PKL --bowtie2db $METAPHLAN_DB --nproc 16 --input_type fastq <(zcat $left_reads $right_reads ) > $metaphlan_out

   """

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}




//#############################################################################################################
//#############################################################################################################
//
// FUNCTIONS
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}

