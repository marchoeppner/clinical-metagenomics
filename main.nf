// Pipeline variables

OUTDIR = params.outdir 

PATHOSCOPE_INDEX_DIR=file(params.pathoscope_index_dir)

METAPHLAN_PKL=file(params.metaphlan_pkl)
METAPHLAN_DB=params.metaphlan_db

KAIJU_REPORT=file(params.kaiju_report)
KAIJU_DB=params.kaiju_db


ARIBA_DB=file(params.ariba_db)

REF = file(params.ref)

inputFile=file(params.samples)

params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2


// Logging and reporting

logParams(params, "nextflow_parameters.txt")

params.version = "0.1" 
// Header log info 

log.info "=========================================" 
log.info "IKMB pipeline version v${params.version}" 
log.info "Nextflow Version: $workflow.nextflow.version" 
log.info "Command Line: $workflow.commandLine" 
log.info "=========================================" 


// Starting the workflow
Channel.from(inputFile)
       	.splitCsv(sep: ';', header: true)
       	.set { inputTrimgalore }

process runTrimgalore {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/trimgalore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else params.saveTrimmed ? filename : null
        }

   input:
   set val(patientID),val(sampleID),val(sampleType),val(readType),val(platform),left,right from inputTrimgalore

   output:
   set val(patientID),val(sampleID),file("*val_1.fq.gz"),file("*val_2.fq.gz") into trimgaloreOutput
   file "*trimming_report.txt" into trimgalore_results, trimgalore_logs 
   file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
   file "v_trimgalore.txt" into version_trimgalore

   script:

    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''

    """
    trim_galore --paired --fastqc --length 35 --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $left $right
    trim_galore --version &> v_trimgalore.txt
    """

}

inputMerge = trimgaloreOutput.groupTuple(by: [0,1])

process Merge {

        tag "${patientID}|${sampleID}"

        input:
        set patientID,sampleID,forward_reads,reverse_reads from inputMerge

        scratch true

        output:
        set patientID,sampleID,file(left_merged),file(right_merged) into inputPathoscopeMap,inputBwa

        script:
        left_merged = sampleID + "_R1.fastq.gz"
        right_merged = sampleID + "_R2.fastq.gz"

        """
                zcat ${forward_reads.join(" ")} | gzip > $left_merged
                zcat ${reverse_reads.join(" ")} | gzip > $right_merged
        """
}

process runBwa {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Host", mode: 'copy'

   input:
   set patientID,sampleID,file(left),file(right) from inputBwa

   output:
   set patientID,sampleID,file(bam) into alignedBam
   file(stats) into BamStats

   file(samtools_version) into version_samtools

   script:

   bam = sampleID + ".bam"
   stats = sampleID + "_bwa_stats.txt"

   samtools_version = "v_samtools.txt"

   """
        samtools --version &> $samtools_version
	bwa mem -M -t ${task.cpus} ${REF} $left $right | samtools sort -O BAM - > $bam
	samtools stats $bam > $stats
	
   """

}

// We extract the reads not mapping to the host genome
process extractUnmapped {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Host", mode: 'copy'

   input:
   set patientID,sampleID,file(bam) from alignedBam

   output:
   set patientID,sampleID,file(left),file(right) into inputMetaphlan,inputKaiju,inputAriba
  
   script:
   left = sampleID + "_R1.fastq.gz"
   right = sampleID + "_R2.fastq.gz"

   """
	samtools fastq -f 4 -1 $left -2 $right $bam
   """

}

process runAriba {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Ariba", mode: 'copy'

   input:
   set patientID,sampleID,file(left),file(right) from inputAriba

   output:
   file(report) into AribaReport

   script:

   report = sampleID + "_report.txt"

   """
        ariba run $ARIBA_DB $left $right out.${sampleID}.run && cp out.${sampleID}.run/report.tsv report.${sampleID}.tsv
   """  

}

process runPathoscopeMap {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Pathoscope", mode: 'copy'

   input:
   set patientID,sampleID,file(left_reads),file(right_reads) from inputPathoscopeMap

   output:
   set patientID,sampleID,file(pathoscope_sam) into inputPathoscopeId
   file "v_pathoscope.txt" into version_pathoscope

   script:
   pathoscope_sam = sampleID + ".sam"

   """
	pathoscope MAP -1 $left_reads -2 $right_reads -indexDir $PATHOSCOPE_INDEX_DIR -filterIndexPrefixes hg19_rRNA \
	-targetIndexPrefix A-Lbacteria.fa,M-Zbacteria.fa,virus.fa -outAlign $pathoscope_sam -expTag $sampleID -numThreads 8
	pathoscope --version &> v_pathoscope.txt
   """

}

process runPathoscopeId {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Pathoscope", mode: 'copy'

   input:
   set patientID,sampleID,file(samfile) from inputPathoscopeId

   output:
   set patientID,sampleID,file(pathoscope_tsv) into outputPathoscopeId

   script:

   //pathoscope_sam = "updated_" + samfile
   pathoscope_tsv = sampleID + "-sam-report.tsv"

   """
	pathoscope ID -alignFile $samfile -fileType sam -expTag $sampleID
   """

}

process runMetaphlan {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Metaphlan2", mode: 'copy'

   input:
   set patientID,sampleID,file(left_reads),file(right_reads) from inputMetaphlan

   output:
   file(metaphlan_out) into outputMetaphlan
   file "v_metaphlan.txt" into version_metaphlan

   script:

   metaphlan_out = sampleID + "_metaphlan_report.txt"

   """
     metaphlan2.py --version &> v_metaphlan.txt
     metaphlan2.py --bowtie2db $METAPHLAN_DB --nproc ${task.cpus} --input_type fastq <(zcat $left_reads $right_reads ) > $metaphlan_out

   """

}

process runKaiju {

   tag "${patientID}|${sampleID}"

   input:
   set patientID,sampleID,file(left_reads),file(right_reads) from inputKaiju

   output:
   set patientID,sampleID,file(kaiju_out) into inputKaijuReport

   script:

   kaiju_out = sampleID + "_kaiju.out"

   """
	kaiju -z 16 -t $KAIJU_DB/nodes.dmp -f $KAIJU_DB/kaiju_db.fmi -i <(gunzip -c $left_reads) -j <(gunzip -c $right_reads) -o $kaiju_out
   """


}

process runKaijuReport {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Kaiju", mode: 'copy'

   input:
   set patientID,sampleID,file(kaiju_out) from inputKaijuReport

   output:
   file(kaiju_report) into outputKaijuReport

   script:
   kaiju_report = sampleID + "_kaiju_report.txt"

   """
	kaijuReport -t $KAIJU_DB/nodes.dmp -n $KAIJU_DB/names.dmp -i $kaiju_out -r species -o $kaiju_report
   """

   
}

process makeReport {

	tag "Generating Report|${patientID}|${sampleID}"
	publishDir "${OUTDIR}/Reports"

	input:
	set patientID,sampleID,pathoscope from outputPathoscopeId
	file(kaiju) from outputKaijuReport
	file(metaphlan) from outputMetaphlan
	file(ariba) from AribaReport
	file(bwa) from BamStats

	output:
	file(report) into Report

	script:
	
	json = patientID + "_" + sampleID + ".json"
	report = patientID + "_" + sampleID + ".pdf"

	"""
		cp $baseDir/assets/*.tlf . 
		ruby $baseDir/bin/pipeline2json.rb \
			--ariba $ariba \
			--metaphlan $metaphlan \
			--pathoscope $pathoscope \
			--kaiju $kaiju \
			--patient_id $patientID \
			--sample_id $sampleID \
			--samplesheet $inputFile \
			--bam-stats $bwa \
			-o $json
		ruby $baseDir/bin/json2report.rb -i $json -o $report
	"""

}

process runMultiQCFastq {

    tag "Generating fastq level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Fastqc", mode: 'copy'

    input:
    file('*') from trimgalore_fastqc_reports.flatten().toList()

    output:
    file("fastq_multiqc*") into runMultiQCFastqOutput

    script:

    """
    multiqc -n fastq_multiqc *.zip *.html
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

