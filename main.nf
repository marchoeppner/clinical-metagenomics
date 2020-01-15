// Pipeline variables

OUTDIR = params.outdir 

PATHOSCOPE_DB=file(params.pathoscope_db)

METAPHLAN_DB=file(params.metaphlan_db)
METAPHLAN_PKL = file(params.metaphlan_pkl)

REF = params.ref

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
       	.set { inputTrim  }

process runFastp {

	input:
        set patientID,sampleID,Platform,R1,R2 from inputTrim

        output:
        set patientID,sampleID,file("${left}"),file("${right}") into inputBwa

	script:
	
	left = file(R1).getBaseName() + "_trimmed.fastq.gz"
        right = file(R2).getBaseName() + "_trimmed.fastq.gz"
        json = file(R1).getBaseName() + ".fastp.json"
        html = file(R1).getBaseName() + ".fastp.html"

	"""
                fastp --in1 $R1 --in2 $R2 --out1 $left --out2 $right \
		--unpaired1 unpaired.fastq.gz --unpaired2 unpaired.fastq.gz \
		--detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
        """
	
}

process runBwa {

   publishDir "${OUTDIR}/${patientID}/${sampleID}/Host", mode: 'copy'


   input:
   set patientID,sampleID,file(left),file(right) from inputBwa

   output:
   set patientID,sampleID,file(bam) into alignedBam
   file(stats) into BamStats

   script:

   bam = sampleID + ".bam"
   stats = sampleID + "_bwa_stats.txt"

   samtools_version = "v_samtools.txt"

   """
	bwa mem -M -t ${task.cpus} ${REF} $left $right | /opt/samtools/1.9/bin/samtools sort -m 16G -@4 -O BAM - > $bam
	/opt/samtools/1.9/bin/samtools stats $bam > $stats
	
   """

}

// We extract the reads not mapping to the host genome
process extractUnmapped {

   publishDir "${OUTDIR}/${patientID}/${sampleID}/Host", mode: 'copy'

   input:
   set patientID,sampleID,file(bam) from alignedBam

   output:
   set patientID,sampleID,file(left),file(right) into inputMetaphlan
   set patientID,sampleID,file(left),file(right) into inputBBmap
     
   script:
   left = sampleID + "_R1.fastq.gz"
   right = sampleID + "_R2.fastq.gz"

   """
	/opt/samtools/1.9/bin/samtools fastq -f 4 -1 $left -2 $right $bam
   """

}

process runBBmapFix {


   input:
   set patientID,sampleID,file(left_reads),file(right_reads) from inputBBmap

   output:
   set patientID,sampleID,file(left_fixed),file(right_fixed) into inputPathoscopeMap

   script:

   left_fixed = left_reads.getBaseName() + ".fixed.fastq.gz"
   right_fixed = right_reads.getBaseName() + ".fixed.fastq.gz"
   
   """
	BBmap repair.sh in1=$left_reads in2=$right_reads out1=$left_fixed out2=$right_fixed outsingle=singletons.fastq.gz
   """

}

process runPathoscopeMap {

   publishDir "${OUTDIR}/${patientID}/${sampleID}/Pathoscope", mode: 'copy'

   input:
   set patientID,sampleID,file(left_reads),file(right_reads) from inputPathoscopeMap

   output:
   set patientID,sampleID,file(pathoscope_sam) into inputPathoscopeId

   script:
   pathoscope_sam = sampleID + ".sam"

   """
	pathoscope MAP -1 $left_reads -2 $right_reads -indexDir $PATHOSCOPE_DB -filterIndexPrefixes hg19_rRNA \
	-targetIndexPrefix A-Lbacteria.fa,M-Zbacteria.fa,virus.fa -outAlign $pathoscope_sam -expTag $sampleID -numThreads 8
   """

}

process runPathoscopeId {

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

   publishDir "${OUTDIR}/${patientID}/${sampleID}/Metaphlan2", mode: 'copy'

   input:
   set patientID,sampleID,file(left_reads),file(right_reads) from inputMetaphlan

   output:
   file(metaphlan_out) into outputMetaphlan
   file(sam_out)
   file(bowtie_out)
   file "v_metaphlan.txt" into version_metaphlan

   script:

   metaphlan_out = sampleID + "_metaphlan_report.txt"
   sam_out = sampleID + "_metaphlan.sam.bz2"
   bowtie_out = sampleID + "_bowtie.out"

   """
     metaphlan2.py --version &> v_metaphlan.txt
     metaphlan2.py --mpa_pkl  $METAPHLAN_PKL --bowtie2db $METAPHLAN_DB --samout $sam_out --bowtie2out $bowtie_out --nproc ${task.cpus} --input_type fastq <(zcat $left_reads $right_reads ) > $metaphlan_out

   """

}

process makeReport {

	publishDir "${OUTDIR}/Reports"

	input:
	set patientID,sampleID,pathoscope from outputPathoscopeId
	file(metaphlan) from outputMetaphlan
	file(bwa) from BamStats

	output:
	file(report) into Report

	script:
	
	json = patientID + "_" + sampleID + ".json"
	report = patientID + "_" + sampleID + ".pdf"

	"""
		cp $baseDir/assets/*.tlf . 
		ruby $baseDir/bin/pipeline2json.rb \
			--metaphlan $metaphlan \
			--pathoscope $pathoscope \
			--patient_id $patientID \
			--sample_id $sampleID \
			--samplesheet $inputFile \
			--bam-stats $bwa \
			-o $json
		ruby $baseDir/bin/json2report.rb -i $json -o $report
	"""

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}


