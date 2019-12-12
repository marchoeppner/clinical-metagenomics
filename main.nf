// Pipeline variables

OUTDIR = params.outdir 

PATHOSCOPE_DB=file(params.pathoscope_db)

METAPHLAN_DB=file(params.metaphlan_db)

KAIJU_DB=file(params.kaiju_db)

REF = params.ref

inputFile=file(params.samples)

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
       	.into { inputTrim  }

process runFastp {

	input:
        set patientID,sampleID,Platform,R1,R2 from inputTrim

        output:
        set patientID,sampleID,file("${left}"),file("${right}") into fastpOutput

	script:
	
	left = file(R1).getBaseName() + "_trimmed.fastq.gz"
        right = file(R2).getBaseName() + "_trimmed.fastq.gz"
        json = file(R1).getBaseName() + ".fastp.json"
        html = file(R1).getBaseName() + ".fastp.html"

	"""
                fastp --in1 $R1 --in2 $R2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
        """
	
}

inputMerge = fastpOutput.groupTuple(by: [0,1])

process Merge {

        tag "${patientID}|${sampleID}"

        input:
        set patientID,sampleID,forward_reads,reverse_reads from inputMerge

        scratch true

        output:
        set patientID,sampleID,file(left_merged),file(right_merged) into inputBwa

        script:
        left_merged = sampleID + "_R1.fastq.gz"
        right_merged = sampleID + "_R2.fastq.gz"

        """
                zcat ${forward_reads.join(" ")} | gzip > $left_merged
                zcat ${reverse_reads.join(" ")} | gzip > $right_merged
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
	bwa mem -M -t ${task.cpus} ${REF} $left $right | samtools sort -O BAM - > $bam
	samtools stats $bam > $stats
	
   """

}

// We extract the reads not mapping to the host genome
process extractUnmapped {

   publishDir "${OUTDIR}/${patientID}/${sampleID}/Host", mode: 'copy'

   input:
   set patientID,sampleID,file(bam) from alignedBam

   output:
   set patientID,sampleID,file(left),file(right) into (inputPathoscopeMap,inputMetaphlan,inputKaiju)
  
   script:
   left = sampleID + "_R1.fastq.gz"
   right = sampleID + "_R2.fastq.gz"

   """
	samtools fastq -f 4 -1 $left -2 $right $bam
   """

}

process runPathoscopeMap {

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
   file "v_metaphlan.txt" into version_metaphlan

   script:

   metaphlan_out = sampleID + "_metaphlan_report.txt"
   sam_out = sampleID + "_metaphlan.sam.bz2"

   """
     metaphlan2.py --version &> v_metaphlan.txt
     metaphlan2.py $left_reads,$right_reads --bowtie2db $METAPHLAN_DB --samout $sam_out --bowtie2out $bowtie_out --nproc ${task.cpus} --input_type fastq -o $metaphlan_out

   """

}

process runKaiju {

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

	publishDir "${OUTDIR}/Reports"

	input:
	set patientID,sampleID,pathoscope from outputPathoscopeId
	file(kaiju) from outputKaijuReport
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
			--kaiju $kaiju \
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

