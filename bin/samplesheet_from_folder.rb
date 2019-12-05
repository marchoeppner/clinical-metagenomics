#!/bin/env ruby

require 'optparse'
require 'ostruct'

### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "Reads Fastq files from a folder and writes a sample sheet to STDOUT"
opts.separator ""
opts.on("-f","--folder", "=FOLDER","Folder to scan") {|argument| options.folder = argument }
opts.on("-t","--type","=SAMPLE_TYPE","Sample type") {|argument| options.source = argument }
opts.on("-r","-read_type","=READ_TYPE","Read type") {|argument| options.read_type = argument }
opts.on("-s","--sanity", "Perform sanity check of md5 sums") { options.sanity = true }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

abort "Folder not found (#{options.folder})" unless File.directory?(options.folder)

date = Time.now.strftime("%Y-%m-%d")

fastq_files = Dir["#{options.folder}/*_R*.fastq.gz"]

groups = fastq_files.group_by{|f| f.split("/")[-1].split(/_L0/)[0] }

warn "Building input sample sheet from FASTQ folder"
warn "Performing sanity check on md5sums" if options.sanity

puts "PatientID;SampleID;SampleType;ReadType;Platform;R1;R2"

#G00076-L2_S19_L003_R1_001.fastq.gz

# group = the library id, may be split across lanes
groups.each do |group, files|

	warn "...processing library #{group}"

	pairs = files.group_by{|f| f.split("/")[-1].split(/_R[1,2]/)[0] }

	pairs.each do |p,reads|

        	left,right = reads.sort.collect{|f| File.absolute_path(f)}

		if options.sanity
			Dir.chdir(options.folder) {
				[left,right].each do |fastq|
					fastq_simple = fastq.split("/")[-1].strip
					raise "Aborting - no md5sum found for fastq file #{fastq}" unless File.exists?(fastq_simple + ".md5")
					status = `md5sum -c #{fastq_simple}.md5`
					raise "Aborting - failed md5sum check for #{fastq}" unless status.strip.include?("OK")
				end
			}
		end
                sample = group.split("-")[0]

        	puts "Indiv_#{sample};Sample_#{sample};#{options.source};#{options.read_type};Illumina;#{left};#{right}"
	end
end


