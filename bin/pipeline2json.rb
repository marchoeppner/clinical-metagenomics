#!/bin/env ruby
# == NAME
# script_skeleton.rb
#
# == USAGE
# ./this_script.rb [ -h | --help ]
#[ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A skeleton script for Ruby
#
# == OPTIONS
# -h,--help Show help
# -i,--infile=INFILE input file
# -o,--outfile=OUTFILE : output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'
require 'json'

BUCKET = []
Entry = Struct.new(:level, :parent, :amount )

### Define modules and classes here

def children_nodes(node,bucket)
  
  answer = []
  children = bucket.select{|e| e.parent == node.level }
  
  if children.empty?
    BUCKET << node
      #puts "\t#{node.level}\t#{node.amount}"
  else
    children.each do |child|
        children_nodes(child,bucket)
    end
  end
    
end

def parse_pathoscope(file_name)
  
  warn "File not found (#{file_name})" unless File.exists?(file_name)
  
  lines = IO.readlines(file_name)
  
  stats = lines.shift
  header = lines.shift
  
  data = {}
  
  lines.each do |line|
    e = line.strip.split("\t")
    species = e[0].split("|")[-1].split(",")[0].gsub(/_/, ' ').strip.split(/\s/)[0..1].join(" ")
    abundance = ("%.4f" % e[1]).to_f
    next if abundance < 0.0001
    
    data[species] = abundance
    
  end
  
  return data
  
end

def parse_kaiju(file_name)
  
  warn "File not found (#{file_name})" unless File.exists?(file_name)
  
  lines = IO.readlines(file_name)
  
  data = {}
  
  lines.each do |line|
  
    next unless line.strip.match(/^\d+.*/)
    
    e = line.strip.split(/\s+/)
    data[e[2..-1].join(" ").gsub(/_/, ' ')] = e[0]
    
  end
  
  return data
  
end

def parse_metaphlan(file_name)
  
  warn "File not found (#{file_name})" unless File.exists?(file_name)
  
  bin = [] 
  data = {}
  
  IO.readlines(file_name).each do |line|
  
    # k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Potyviridae|g__Potyvirus|s__Dasheen_mosaic_virus    23.65111
  
    next unless line.match(/^k_.*/)
  
    levels = line.strip.split("|")
    
    # Exclude subspecies and strain level classifications for readability
    next if levels[-1].match(/^t__.*/)
  
    deepest_level,amount = levels[-1].split(" ")
    parent_level = levels[-2]
  
    bin << Entry.new(deepest_level, parent_level, amount)

  end
  
  roots = bin.select{|e| e.parent.nil? }

  roots.each do |root|
      
    children_nodes(root,bin)
  
    BUCKET.each do |entry|
      data[entry.level[3..-1]] = entry.amount
    end
   
    BUCKET.clear
  
  end
  
  return data
  
end
  
def parse_bwa(file_name)
 
  warn "File not found (#{file_name})" unless File.exists?(file_name)
  
  lines = IO.readlines(file_name)
  data = {}
  lines.each do |line|
    if line.match(/^SN.*/)
      e = line.strip.split(/\t/)
      case e[1].strip
        when "raw total sequences:"
          data["read_count"] = e[2]
        when "reads mapped:"
          data["read_mapped"] = e[2]
        when "reads unmapped:"
          data["reads_unmapped"] = e[2]
      end
    end
  end
  
  return data
  
end

def parse_ariba(file_name)

	data = []
	lines = IO.readlines(file_name)
	header = lines.shift
	
	lines.each do |line|
		e = line.strip.split("\t")
		description = e[-1]
		data << description unless data.include?(description)
	end
	return data
end

def parse_nextflow_version(file_name)

	lines = IO.readlines(file_name)
	return lines.shift.strip
end

def parse_pipeline_version(file_name)

	lines = IO.readlines(file_name)
	return lines.shift.strip
	
end

def parse_logs(files)
  
  data = {}

  files.each do |file|
    
    lines = IO.readlines(file)
    
    if file.include?("bwa")
      software,version = lines[2].strip.split(/\s+/)
      data["BWA"] = version
    elsif file.include?("kaiju")
      software,version = lines[2].strip.split(/\s+/)
      data[software] = version
    elsif file.include?("pathoscope")
      software,version = lines[1].strip.split(/\s+/)
      data[software] = version
    elsif file.include?("metaphlan")
      e = lines[1].strip.split(/\s+/)
      data[e[0]] = e[2]
    else
      software,version = lines.shift.strip.split(/\s+/)
      data[software] = version   
    end 

  end
  
  return data
  
end
### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "A script description here"
opts.separator ""
opts.on("-f","--folder", "=Folder","Input folder") {|argument| options.folder = argument }
opts.on("-p","--patient_id", "=PATIENT","Patient ID") {|argument| options.patient = argument }
opts.on("-q","--sample_id", "=SAMPLE","SAMPLE ID") {|argument| options.sample = argument }
opts.on("-s","--samplesheet", "=Samples","Samplesheet") {|argument| options.samplesheet = argument }
opts.on("-t","--pathoscope", "=PATHOSCOPE","Pathoscope report") {|argument| options.pathoscope = argument }
opts.on("-m","--metaphlan", "=METAPHLAN","Metaphlan") {|argument| options.metaphlan = argument }
opts.on("-k","--kaiju", "=KAIJU","Kaiju") {|argument| options.kaiju = argument }
opts.on("-b","--bam-stats", "=BWA","BWA stats") {|argument| options.bwa = argument }
opts.on("-o","--outfile", "=OUTFILE","Output file") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

options.outfile ? output_stream = File.new(options.outfile,'w') : output_stream = $stdout

# We assume the follwing structure:
# Patient ID
# => sample ID
# => => Pathoscope/
# => => Kauju/
# => => Metaphlan/
# => => Host/
# => => trimgalore/fastqc

# Parse metadata from the sample sheet
metadata = {}

lines = IO.readlines(options.samplesheet)

header = lines.shift

lines.each do |line|
  patient,sample,sample_type,read_type,platform,left,right = line.strip.split(";")
  data = { "sample_type" => sample_type, "read_type" => read_type, "platform" => platform, "left" => left, "right" => right }
  if metadata.has_key?(patient)
    metadata[patient].has_key?(sample) ? metadata[patient][sample] << data :  metadata[patient][sample] = [ data ]
  else
    metadata[patient] = {}
    metadata[patient][sample] = [ data ]
  end
end  

# Traverse the output directory

            
sbucket = { "patient" => options.patient, "sample" => options.sample, "metadata" => metadata[options.patient][options.sample] }
        
# parse run info  
logs = Dir["*.txt"]
run_info = parse_logs(logs)
  
sbucket["run_info"] = run_info
  
# parse Pathoscope
if options.pathoscope 
  pathoscope = parse_pathoscope(options.pathoscope)
        
  sbucket["pathoscope"] = pathoscope
end
        
# parse Kaiju
if options.kaiju
  kaiju = parse_kaiju(options.kaiju)

  sbucket["kaiju"] = kaiju
end

if options.metaphlan      
  # Parse Metaphlan
  metaphlan = parse_metaphlan(options.metaphlan)
  sbucket["metaphlan"] = metaphlan
end
        
if options.bwa
  bwa_stats = parse_bwa(options.bwa)
                
  sbucket["bwa"] = bwa_stats
end

if options.ariba		

  ariba = parse_ariba(options.ariba)
	unless ariba.empty?
	  sbucket["ariba"] = ariba
	end
  
end
        
output_stream.puts sbucket.to_json
      
output_stream.close
