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
require 'thinreports'
require 'gruff'
require 'json'

### Define modules and classes here


### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "A script description here"
opts.separator ""
opts.on("-i","--infile", "=INFILE","Input file") {|argument| options.infile = argument }
opts.on("-o","--outfile", "=OUTFILE","Output file") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

json = JSON.parse(IO.readlines(options.infile).join)
patient_id = json["patient"]
sample_id = json["sample"]

# Generate empty report
report = Thinreports::Report.new

report.use_layout 'metagenomics_default.tlf', :default => true
report.use_layout 'metagenomics_cover.tlf', id: :cover
report.use_layout 'metagenomics_tool_output.tlf', id: :tool
report.use_layout 'metagenomics_ariba.tlf', id: :ariba

# Cover page
report.start_new_page layout: :cover do |page|

  page.item(:header).value("Patient: " + json["patient"])

  page.item(:analysis_date).value(Time.now.strftime("%d/%m/%Y %H:%M"))
  page.item(:patient_name).value(json["patient"])
  page.item(:ordered_by).value("Dr. Leypoldt ")
  page.item(:order_number).value("SF_Liquor_Metagenome_Leypoldt_01")
  
  page.item(:sample).value(json["metadata"][0]["sample_type"])
  page.item(:sequencer).value(json["metadata"][0]["platform"])
  page.item(:read_type).value(json["metadata"][0]["read_type"])
  page.item(:readcount).value(json["bwa"]["read_count"])
  percentage_host = ( json["bwa"]["read_mapped"].to_f / json["bwa"]["read_count"].to_f ).round(4)
  host_reads = json["bwa"]["read_mapped"] + " (#{percentage_host.to_s} perc.)"
  page.item(:host_reads).value(host_reads)
  page.item(:pipeline_version).value("0.1")
  
  report.list.add_row do |row|
	  row.item(:table_text).value("Automatischer Report der IKMB Metagenom Pipeline.")
  end
  
end

# Pathoscope results
report.start_new_page layout: :tool do |page|
     
  pathoscope = json["pathoscope"]
  
  page.item(:header).value("Patient: " + json["patient"])
  
  page.item(:title).value("Pathoscope")
  
  g = Gruff::Pie.new
  g.title = "Prädizierte Bakterien/Viren"
  g.theme_pastel
  #g.LABEL_MARGIN = 5
  
  g.legend_font_size = 14
  g.marker_font_size = 14
  g.title_font_size = 24
  
  pathoscope.each do |species,abundance|
    g.data(species,abundance.to_f)
  end
  
  g.write("pie_#{sample_id}_pathoscope.png")
  
  page.item(:chart).src("pie_#{sample_id}_pathoscope.png")
  
  page.item(:reference).value("https://github.com/PathoScope/PathoScope\nPMID: 23843222")
  
  pathoscope.each do |species,abundance|
    report.list.add_row do |row|
  	  row.item(:species).value("#{species}")
  	  row.item(:abundance).value("#{abundance}")
    end
  end
  
end

# Metaphlan results
report.start_new_page layout: :tool do |page|
  
  metaphlan = json["metaphlan"]
  
  page.item(:header).value("Patient: " + json["patient"])
  
  page.item(:title).value("Metaphlan2")  
    
  g = Gruff::Pie.new
  g.title = "Prädizierte Bakterien/Viren"
  g.theme_pastel
  #g.LABEL_MARGIN = 5
  
  g.legend_font_size = 14
  g.marker_font_size = 14
  g.title_font_size = 24
    
  metaphlan.each do |species,abundance|
    g.data(species,abundance.to_f)
  end
  
  g.write("pie_#{sample_id}_metaphlan.png")
  
  page.item(:chart).src("pie_#{sample_id}_metaphlan.png")
  
  page.item(:reference).value("https://bitbucket.org/biobakery/metaphlan2\nPMID: 26418763")
    
  metaphlan.each do |species,abundance|
    report.list.add_row do |row|
  	  row.item(:species).value("#{species}")
  	  row.item(:abundance).value("#{abundance}")
    end
  end
  
end

if json.has_key?("ariba")

	report.start_new_page layout: :ariba do |page|
	
    page.item(:header).value("Patient: " + json["patient"])
  
		ariba = json["ariba"]
		ariba.each do |a|
			 report.list.add_row do |row|
				row.item(:ariba_result).value(a)
			end
		end
		
	end
	
end
report.generate(filename: "#{options.outfile}")

system("rm *_#{sample_id}*.png")

