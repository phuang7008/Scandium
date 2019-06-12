#!/hgsc_software/python/latest/bin/python
# Python 3.6.5

# this script is used to parse Scandium Summary Report just like we did for alignstats!
#
import argparse
import os

def main(args_in):
	# since we have two types of summary reports, one is for WGS, while another one is for Capture
	# we need to parse them separately
	#
	type_of_reports = ["Capture", "WGS"]

	# In order to write the header only once, First we need to remove output files if they exist
	#
	for tmp_type in type_of_reports:
		file_exists = os.path.isfile(args_in.output+"-"+tmp_type+"."+"csv")
		if file_exists:
			os.remove(args_in.output+"-"+tmp_type+"."+"csv")

	# open input file for reading
	#
	with open (args_in.input, 'r') as fh:
		for cur_dir in fh:
			# one line is one directory, need to check for all the summary files in this directory
			#
			cur_dir = cur_dir.rstrip("\n")
			for tmp_type in type_of_reports:
				for tmp_file in os.listdir(cur_dir):
					if tmp_file.endswith(tmp_type + "_Coverage_Summary_Report.txt"):
						tmp_input_file = os.path.abspath(os.path.join(cur_dir, tmp_file))
						#print(tmp_input_file)
				
						# now we need to process one file at a time
						#
						if not os.path.isfile(args_in.output+"-"+tmp_type+"."+"csv"):
							output_header(tmp_input_file, tmp_type, args_in)
							
						process_summary_report_file(tmp_input_file, args_in, tmp_input_file, tmp_type)

# process the summary file passed in and append everything to the keys
#
def process_summary_report_file(file_in, args_in, file_name_in, type_in):
	# open files one for reading and one for writing
	#
	#fout = open(os.path.abspath(os.path.join(args_in.outDir, args_in.output+"."+"csv")), 'a')
	fout = open (args_in.output+"-"+type_in+"."+"csv", 'a')
	fout.write(file_name_in)

	with open (file_in, 'r') as fin:

		for line in fin:
			if not line.startswith('#'):
				line = line.rstrip("\n")
				tmp_list = line.split("\t")

				if len(tmp_list) == 2:
					fout.write(",")
					fout.write(tmp_list[1])

		fout.write("\n")
	fout.close()


# write header here
#
def output_header(file_in, type_in, args_in):
	# open file for writing
	#
	#fout = open (os.path.abspath(os.path.join(args_in.outDir, args_in.output+"."+"csv")), 'w')
	fout = open (args_in.output+"-"+type_in+"."+"csv", 'w')
	fout.write("Samples")

	with open (file_in, 'r') as fin:

		for line in fin:
			if not line.startswith('#'):
				line = line.rstrip("\n")
				tmp_list = line.split("\t")

				if len(tmp_list) == 2:
					fout.write(",")
					fout.write(tmp_list[0])

		fout.write("\n")

	fout.close()


#########################################################################################
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Process Scandium Summary Reports")
	parser.add_argument('-i', '--input', help="Input file that contains all the directories to the Scandium Summary Report [MANDATORY]", required=True)
	parser.add_argument('-o', '--output', help="Output file name's prefix to be written in .csv format [MANDATORY]", required=True)
	#parser.add_argument('-d', '--outDir', help="Output directory", required=True)
	
	# grap all the arguments
	#
	args = parser.parse_args()

	main(args)
