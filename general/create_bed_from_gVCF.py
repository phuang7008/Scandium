#!/hgsc_software/python/latest/bin/python
# Python 3.6.5

# This script is used to process gVCF file and create a bed formatted file
#
import os
import argparse
import re

def main(args_in):
	# create the name of output file based on input file and open output file
	#
	fout = open(args_in.input + ".bed", "w")

	# open input file for reading
	#
	with open (args_in.input, 'r') as fin:
		for line in fin:
			if not line.startswith("#"):
				line = line.rstrip("\n")
				tmp_list = line.split("\t")
				patt = re.compile(r'END=(\d+);.*', re.M|re.I)
				#patt = re.compile(r'END=([0..9]+);.*', re.M|re.I)
				matchObj = re.search(patt, str(tmp_list[7]))
				if not matchObj:
					tmp_num = int(tmp_list[1]) - 1
					tmp_list.insert(1, str(tmp_num))
				else:
					tmp_list[1] = str(int(tmp_list[1]) - 1)
					tmp_list.insert(2, matchObj.group(1))

				out_string = ""
				for item in tmp_list:
					if (out_string == ""):
						out_string += item
					else:
						out_string += "\t"+item

				fout.write(out_string)
				fout.write("\n")


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Process gVCF file to create a bed file")
	parser.add_argument('-i', '--input', help="input gVCF file [MANDATORY]", required=True)

	args = parser.parse_args()

	main(args)
