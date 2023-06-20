#!/hgsc_software/python/anaconda3/bin/python

import sys, getopt

def main(argv):
	in_file = ''
	try:
		opts, args = getopt.getopt(argv, "i:", ["ifile="])
	except getopt.GetoptError:
		print("rearrange_miRNA_bed.py -i <inputfile>" )
	
	for opt, arg in opts:
		if opt in ("-i", "--ifile"):
			in_file = arg
			#print("in\n")
	#print("The input file is ", in_file)
	processFile(in_file)

def processFile(file_in):

	with open(file_in, 'r') as rfh:
		for line in rfh:
			if "chrom" not in line:
				items  = line.rstrip("\n").split()
				cds_length = int(items[2]) - int(items[1])

				print("%s\t%d\t%d\t%s\t%s" % (items[0], int(items[1]), int(items[2]), items[3]+"_exon_0_1="+items[4]+"="+str(items[1])+"="+str(items[2])+"="+str(cds_length), items[4]))


if __name__ == "__main__":
	main(sys.argv[1:])

