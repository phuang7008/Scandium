#!/users/phuang/anaconda3/bin/python

import sys, getopt

def main(argv):
	in_file = ''
	try:
		opts, args = getopt.getopt(argv, "i:", ["ifile="])
	except getopt.GetoptError:
		print("generate_bedfile.py -i <inputfile>" )
	
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
				starts = items[6].split(",")
				ends   = items[7].split(",")
				for idx in range(0,int(items[5])):
					print("%s\t%d\t%d\t%s\t%s\n" % (items[1], int(starts[idx]), int(ends[idx]), items[0]+"_"+str(idx)+"="+items[8]+"="+items[3]+"="+items[4], items[8]))


if __name__ == "__main__":
	main(sys.argv[1:])

