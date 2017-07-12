#!/users/phuang/anaconda3/bin/python

import sys, getopt

def main(argv):
	in_file = ''
	try:
		opts, args = getopt.getopt(argv, "i:", ["ifile="])
	except getopt.GetoptError:
		print("getTranscripts.py -i <inputfile>" )
	
	for opt, arg in opts:
		if opt in ("-i", "--ifile"):
			in_file = arg
			#print("in\n")
	print("The input file is ", in_file)
	processFile(in_file)

def processFile(file_in):
	wfh = open(file_in + ".transcript", "w")

	with open(file_in, 'r') as rfh:
		for line in rfh:
			if "chrom" not in line:
				items  = line.rstrip("\n").split()
				wfh.write("%s\t%d\t%d\n" % (items[0], int(items[1]), int(items[2])))
				#wfh.write("%s\t%d\t%d\n" % (items[1], int(items[3]), int(items[4])))


if __name__ == "__main__":
	main(sys.argv[1:])

