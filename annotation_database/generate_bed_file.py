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
	print("The input file is ", in_file)
	processFile(in_file)

def processFile(file_in):
	wfh = open(file_in + ".bed", "w")

	with open(file_in, 'r') as rfh:
		for line in rfh:
			if "chrom" not in line:
				items  = line.rstrip("\n").split()
				#starts = items[6].split("|")
				starts = items[6].split(",")
				#ends   = items[7].split("|")
				ends   = items[7].split(",")
				for idx in range(0,int(items[5])):
					#wfh.write("%s\t%d\t%d\t%s\t%s\t%d\n" % (items[1], int(starts[idx]), int(ends[idx]), items[0]+"_"+str(idx), '', idx))
					wfh.write("%s\t%d\t%d\t%s\t%s\n" % (items[1], int(starts[idx]), int(ends[idx]), items[0]+"_"+str(idx)+"="+items[8], items[8]))


if __name__ == "__main__":
	main(sys.argv[1:])

