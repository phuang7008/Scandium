#!/hgsc_software/python/anaconda3/bin/python

import sys, getopt

def main(argv):
	in_file = ''
	try:
		opts, args = getopt.getopt(argv, "i:", ["ifile="])
	except getopt.GetoptError:
		print("extractIntrons.py -i <inputfile>" )
	
	for opt, arg in opts:
		if opt in ("-i", "--ifile"):
			in_file = arg
	#print("The input file is ", in_file)
	processFile(in_file)

def processFile(file_in):
	with open(file_in, 'r') as rfh:
		for line in rfh:
			if "chrom" not in line:
				items  = line.rstrip("\n").split()
				starts = items[6].split(",")
				ends   = items[7].split(",")
				prev_start = int(starts[0])
				prev_end   = int(ends[0])

				# some version of gene annotation has 'chr' in front of chromosome id, so let's remove them
				items[1] = items[1].replace("chr", "")
				items[1] = items[1].replace("Chr", "")
				items[1] = items[1].replace("CHR", "")

				for idx in range(1,int(items[5])):
					if (len(items) == 9):
						print("%s\t%d\t%d\t%s\t%s" % (items[1], prev_end, int(starts[idx]), items[0]+"="+items[8], items[8]))
					else:
						print("%s\t%d\t%d\t%s\t%s" % (items[1], prev_end, int(starts[idx]), items[0], '.'))

					prev_start = int(starts[idx])
					prev_end = int(ends[idx])


if __name__ == "__main__":
	main(sys.argv[1:])

