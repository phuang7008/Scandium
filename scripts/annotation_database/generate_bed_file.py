#!/hgsc_software/python/anaconda3/bin/python

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

				# some version of gene annotation has 'chr' in front of chromosome id, so let's remove them
				items[1] = items[1].replace("chr", "");
				items[1] = items[1].replace("Chr", "");
				items[1] = items[1].replace("CHR", "");

				for idx in range(0,int(items[5])):
					# take care of the reverse strand
					exon_id = idx
					if (items[2] == "-"):
						exon_id = int(items[5]) - idx - 1

					if (len(items) == 9): 
						print("%s\t%d\t%d\t%s\t%s" % (items[1], int(starts[idx]), int(ends[idx]), items[0]+"_exon_"+str(exon_id)+"_"+items[5]+"="+items[8], items[8]))
					else:
						# for CCDS as it doesn't have name2 field
						print("%s\t%d\t%d\t%s\t%s" % (items[1], int(starts[idx]), int(ends[idx]), items[0]+"_exon_"+str(exon_id)+"_"+items[5], "."))



if __name__ == "__main__":
	main(sys.argv[1:])
