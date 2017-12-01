#!/users/phuang/anaconda3/bin/python

import sys, getopt

def main(argv):
	in_file = ''
	try:
		opts, args = getopt.getopt(argv, "i:", ["ifile="])
	except getopt.GetoptError:
		print("./generate_refseq_exon_bed.py -i <inputfile>" )
	
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
				# the order of item list is: name, chrom_id, strand, gene_start, gene_end, exon_count, exon_starts, exon_ends, gene_symbol
				items  = line.rstrip("\n").split()
				cds_start = int(items[5])
				cds_end   = int(items[6])
				starts = items[8].split(",")
				ends   = items[9].split(",")
				
				# some version of gene annotation has 'chr' in front of chromosome id, so let's remove them
				items[1] = items[1].replace("chr", "");
				items[1] = items[1].replace("Chr", "");
				items[1] = items[1].replace("CHR", "");

				# we need to loop through the exon regions twice. The first time to calculate the cds size
				# while the second time will be to print the combined results out to a file
				# for 'NR' genes, the cds start = cds end, so we will use transcript start and end!
				#
				if (cds_start == cds_end):
					cds_start = int(items[3])
					cds_end   = int(items[4])

				cds_length = 0

				for idx in range(0, int(items[7])):
					# exon_start ========= exon end
					#                     cds_start ---------- cds_end
					#
					if (int(ends[idx]) < cds_start):
						continue

					#               exon_start ========== exon_end
					# cds_start -------- cds_end
					#
					if (int(starts[idx]) > cds_end):
						continue

					# exon_start ======= exon_end
					#     cds_start -------- cds_end
					#
					if ( int(starts[idx]) < cds_start and int(ends[idx]) < cds_end):
						cds_length += int(ends[idx]) - cds_start
						continue

					#      exon_start ======== exon_end
					# cds_start --------- cds_end
					#
					if (cds_start < int(starts[idx]) and cds_end < int(ends[idx])):
						cds_length += cds_end - int(starts[idx])
						continue

					cds_length += int(ends[idx]) - int(starts[idx])

				# second loop will dump all the results out to an output file
				#
				for idx in range(0,int(items[7])):
					# take care the reverse strand issue
					exon_id = idx
					if (items[2] == "-"): 
						exon_id = int(items[7]) - idx - 1

					print("%s\t%d\t%d\t%s\t%s" % (items[1], int(starts[idx]), int(ends[idx]), items[0]+"_"+str(exon_id)+"_"+items[7]+"="+items[10]+"="+str(cds_start)+"="+str(cds_end)+"="+str(cds_length), items[10]))


if __name__ == "__main__":
	main(sys.argv[1:])

