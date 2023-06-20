#!/usr/bin/python

import argparse
import sys
import mysql.connector


def main(argv):
    parser = argparse.ArgumentParser(description='From Exons to CDSs')
    parser.add_argument('-i', '--input', help='The input file [MANDATORY]', required=True)
    parser.add_argument('-t', '--type',  help='The type of annotation (refseq, ccds or gencode etc) [MANDATORY]', required=True)
    parser.add_argument('-v', '--version', help='The version of reference genome [MANDATORY]', required=True)
    args = parser.parse_args()
	
    processFile(args)

def processFile(args):

    con = ''
    hgnc = ''
    if (args.type == 'ccds' or args.type == 'CCDS'):
        if args.version == 'hg38':
            hgnc = 'HGNC38'
        else:
            hgnc = 'HGNC37'
        con = mysql.connector.connect(user='phuang468', password='phuang468', host='sug-esxa-db1', database='GeneAnnotations')

    with open(args.input, 'r') as rfh:
        for line in rfh:
            if "chrom" not in line:
                # the order of item list is: name, chrom_id, strand, gene_start, gene_end, exon_count, exon_starts, exon_ends, gene_symbol
                items  = line.rstrip("\n").split()
                if 'XM_' in items[0]:
                    continue

                cds_start = int(items[5])
                cds_end   = int(items[6])
                starts = items[8].split(",")
                ends   = items[9].split(",")
				
                # some version of gene annotation has 'chr' in front of chromosome id, so let's remove them
                if (args.version != 'hg38'):
                    items[1] = items[1].replace("chr", "");
                    items[1] = items[1].replace("Chr", "");
                    items[1] = items[1].replace("CHR", "");

                # for CCDS, we don't have gene symbol available, we need to query HGNC table to get it!
                #
                if (args.type == 'ccds' or args.type == 'CCDS'):
                    c = con.cursor()
                    query  = 'SELECT symbol FROM ' + hgnc + ' WHERE (ccds_id IS NOT NULL AND find_in_set('
                    query += items[0] +', replace(ccds_id, '|', ',')))'
                    c.execute(query)
                    for row in c:
                        items[10] = row[0]

                # we need to loop through the exon regions twice. The first time to calculate the cds size
                # while the second time will be to print the combined results out to a file
                # for 'NR' genes, the cds start = cds end, so we will use transcript start and end!
                # we should skip this one 
                #
                if (cds_start == cds_end):
                    continue;
                    #cds_start = int(items[3])
                    #cds_end   = int(items[4])

                cds_length = 0
                excluded_exon_count = 0

                for idx in range(0, int(items[7])):
                    # exon_start ========= exon end
                    #                     cds_start ---------- cds_end
                    #
                    if (int(ends[idx]) < cds_start):
                        starts[idx] = 0;
                        ends[idx] = 0;
                        excluded_exon_count += 1
                        continue

                    #               exon_start ========== exon_end
                    # cds_start -------- cds_end
                    #
                    if (int(starts[idx]) > cds_end):
                        starts[idx] = 0;
                        ends[idx] = 0;
                        excluded_exon_count += 1
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
                    if starts[idx] == 0:
                        continue

                    if (items[2] == "-"): 
                        exon_id = int(items[7]) - idx - 1

                    exon_count = int(items[7]) - excluded_exon_count
                    print("%s\t%d\t%d\t%s\t%s" % (items[1], int(starts[idx]), int(ends[idx]), items[0]+"_"+str(exon_id)+"_"+str(exon_count)+"="+items[10]+"="+str(cds_start)+"="+str(cds_end)+"="+str(cds_length), items[10]))


if __name__ == "__main__":
    main(sys.argv[1:])

