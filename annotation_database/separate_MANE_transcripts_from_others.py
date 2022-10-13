#!/usr/bin/env python3

import argparse
import gzip

def parseArgument():
    parser = argparse.ArgumentParser()
    parser.add_argument('--MANE_file', required=True)
    parser.add_argument('--annotation_file', required=True)
    parser.add_argument('--type', required=True)

    return parser.parse_args()


def main():
    args = parseArgument()

    mainGeneTranscript  = {}
    with open(args.MANE_file, 'r') as fh:
        for line in fh:
            line = line.rstrip()

            (lid, gencode_gene, gene, refseq, gencode) = line.split('\t')
            (gencode_gene_name, tmp) = gencode_gene.split('.')
            (refseq_name, tmp) = refseq.split('.')
            (gencode_transcript, tmp) = gencode.split('.')

            if args.type == 'refseq':
                mainGeneTranscript[gene] = refseq_name
            else:
                mainGeneTranscript[gene] = gencode_transcript
    
    # here is the refseq info
    # 2085 NR_046630 chr3 + 19666 19669 19666 19666 3 19666,196667,19692, 19669,19613,19405, 0 NCBP2-AS1 unk unk -1,-1,-1,
    #

    # here is the gencode info
    # 585 ENST00000473358.1_4 chr1 + 29553 31097 29553 29553 3 29553,30563,30975, 30039,30667,31097, 0 MIR1302-2HG none none -1,-1,-1,

    outfile_MANE=''
    outfile_Non_MANE=''
    if (args.type == 'refseq'):
        outfile_MANE = 'MANE_refseq.txt'
        outfile_Non_MANE = 'non_MANE_refseq.txt'
    else:
        outfile_MANE = 'MANE_gencode.txt'
        outfile_Non_MANE = 'non_MANE_gencode.txt'
    
    fw = open(outfile_MANE, 'w')
    fq = open(outfile_Non_MANE, 'w')

    with gzip.open(args.annotation_file, 'rt') as hf:
        for line in hf:
            line = line.rstrip()

            info = line.split('\t')
            starts = info[9].split(',')
            ends   = info[10].split(',')
            transcript=''
            if (args.type == 'refseq'):
                transcript = info[1]
            else:
                (transcript, tmp) = info[1].split('.')

            for sid in range(int(info[8])):
                index = sid + 1

                # need to adjust the index if it is on the '-' strand
                #
                if info[3] == '-':
                    index = int(info[8]) - sid

                if info[12] in mainGeneTranscript.keys():
                    if (mainGeneTranscript[info[12]] == transcript):
                        fw.write(info[2] + '\t' + starts[sid] + '\t' + ends[sid] + '\t' + info[12] + '\t' + info[1] + '\t' + str(index) + '\n')
                    else:
                        fq.write(info[2] + '\t' + starts[sid] + '\t' + ends[sid] + '\t' + info[12] + '\t' + info[1] + '\t' + str(index) + '\n')
                else:
                    fq.write(info[2] + '\t' + starts[sid] + '\t' + ends[sid] + '\t' + info[12] + '\t' + info[1] + '\t' + str(index) + '\n')


if __name__=='__main__':
    main()
