/*
 * =====================================================================================
 *
 *       Filename:  stats.c
 *
 *    Description:  the detailed implementation of sequencing statistics
 *
 *        Version:  1.0
 *        Created:  02/22/2017 01:55:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang (phuang@bcm.edu)
 *        Company:  Baylor College of medicine
 *
 * =====================================================================================
 */

#include <time.h>
#include <stdio.h>
#include "stats.h"
#include "utils.h"

void read_buff_init(Read_Buffer *read_buff_in) {
	int i;
	for (i=0; i<read_buff_in->size; i++) {
        read_buff_in->chunk_of_reads[i] = bam_init1();
    }
}

void read_buff_destroy(Read_Buffer *read_buff_in) {
	int i;
	for (i=0; i<read_buff_in->size;i++) {
		bam_destroy1(read_buff_in->chunk_of_reads[i]);
	}
}

uint32_t read_bam(samFile *sfin, bam_hdr_t *header, bool more_to_read, Read_Buffer *read_buff_in) {
    uint32_t record_idx = 0;
    while (record_idx < read_buff_in->size && more_to_read) {
        if (sam_read1(sfin, header, read_buff_in->chunk_of_reads[record_idx]) < 0) {
        	//more_to_read = false;
            break;
        }
        ++record_idx;
    }

	return record_idx;
}

void process_chunk_of_bam(int thread_id, Chromosome_Tracking *chrom_tracking, khash_t(32)  *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in) {
	uint32_t i, current_chrom_size;

	for (i=0; i<read_buff_in->size; i++) {
		// we need to check for various flags here before conduction statistics analysis
		if(_PERCENTAGE < 1.0)
        {
            srand((unsigned int)time(NULL));        // set random seed
            float random_num = (float)rand() / (float)RAND_MAX;
            if(random_num < _PERCENTAGE) continue;
        }

		TOTAL_READS_PRODUCED++;

        // get core alignment info
        bam1_core_t *core = &read_buff_in->chunk_of_reads[i]->core;

		if(core->qual < MIN_MAP_SCORE)
        {
            continue;
        }

        if(core->flag & BAM_FQCFAIL){           /* Read Fails Vendor Quality Check*/
			continue;
        }
        if(core->flag & BAM_FSECONDARY){        /* Read Alignment is not Primary Alignment*/
            continue;
        }
 
        if(core->flag & BAM_FUNMAP)                     /* Read Unmapped*/
        {
            //printf("Unmapped! region.  breaking!!!\n");
             continue;
		}

		TOTAL_READS_ALIGNED++;
        
		if(core->flag & BAM_FPROPER_PAIR)               /* Read is properly paired */
        {
			TOTAL_READS_PAIRED++;
            if(!(core->flag & BAM_FMUNMAP)){   /* Read is Paird with Mapped Mate */
                TOTAL_PAIRED_READS_WITH_MAPPED_MATES++;
            }
        }

		if(core->flag & BAM_FDUP)               /* Read is a Duplicate (either optical or PCR) */
        {
            DUPLICATE_READS++;
            if(REMOVE_DUPLICATES){continue;}
        }
 
        /////////////////Added on Feb 11, 2015/////////////////////////////
        //READLENGTH = read_buff_in[i].getReadLength();
        ///////////////////////////////////////////////////////////////////

		strcpy(chrom_tracking->curr_chromosome, header->target_name[core->tid]);
        if (strcmp(chrom_tracking->curr_chromosome, chrom_tracking->prev_chromosome) != 0) {	// 0 means two stings match
            if (strcmp(chrom_tracking->prev_chromosome, "nothing") == 0) {
				// update the prev_chromosome value to current one
				strcpy(chrom_tracking->prev_chromosome, chrom_tracking->curr_chromosome);

				// now need to initialize the COVERAGE array for current chromosome if it hasn't been initialized by other threads
				printf("Chromosome %s\n", chrom_tracking->curr_chromosome);
            	current_chrom_size = header->target_len[core->tid];         /* get the chromosome sequence length */
            	//COVERAGE = new short[size];
            }
		} else {
			printf("They are the same %s and %s in thread %d\n", chrom_tracking->curr_chromosome, chrom_tracking->prev_chromosome, thread_id);
		}
		
        processRecord(coverage_hash, read_buff_in->chunk_of_reads[i]);
    }

    printf("Done read bam\n");
    //getTargetsAndWriteCoverage(prev_chr,  COVERAGE, covFasta,  missTraget,  wgCoverage);
    //findWhereReadsHit(prev_chr,  COVERAGE,wig);
    //COVERAGE = null;
}

void processRecord(khash_t(32) *coverage_hash, bam1_t *rec) {
	bam1_core_t *core = &rec->core;
	int i, ret;
	long start = rec->core.pos + 1;		// left most position of alignment in zero based coordinates (as BAM is 0-based)
	khiter_t k_iter;

	for (i=0; i<core->l_qseq; i++) {
		TOTAL_ALIGNED_BASES++;

		int pos = start + i;
		k_iter = kh_put(32, coverage_hash, pos, &ret);		// set/get the key iter
		uint32_t val = kh_value(coverage_hash, k_iter);		// get the value
		kh_value(coverage_hash, k_iter) = val + 1;			// set the value
	}
}

void combine_thread_results(short *one_chromosome_array, khash_t(32) *coverage_hash) {
	int ret;
	khint_t k;
	khiter_t k_iter;
	for (k = kh_begin(coverage_hash); k != kh_end(coverage_hash); ++k) {
		k_iter = kh_put(32, coverage_hash, k, &ret);				// get the key iter
		one_chromosome_array[k] = kh_value(coverage_hash, k_iter);	// get the value based on the key iter
	}
}
