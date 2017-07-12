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

#include <omp.h>
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

void process_chunk_of_bam(int thread_id, khash_t(str) *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in) {
	int not_added = true;	// it is a flag that is used to indicate if we have added the khash into the coverage_hash
	int i = 0;
	char *cur_chr  = malloc(15 * sizeof(char));
	char *last_chr = malloc(30 * sizeof(char));
	strcpy(last_chr, "SOMEVERYFAKETHINGGOESHERE");

	for (i=0; i<read_buff_in->size; i++) {
		if(_PERCENTAGE < 1.0) {
            srand((unsigned int)time(NULL));	// set random seed
            float random_num = (float)rand() / (float)RAND_MAX;
            if(random_num < _PERCENTAGE) continue;
        }

		TOTAL_READS_PRODUCED++;

        // get core alignment info
        bam1_core_t *core = &read_buff_in->chunk_of_reads[i]->core;

		// Need to check various 'READ' flags regarding the current read before doing statistics analysis
		if(core->qual < MIN_MAP_SCORE) {
            continue;
        }

        if(core->flag & BAM_FQCFAIL){           /* Read Fails Vendor Quality Check*/
			continue;
        }
        if(core->flag & BAM_FSECONDARY){        /* Read Alignment is not Primary Alignment*/
            continue;
        }
 
        if(core->flag & BAM_FUNMAP)	{			/* Read Unmapped*/
            //printf("Unmapped! region.  breaking!!!\n");
             continue;
		}

		TOTAL_READS_ALIGNED++;
        
		if(core->flag & BAM_FPROPER_PAIR) {		/* Read is properly paired */
			TOTAL_READS_PAIRED++;
            if(!(core->flag & BAM_FMUNMAP)) {	/* Read is Paird with Mapped Mate */
                TOTAL_PAIRED_READS_WITH_MAPPED_MATES++;
            }
        }

		if(core->flag & BAM_FDUP) {				/* Read is a Duplicate (either optical or PCR) */
            TOTAL_DUPLICATE_READS++;
            if(REMOVE_DUPLICATES){continue;}
        }

		if (core->flag & BAM_FSUPPLEMENTARY) {
			TOTAL_SUPPLEMENTARY_ALIGNMENTS++;
			if (REMOVE_SUPPLEMENTARY_ALIGNMENTS)
				continue;
		}
 
		// check to see if we have changed the chromosome
		strcpy(cur_chr, header->target_name[core->tid]);
		if ( strcmp(last_chr, cur_chr) != 0) {

			// update the last_chr value, as we are going to handle new chromosome 
			strcpy(last_chr, cur_chr);

			printf("Start processing chromosome id %s\n", cur_chr);

			not_added = true;
		}

		// Add the instance of khash_t(m32) to the khash_t(str) object (ie coverage_hash) based on the chrom_id
		if (not_added) {
			int absent;
			// Create an instance of khash_t(m32) and initialize it, Note each thread should have its own khash instance
			khash_t(m32) *tmp_hash = kh_init(m32);
			khiter_t k_iter = kh_put(str, coverage_hash, cur_chr, &absent);    // get key iterator for chrom_id
			if (absent) 
				kh_key(coverage_hash, k_iter) = strdup(cur_chr); 
			kh_value(coverage_hash, k_iter) = tmp_hash;
			not_added = false;

			/////////////////Added on Feb 11, 2015/////////////////////////////
        	if (READLENGTH <= 0)
				READLENGTH = read_buff_in->chunk_of_reads[i]->core.l_qseq;
        	///////////////////////////////////////////////////////////////////
		}

        processRecord(coverage_hash, header->target_name[core->tid], read_buff_in->chunk_of_reads[i]);
    }

	free(cur_chr);
	free(last_chr);
    printf("Done read bam\n");
}

void processRecord(khash_t(str) *coverage_hash, char * chrom_id, bam1_t *rec) {
	int i, ret;
	uint32_t start = rec->core.pos + 1;		// left most position of alignment in zero based coordinates (as BAM is 0-based)
	khiter_t outer_iter, inner_iter;

	// get key iterator for chrom_id (outer hash table: khash_t(str))
	outer_iter = kh_put(str, coverage_hash, chrom_id, &ret);	
	//if (ret)
	//	fprintf(stderr, "Something went wrong, has key for %s is not added", chrom_id);

	for (i=0; i<rec->core.l_qseq; i++) {
		TOTAL_ALIGNED_BASES++;

		int pos = start + i;
		inner_iter = kh_put(m32, kh_value(coverage_hash, outer_iter), pos, &ret);			// iterator for inner hash table khash_t(m32)

		if (ret == -1) {	// failed!
			fprintf(stderr, "Can't insert the pos key into the hash table at pos %"PRIu32"\n", pos);
			exit(1);
		}

		if (ret == 1)		// 1 if the bucket is empty (never used)
			// initialize the value to 0
			kh_value(kh_value(coverage_hash, outer_iter), inner_iter) = 0;
		
		unsigned short val = kh_value(kh_value(coverage_hash, outer_iter), inner_iter);		// fetch the current value for 'pos' key
		kh_value(kh_value(coverage_hash, outer_iter), inner_iter) = val + 1;				// increment value by 1

		//if (pos == 8390000) {
		//	printf("Before tracking position 8390000 value is %d\n", val);
		//	printf("After tracking position 8390000 value is %d\n", kh_value(kh_value(coverage_hash, outer_iter), inner_iter));
		//}
	}
}

//Note: this function needs to be run unde the critical condition, that is only one thread will run this function
//
void combine_thread_results(Chromosome_Tracking *chrom_tracking, khash_t(str) *coverage_hash, bam_hdr_t *header) {
	khiter_t outer_iter, inner_iter;		// Note: khint_t and khiter_t are the same thing!
	char *chrom_id = malloc(15 * sizeof(char));
	int i;

	// First we need to loop through the coverage_hash to find out how many chromosomes it tracks
	// if it tracks new chromosome, we need to allocate memory for the coverage array of the new chromosome
	//
	for (outer_iter = kh_begin(coverage_hash); outer_iter != kh_end(coverage_hash); ++outer_iter) {
		if (kh_exist(coverage_hash, outer_iter)) {
			strcpy(chrom_id, kh_key(coverage_hash, outer_iter));
			bool need_to_allocate = true;

			// check to see if the memories have been allocated for the chrom_id
			for(i=0; i<chrom_tracking->number_tracked; i++) {
				if (strcmp(chrom_tracking->chromosome_ids[i], chrom_id) == 0) {
					// it is already tracking the indicated chromosome 
					need_to_allocate = false;
					break;
				} else {
					// if it goes to the next chromosome, the previous chromosome should be done processing
					// Thus, we need to update the previous chromosome status
					chrom_tracking->chromosome_status[i] = 2;
                	printf("The previous chromosome is done %s\n", chrom_id);
				}
			}

			printf("processing chromosome id %s\n", chrom_id);

            // Process the current chromosome for the first time!
			if (need_to_allocate) {
				uint32_t chrom_len = header->target_len[get_chrom_index_from_id(header, chrom_id)];
				chromosome_tracking_init(chrom_tracking, chrom_id, chrom_len, chrom_tracking->number_tracked);
			}
			printf("After memory allocation\n");
		}
	}

	//printf("chr\tposition\tcounts\n");

	// start with the current thread hash table, and find what chromosome ids they are working at
	// by looking through the outer key!
	//
	for (outer_iter = kh_begin(coverage_hash); outer_iter != kh_end(coverage_hash); ++outer_iter) {
		printf("Inside the outer key\n");
		if (kh_exist(coverage_hash, outer_iter)) {
			strcpy(chrom_id, kh_key(coverage_hash, outer_iter));
			printf("At combining stage: Current chromosome ID is %s\n", chrom_id);

			// locate the index for the chromosome the current thread (hash) tracks!
			//
			int i, idx;
			for (i=0; i<chrom_tracking->number_tracked; i++) {
				printf("chrom_tracking->chromosome_ids[i] is %s\n", chrom_tracking->chromosome_ids[i]);
				if (strcmp(chrom_id, chrom_tracking->chromosome_ids[i]) == 0) {
					idx = i;
					break;
				}
			}
			printf("The tracking index is %d\n", idx);

			// update the coverage information for each position at current chromosome
			//
			for (inner_iter=kh_begin(kh_value(coverage_hash, outer_iter));
					inner_iter!=kh_end(kh_value(coverage_hash, outer_iter));
						inner_iter++) {
				if (kh_exist(kh_value(coverage_hash, outer_iter), inner_iter)) {
					unsigned short val = kh_value(kh_value(coverage_hash, outer_iter), inner_iter);
					uint32_t pos = kh_key(kh_value(coverage_hash, outer_iter), inner_iter);
					chrom_tracking->coverage[idx][pos] = chrom_tracking->coverage[idx][pos] + val;

					//printf("inside the inner loop with val %d, and pos %"PRIu32"\n", val, pos);
					//if ((key % 1000000) == 0) {
					//	printf("%s\t %"PRIu32"\t %d\n",chrom_tracking->chromosome_ids[idx], key, chrom_tracking->coverage[idx][key]);
					//}
				}
			}
		}
	}
	free(chrom_id);
}
