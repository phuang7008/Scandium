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
#include <math.h>
#include <stdlib.h>
#include "stats.h"
#include "utils.h"

void readBufferInit(Read_Buffer *read_buff_in) {
	uint32_t i=0;
	for (i=0; i<read_buff_in->size; i++) {
        read_buff_in->chunk_of_reads[i] = bam_init1();
    }
}

void readBufferDestroy(Read_Buffer *read_buff_in) {
	uint32_t i=0;
	for (i=0; i<read_buff_in->size;i++) {
		if (read_buff_in->chunk_of_reads[i] != NULL) {
			bam_destroy1(read_buff_in->chunk_of_reads[i]);
			read_buff_in->chunk_of_reads[i]=NULL;
		}
		//} else {
		//	printf("Something is wrong with the bam_init1 destroy\n");
		//fprintf(stderr, "at position %d\n", i);
	}
	//read_buff_in->size = 0;

	//if (read_buff_in) free(read_buff_in);
}

uint32_t readBam(samFile *sfin, bam_hdr_t *header, Chromosome_Tracking *chrom_tracking, Read_Buffer *read_buff_in) {
    uint32_t record_idx = 0;
    while (record_idx < read_buff_in->size && chrom_tracking->more_to_read) {
        if (sam_read1(sfin, header, read_buff_in->chunk_of_reads[record_idx]) < 0) {
        	chrom_tracking->more_to_read = false;
			//fprintf(stderr, "Reading Bam has encountered some problem\n");
            break;
        }
        ++record_idx;
    }

	return record_idx;
}

void processBamChunk(User_Input *user_inputs, Coverage_Stats *cov_stats, khash_t(str) *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in, Target_Buffer_Status *target_buffer_status, int thread_id) {
	// it is the flag that is used to indicate if we need to add the khash into the coverage_hash
	bool not_added = true;	
	uint32_t i = 0;
	char cur_chr[50];
	strcpy(cur_chr, "NOTTHEREALONE");

	for (i=0; i<read_buff_in->size; i++) {

		if(user_inputs->percentage < 1.0) {
            //srand((uint32_t)time(NULL));	// set random seed and only need to be set ONCE
            float random_num = (float)rand() / (float)RAND_MAX;
            //if(random_num < user_inputs->percentage) continue;
            if(random_num > user_inputs->percentage) continue;
        }

		cov_stats->total_reads_produced++;

		// Need to check various 'READ' flags regarding the current read before doing statistics analysis
		if(read_buff_in->chunk_of_reads[i]->core.qual < user_inputs->min_map_quality) {
            continue;
        }

        if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FQCFAIL){           // Read Fails Vendor Quality Check
			continue;
        }

        if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FSECONDARY){        // Read Alignment is not Primary Alignment
            continue;
        }
 
        if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FUNMAP)	{			// Read Unmapped
             continue;
		}

		cov_stats->total_reads_aligned++;
        
		if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FPAIRED) {			// Read is properly paired
			cov_stats->total_reads_paired++;
            if(!(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FMUNMAP)) {	// Read is Paird with Mapped Mate 
				cov_stats->total_paired_reads_with_mapped_mates++;
            }
        }

		if(read_buff_in->chunk_of_reads[i]->core.flag & BAM_FDUP) {				// Read is a Duplicate (either optical or PCR)
			cov_stats->total_duplicate_reads++;
            if(user_inputs->remove_duplicate) continue;
        }

		if (read_buff_in->chunk_of_reads[i]->core.flag & BAM_FSUPPLEMENTARY) {
			cov_stats->total_supplementary_reads++;
			if (user_inputs->remove_supplementary_alignments)
				continue;
		}
 
		// check to see if we have changed the chromosome
		if ( strcmp(cur_chr, removeChr(header->target_name[read_buff_in->chunk_of_reads[i]->core.tid])) != 0) {

			// update the last_chr value, as we are going to handle new chromosome 
			strcpy(cur_chr, removeChr(header->target_name[read_buff_in->chunk_of_reads[i]->core.tid]));

			printf("Start processing chromosome id %s for thread %d\n", cur_chr, thread_id);

			not_added = true;
		}
		//if (strcmp(cur_chr, "17") != 0) return;		// for debugging
		//if (strcmp(cur_chr, "17") != 0) continue;		// for debugging
		//if (strstr(cur_chr, "GL000") == 0) continue;		// for debugging

		// Add the instance of khash_t(m32) to the khash_t(str) object (ie coverage_hash) based on the chrom_id
		if (not_added) {
			int absent;
			// Create an instance of khash_t(m32) and initialize it, Note each thread should have its own khash instance
			Temp_Coverage_Array *temp_cov_array = calloc(1, sizeof(Temp_Coverage_Array));
			temp_cov_array->cov_array = calloc(header->target_len[read_buff_in->chunk_of_reads[i]->core.tid]+1, sizeof(uint32_t));
			temp_cov_array->size = header->target_len[read_buff_in->chunk_of_reads[i]->core.tid]+1;

			khiter_t k_iter = kh_put(str, coverage_hash, strdup(cur_chr), &absent);    // get key iterator for chrom_id
			kh_value(coverage_hash, k_iter) = temp_cov_array;
			not_added = false;

			/////////////////Added on Feb 11, 2015/////////////////////////////
        	if (cov_stats->read_length <= 0)
				cov_stats->read_length = read_buff_in->chunk_of_reads[i]->core.l_qseq;
        	///////////////////////////////////////////////////////////////////
		}

        processRecord(user_inputs, cov_stats, coverage_hash, removeChr(header->target_name[read_buff_in->chunk_of_reads[i]->core.tid]), read_buff_in->chunk_of_reads[i], target_buffer_status);
    }

    printf("Done read bam for thread %d\n", thread_id);
}

void processRecord(User_Input *user_inputs, Coverage_Stats *cov_stats, khash_t(str) *coverage_hash, char * chrom_id, bam1_t *rec, Target_Buffer_Status *target_buffer_status) {
	uint32_t i=0, j=0, chrom_len=0;
	bool on_target=false, in_buffer=false;

	// get the target_buffer_status index,
	// need to do the following no matter target file or Ns regions specified or not
	// because we need chromosome length information here
	//
	int32_t idx = -1;
	for (i=0; i<target_buffer_status[0].num_of_chromosomes; i++) {
		if (strcmp(chrom_id, target_buffer_status[i].chrom_id) == 0) {
			chrom_len = target_buffer_status[i].size;
			idx = i;
			break;
		}
	}

	khiter_t outer_iter;
    int	ret=0;

	// get key iterator for chrom_id (outer hash table: khash_t(str))
	//
	outer_iter = kh_put(str, coverage_hash, chrom_id, &ret);	
	if (ret)
		fprintf(stderr, "Something went wrong, hash key for %s is not added\n", chrom_id);

	// get the quality information
	uint8_t *qual = bam_get_qual(rec);

	// Need to take care of soft-clip here as it will be part of the length, especially the soft-clip after a match string
	//
    uint32_t * cigar = bam_get_cigar(rec);		// get cigar info
	uint32_t start = rec->core.pos;
	uint32_t qual_pos = 0;

	for (i=0; i<rec->core.n_cigar; ++i) {
        int cop = cigar[i] & BAM_CIGAR_MASK;    // operation
        int cln = cigar[i] >> BAM_CIGAR_SHIFT;  // length
        if (cop == BAM_CMATCH) {

			for (j=0; j<cln; ++j) {
				cov_stats->total_aligned_bases++;
				uint32_t pos = start + qual_pos + 1;

				if (pos < 0) continue;
				if (pos >= chrom_len) break;

				if (TARGET_FILE_PROVIDED && (idx >=0)) {
					if (!on_target) {
						if (target_buffer_status[idx].status_array[pos] == 1 || target_buffer_status[idx].status_array[pos] == 4)
							on_target = true;
					}

					if (!in_buffer) {
						if (target_buffer_status[idx].status_array[pos] == 2 || target_buffer_status[idx].status_array[pos] == 5)
							in_buffer = true;
					}
				}

				if ( (user_inputs->min_base_quality) > 0 && (qual[qual_pos] < user_inputs->min_base_quality) ) {	
					qual_pos++;
					// filter based on the MIN_BASE_SCORE
					continue;
				}

				kh_value(coverage_hash, outer_iter)->cov_array[pos]++;
				qual_pos++;
			}
		} else if (cop == BAM_CREF_SKIP || cop == BAM_CDEL || cop == BAM_CPAD) {
			qual_pos += cln;
		}
	}

	if (TARGET_FILE_PROVIDED) {
		if (on_target) {
			cov_stats->on_target_read_hit_count += 1;
		} else if (in_buffer) {
			cov_stats->in_buffer_read_hit_count += 1;
		} else {
			cov_stats->off_target_read_hit_count += 1;
		}
	}
}

//Note: this function needs to be run unde the critical condition, that is only one thread will run this function
//
void combineThreadResults(Chromosome_Tracking *chrom_tracking, khash_t(str) *coverage_hash, bam_hdr_t *header) {
	khiter_t outer_iter;		// Note: khint_t and khiter_t are the same thing!
	int32_t i=0, j=0;

	/**
		First we need to loop through the coverage_hash to find out how many chromosomes it tracks
		As the keys in hash is not in order, I need to order them here!!!
		The chromosome id list in the header is ordered, so I am going to use it!
		Note: not all of the chromosomes will be presented in the alignments. 
		In this case, the bucket will be empty.  But I will still keep it as empty. 
		so everything we do something after this point, we will still have to check for the existence of the chromosome id
	*/

	for (i = 0; i < header->n_targets; i++) {
		// Now go through the coverage_hash map to see if it is currently tracking it
		// if it tracks new chromosome, we need to allocate memory for the coverage array of the new chromosome
		//
		if (chrom_tracking->chromosome_status[i] > 2) continue;

		for (outer_iter = kh_begin(coverage_hash); outer_iter != kh_end(coverage_hash); ++outer_iter) {
			if (kh_exist(coverage_hash, outer_iter)) {
				// compare with ordered chromosome ID in the header, if it is matched, then process it!
				//
				if (strcmp(removeChr(header->target_name[i]), kh_key(coverage_hash, outer_iter)) == 0) {
					bool need_to_allocate = true;

					if (chrom_tracking->chromosome_lengths[i] > 0)
						// it is already tracked and memory allocated!
						need_to_allocate = false;

		            // Process the current chromosome for the first time!
					if (need_to_allocate) {
						uint32_t chrom_len = header->target_len[i];
						chromosomeTrackingUpdate(chrom_tracking, removeChr(header->target_name[i]), chrom_len, i);

						// if it goes to the next chromosome, the previous chromosome should be done processing
				        // Thus, we need to update the previous chromosome status
						// Here we will have to make sure that j is signed. Otherwise, you won't get negative value
						for(j=i-1; j >= 0; j--) {
							if (chrom_tracking->chromosome_status[j] == 1)
								chrom_tracking->chromosome_status[j] = 2;
						}

					}

					// update the coverage information for each position at current chromosome
					//
					for(j=0; j<kh_value(coverage_hash, outer_iter)->size; j++) {
						chrom_tracking->coverage[i][j] += kh_value(coverage_hash, outer_iter)->cov_array[j];
					}
				}
			}
		}
	}
}
