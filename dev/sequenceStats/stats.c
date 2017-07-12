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
		if (read_buff_in->chunk_of_reads[i]) {
			bam_destroy1(read_buff_in->chunk_of_reads[i]);
		} else {
			printf("Something is wrong with the bam_init1 destroy\n");
		}
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
            srand((uint32_t)time(NULL));	// set random seed
            float random_num = (float)rand() / (float)RAND_MAX;
            if(random_num < user_inputs->percentage) continue;
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
		if ( strcmp(cur_chr, header->target_name[read_buff_in->chunk_of_reads[i]->core.tid]) != 0) {

			// update the last_chr value, as we are going to handle new chromosome 
			strcpy(cur_chr, header->target_name[read_buff_in->chunk_of_reads[i]->core.tid]);

			printf("Start processing chromosome id %s for thread %d\n", cur_chr, thread_id);

			not_added = true;
		}
		//if (strcmp(cur_chr, "1") == 0) return;		// for debugging
		//if (strcmp(cur_chr, "1") == 0) continue;		// for debugging

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

        processRecord(user_inputs, cov_stats, coverage_hash, header->target_name[read_buff_in->chunk_of_reads[i]->core.tid], read_buff_in->chunk_of_reads[i], target_buffer_status);
    }

    printf("Done read bam for thread %d\n", thread_id);
}

void processRecord(User_Input *user_inputs, Coverage_Stats *cov_stats, khash_t(str) *coverage_hash, char * chrom_id, bam1_t *rec, Target_Buffer_Status *target_buffer_status) {
	uint32_t i=0, j=0, chrom_len=0;
	bool on_target=false, in_buffer=false;

	// get the target_buffer_status index
	uint32_t idx = 0;
	if (TARGET_FILE_PROVIDED) {
		for (i=0; i<35; i++) {
			if (strcmp(chrom_id, target_buffer_status[i].chrom_id) == 0) {
				chrom_len = target_buffer_status[i].size;
				idx = i;
				break;
			}
		}
	}

	khiter_t outer_iter;
    int	ret=0;

	// get key iterator for chrom_id (outer hash table: khash_t(str))
	outer_iter = kh_put(str, coverage_hash, chrom_id, &ret);	
	if (ret)
		fprintf(stderr, "Something went wrong, hash key for %s is not added\n", chrom_id);

	// get the quality information
	uint8_t *qual = bam_get_qual(rec);

	// Need to take care of soft-clip here as it will be part of the length, especially the soft-clip after a match string
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

				if (TARGET_FILE_PROVIDED) {
					if (!on_target) {
						if (target_buffer_status[idx].status_array[pos] == 1)
							on_target = true;
					}

					if (!in_buffer) {
						if (target_buffer_status[idx].status_array[pos] == 2)
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
			//pos += cln;
		} else if (cop == BAM_CREF_SKIP || cop == BAM_CDEL || cop == BAM_CPAD) {
			qual_pos += cln;
			//pos += cln;
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
		if (chrom_tracking->chromosome_status[i] == 3) continue;

		for (outer_iter = kh_begin(coverage_hash); outer_iter != kh_end(coverage_hash); ++outer_iter) {
			if (kh_exist(coverage_hash, outer_iter)) {
				// compare with ordered chromosome ID in the header, if it is matched, then process it!
				//
				if (strcmp(header->target_name[i], kh_key(coverage_hash, outer_iter)) == 0) {
					bool need_to_allocate = true;

					if (chrom_tracking->chromosome_lengths[i] > 0)
						// it is already tracked and memory allocated!
						need_to_allocate = false;

		            // Process the current chromosome for the first time!
					if (need_to_allocate) {
						uint32_t chrom_len = header->target_len[i];
						chromosomeTrackingUpdate(chrom_tracking, header->target_name[i], chrom_len, i);

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

void writeCoverage(char * chrom_id, Bed_Info *Ns_bed_info, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, MYSQL *con) {

    // First, we need to find the index that is used to track current chromosome chrom_id
    uint32_t idx = locateChromosomeIndex(chrom_id, chrom_tracking);

    // for the whole genome, we need to use the file that contains regions of all Ns in the reference
    // As they will be not used, so we are going to set the count info in these regions to 0
    if (N_FILE_PROVIDED)
        zeroAllNsRegions(chrom_id, Ns_bed_info, chrom_tracking);

    // write to the file that contains whole genome info
	uint32_t i=0;
    if(user_inputs->wgs_coverage) {

		if (user_inputs->Write_WGS) {
			FILE *wgs_fp = fopen(user_inputs->wgs_file, "a");
			printf("Whole Genome output for chrom id %s is on\n", chrom_id);

			// no need to add newline as the next line will take care of it!
			fprintf(wgs_fp, ">%s 1 %"PRIu32, chrom_id, chrom_tracking->chromosome_lengths[idx]-1);

	        for (i = 0; i < chrom_tracking->chromosome_lengths[idx]; i++) {
				if(i%100==0) fputc('\n', wgs_fp);
				fprintf(wgs_fp, "%d ", chrom_tracking->coverage[idx][i]);

				addBaseStats(stats_info, (uint32_t) chrom_tracking->coverage[idx][i], 0, 1);
			}

			fputc('\n', wgs_fp);
			fclose(wgs_fp);
		} else {
			for (i = 0; i < chrom_tracking->chromosome_lengths[idx]; i++)
				addBaseStats(stats_info, (uint32_t) chrom_tracking->coverage[idx][i], 0, 1);
		}

		// now need to report those regions with low or too high coverages
		// NOTE: the bed format is different here, the end position is included!
		//
		FILE *wgs_low_x_fp  = fopen(user_inputs->wgs_low_cov_file, "a");
		FILE *wgs_high_x_fp = fopen(user_inputs->wgs_high_cov_file, "a");

		//fprintf(wgs_low_x_fp,  "Chr\tStart\tStop\tLength\tCoverage\n");
		//fprintf(wgs_high_x_fp, "Chr\tStart\tStop\tLength\tCoverage\n");
		writeLow_HighCoverageReport(0, chrom_tracking->chromosome_lengths[idx], chrom_tracking, idx, user_inputs, wgs_low_x_fp, wgs_high_x_fp, con);

		fclose(wgs_low_x_fp);
		fclose(wgs_high_x_fp);
    }

	// if the target bed file is available, we will need to handle it here and write the results to cov.fasta file
	if (TARGET_FILE_PROVIDED) {
		FILE * cov_fp = fopen(user_inputs->cov_file, "a");
		FILE * missed_target_fp = fopen(user_inputs->missed_targets_file, "a");

		for(i = 0; i < target_info->size; i++) {
            if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
                continue;

            stats_info->cov_stats->total_targets++;
            uint32_t start = target_info->coords[i].start;
            uint32_t end = target_info->coords[i].end;
            int length = end - start + 1;
            bool collect_target_cov = length > 99 ? true : false ;  // TODO: is it the min length of the exon? Check.
			//if (collect_target_cov) 
			//	fprintf(stderr, "%s %"PRIu32" %"PRIu32" %d\n", chrom_id, start, end, length);

            uint32_t j=0;
            if (collect_target_cov) {
                for(j = 0; j < PRIMER_SIZE; j++) {
                    if ( (start - j) < 0 || (end + j) >= chrom_tracking->chromosome_lengths[idx])
                        continue;

					stats_info->five_prime[j]  += chrom_tracking->coverage[idx][start-j];
					stats_info->three_prime[j] += chrom_tracking->coverage[idx][end+j];
                }
			}

			//In the original java code, the following is used for the pc and pc2 definition
            //short pc[101], pc2[101];
			//After discussion with Divya and Qiaoyan, here is what we understand.
			//pc and pc2 is used to setup the 'target_coverage' variable which is defined as an array of 101 in the origain java code
			//here I will use percentage_bin and percentage_count_per_bin
			uint32_t percentage_bin[101], percentage_count_per_bin[101];
			for(j=0; j<101; j++) {
				percentage_bin[j] = 0;
				percentage_count_per_bin[j] = 0;
			}

            bool target_hit = false;
            fprintf(cov_fp, ">%s %"PRIu32" %"PRIu32"\n", chrom_id, start, end);

            bool space_it = false;
            if(end - start > 10000) space_it = true;

            for(j = 0; j < length; j++) {
                if (j+start >= chrom_tracking->chromosome_lengths[idx])
                    continue;

                if (space_it && j%100 == 0) fputc('\n', cov_fp);    // enter a new line after every 100 bases

				uint32_t cov = chrom_tracking->coverage[idx][j+start];
				addBaseStats(stats_info, cov, 1, 0);

				if (!target_hit && (cov > 0))
					target_hit = true;

                // output to the cov.fasta file
                fprintf(cov_fp, "%d ", cov);

                if (collect_target_cov) {
					float num, den;
					memcpy(&num, &j, sizeof(num));
					memcpy(&den, &length, sizeof(den));
                    int percentage_pos = (int)((num*100)/den + 0.5);
                    //int percentage_pos = lround((num/den)*100);
                    percentage_bin[percentage_pos] += cov;
                    percentage_count_per_bin[percentage_pos] += 1;
					//fprintf(stderr, "%"PRIu32" %d %d\n", j, length, percentage_pos);
                }
            }

            // output a newline char to the cov.fasta file 
            fputc('\n', cov_fp);

			if (collect_target_cov) {
				for (j = 0; j < 101; j++) {
					if(percentage_count_per_bin[j] != 0) {
						float num, den;
				        memcpy(&num, &percentage_bin[j], sizeof(num));
					    memcpy(&den, &percentage_count_per_bin[j], sizeof(den));
						int d = (int) (num/den+0.5);
	                    //int d = lround(num/den);
		                percentage_bin[j] = d;	
					}
                }

				for(j = 0; j < 101; j++)
					stats_info->target_coverage[j] += percentage_bin[j];
            }

			if (target_hit) {
				stats_info->cov_stats->hit_target_count += 1;
			} else {
				// need to write to the missed target file
				fprintf(missed_target_fp, "%s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, start, end);
				bool hit = false;
				for (j = start - BUFFER; j < start && !hit; j++) {
					if (j >= chrom_tracking->chromosome_lengths[idx])
						continue;

					if ( chrom_tracking->coverage[idx][j] > 0 )
						hit = true;
				}

				if (hit)
					stats_info->cov_stats->hit_target_buffer_only_count += 1;
			}

			// now need to report those Capture regions with low or too high coverages
	        // NOTE: the bed format is different here, the end position is included!
    	    //
        	FILE *capture_low_x_fp  = fopen(user_inputs->capture_low_cov_file, "a");
	        FILE *capture_high_x_fp = fopen(user_inputs->capture_high_cov_file, "a");

    	    //fprintf(capture_low_x_fp,  "Chr\tStart\tStop\tLength\tCoverage\n");
        	//fprintf(capture_high_x_fp, "Chr\tStart\tStop\tLength\tCoverage\n");
			writeLow_HighCoverageReport(start, length, chrom_tracking, idx, user_inputs, capture_low_x_fp, capture_high_x_fp, con);

	        fclose(capture_low_x_fp);
    	    fclose(capture_high_x_fp);
		}
		fclose(cov_fp);
		fclose(missed_target_fp);
	}
}

uint32_t writeLow_HighCoverageReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, uint16_t chrom_idx, User_Input *user_inputs,FILE *fh_low, FILE *fh_high, MYSQL *con) {
	uint32_t i=0;
	// for debugging
	//if (strcmp(chrom_tracking->chromosome_ids[chrom_idx], "1") == 0) return i;

	for (i = begin; i < begin+length; i++) {
        uint32_t start=0, end=0;
        uint64_t cov_total=0;

        // for low coverage
        if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] < user_inputs->low_coverage_to_report) {
            start = i;

            while(i < begin+length && chrom_tracking->coverage[chrom_idx][i] < user_inputs->low_coverage_to_report) {
                cov_total += chrom_tracking->coverage[chrom_idx][i];
                i++;
            }
            end = i;
			if (start == end) continue;

            //float ave_coverage = (float)cov_total / (float)(end - start);
            uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);

			// generate the annotation here
			if (user_inputs->annotation_on) {
				char *annotation = produceGeneAnnotations(start, end-1, chrom_tracking->chromosome_ids[chrom_idx], con);
				fprintf(fh_low, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\t%s\n", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage, annotation);
				if (annotation) free(annotation);
			} else {
		        fprintf(fh_low, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage);
			}
        }
		//fflush(fh_low);

        // for high coverage
		start = 0;
		end = 0;
		cov_total = 0;
        if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] >= user_inputs->high_coverage_to_report) {
            start = i;

            while( i < begin+length && chrom_tracking->coverage[chrom_idx][i] >= user_inputs->high_coverage_to_report) {
				if ( (user_inputs->upper_bound_to_report == -1 ) || 
						(user_inputs->upper_bound_to_report > -1 && chrom_tracking->coverage[chrom_idx][i] < user_inputs->upper_bound_to_report)) {
					cov_total += chrom_tracking->coverage[chrom_idx][i];
					i++;
				} else {
					break;
				}
            }
            end = i;
			if (start == end) continue;

            //float ave_coverage = (float)cov_total / (float)(end - start);
            uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);

            //fprintf(fh_low, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%.2f\n", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage);
            fprintf(fh_high, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage);
        }
		//fflush(fh_high);
    }

	return i;
}

void addBaseStats(Stats_Info *stats_info, uint32_t cov_val, uint8_t target, uint8_t wgs) {
	// need to check and update max coverage information
	if (cov_val > stats_info->cov_stats->max_coverage) {
		stats_info->cov_stats->max_coverage = cov_val;
		stats_info->cov_stats->base_with_max_coverage = 1;
	} else if (cov_val == stats_info->cov_stats->max_coverage) {
		stats_info->cov_stats->base_with_max_coverage += 1;
	}

	if (cov_val < 0) {
        fprintf(stderr, "Coverage less than 0!!!!!!!\n");
    }

	// for histogram only
	uint32_t tmp_val = cov_val;
	if (tmp_val > 1000) tmp_val = 1000;
	if (target == 1) addValueToKhashBucket32(stats_info->target_cov_histogram, tmp_val, 1);
	if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_cov_histogram, tmp_val, 1);

    if (target == 1) stats_info->cov_stats->total_target_coverage += (uint64_t) cov_val;
    if (wgs == 1)    stats_info->cov_stats->total_genome_coverage += (uint64_t) cov_val;

    if (cov_val > 0) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 1, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 1, 1);
    } else if (cov_val == 0) {
		if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 0, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 0, 1);
	}

	if (cov_val >= 5) {
		if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 5, 1);
		if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 5, 1);
	}

	if (cov_val >= 10) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 10, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 10, 1);
    }

	if (cov_val >= 15) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 15, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 15, 1);
    }

	if (cov_val >= 20) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 20, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 20, 1);
    }

	if (cov_val >= 30) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 30, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 30, 1);
    }

	if (cov_val >= 40) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 40, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 40, 1);
    }

	if (cov_val >= 50) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 50, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 50, 1);
    }

	if (cov_val >= 60) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 60, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 60, 1);
    }

	if (cov_val >= 100) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 100, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 100, 1);
    }

	if (cov_val >= 500) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 500, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 500, 1);
    }

	if (cov_val >= 1000) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 1000, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 1000, 1);
    }

	if (target == 1) addValueToKhashBucket32(stats_info->target_coverage_for_median, cov_val, 1);
	if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_coverage_for_median, cov_val, 1);
}

void writeReport(Stats_Info *stats_info, User_Input *user_inputs) {
	if (TARGET_FILE_PROVIDED && stats_info->cov_stats->total_targeted_bases == 0) {
		fprintf(stderr, "Total targeted bases is zero.  This means that no read has aligned to a chromosome that contains a target.");
		fprintf(stderr, "No target matches a chromosome in the BAM, or something else went wrong.  Aborting.\n");
		exit(1);
	}

	if (stats_info->cov_stats->total_reads_aligned == 0) {
		fprintf(stderr, "No reads aligned. Aborting.\n");
		exit(1);
	}

	uint32_t non_duplicate_reads = stats_info->cov_stats->total_reads_aligned - stats_info->cov_stats->total_duplicate_reads;
	if (non_duplicate_reads == 0) {
        fprintf(stderr, "All reads are duplicates. Aborting.\n");
        exit(1);
    }

	if (TARGET_FILE_PROVIDED && stats_info->cov_stats->total_targets == 0) {
        //I don't think we should ever see this error, as its dealt with above.
        fprintf(stderr, "No target regions given.  Aborting.\n");
        exit(1);
    }

	//printf("Before Median calculation\n");
   	uint64_t sum=0;
	int32_t i=0;
	double average_coverage;
	int ret;
	uint16_t bins[13] = { 0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 100, 500, 1000 };
	khiter_t k_iter;
	char *file_name;

	if (user_inputs->wgs_coverage) {
		//Do not consider the Ns for Median calculation.
		uint64_t total_genome_non_Ns_bases = stats_info->cov_stats->total_genome_bases - stats_info->cov_stats->total_Ns_bases;
		for (k_iter=0; k_iter!=kh_end(stats_info->genome_coverage_for_median); k_iter++) {
			if (kh_exist(stats_info->genome_coverage_for_median, k_iter)) {

				if (sum >= (total_genome_non_Ns_bases/2)) {
					stats_info->cov_stats->median_genome_coverage = k_iter--;
					break;
				}else{
					sum += kh_value(stats_info->genome_coverage_for_median, k_iter);;
				}
			}
        }

		file_name = calloc(strlen(user_inputs->out_file) + 30, sizeof(char));
		strcpy(file_name, user_inputs->out_file);
		strcat(file_name, ".WGCoverageReport.csv");
		FILE *out_fp = fopen(file_name, "w");

		average_coverage = (double) stats_info->cov_stats->total_genome_coverage/ (double) total_genome_non_Ns_bases;
		outputGeneralInfo(out_fp, stats_info, average_coverage, 1);

		fprintf(out_fp, "Base Stats\n");
	    fprintf(out_fp, "Bases Targeted:,%"PRIu64"\n", total_genome_non_Ns_bases);

		for(i=0; i<13; i++) {
			int32_t val = getValueFromKhash32(stats_info->genome_base_with_N_coverage, bins[i]);
			if (i==0) { val -= stats_info->cov_stats->total_Ns_bases; }		// need to remove all Ns

			float percent = calculatePercentage(val, total_genome_non_Ns_bases);
			if (i==0) fprintf(out_fp, "Bases with %d coverage:,%"PRIu32",(%.2f%%)\n", bins[i], val, percent);
			if (i>0)  fprintf(out_fp, "Bases with %d+ coverage:,%"PRIu32",(%.2f%%)\n", bins[i], val, percent);
		}

		fprintf(out_fp, "Max Coverage:,%"PRIu32",(%"PRIu32")\n", stats_info->cov_stats->max_coverage, stats_info->cov_stats->base_with_max_coverage);
	    fprintf(out_fp, "\n");
		fprintf(out_fp, "Coverage Histogram for Whole Genome (may look weird if target regions overlap...)\n");

		for (i=1; i<=1000; i++) {
            k_iter = kh_put(m32, stats_info->genome_cov_histogram, i, &ret);
            if (ret == 1)
                kh_value(stats_info->genome_cov_histogram, k_iter) = 0;

			fprintf(out_fp, "%"PRIu32",", kh_key(stats_info->genome_cov_histogram, k_iter)-1);
        }
		fprintf(out_fp, "1000,\n");

		for (k_iter=kh_begin(stats_info->genome_cov_histogram); k_iter!=kh_end(stats_info->genome_cov_histogram); k_iter++) {
			if (kh_exist(stats_info->genome_cov_histogram, k_iter)) {
				fprintf(out_fp, "%"PRIu32",", kh_value(stats_info->genome_cov_histogram, k_iter));
			}
		}
	    fprintf(out_fp, "\n");

		fclose(out_fp);
		free(file_name);
	}

	// Now we need to process target information if target bed file is provided
	if (TARGET_FILE_PROVIDED) {
		// First we need to calculate the coverage for median, this is like N50 for sequencing
		sum = 0;
		for (k_iter=0; k_iter!=kh_end(stats_info->target_coverage_for_median); k_iter++) {
			if (kh_exist(stats_info->target_coverage_for_median, k_iter)) {

				// Divya: 29826 number of Ns in our VCrome + PKv2 design; will need to change if design changes
				//if (sum >= (stats_info->cov_stats->total_targeted_bases - 29826)/2) {
				if (sum >= stats_info->cov_stats->total_targeted_bases/2) {
					stats_info->cov_stats->median_target_coverage = k_iter--;
					break;
				} else {
					sum += kh_value(stats_info->target_coverage_for_median, k_iter);
				}
			}
        }

		average_coverage = (double)stats_info->cov_stats->total_target_coverage/(double)stats_info->cov_stats->total_targeted_bases;
		file_name = calloc(strlen(user_inputs->out_file) + 30, sizeof(char));
		strcpy(file_name, user_inputs->out_file);
		strcat(file_name, ".CoverageReport.csv");
		FILE *trt_fp = fopen(file_name, "w");	// trt_fp: target_fp

		//printf("Before output general information\n");
		outputGeneralInfo(trt_fp, stats_info, average_coverage, 2);
		//printf("After output general information\n");

		fprintf(trt_fp, "Target Stats\n");

		float percent = calculatePercentage(stats_info->cov_stats->hit_target_count, stats_info->cov_stats->total_targets);
        fprintf(trt_fp, "Targets Hit:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->hit_target_count, percent);

		percent = calculatePercentage(stats_info->cov_stats->hit_target_buffer_only_count, stats_info->cov_stats->total_targets);
        fprintf(trt_fp, "Target Buffers Hit:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->hit_target_buffer_only_count, percent);
        fprintf(trt_fp, "Total Targets:,%"PRIu32"\n", stats_info->cov_stats->total_targets);
        fprintf(trt_fp, "Non target regions with high coverage:,%"PRIu32"\n", stats_info->cov_stats->non_target_good_hits);
        fprintf(trt_fp, "Base Stats\n");
        fprintf(trt_fp, "Bases Targeted:,%"PRIu32"\n", stats_info->cov_stats->total_targeted_bases);
        fprintf(trt_fp, "Buffer Bases:,%"PRIu32"\n", stats_info->cov_stats->total_buffer_bases);

		for(i=0; i<13; i++) {
            uint32_t val = getValueFromKhash32(stats_info->targeted_base_with_N_coverage, bins[i]);

            float percent = calculatePercentage(val, stats_info->cov_stats->total_targeted_bases);
            if (i==0) fprintf(trt_fp, "Bases with %d coverage:,%"PRIu32",(%.2f%%)\n", bins[i], val, percent);
            if (i>0)  fprintf(trt_fp, "Bases with %d+ coverage:,%"PRIu32",(%.2f%%)\n", bins[i], val, percent);
        }

		//printf("After the base count for 1000\n");

		fprintf(trt_fp, "\n");
		fprintf(trt_fp, "Coverage Histogram (may look weird if target regions overlap...)\n");

		for (i=1; i<=1000; i++) {
			k_iter = kh_put(m32, stats_info->target_cov_histogram, i, &ret);
			if (ret == 1)
				kh_value(stats_info->target_cov_histogram, k_iter) = 0;

			fprintf(trt_fp, "%"PRIu32",", kh_key(stats_info->target_cov_histogram, k_iter)-1);
		}
    	fprintf(trt_fp, "1000,\n");

		for (k_iter=kh_begin(stats_info->target_cov_histogram); k_iter!=kh_end(stats_info->target_cov_histogram); k_iter++) {
			if (kh_exist(stats_info->target_cov_histogram, k_iter)) {
				fprintf(trt_fp, "%"PRIu32",", kh_value(stats_info->target_cov_histogram, k_iter));
			}
		}
		fprintf(trt_fp, "\n");
		printf("After histogram plotting\n");

		fprintf(trt_fp, "Target and region coverage plot\n");
		fprintf(trt_fp, "Position,5'count,3'count\n");
		for(i = 20; i <= PRIMER_SIZE; i+=20) {
        	fprintf(trt_fp, "%d,%"PRIu32",%"PRIu32"\n", i, stats_info->five_prime[PRIMER_SIZE-(i-1)-1], stats_info->three_prime[i-1]);
        }

		//fputs("%tar-Pos,count\n", trt_fp);
		fprintf(trt_fp, "target_Pos,count\n");
        for(i = 0; i < 101; i+=2) {
            fprintf(trt_fp, "%d,%"PRIu32"\n", i, stats_info->target_coverage[i]);
        }

		//printf("Before free the file name second time\n");
		fclose(trt_fp);
		free(file_name);
	}
}

void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage, uint8_t type) {
	fprintf(fp, "Version: %s\n", VERSION);
    if (type == 2) fprintf(fp, "BUFFER size:,%d\n", BUFFER);
    fprintf(fp, "Read Stats\n");
     
    fprintf(fp, "Total Reads Produced:,%"PRIu32"\n", stats_info->cov_stats->total_reads_produced);

    uint64_t yield = stats_info->cov_stats->read_length * (uint64_t) stats_info->cov_stats->total_reads_produced;
    fprintf(fp, "Total Yield Produced:,%d,%"PRIu64"\n", stats_info->cov_stats->read_length, yield);

    yield = stats_info->cov_stats->read_length * (uint64_t) (stats_info->cov_stats->total_reads_aligned - stats_info->cov_stats->total_duplicate_reads);
    fprintf(fp, "Total Unique Yield Produced:,%"PRIu64"\n", yield);

	float percent = calculatePercentage(stats_info->cov_stats->total_duplicate_reads,stats_info->cov_stats->total_reads_aligned);
    fprintf(fp, "Duplicate Reads:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->total_duplicate_reads, percent);

    percent = calculatePercentage(stats_info->cov_stats->total_reads_aligned,stats_info->cov_stats->total_reads_produced);
    fprintf(fp, "Total Reads Aligned:,%"PRIu32",(%.2f%%)", stats_info->cov_stats->total_reads_aligned, percent);

	fprintf(fp, ",reads paired:,%"PRIu32, stats_info->cov_stats->total_reads_paired);
    fprintf(fp, ",reads paired with mapped mates:,%"PRIu32"\n", stats_info->cov_stats->total_paired_reads_with_mapped_mates);

	if (type == 2) {
		percent = calculatePercentage(stats_info->cov_stats->in_buffer_read_hit_count, stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Aligned Reads On-Buffer:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->in_buffer_read_hit_count, percent);

		percent = calculatePercentage(stats_info->cov_stats->on_target_read_hit_count, stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Aligned Reads On-Target:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->on_target_read_hit_count, percent);
	}

    fprintf(fp, "Average Coverage:,-,(%.2f)\n", average_coverage);
    if (type == 1)
		fprintf(fp, "Median Coverage:,-,(%d)\n", stats_info->cov_stats->median_genome_coverage);

    if (type == 2)
		fprintf(fp, "Median Coverage:,-,(%d)\n", stats_info->cov_stats->median_target_coverage);

	if (type == 2) {
		uint32_t reads_hit_target_or_buffer = stats_info->cov_stats->on_target_read_hit_count + stats_info->cov_stats->in_buffer_read_hit_count;
		percent = calculatePercentage(reads_hit_target_or_buffer, stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Reads that hit target or buffer:,%"PRIu32",(%.2f%%)\n", reads_hit_target_or_buffer, percent);

		fprintf(fp, "Total Aligned Reads (expected):,%"PRIu32"\n", stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Total Aligned Reads (calculated):,%"PRIu32"\n", stats_info->cov_stats->on_target_read_hit_count + stats_info->cov_stats->in_buffer_read_hit_count + stats_info->cov_stats->off_target_read_hit_count);
	}
}

void produceOffTargetWigFile(Chromosome_Tracking *chrom_tracking, char *chrom_id, Bed_Info *target_bed_info, User_Input *user_inputs, Stats_Info *stats_info) {
	int32_t i, j;
	uint16_t idx = locateChromosomeIndex(chrom_id, chrom_tracking);
	FILE *wp = fopen(user_inputs->wig_file, "a");
	fprintf(wp, "track type=wiggle_0 name=%s\n", user_inputs->bam_file);

	for(i=0; i<target_bed_info->size; i++) {

		// here we are going to erase everything that is on-target and in-buffer and set them to 0
		if (strcmp(chrom_id, target_bed_info->coords[i].chrom_id) == 0) {
			uint32_t start = target_bed_info->coords[i].start;
			uint32_t stop  = target_bed_info->coords[i].end;
			//for(j=start-BUFFER; j<=stop+BUFFER; j++) {
			for(j=start-500; j<=stop+500; j++) {
				if (j < 0 || j > chrom_tracking->chromosome_lengths[idx])
					continue;

				chrom_tracking->coverage[idx][j] = 0;
			}
		}
	}

	// once captured area + BUFFER is initialized to 0, we want to check if there is any coverage > 20 in off-target regions
	for(i=0; i<chrom_tracking->chromosome_lengths[idx]; i++) {
		if (chrom_tracking->coverage[idx][i] > 20) {
			j = i;
			stats_info->cov_stats->non_target_good_hits += 1;

			// right side 
			while(i < chrom_tracking->chromosome_lengths[idx] && chrom_tracking->coverage[idx][i] > 0)
				i++;

			// left side
			while(j > 0 && chrom_tracking->coverage[idx][i] > 0)
				j--;

			// now write to the off target wig file
			fprintf(wp, "fixedStep chrom_id=%s start=%"PRId32" step=1\n", chrom_id, j);
			uint32_t h = j;
			for(h=j; h<i; h++)
				fprintf(wp,"%"PRIu32"\n", chrom_tracking->coverage[idx][h]);
		}
	}
	fclose(wp);
}

// The followings are codes that related to gene annotation
void finish_with_error(MYSQL *con)
{
    fprintf(stderr, "%s\n", mysql_error(con));
    mysql_close(con);
    exit(1);
}

void fromStringToIntArray(char *str_in, uint32_t *array_in) {
    uint16_t counter = 0;
    char *pch;
    pch = strtok(str_in, ",|\"");
    while (pch != NULL) {
        if (strcmp(pch, "\"") != 0) {
            array_in[counter] = (uint32_t) atol(pch);
            //printf("value is %s and counter is %d\n", pch, counter);
            counter++;
        }
        pch = strtok(NULL, ",|");
    }
}

void processExonArrays(uint16_t exon_count, uint32_t *exon_starts, uint32_t *exon_ends, char *gene_name, uint32_t pos, khash_t(str) *ret_hash) {
    char string_to_add[100]="";
    uint16_t i = 0;
	int ret;
	khiter_t k_iter;

	// Create an instance of Temp_Coverage_Array and initialize it, 
	// Here we have to use Temp_Coverage_Array instead of uint32_t because it seems that I could have only one khash_t(str) type declared
	//
    //Temp_Coverage_Array *temp_cov_array = calloc(1, sizeof(Temp_Coverage_Array));
    //temp_cov_array->size = 0;

    for (i=0; i<exon_count; i++) {
        if (exon_starts[i] <= pos && pos <= exon_ends[i]) {
            //printf("Found with exon start at %d and end at %d with iteration of %d\n", exonStarts[i], exonEnds[i], i);
            sprintf(string_to_add, "%s_exon_%d", gene_name, i);

			// create the key
			k_iter = kh_put(str, ret_hash, string_to_add, &ret);
			if (ret)
				kh_key(ret_hash, k_iter) = strdup(string_to_add);
            //kh_value(coverage_hash, k_iter) = temp_cov_array;
        }
		strcpy(string_to_add, "");
    }
}

char* combinedEachAnnotation(khash_t(str) *hash_in) {
    khiter_t k_iter;
    int flag = 0, str_total_size=200;
	char *ret_string;
	ret_string = calloc(str_total_size, sizeof(char));
    strcpy(ret_string, ".");

    for (k_iter = kh_begin(hash_in); k_iter != kh_end(hash_in); ++k_iter) {
        if (kh_exist(hash_in, k_iter) && strlen(kh_key(hash_in, k_iter)) > 0) {
            if (flag == 0) {
                strcpy(ret_string, kh_key(hash_in, k_iter));
                flag++;
            } else {
				// check to see if the ret_string can hold the s_length
            	if (strlen(ret_string) >= (str_total_size - 30)) {
                	// need to dynamically allocate the memory
	                char *tmp=NULL;
					str_total_size *= 2;
					//printf("Dynamic allocation memory with size %d and string is %s\n", str_total_size, ret_string);
    	            tmp = realloc(ret_string, str_total_size);
        	        if (!tmp) {
						free(ret_string);
						ret_string = NULL;
                    	fprintf(stderr, "String realloc() failed! Exiting...");
	                    //exit(1);
    	            }
            	    ret_string = tmp;
				}
                strcat(ret_string, ";");
                strcat(ret_string, kh_key(hash_in, k_iter));
            }
        }
    }

	return ret_string;
}

// here hash_in refers to the refseq_hash, or ccds_hash or vega_hash or miRNA_hash 
void processingMySQL(MYSQL *con, char *sql, uint32_t pos_start, uint32_t pos_end, char *gene, khash_t(str) *prev_gene, khash_t(str) *Synonymous, khash_t(str) *hash_in) {
    if (mysql_query(con,sql))
        finish_with_error(con);

    MYSQL_RES *result = mysql_store_result(con);
    if (result == NULL)
        finish_with_error(con);

    //int num_fields = mysql_num_fields(result);
    //printf("Total number of fields is %d\n", num_fields);

    MYSQL_ROW row;
    int absent;
    uint16_t exon_count=0;
    khiter_t Synonymous_iter, prev_gene_iter;

    while ((row = mysql_fetch_row(result))) {
        // here I need to locate which exon it is part of and process them accordingly
        exon_count = (uint16_t) atoi(row[3]);
        uint32_t exon_starts[exon_count], exon_ends[exon_count];
        fromStringToIntArray(row[4], exon_starts);
        fromStringToIntArray(row[5], exon_ends);

        processExonArrays(exon_count, exon_starts, exon_ends, row[0], pos_start, hash_in);
		if (pos_start != pos_end)
	        processExonArrays(exon_count, exon_starts, exon_ends, row[0], pos_end, hash_in);

		//id_iter_hash = kh_put(m32, id_list, atoi(row[1]), &absent);

        // for gene row[6], Synonymous row[7], and prev_gene row[8]
        if (row[6] && strlen(row[6]) > 0 && strcmp(row[6], "NULL") != 0) {
            if (strlen(gene) == 0) {
                strcpy(gene, row[6]);
            } else {
                if (strcmp(gene, row[6]) != 0) {
                    Synonymous_iter = kh_put(str, Synonymous, row[6], &absent);
					if (absent)
						kh_key(Synonymous, Synonymous_iter) = strdup(row[6]);
                }
            }
        }

		/*
		if (row[7] && strlen(row[7]) > 0 && strcmp(row[7], "NULL") != 0) {
            Synonymous_iter = kh_put(str, Synonymous, row[7], &absent);
            if (absent)
				kh_key(Synonymous, Synonymous_iter) = strdup(row[7]);
		}*/

        if (row[8] && strlen(row[8]) > 0 && strcmp(row[8], "NULL") != 0) {
            prev_gene_iter = kh_put(str, prev_gene, row[8], &absent);
			if (absent)
				kh_key(prev_gene, prev_gene_iter) = strdup(row[8]);
        }
	}

	// Need to clean it here otherwise, valgrind will give possible memory leak error
	if (result) mysql_free_result(result);
}

char * produceGeneAnnotations(uint32_t start_in, uint32_t stop_in, char *chrom_id, MYSQL *con) {
    char *prev_sql  = " start, end, exon_count, exon_starts, exon_ends, gene_symbol, alias_gene_symbol, prev_gene_symbol ";
	char *mid_sql   = calloc(150, sizeof(char)); 
    sprintf(mid_sql,  " chrom='%s' AND ((start <= %"PRIu32" AND %"PRIu32" <= end) OR (start <= %"PRIu32" AND %"PRIu32" <= end))", chrom_id, start_in, start_in, stop_in, stop_in);

	// for debugging
	/*
	//if (strcmp(chrom_id, "1") == 0)
	//	return;

	//if (start_in < 70276350) {
	if (start_in < 53676957) {
		//printf("%s\n", sql);
		return;
	} else {
		printf("%s\n", mid_sql);
	}
	if (start_in == 1246709) {
		printf("%s\n", mid_sql);
	}*/

    char gene[100]="";
	khash_t(str) *refseq     = kh_init(str);     // hash_table using string as key
	khash_t(str) *ccds       = kh_init(str);     // hash_table using string as key
	khash_t(str) *vega       = kh_init(str);     // hash_table using string as key
	khash_t(str) *Synonymous = kh_init(str);     // hash_table using string as key
	khash_t(str) *prev_gene  = kh_init(str);     // hash_table using string as key
	khash_t(str) *miRNA      = kh_init(str);     // hash_table using string as key

	// for refseq
	char *sql = calloc(strlen(prev_sql) + strlen(mid_sql) + 75, sizeof(char));
    sprintf(sql, "SELECT DISTINCT refseq_name, %s FROM RefSeq_annotation WHERE %s ", prev_sql, mid_sql);
    //printf("%s\n", sql);
    processingMySQL(con, sql, start_in, stop_in, gene, prev_gene, Synonymous, refseq);
	memset(sql,0,strlen(sql));	 

    // for ccds
    sprintf(sql, "SELECT DISTINCT ccds_name, %s FROM CCDS_annotation WHERE %s ", prev_sql, mid_sql);
    //printf("%s\n", sql);
    processingMySQL(con, sql, start_in, stop_in, gene, prev_gene, Synonymous, ccds);
	memset(sql,0,strlen(sql));	

    // for vega
    sprintf(sql, "SELECT DISTINCT vega_name, %s FROM VEGA_annotation WHERE %s ", prev_sql, mid_sql);
    //printf("%s\n", sql);
    processingMySQL(con, sql, start_in, stop_in, gene, prev_gene, Synonymous, vega);
	memset(sql,0,strlen(sql));	 

    // now for miRNA
    sprintf(sql, "SELECT DISTINCT miRNA_name, %s FROM miRNA_annotation WHERE %s ", prev_sql, mid_sql);
    //printf("%s\n", sql);
    processingMySQL(con, sql, start_in, stop_in, gene, prev_gene, Synonymous, miRNA);
	memset(sql,0,strlen(sql));	 

    // now we need to combine everything together
	char *refseq_str = combinedEachAnnotation(refseq);
	char *ccds_str   = combinedEachAnnotation(ccds);
	char *vega_str   = combinedEachAnnotation(vega);
	char *miRNA_str  = combinedEachAnnotation(miRNA);
	char *prev_gene_str  = combinedEachAnnotation(prev_gene);
	char *Synonymous_str = combinedEachAnnotation(Synonymous);

	uint16_t annotation_size = strlen(gene) + strlen(prev_gene_str) + strlen(Synonymous_str) + strlen(refseq_str) + strlen(ccds_str) 
								+ strlen(vega_str) + strlen(miRNA_str);
	char *annotation = calloc(annotation_size + 50, sizeof(char));

	sprintf(annotation, "%s\t%s\t%s\t%s\t%s\t%s\t%s", gene, prev_gene_str, Synonymous_str, refseq_str, ccds_str, vega_str, miRNA_str);

	cleanKhashStr(refseq, 2);
	cleanKhashStr(ccds, 2);
	cleanKhashStr(vega, 2);
	cleanKhashStr(miRNA, 2);
	cleanKhashStr(Synonymous, 2);
	cleanKhashStr(prev_gene, 2);

	//printf("%s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, start_in, stop_in);
	//fflush(stdout);

	// properly empty/reset the strings I declared
	if (sql) free(sql);
	//free(prev_sql);
	if (mid_sql) free(mid_sql);

	memset(gene,0,sizeof(gene));
	if (refseq_str) free(refseq_str);
	if (prev_gene_str) free(prev_gene_str);
	if (Synonymous_str) free(Synonymous_str);
	if (ccds_str) free(ccds_str);
	if (vega_str) free(vega_str);
	if (miRNA_str) free(miRNA_str);
	//memset(,0,sizeof());
	//printf("return produceGeneAnnations %d\n", start_in);
	
	return annotation;
}
