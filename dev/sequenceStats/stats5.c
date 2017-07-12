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

void readBufferInit(Read_Buffer *read_buff_in) {
	int i;
	for (i=0; i<read_buff_in->size; i++) {
        read_buff_in->chunk_of_reads[i] = bam_init1();
    }
}

void readBufferDestroy(Read_Buffer *read_buff_in) {
	int i;
	for (i=0; i<read_buff_in->size;i++) {
		bam_destroy1(read_buff_in->chunk_of_reads[i]);
	}
}

uint32_t readBam(samFile *sfin, bam_hdr_t *header, bool more_to_read, Read_Buffer *read_buff_in) {
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

void processBamChunk(User_Input *user_inputs, Stats_Info *stats_info, khash_t(str) *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in) {
	// it is the flag that is used to indicate if we need to add the khash into the coverage_hash
	int not_added = true;	
	int i = 0;
	char *cur_chr  = malloc(15 * sizeof(char));
	char *last_chr = malloc(30 * sizeof(char));
	strcpy(last_chr, "SOMEVERYFAKETHINGGOESHERE");

	for (i=0; i<read_buff_in->size; i++) {
		if(user_inputs->percentage < 1.0) {
            srand((unsigned int)time(NULL));	// set random seed
            float random_num = (float)rand() / (float)RAND_MAX;
            if(random_num < user_inputs->percentage) continue;
        }

		stats_info->total_reads_produced++;

        // get core alignment info
        bam1_core_t *core = &read_buff_in->chunk_of_reads[i]->core;

		// Need to check various 'READ' flags regarding the current read before doing statistics analysis
		if(core->qual < user_inputs->min_map_quality) {
            continue;
        }

        if(core->flag & BAM_FQCFAIL){           // Read Fails Vendor Quality Check
			continue;
        }
        if(core->flag & BAM_FSECONDARY){        // Read Alignment is not Primary Alignment
            continue;
        }
 
        if(core->flag & BAM_FUNMAP)	{			// Read Unmapped
            //printf("Unmapped! region.  breaking!!!\n");
             continue;
		}

		stats_info->total_reads_aligned++;
        
		if(core->flag & BAM_FPROPER_PAIR) {		// Read is properly paired
			stats_info->total_reads_paired++;
            if(!(core->flag & BAM_FMUNMAP)) {	// Read is Paird with Mapped Mate 
				stats_info->total_paired_reads_with_mapped_mates++;
            }
        }

		if(core->flag & BAM_FDUP) {				// Read is a Duplicate (either optical or PCR)
			stats_info->total_duplicate_reads++;
            if(user_inputs->remove_duplicate) {continue;}
        }

		if (core->flag & BAM_FSUPPLEMENTARY) {
			stats_info->total_supplementary_reads++;
			if (user_inputs->remove_supplementary_alignments)
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
        	if (stats_info->read_length <= 0)
				stats_info->read_length = read_buff_in->chunk_of_reads[i]->core.l_qseq;
        	///////////////////////////////////////////////////////////////////
		}

        processRecord(user_inputs, stats_info, coverage_hash, header->target_name[core->tid], read_buff_in->chunk_of_reads[i]);
    }

	free(cur_chr);
	free(last_chr);
    printf("Done read bam\n");
}

void processRecord(User_Input *user_inputs, Stats_Info *stats_info, khash_t(str) *coverage_hash, char * chrom_id, bam1_t *rec) {
	int i, ret;
	khiter_t outer_iter, inner_iter;

	// left most position of alignment in zero based coordinates (as BAM is 0-based)
	uint32_t start = rec->core.pos + 1;		

	// get key iterator for chrom_id (outer hash table: khash_t(str))
	outer_iter = kh_put(str, coverage_hash, chrom_id, &ret);	
	if (ret)
		fprintf(stderr, "Something went wrong, has key for %s is not added\n", chrom_id);

	// get the quality information
	uint8_t *qual = bam_get_qual(rec);

	for (i=0; i<rec->core.l_qseq; i++) {
		stats_info->total_aligned_bases++;

		if (qual[i] < user_inputs->min_base_quality)	// filter based on the MIN_BASE_SCORE
			continue;

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
void combineThreadResults(Chromosome_Tracking *chrom_tracking, khash_t(str) *coverage_hash, bam_hdr_t *header) {
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
				uint32_t chrom_len = header->target_len[getChromIndexFromID(header, chrom_id)];
				chromosomeTrackingInit(chrom_tracking, chrom_id, chrom_len, chrom_tracking->number_tracked);
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

void writeCoverage(char * chrom_id, khash_t(str) *Ns_buffer_hash, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, FILE *cov_fp, FILE *wig_fp, FILE *wgs_fp) {
    // First, we need to find the index that is used to track current chromosome chrom_id
    uint8_t idx = locateChromosomeIndex(chrom_id, chrom_tracking);

    // for the whole genome, we need to use the file that contains regions of all Ns in the reference
    // As they will be not used, so we are going to set the count info in these regions to 0
    if (N_FILE_PROVIDED)
        zeroAllNsRegions(chrom_id, Ns_buffer_hash, chrom_tracking);

    // write to the file that contains whole genome info
	int i;
    if(user_inputs->wgs_coverage) {

        char header_line[30];
        strcpy(header_line, ">chromosome_");
        strcat(header_line, chrom_id);
        fputs(header_line, wgs_fp);     // no need to add newline as the next line will take care of it!

        for (i = 0; i < chrom_tracking->chromosome_lengths[i]; i++) {
  			if(i%100==0) fputc('\n', wgs_fp);

			addBaseStats(stats_info, chrom_tracking->coverage[idx][i], 0, 1);
            fprintf(wgs_fp, "%d ", chrom_tracking->coverage[idx][i]);
        }
        fputc('\n', wgs_fp);
    }

	// if the target bed file is available, we will need to handle it here and write the results to cov.fasta file
	if (TARGET_FILE_PROVIDED) {
		for(i = 0; i < target_info->size; i++) {
            if ( strcmp(target_info->coords[i].chr, chrom_id) != 0)
                continue;

            //stats_info->total_targets++;
            int j;
            uint32_t start = target_info->coords[i].start;
            uint32_t end = target_info->coords[i].end;
            int length = end - start + 1;
            bool collect_target_cov = length > 99 ? true : false ;  // TODO: is it the min length of the exon? Check.

            if (collect_target_cov) {
                for(j = 0; j < PRIMER_SIZE; j++) {
                    if ( (start - j) < 0 || (end + j) >= chrom_tracking->chromosome_lengths[idx])
                        continue;

					//addValueToKhashBucket32(stats_info->five_prime,  j, chrom_tracking->coverage[idx][start-j]);
					//addValueToKhashBucket32(stats_info->three_prime, j, chrom_tracking->coverage[idx][end+j]);
					stats_info->five_prime[j]  += chrom_tracking->coverage[idx][start-j];
					stats_info->three_prime[j] += chrom_tracking->coverage[idx][end+j];
                }
			}

			//In the original java code, the following is used for the pc and pc2 definition
            //short pc[101], pc2[101];
			//After discussion with Divya and Qiaoyan, here is what we understand.
			//pc and pc2 is used to setup the 'target_coverage' variable which is defined as an array of 101 in the origain java code
			//here I will use percentage_bin and percentage_count_per_bin
			short percentage_bin[101], percentage_count_per_bin[101];


            bool target_hit = false;
            fprintf(cov_fp, ">%s %"PRIu32" %"PRIu32, chrom_id, start, end);

            bool space_it = false;
            if(end - start > 10000) space_it = true;

            for(j = 0; j < length; j++) {
                if (j > chrom_tracking->chromosome_lengths[idx])
                    continue;

                if (space_it && j%100 == 0) fputc('\n', cov_fp);    // enter a new line after every 100 bases

				short cov = chrom_tracking->coverage[idx][j+start];
                target_hit = addBaseStats(stats_info, cov, 1, 0);

                // output to the cov.fasta file
                fprintf(cov_fp, "%d ", cov);

                if (collect_target_cov) {
                    short percentage_pos = (short)((float)(j)/(float)length*100+0.5);
                    percentage_bin[percentage_pos] += cov;
                    percentage_count_per_bin[percentage_pos]++;
                }
            }

            // output a newline char to the cov.fasta file 
            fputc('\n', cov_fp);

			for (j = 0; j < 101; j++) {
				if(percentage_count_per_bin[j] != 0) {
                    int d = (int) (((double)percentage_bin[j]/(double)percentage_count_per_bin[j])+0.5);
                    percentage_bin[j] = (short) d;	
                }
            }

            for(j = 0; j < 101; j++)
                //addValueToKhashBucket32(stats_info->target_coverage, i, percentage_bin[i]);
				stats_info->target_coverage[j] +=percentage_bin[j];

			if (target_hit) {
				stats_info->hit_target_count += 1;
			} else {
				// need to write to the missed target file
				fprintf(wig_fp, "%s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, start, end);
				bool hit = false;
				for (j = start - BUFFER; j < start && !hit; j++) {
					if (j >= chrom_tracking->chromosome_lengths[idx])
						continue;

					if ( chrom_tracking->coverage[idx][j] > 0 )
						hit = true;
				}

				if (hit)
					stats_info->hit_buffer_only_count += 1;
			}
		}
	}
}

bool addBaseStats(Stats_Info *stats_info, uint16_t cov_val, uint8_t target, uint8_t wgs) {
	// need to check and update max coverage information
	if (cov_val > stats_info->max_coverage) {
		stats_info->max_coverage = cov_val;
		stats_info->base_with_max_coverage = 1;
	} else if (cov_val == stats_info->max_coverage) {
		stats_info->base_with_max_coverage += 1;
	}

	if (cov_val < 0) {
        fprintf(stderr, "Coverage less than 0!!!!!!!\n");
    }

	if (target == 1) addValueToKhashBucket(stats_info->target_cov_histogram, cov_val, 1);
	if (wgs == 1)    addValueToKhashBucket(stats_info->genome_cov_histogram, cov_val, 1);

    if (target == 1) stats_info->total_target_coverage += cov_val;
    if (wgs == 1)    stats_info->total_genome_coverage += cov_val;

	bool target_hit = false;

    if (cov_val > 0) {
        target_hit=true;
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 1, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 1, 1);
    }

	if (cov_val >= 5) {
		if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 5, 1);
		if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 5, 1);
	}

	if (cov_val >= 10) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 10, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 10, 1);
    }

	if (cov_val >= 15) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 15, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 15, 1);
    }

	if (cov_val >= 20) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 20, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 20, 1);
    }

	if (cov_val >= 30) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 30, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 30, 1);
    }

	if (cov_val >= 40) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 40, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 40, 1);
    }

	if (cov_val >= 50) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 50, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 50, 1);
    }

	if (cov_val >= 60) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 50, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 50, 1);
    }

	if (cov_val >= 100) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 100, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 100, 1);
    }

	if (cov_val >= 500) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 500, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 500, 1);
    }

	if (cov_val >= 1000) {
        if (target == 1) addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 1000, 1);
        if (wgs == 1)    addValueToKhashBucket(stats_info->genome_base_with_N_coverage, 1000, 1);
    }

	if (target == 1) addValueToKhashBucket(stats_info->target_coverage_for_median, cov_val, 1);
	if (wgs == 1)    addValueToKhashBucket(stats_info->genome_coverage_for_median, cov_val, 1);

	return target_hit;
}

void writeReport(Stats_Info *stats_info, char *output_file) {
	uint32_t non_duplicate_reads = stats_info->total_reads_aligned - stats_info->total_duplicate_reads;

	if (TARGET_FILE_PROVIDED && stats_info->total_targeted_bases == 0) {
		fprintf(stderr, "Total targeted bases is zero.  This means that no read has aligned to a chromosome that contains a target. No target matches a chromosome in the BAM, or something else very weird.  Aborting.");
		exit(1);
	}

	if (stats_info->total_reads_aligned == 0) {
		fprintf(stderr, "No reads aligned. Aborting.");
		exit(1);
	}

	if (non_duplicate_reads == 0) {
        fprintf(stderr, "All reads are duplicates. Aborting.");
        exit(1);
    }

	if (TARGET_FILE_PROVIDED && stats_info->total_targets == 0) {
        //I don't think we should ever see this error, as its dealt with above.
        fprintf(stderr, "No target regions given.  Aborting.");
        exit(1);
    }

	//Do not consider the Ns for Median calculation.
	uint64_t total_genome_non_Ns_bases = stats_info->total_genome_bases- stats_info->total_Ns_bases;
	int i=0, sum=0;
	khiter_t k_iter;
	for (k_iter=kh_begin(stats_info->genome_coverage_for_median); 
			k_iter!=kh_end(stats_info->genome_coverage_for_median); k_iter++) {

		int tmp_val = kh_value(stats_info->genome_coverage_for_median, k_iter);
		i++;

        if (tmp_val >= (total_genome_non_Ns_bases/2)) {
            stats_info->median_genome_coverage = i--;
            break;
        }else{
            sum += tmp_val;
        }
    }

	int bins[15] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 100, 500, 1000};
	
	int file_len = strlen(output_file);
	char *file_name = calloc(file_len + 30, sizeof(char));
	strcpy(file_name, output_file);
	strcat(file_name, ".WGCoverageReport.csv");
	FILE *out_fp = fopen(file_name, "w");

	double average_coverage = (double)stats_info->total_genome_coverage/(double)total_genome_non_Ns_bases;

    fprintf(out_fp, "Version: %s\n", VERSION);
    fprintf(out_fp, "Read Stats\n");
    fprintf(out_fp, "Total Reads Produced:,%"PRIu32"\n", stats_info->total_reads_produced);
    fprintf(out_fp, "Total Yield Produced:,%d,%"PRIu32"\n", stats_info->read_length, stats_info->read_length * stats_info->total_reads_produced);
    fprintf(out_fp, "Total Unique Yield Produced:,%"PRIu32"\n", stats_info->read_length * (stats_info->total_reads_aligned-stats_info->total_duplicate_reads));
    fprintf(out_fp, "Duplicate Reads:,%"PRIu32",(%.2f %%)\n", stats_info->total_duplicate_reads,calculatePercentage(stats_info->total_duplicate_reads,stats_info->total_reads_aligned));
    fprintf(out_fp, "Total Reads Aligned:,%"PRIu32",(%.2f%%)", stats_info->total_reads_aligned, calculatePercentage(stats_info->total_reads_aligned,stats_info->total_reads_produced));
    fprintf(out_fp, ",reads paired:,%"PRIu32, stats_info->total_reads_paired);
    fprintf(out_fp, ",reads paired with mapped mates:,%"PRIu32"\n", stats_info->total_paired_reads_with_mapped_mates);
    fprintf(out_fp, "Average Coverage:,-,(%.2f)\n", average_coverage);
    fprintf(out_fp, "Median Coverage:,-,(%d)\n", stats_info->median_genome_coverage);           
    fprintf(out_fp, "Base Stats\n");
    fprintf(out_fp, "Bases Targeted:,%"PRIu64"\n", total_genome_non_Ns_bases);
    fprintf(out_fp, "Bases with 0 coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 0) - stats_info->total_Ns_bases, calculatePercentage((getValueFromKhash(stats_info->genome_cov_histogram, NULL, 0) - stats_info->total_Ns_bases), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 1+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 1), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 1), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 5+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 5), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 5), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 10+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 10), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 10), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 15+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 15), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 15), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 20+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 20), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 20), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 30+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 30), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 30), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 40+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 40), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 40), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 50+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 50), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 50), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 60+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 60), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 60), total_genome_non_Ns_bases));
    fprintf(out_fp, "Bases with 100+ coverage:,%"PRIu16",(%.2f%%)\n", getValueFromKhash(stats_info->genome_cov_histogram, NULL, 100), calculatePercentage(getValueFromKhash(stats_info->genome_cov_histogram, NULL, 100), total_genome_non_Ns_bases));
    fprintf(out_fp, "Max Coverage:,%"PRIu16",(%"PRIu16")\n", stats_info->max_coverage, stats_info->base_with_max_coverage);
    fprintf(out_fp, "\n");
    fprintf(out_fp, "Coverage Histogram for Whole Genome (may look weird if target regions overlap...)\n");

	for (k_iter=kh_begin(stats_info->genome_cov_histogram); k_iter!=kh_end(stats_info->genome_cov_histogram); k_iter++)
        fprintf(out_fp, "%"PRIu16",", kh_key(stats_info->genome_cov_histogram, k_iter));
    fprintf(out_fp, "\n");

	for (k_iter=kh_begin(stats_info->genome_cov_histogram); k_iter!=kh_end(stats_info->genome_cov_histogram); k_iter++)
        fprintf(out_fp, "%"PRIu16",", kh_value(stats_info->genome_cov_histogram, k_iter));
    fprintf(out_fp, "\n");

	free(file_name);
	fclose(out_fp);

	// Now we need to process target information if target bed file is provided
	if (TARGET_FILE_PROVIDED) {

	}
}
