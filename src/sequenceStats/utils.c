/*
 * =====================================================================================
 *
 *		Filename:		utils.c
 *
 *		Description:	The implementation file for the utility functions
 *
 *      Version:		1.0
 *      Created:		02/06/2017 04:45:04 PM
 *      Revision:		none
 *      Compiler:		gcc
 *
 *      Author:			Peiming (Peter) Huang, phuang@bcm.edu
 *      Company:		Baylor College of Medicine
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>		// for file access() and getopt()
#include <dirent.h>		// for checking output directory
#include <libgen.h>		// for function basename()
#include "terms.h"
#include "utils.h"
#include "user_inputs.h"
#include "annotation.h"


// split a string into an array of words
// 0: RefSeq	1:CCDS		2:VEGA/Gencode		3:miRNA		4:everything else (such as SNP, Pseudo gene etc)
// As the string size varies, it might be a good idea to set them to a pre-defined size
// The reason is that when we need to remove duplicates by shifting, the size using strlen
// for the previous one might not be big enough for the later one
// So I will set it here as 50. The longest one seems to be 26
//
void splitStringToKhash(char *stringPtr, khash_t(khStrInt) **khashArrayPtr, uint8_t index) {
	// make a copy of the input string and only work on the copied one
	//
	char *copy_info = calloc(strlen(stringPtr)+2, sizeof(char));
	strcpy(copy_info, stringPtr);

	//
	char *savePtr = copy_info;
	char *tokPtr;
	int absent;
	khiter_t iter;

	// process the copied string and fetch all of the remaining elements
	//
	uint32_t j=0;
	while ((tokPtr = strtok_r(savePtr, ";", &savePtr))) {
		iter = kh_put(khStrInt, khashArrayPtr[index], tokPtr, &absent);
		if (absent) {
			kh_key(khashArrayPtr[index], iter)   = strdup(tokPtr);
			kh_value(khashArrayPtr[index], iter) = 0;
		}
		kh_value(khashArrayPtr[index], iter)++;
		j++;
    }

	if (copy_info != NULL) {
		free(copy_info);
		copy_info=NULL;
	}
}

// to process StrInt Hash variable to get the low coverage regions base count
//
uint32_t processLowCovRegionFromKhash(khash_t(khStrInt) *low_cov_regions, char ** output) {
	uint32_t low_cov_region_size=0, stringLength=0;
	uint32_t *start_L, *end_L;
	uint32_t size_L=0;
	uint32_t capacity = 2;
	start_L = calloc(capacity, sizeof(uint32_t));
	end_L   = calloc(capacity, sizeof(uint32_t));

	khiter_t iter;
	for (iter = kh_begin(low_cov_regions); iter != kh_end(low_cov_regions); ++iter) {
		if (kh_exist(low_cov_regions, iter)) {
			// get the key information which contains the start and end of low coverage region
			// Here is an example:	[235856749-235856786] = 1
			//						[235878500-235878529] = 1
			// kh_key() returns a const char*, so we need to make a copy first 
			// because strtok_r() will destroy the string passed in
			//
			char key[80];
			strcpy(key, kh_key(low_cov_regions,iter));
			stringLength += strlen(key) + 1;

			char *savePtr = key;
			char *tokPtr;
			uint32_t start=0, end=0;
			int idx = 0;
			while ((tokPtr = strtok_r(savePtr, "-", &savePtr))) {
				if (idx == 0) {
					// check if we need more space for memory allocation
					//
					if (capacity == size_L) {
						capacity *= 2;
						start_L = realloc(start_L, capacity * sizeof(uint32_t));
						end_L   = realloc(end_L, capacity * sizeof(uint32_t));
						if (start_L == NULL || end_L == NULL) {
							fprintf(stderr, "ERROR: Memory re-allocation failed at processLowCovRegionFromKhash()!\n");
							exit(EXIT_FAILURE);
						}
					}
					start = (uint32_t) strtol(tokPtr, NULL, 10);
					start_L[size_L] = start;
				}

				if (idx == 1) {
					end = (uint32_t) strtol(tokPtr, NULL, 10);
					end_L[size_L] = end;
				}

				idx++;
			}
			size_L++;

			low_cov_region_size += end - start;
		}
	}

	// Now need to sort the uint32_t array
	//
	qsort(start_L, size_L, sizeof(uint32_t), compare);
	qsort(end_L, size_L, sizeof(uint32_t), compare);

	// Finally, combine them together
	//
	int i;
	char tmp_str[30];

	// expand the memory for output
	//
	if (stringLength > 0) {
		if (strlen(*output) <= stringLength+50) {
			*output = realloc(*output, stringLength+50*sizeof(char));

			if (*output == NULL) {
				fprintf(stderr, "ERROR: Memory re-allocation failed at the processLowCovRegionFromKhash()!\n");
				exit(EXIT_FAILURE);
			}
		}

		*output[0] = '\0';	// set to the null terminator, so we could use strcat() all the way

		for (i=0; (unsigned int) i<size_L; i++) {
			if (i > 0) 
				strcat(*output, ";");

			sprintf(tmp_str, "%"PRIu32, start_L[i]);
			strcat(*output, tmp_str);

			sprintf(tmp_str, "%"PRIu32, end_L[i]);
			strcat(*output, "-");
			strcat(*output, tmp_str);
		}
	}

	// clean-up
	//
	if (start_L != NULL) {
		free(start_L);
		start_L=NULL;
	}

	if (end_L != NULL) {
		free(end_L);
		end_L=NULL;
	}

	return low_cov_region_size;
}

uint32_t processLowCovRegionFromStrArray(StringArray *merged_low_cov_regions, char **output) {
	uint32_t low_cov_region_size=0, stringLength=0;
	uint32_t i;

	for (i=0; i<merged_low_cov_regions->size; i++) {
		stringLength += strlen(merged_low_cov_regions->theArray[i]) + 1;
	}

	*output = realloc(*output, (stringLength + 50) * sizeof(char));
	if (*output == NULL) {
		fprintf(stderr, "ERROR: Memory re-allocation failed at processLowCovRegionFromStrArray()\n");
		exit(EXIT_FAILURE);
	}

	*output[0] = '\0';  // set to the null terminator, so we could use strcat() all the way

	for (i=0; i<merged_low_cov_regions->size; i++) {
		if (i > 0)
			strcat(*output, ";");

		strcat(*output, merged_low_cov_regions->theArray[i]);

		// find the size of low coverage region
		//
		char *savePtr = merged_low_cov_regions->theArray[i];
		char *tokPtr;
		uint32_t start=0, end=0, k=0;

		while ((tokPtr = strtok_r(savePtr, "-", &savePtr))) {
			if (k==0)
				start = (uint32_t) strtol(tokPtr, NULL, 10);

			if (k==1)
				end = (uint32_t) strtol(tokPtr, NULL, 10);

			k++;
		}
		low_cov_region_size += end - start;
	}

	return low_cov_region_size;
}

void addToGeneTranscriptKhashTable(char *gene_symbol, char *transcript_name, khash_t(khStrStrArray) *gene_transcripts, khash_t(khStrInt) *seen_transcript) {
	uint32_t i;
	int absent;

	// First check if the current gene_symbol exists
	//
	khiter_t iter = kh_put(khStrStrArray, gene_transcripts, gene_symbol, &absent);
	if (absent) {
		// key doesn't exist
		//
		kh_key(gene_transcripts, iter) = strdup(gene_symbol);
		kh_value(gene_transcripts, iter) = calloc(1, sizeof(StringArray));
		kh_value(gene_transcripts, iter)->size = 0;
		kh_value(gene_transcripts, iter)->capacity = 3;
		kh_value(gene_transcripts, iter)->theArray = calloc(kh_value(gene_transcripts, iter)->capacity, sizeof(char*));

		// Initialize the theArray
		//
		for (i=0; i<kh_value(gene_transcripts, iter)->capacity; i++)
			kh_value(gene_transcripts, iter)->theArray[i]=NULL;
	}

	// check to see if we have seen the current transcript_name
	//
	khiter_t its = kh_put(khStrInt, seen_transcript, transcript_name, &absent);
	if (absent) {
		// key doesn't exist, which means we haven't seen it
		//
		kh_key(seen_transcript, its) = strdup(transcript_name);
		kh_value(seen_transcript, its) = 1;

		// Need to check if we have allocated enough space for the current gene symbol of the gene_transcripts
		//
		if (kh_value(gene_transcripts, iter)->size == kh_value(gene_transcripts, iter)->capacity) {
			kh_value(gene_transcripts, iter)->capacity *= 2;
			kh_value(gene_transcripts, iter)->theArray =
				realloc(kh_value(gene_transcripts, iter)->theArray, kh_value(gene_transcripts, iter)->capacity * sizeof(char*));

			if (kh_value(gene_transcripts, iter)->theArray == NULL) {
				fprintf(stderr, "ERROR: Memory re-allocation failed at addToGeneTranscriptKhashTable()!\n");
				exit(EXIT_FAILURE);
			}

			// Initialize the theArray
			//
			for (i=kh_value(gene_transcripts, iter)->size; i<kh_value(gene_transcripts, iter)->capacity; i++)
				kh_value(gene_transcripts, iter)->theArray[i]=NULL;
		}

		// add the transcript_name to the gene_transcripts array for the current gene_symbol key
		//
		uint16_t idx = kh_value(gene_transcripts, iter)->size;
		kh_value(gene_transcripts, iter)->theArray[idx] = calloc(strlen(transcript_name)+1, sizeof(char));
		strcpy(kh_value(gene_transcripts, iter)->theArray[idx], transcript_name);
		kh_value(gene_transcripts, iter)->size++;
	}
}

void cleanKhashStr(khash_t(str) *hash_to_clean, uint8_t type) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
		if (kh_exist(hash_to_clean, k)) {
			// clean key if the key exist
			if (kh_key(hash_to_clean, k)) free((char *) kh_key(hash_to_clean, k));

			if (type == 1) {
				// clean Temp_Coverage_Array
				free(kh_value(hash_to_clean, k)->cov_array);
				free(kh_value(hash_to_clean, k));
				//kh_del(str, hash_to_clean, k);
			}
		}
	}
	//printf("before clean hash string\n");

	if (hash_to_clean) kh_destroy(str, hash_to_clean);
	//printf("after clean hash string\n");
}

void cleanKhashStrStrArray(khash_t(khStrStrArray) * hash_to_clean) {
	khint_t k;
	int i;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
		if (kh_exist(hash_to_clean, k)) {
			// clean the value first
			//
			for (i=0; i<kh_value(hash_to_clean, k)->size; i++) {
				if (kh_value(hash_to_clean, k)->theArray[i] != NULL)
					free(kh_value(hash_to_clean, k)->theArray[i]);
			}

			// clean theArray pointer
			//
			if (kh_value(hash_to_clean, k)->theArray)
				free(kh_value(hash_to_clean, k)->theArray);

			// clean key if the key exist
			//
			if (kh_key(hash_to_clean, k)) free((char *) kh_key(hash_to_clean, k));

			// clean value if it exists
			//
			if (kh_value(hash_to_clean, k)) free (kh_value(hash_to_clean, k));
		}
	}

	if (hash_to_clean) kh_destroy(khStrStrArray, hash_to_clean);
}

void cleanKhashStrStr(khash_t(khStrStr) * hash_to_clean) {
    khint_t k;
    for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
        if (kh_exist(hash_to_clean, k)) {
            // clean key if the key exist
            //
            if (kh_key(hash_to_clean, k)) 
                free((char *) kh_key(hash_to_clean, k));

            // clean value if it exists
            //
            if (kh_value(hash_to_clean, k)) 
                free ((char *) kh_value(hash_to_clean, k));
        }
    }

    if (hash_to_clean) kh_destroy(khStrStr, hash_to_clean);
}

void cleanGeneTranscriptPercentage(khash_t(khStrGTP) *gene_transcript_percentage_hash) {
	khint_t k;
	for (k = kh_begin(gene_transcript_percentage_hash); k != kh_end(gene_transcript_percentage_hash); ++k) {
		if (kh_exist(gene_transcript_percentage_hash, k)) {
			// clean transcript_percentage.transcript_name value first
			//
			int i;
			for (i=0; i<kh_value(gene_transcript_percentage_hash, k)->size; i++) {
				if (kh_value(gene_transcript_percentage_hash, k)->transcript_percentage[i].transcript_name)
					free(kh_value(gene_transcript_percentage_hash, k)->transcript_percentage[i].transcript_name);
			}

			// clean transcript_percentage array variable
			//
			if (kh_value(gene_transcript_percentage_hash, k)->transcript_percentage)
				free(kh_value(gene_transcript_percentage_hash, k)->transcript_percentage);

			// clean Gene_Transcript_Percentage variable
			//
			if (kh_value(gene_transcript_percentage_hash, k))
				free(kh_value(gene_transcript_percentage_hash, k));

			// clean gene_symbol key
			//
			if (kh_key(gene_transcript_percentage_hash, k))
				free((char*) kh_key(gene_transcript_percentage_hash, k));
		}
	}

	if (gene_transcript_percentage_hash) kh_destroy(khStrGTP, gene_transcript_percentage_hash);
}

// This function is used to dynamically allocate the memory and then copy everything in
//
void dynamicStringAllocation(char *str_in, char **storage_str) {
    if (str_in == NULL) return;     // nothing to add

	if (*storage_str) {
		if (strlen(str_in) > strlen(*storage_str)) {
			*storage_str = realloc(*storage_str, (strlen(*storage_str) + strlen(str_in) + 2)*sizeof(char));
			if (*storage_str == NULL) {
				fprintf(stderr, "ERROR: Dynamic Memory allocation failed\n");
				exit(EXIT_FAILURE);
			}
		}
	} else {
		*storage_str = calloc(strlen(str_in) + 1, sizeof(char));
	}

	if (str_in == NULL || strlen(str_in) == 0) {
		strcpy(*storage_str, ".");
	} else {
		strcpy(*storage_str, str_in);
	}
}

int32_t locateChromosomeIndexForRegionSkipMySQL(char *chrom_id, Regions_Skip_MySQL *regions_in) {
	int32_t i=0;
    for (i = 0; i < regions_in->chrom_list_size; i++) {
		if (regions_in->chromosome_ids[i]) {
			if (strcmp(chrom_id, regions_in->chromosome_ids[i]) == 0) {
				return i;
			}
        }
    }

	//fprintf(stderr, "Something is wrong because the chromosome %s couldn't be found\n", chrom_id);
	return -1;
}

int32_t findChromsomeIndex(Chromosome_Tracking *chrom_tracking, bam_hdr_t *header, int32_t index) {
    int32_t chrom_index;
    for (chrom_index=0; chrom_index<(int32_t)chrom_tracking->number_of_chromosomes; chrom_index++) {
        if (strcmp(header->target_name[index], chrom_tracking->chromosome_ids[chrom_index]) == 0)
            return chrom_index;
    }

    fprintf(stderr, "Can't find the corresponding chromosome tracking index for header's index %d\n", index);
    return -1;      // something might be wrong here
}

int32_t findTargetBufferIndex(Target_Buffer_Status *target_buffer_status, int32_t number_of_chromosomes, char* chrom_id) {
    int32_t i;
    for (i=0; i<number_of_chromosomes; i++) {
        if (strcmp(chrom_id, target_buffer_status[i].chrom_id) == 0) {
            return i;
        }
    }

    return -1;
}

char* createTmpFileName(char* base_file_name, char* string_to_add) {
    uint16_t s_len = strlen(base_file_name) + strlen(string_to_add) + 5;
    char *tmp_file_name = calloc(s_len, sizeof(char));
    strcpy(tmp_file_name, base_file_name);
    strcat(tmp_file_name, ".");
    strcat(tmp_file_name, string_to_add);

    return tmp_file_name;
}

void statsInfoInit(Stats_Info *stats_info, User_Input *user_inputs) {
	if (!stats_info) {
		fprintf(stderr, "ERROR: Memory allocation failed for Stats_Info\n");
		exit(EXIT_FAILURE);
	}

    stats_info->read_cov_stats = calloc(1, sizeof(Read_Coverage_Stats));
    stats_info->wgs_cov_stats  = calloc(1, sizeof(WGS_Coverage_Stats));
    stats_info->capture_cov_stats = calloc(user_inputs->num_of_target_files, sizeof(Capture_Coverage_Stats*));

    int i, j;
    for (i=0; i<=1000; i++) {
        stats_info->wgs_cov_stats->genome_cov_histogram[i]=0;
    }

    for (i=0; i<user_inputs->num_of_target_files; i++) {
        stats_info->capture_cov_stats[i] = calloc(1, sizeof(Capture_Coverage_Stats));
        stats_info->capture_cov_stats[i]->target_base_with_N_coverage = kh_init(m32);
        stats_info->capture_cov_stats[i]->target_coverage_for_median  = kh_init(m32);
        captureCoverageStatsInit(stats_info->capture_cov_stats[i]);

        for (j=0;j<=1000; j++) {
            stats_info->capture_cov_stats[i]->target_cov_histogram[j]=0;
        }
    }

    stats_info->wgs_cov_stats->genome_base_with_N_coverage = kh_init(m32);
    stats_info->wgs_cov_stats->genome_coverage_for_median  = kh_init(m32);

    readCoverageStatsInit(stats_info->read_cov_stats);
    WGSCoverageStatsInit(stats_info->wgs_cov_stats);   
}

void readCoverageStatsInit(Read_Coverage_Stats * read_cov_stats) {
    if (!read_cov_stats) {
        fprintf(stderr, "ERROR: Memory allocation failed for Read_Coverage_Stats\n");
        exit(EXIT_FAILURE);
    }

    read_cov_stats->total_reads_produced = 0;
    read_cov_stats->total_reads_aligned = 0;
    read_cov_stats->total_reads_paired = 0;
    read_cov_stats->total_reads_proper_paired = 0;
    read_cov_stats->total_duplicate_reads = 0;
    read_cov_stats->total_chimeric_reads = 0;
    read_cov_stats->total_supplementary_reads = 0;
    read_cov_stats->total_paired_reads_with_mapped_mates = 0;
    read_cov_stats->read_length = 0;
}

void WGSCoverageStatsInit(WGS_Coverage_Stats * wgs_cov_stats) {
    if (!wgs_cov_stats) {
        fprintf(stderr, "ERROR: Memory allocation failed for WGS Coverage Stats\n");
        exit(EXIT_FAILURE);
    }

    wgs_cov_stats->total_genome_bases = 0;
    wgs_cov_stats->total_Ns_bases = 0;
    wgs_cov_stats->total_Ns_bases_on_chrX = 0;
    wgs_cov_stats->total_Ns_bases_on_chrY = 0;
    wgs_cov_stats->total_mapped_bases = 0;
    wgs_cov_stats->total_uniquely_aligned_bases = 0;
    wgs_cov_stats->total_genome_coverage = 0;
    wgs_cov_stats->base_quality_20 = 0;
    wgs_cov_stats->base_quality_30 = 0;
    wgs_cov_stats->total_overlapped_bases = 0;

    wgs_cov_stats->wgs_max_coverage = 0;
    wgs_cov_stats->base_with_wgs_max_coverage = 0;
    wgs_cov_stats->median_genome_coverage = 0;

    wgs_cov_stats->mode = 0;
    wgs_cov_stats->uniformity_metric_all = 0.0;
    wgs_cov_stats->uniformity_metric_all_primary = 0.0;
    wgs_cov_stats->uniformity_metric_autosome_only = 0.0;
    wgs_cov_stats->uniformity_metric_primary_autosome_only = 0.0;
}

void captureCoverageStatsInit(Capture_Coverage_Stats * capture_cov_stats) {
    if (!capture_cov_stats) {
        fprintf(stderr, "ERROR: Memory allocation failed for Capture Coverage Stats\n");
        exit(EXIT_FAILURE);
    }

    capture_cov_stats->total_targeted_bases = 0;
    capture_cov_stats->total_target_coverage = 0;
    capture_cov_stats->total_buffer_bases = 0;
    capture_cov_stats->on_target_read_hit_count = 0;
    capture_cov_stats->off_target_read_hit_count = 0;
    capture_cov_stats->in_buffer_read_hit_count = 0;
    capture_cov_stats->hit_target_count = 0;
    capture_cov_stats->hit_target_buffer_only_count = 0;
    capture_cov_stats->non_target_good_hits = 0;

    capture_cov_stats->total_targets = 0;
    capture_cov_stats->target_max_coverage = 0;
    capture_cov_stats->base_with_target_max_coverage = 0;
    capture_cov_stats->median_target_coverage = 0;

    // initializing all numbers to 0
    //
    int i;
    for (i=0; i<PRIMER_SIZE; i++) {
        capture_cov_stats->five_prime[i] = 0;
        capture_cov_stats->three_prime[i] = 0;
    }

    for (i=0; i<101; i++)
        capture_cov_stats->target_coverage[i] = 0;
}

void statsInfoDestroy(Stats_Info *stats_info, User_Input *user_inputs) {
    if (stats_info->read_cov_stats) {
        free(stats_info->read_cov_stats);
        stats_info->read_cov_stats = NULL;
    }

    if (stats_info->wgs_cov_stats) {
        //kh_destroy(m32, stats_info->wgs_cov_stats->genome_cov_histogram);
        kh_destroy(m32, stats_info->wgs_cov_stats->genome_base_with_N_coverage);
        kh_destroy(m32, stats_info->wgs_cov_stats->genome_coverage_for_median);
        free(stats_info->wgs_cov_stats);
        stats_info->wgs_cov_stats = NULL;
    }

    if (stats_info->capture_cov_stats) {
        int i;
        for (i=0; i<user_inputs->num_of_target_files; i++) {
            kh_destroy(m32, stats_info->capture_cov_stats[i]->target_coverage_for_median);
            kh_destroy(m32, stats_info->capture_cov_stats[i]->target_base_with_N_coverage);
            free(stats_info->capture_cov_stats[i]);
        }       
        free(stats_info->capture_cov_stats);
        stats_info->capture_cov_stats=NULL;
    }   
    
    if (stats_info) { free(stats_info); stats_info=NULL; }
}

void addValueToKhashBucket16(khash_t(m16) *hash_in, uint16_t pos_key, uint16_t val) {
    int ret;
    khiter_t k_iter = kh_put(m16, hash_in, pos_key, &ret);
    if (ret == 1) {
        kh_value(hash_in, k_iter) = 0;
    } else if (ret == -1) {
        fprintf(stderr, "ERROR: can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        exit(EXIT_FAILURE);
    }

    kh_value(hash_in, k_iter) += val;

    return;
}

void addValueToKhashBucket32(khash_t(m32) *hash_in, uint32_t pos_key, uint32_t val) {
    int ret;
    khiter_t k_iter = kh_put(m32, hash_in, pos_key, &ret);
    if (ret == 1) {
        kh_value(hash_in, k_iter) = 0;
		//printf("add value is 0 %d ret, with key %"PRIu32"\n", ret, pos_key);
    } else if (ret == -1) {
        fprintf(stderr, "ERROR: can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        exit(EXIT_FAILURE);
    }

    kh_value(hash_in, k_iter) += val;
	//if (pos_key == 1 && kh_value(hash_in, k_iter)%1000000 == 0) printf("key 1 value is %"PRIu32"\n", kh_value(hash_in, k_iter));
	if (kh_value(hash_in, k_iter) > 4294967290) printf("larger value %"PRIu32"\n", kh_value(hash_in, k_iter));

    return;
}

uint32_t getValueFromKhash32(khash_t(m32) *hash32, uint32_t pos_key) {
    khiter_t k_iter;
	if (hash32 != NULL) {
		k_iter = kh_get(m32, hash32, pos_key);

		if (k_iter == kh_end(hash32))
			// this is needed as if the bucket is never used, the value will be whatever left there, and the value is undefined!
			//
			return 0;

        return kh_value(hash32, k_iter);
    }

	return 0;
}

void addValueToKhashBucketStrStr(khash_t(khStrStr) *hash_in, char *key, char * val) {
    int ret;
    khiter_t k_iter = kh_put(khStrStr, hash_in, key, &ret);
    if (ret == 1) {
        kh_key(hash_in, k_iter) = strdup(key);
        kh_value(hash_in, k_iter) = strdup(val);
        //strcpy(kh_value(hash_in, k_iter), val);
    } else if (ret == -1) {
        fprintf(stderr, "ERROR: can't find the key  %s\n", key);
        exit(EXIT_FAILURE);
    }
}

char * getValueFromKhashStrStr(khash_t(khStrStr) *hash_in, char* key) {
    khiter_t k_iter;
    if (hash_in != NULL) {
        k_iter = kh_get(khStrStr, hash_in, key);
        
        if (k_iter == kh_end(hash_in))
            // this is needed as if the bucket is never used, the value will be whatever left there, and the value is undefined!
            //
            return NULL;
        
        return kh_value(hash_in, k_iter);
    }
    
    return NULL;
}

uint16_t getValueFromKhash(khash_t(m16) *hash16, uint32_t pos_key) {
    khiter_t k_iter;
    if (hash16 != NULL) {
		k_iter = kh_get(m16, hash16, pos_key);

		if (k_iter == kh_end(hash16))
			return 0;

        return kh_value(hash16, k_iter);
    }

    return 0;
}

float calculatePercentage32(uint32_t num, uint32_t dom) {
	double val = (double)num/(double)dom;
	int i_val = val * 10000.0 + 0.5;
	return (float)i_val/100.0;
}

float calculatePercentage32_64(uint32_t num, uint64_t dom) {
	double val = (double)num/(double)dom;
	int i_val = val * 10000.0 + 0.5;
	return (float)i_val/100.0;
}

float calculatePercentage64(uint64_t num, uint64_t dom) {
	double val = (double)num/(double)dom;
	int i_val = val * 10000.0 + 0.5;
	return (float)i_val/100.0;
}

void combineCoverageStats(Stats_Info *stats_info, Stats_Info *tmp_stats_info, User_Input *user_inputs) {
	// For read length
    //
    if (stats_info->read_cov_stats->read_length == 0)
        stats_info->read_cov_stats->read_length = tmp_stats_info->read_cov_stats->read_length;

    copyReadCoverageStats(stats_info->read_cov_stats, tmp_stats_info->read_cov_stats);
    copyWGSCoverageStats(stats_info->wgs_cov_stats, tmp_stats_info->wgs_cov_stats);
    copyCaptureCoverageStats(stats_info->capture_cov_stats, tmp_stats_info->capture_cov_stats, user_inputs);

}

void copyReadCoverageStats(Read_Coverage_Stats *read_cov_stats, Read_Coverage_Stats *tmp_read_cov_stats) {
    fprintf(stderr, "%"PRIu64"\n", tmp_read_cov_stats->total_reads_produced);
	read_cov_stats->total_reads_produced  += tmp_read_cov_stats->total_reads_produced;
	read_cov_stats->total_reads_aligned   += tmp_read_cov_stats->total_reads_aligned;
	read_cov_stats->total_chimeric_reads  += tmp_read_cov_stats->total_chimeric_reads;
	read_cov_stats->total_duplicate_reads += tmp_read_cov_stats->total_duplicate_reads;
	read_cov_stats->total_reads_paired    += tmp_read_cov_stats->total_reads_paired;
	read_cov_stats->total_reads_proper_paired += tmp_read_cov_stats->total_reads_proper_paired;
	read_cov_stats->total_supplementary_reads += tmp_read_cov_stats->total_supplementary_reads;
	read_cov_stats->total_paired_reads_with_mapped_mates += tmp_read_cov_stats->total_paired_reads_with_mapped_mates;
}

void copyWGSCoverageStats(WGS_Coverage_Stats *wgs_cov_stats, WGS_Coverage_Stats *tmp_wgs_cov_stats) {
	// base related stats
	//
	wgs_cov_stats->base_quality_20        += tmp_wgs_cov_stats->base_quality_20;
	wgs_cov_stats->base_quality_30        += tmp_wgs_cov_stats->base_quality_30;
	wgs_cov_stats->total_mapped_bases     += tmp_wgs_cov_stats->total_mapped_bases;
	wgs_cov_stats->total_overlapped_bases += tmp_wgs_cov_stats->total_overlapped_bases;
	wgs_cov_stats->total_uniquely_aligned_bases += tmp_wgs_cov_stats->total_uniquely_aligned_bases;
}

void copyCaptureCoverageStats(Capture_Coverage_Stats **capture_cov_stats, Capture_Coverage_Stats **tmp_capture_cov_stats, User_Input *user_inputs) {
	// target (capture) related stats
	//
    int i;
    for (i=0; i<user_inputs->num_of_target_files; i++) {
	    capture_cov_stats[i]->on_target_read_hit_count  += tmp_capture_cov_stats[i]->on_target_read_hit_count;
    	capture_cov_stats[i]->in_buffer_read_hit_count  += tmp_capture_cov_stats[i]->in_buffer_read_hit_count;
	    capture_cov_stats[i]->off_target_read_hit_count += tmp_capture_cov_stats[i]->off_target_read_hit_count;
    }
}

// get the key for hash table for everything 1000 bases
//
uint32_t getHashKey(uint32_t position_in) {
	uint32_t tmp_key = (uint32_t) position_in / 1000;
	return tmp_key * 1000;
}

void lowCoverageGeneHashBucketKeyInit(khash_t(khStrLCG) *low_cov_gene_hash, char *key_in) {
	uint32_t i;
	int absent = 0;

	// need to find out if the key exists
	//
	khiter_t iter = kh_put(khStrLCG, low_cov_gene_hash, key_in, &absent);

	if (absent) {
		// key doesn't exists!
		//
		kh_key(low_cov_gene_hash, iter) = strdup(key_in);
		kh_value(low_cov_gene_hash, iter) = calloc(1, sizeof(Low_Coverage_Genes));
		kh_value(low_cov_gene_hash, iter)->total_size=0;
		kh_value(low_cov_gene_hash, iter)->capacity=5;
		kh_value(low_cov_gene_hash, iter)->gene_coverage = calloc(kh_value(low_cov_gene_hash, iter)->capacity, sizeof(Gene_Coverage));

		// initialize Gene_Coverage variables to NULL or 0
		//
		for (i=0; i<kh_value(low_cov_gene_hash, iter)->capacity; i++)
			geneCoverageInit(&kh_value(low_cov_gene_hash, iter)->gene_coverage[i]);
	}

	// need to expand the memory size for current bucket array if needed
	//
	if (kh_value(low_cov_gene_hash, iter)->capacity == kh_value(low_cov_gene_hash, iter)->total_size) {
		kh_value(low_cov_gene_hash, iter)->capacity = kh_value(low_cov_gene_hash, iter)->capacity * 2;
		kh_value(low_cov_gene_hash, iter)->gene_coverage =
			realloc(kh_value(low_cov_gene_hash, iter)->gene_coverage, kh_value(low_cov_gene_hash, iter)->capacity * sizeof(Gene_Coverage));

		if (kh_value(low_cov_gene_hash, iter)->gene_coverage == NULL) {
			fprintf(stderr, "ERROR: Memory re-allocation failed at lowCoverageGeneHashBucketKeyInit()!\n");
			exit(EXIT_FAILURE);
		}

		// initialize Gene_Coverage variables to NULL or 0
		//
		for (i=kh_value(low_cov_gene_hash, iter)->total_size; i<kh_value(low_cov_gene_hash, iter)->capacity; i++)
			geneCoverageInit(&kh_value(low_cov_gene_hash, iter)->gene_coverage[i]);
	}
}

void geneCoverageInit(Gene_Coverage *gc) {
	gc->gene_symbol = NULL;
	gc->transcript_name   = NULL;
	gc->low_cov_regions = NULL;

	gc->targeted = false;
	gc->num_of_low_cov_bases = 0;

	gc->cds_target_start = 0;
	gc->cds_target_end   = 0;
	gc->cds_length = 0;
	gc->cds_start  = 0;
	gc->cds_end = 0;
	gc->exon_id = 0;
	gc->exon_count = 0;
}

// copy everything from Gene_Coverage *gc1 to Gene_Coverage *gc2. 
// When first passed in, both gc1 and gc2 should be defined
//
void copyGeneCoverageLowCovRegions(Gene_Coverage* gc1, Gene_Coverage* gc2, bool copy_gene) {
	if (gc1 == NULL || gc2 == NULL) {
		fprintf(stderr, "ERROR: gc1 or gc2 shouldn't be NULL\n");
		exit(EXIT_FAILURE);
	}

	// now copy gene symbol and gene name
	//
	if (copy_gene) {
		// need to copy gene symbol and transcript_name (ie, transcript name) here
		//
		gc2->gene_symbol = calloc(strlen(gc1->gene_symbol)+1, sizeof(char));
		strcpy(gc2->gene_symbol, gc1->gene_symbol);

		gc2->transcript_name = calloc(strlen(gc1->transcript_name)+1, sizeof(char));
		strcpy(gc2->transcript_name, gc1->transcript_name);
	} else {
		gc2->gene_symbol = NULL;
		gc2->transcript_name = NULL;
	}

	// copy coordinates
	//
	gc2->cds_target_start = gc1->cds_target_start;
	gc2->cds_target_end   = gc1->cds_target_end;
	gc2->cds_start  = gc1->cds_start;
	gc2->cds_end    = gc1->cds_end;
	gc2->cds_length = gc1->cds_length;
	gc2->exon_count = gc1->exon_count;
	gc2->exon_id  = gc1->exon_id;
	gc2->targeted = gc1->targeted;

	// copy low_cov_regions from gc1 to gc2
	//
	int i;
	if (gc1->low_cov_regions != NULL) {
		if (gc2->low_cov_regions == NULL && gc1->low_cov_regions->size > 0) {
			gc2->low_cov_regions = calloc(1, sizeof(StringArray));
			gc2->low_cov_regions->size = gc1->low_cov_regions->size;
			gc2->low_cov_regions->capacity = gc1->low_cov_regions->size;
			gc2->low_cov_regions->theArray = calloc(gc2->low_cov_regions->size, sizeof(char*));

			for (i=0; i<gc1->low_cov_regions->size; i++) {
				gc2->low_cov_regions->theArray[i] = calloc(strlen(gc1->low_cov_regions->theArray[i])+1, sizeof(char));
				strcpy(gc2->low_cov_regions->theArray[i], gc1->low_cov_regions->theArray[i]);
			}
		} else {
			// this should never happen!
			//
			fprintf(stderr, "Something is wrong, you need to initialize Gene_Coverage variable first!\n");
			printf("ENTERING .....................................");

			uint16_t start_size = gc2->low_cov_regions->size;
			uint16_t end_size   = gc2->low_cov_regions->size + gc1->low_cov_regions->size;

			// append, but first we need to expand the memory allocation size
			//
			gc2->low_cov_regions->theArray = realloc(gc2->low_cov_regions->theArray, end_size * sizeof(char*));

			if (gc2->low_cov_regions->theArray == NULL) {
				fprintf(stderr, "ERROR: Memory re-allocation failed at copyGeneCoverageLowCovRegions()!\n");
				exit(EXIT_FAILURE);
			}

			int j=0;
			for (i=start_size; i<end_size; i++) {
				gc2->low_cov_regions->theArray[i] = calloc(strlen(gc1->low_cov_regions->theArray[j])+1, sizeof(char));
				strcpy(gc2->low_cov_regions->theArray[i], gc1->low_cov_regions->theArray[j]);
				j++;
			}

			gc2->low_cov_regions->size = end_size;
			gc2->low_cov_regions->capacity = end_size;
		}
	}
}

// if there is only one low coverage region, we just need to make a copy and return the copy
//
void getOneLowCovRegions(khash_t(khStrInt) *low_cov_regions_hash, StringArray *mergedArray) {
	mergedArray->theArray=calloc(1, sizeof(char*));                                                
	mergedArray->capacity = 1;                                                                     
	mergedArray->size = 0;

	khiter_t ph_iter;
	for (ph_iter=kh_begin(low_cov_regions_hash); ph_iter!=kh_end(low_cov_regions_hash); ++ph_iter) {
		if (kh_exist(low_cov_regions_hash, ph_iter)) {
			mergedArray->theArray[0] = calloc(strlen(kh_key(low_cov_regions_hash, ph_iter))+1, sizeof(char));            
			strcpy(mergedArray->theArray[0], kh_key(low_cov_regions_hash, ph_iter));
			mergedArray->size++;
		}
	}
}

void mergeLowCovRegions(khash_t(khStrInt) *low_cov_regions_hash, StringArray *mergedArray, uint32_t size_in,  uint32_t cds_t_start, uint32_t cds_t_end) {
	// first store everything into starts and ends hashs
	//
	uint32_t i, k;
	uint32_t *allNum = calloc(size_in*2, sizeof(uint32_t));
	khash_t(m32) *starts = kh_init(m32);
	khash_t(m32) *ends   = kh_init(m32);
	khiter_t k_iter, ph_iter;
	int ret;

	k=0;
	for (ph_iter=kh_begin(low_cov_regions_hash); ph_iter!=kh_end(low_cov_regions_hash); ++ph_iter) {
		if (kh_exist(low_cov_regions_hash, ph_iter)) {
			// split the string
			// example: "25316999-25317079"
			//
			char *region  = calloc(strlen(kh_key(low_cov_regions_hash, ph_iter))+1, sizeof(char));
			strcpy(region, kh_key(low_cov_regions_hash, ph_iter));
			char *savePtr = region;
			char *tokPtr;

			i=0;
			while ((tokPtr = strtok_r(savePtr, "-", &savePtr))) {
				if (i==0) {
					uint32_t beg = (uint32_t) strtol(tokPtr, NULL, 10);
					if (beg < cds_t_start)
						beg = cds_t_start;

					k_iter = kh_put(m32, starts, beg, &ret);
					if (ret)
						kh_value(starts, k_iter) = 0;	// initialize it to 0

					kh_value(starts, k_iter)++;
				}

				if (i==1) {
					uint32_t end = (uint32_t) strtol(tokPtr, NULL, 10);
					if (end > cds_t_end)
						end = cds_t_end;

					k_iter = kh_put(m32, ends, end, &ret);
					if (ret)
						kh_value(ends, k_iter) = 0;		// initialize it to 0

					kh_value(ends, k_iter)++;
				}

				allNum[k] = (uint32_t) strtol(tokPtr, NULL, 10);
				k++;
				i++;
			}

			if (region != NULL) free(region);
		}
	}

	// sort starts and ends array
	//
	qsort(allNum, k, sizeof(uint32_t), compare);

	// now walk through the allNum array and merge them accordingly
	//  s1      s2    s3               s4       s5    s6       s7  
	//    -------------------------------------------------------------------------------
	//              e2     e3     e1       e4                      e7     e6            e7
	//
	uint16_t flag=0;			// when flag==0, record the region
	uint32_t p_start=0, num_of_items=0;
	khash_t(m32) *low_cov_hash = kh_init(m32);

	for (i=0; i<k; i++) {
		if (flag == 0) {
			// record start position into hash
			//
			p_start = allNum[i];
			flag++;
			continue;
		} 

		// check to see if it is a start or end
		//
		k_iter = kh_get(m32, starts, allNum[i]);
		if (k_iter != kh_end(starts))
			flag++;

		k_iter = kh_get(m32, ends, allNum[i]);
		if (k_iter != kh_end(ends))
			flag--;

		if (flag == 0) {
			// record the region
			//
			k_iter = kh_put(m32, low_cov_hash, p_start, &ret); // add the key
			kh_value(low_cov_hash, k_iter) = allNum[i];
			num_of_items++;
		}
	}

	// now store them into the mergedArray
	//
	mergedArray->theArray=calloc(num_of_items, sizeof(char*));
	mergedArray->capacity = num_of_items;
	mergedArray->size = 0;
	for (i=0; i<num_of_items; i++) {
		mergedArray->theArray[i]=NULL;
	}

	for (k_iter = kh_begin(low_cov_hash); k_iter != kh_end(low_cov_hash); ++k_iter) {
		if (kh_exist(low_cov_hash, k_iter)) {
			char tmp_string[50];
			sprintf(tmp_string, "%"PRIu32"-%"PRIu32, kh_key(low_cov_hash, k_iter), kh_value(low_cov_hash, k_iter)); 
			mergedArray->theArray[mergedArray->size] = calloc(strlen(tmp_string)+1, sizeof(char));
			strcpy(mergedArray->theArray[mergedArray->size], tmp_string);
			mergedArray->size++;
		}
	}

	// clean-up
	//
	cleanKhashInt(starts);
	cleanKhashInt(ends);
	cleanKhashInt(low_cov_hash);
	free(allNum);
}

void calculateUniformityMetrics(Stats_Info *stats_info, User_Input *user_inputs, khash_t(khStrInt) *wanted_chromosome_hash, khash_t(m32) *cov_freq_dist, bool autosome, bool primary_chromosomes_only) {
	// need to set peak_size based on average coverage
	//
	set_peak_size_around_mode(stats_info, user_inputs);

	// initialize the hash table, so that the key will be in order
	//
	int i, absent=0;
	khiter_t iter;
	for (i=0; i<=1000; i++) {
		iter = kh_put(m32, cov_freq_dist, i, &absent);
		if (absent) {
			kh_key(cov_freq_dist, iter) = i;
			kh_value(cov_freq_dist, iter) = 0;
		}
	}

	uint64_t uniformity_total_bases = 0;

	// open uniformity data file for read
	//
	FILE *uniformity_fp = fopen(user_inputs->wgs_uniformity_file, "r");
	char *line = NULL;                                                                                        
	size_t len = 0;                                                                                           
	ssize_t read;                                                                                             
	char *tokPtr;
	char *chrom_id = NULL;

	while ((read = getline(&line, &len, uniformity_fp)) != -1) {
		// skip if it is comment line
		//
		if (strstr(line, "#") != NULL)
			continue;

		char *savePtr = line;                                                                                 
		uint8_t i=0;
		uint32_t tmp_cov=0, tmp_len=0;

		while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
			if (i==0)
                dynamicStringExpansion(tokPtr, &chrom_id);

			if (i==3)
				tmp_len = (uint32_t) strtol(tokPtr, NULL, 10);

			if (i==4) 
				tmp_cov = (uint32_t) strtol(tokPtr, NULL, 10);

			i++;
		}

		// skip sex chromosomes if we are only interested in autosomes
		//
		if (autosome == 1) {
			if ( (strcmp(chrom_id, "chrX") == 0) || (strcmp(chrom_id, "X") == 0) || (strcmp(chrom_id, "CHRX") == 0)
				 || (strcmp(chrom_id, "chrY") == 0) || (strcmp(chrom_id, "Y") == 0) || (strcmp(chrom_id, "CHRY") == 0) ) {
				continue;
			}
		}

		// if we are only handle primary chromosomes for this calculation, we need to check it here
		//
		if (primary_chromosomes_only && wanted_chromosome_hash) {
			khiter_t iter_p = kh_get(khStrInt, wanted_chromosome_hash, chrom_id);
			if (iter_p == kh_end(wanted_chromosome_hash))
				// chrom_id is not one of the primary chromosomes, so skip it!
				//
				continue;
		}

		if (tmp_cov > 1000) tmp_cov = 1000;
		iter = kh_put(m32, cov_freq_dist, tmp_cov, &absent);
		if (absent)
			fprintf(stderr, "Hash table initialization failed for current coverage: %d\n", tmp_cov);

		kh_value(cov_freq_dist, iter) += tmp_len;

		uniformity_total_bases += tmp_len;
	}

	if (line != NULL) free(line);
	fclose(uniformity_fp);
	free(chrom_id);

	// need to remove Ns_bases, (the number of Ns bases need to be calculated based on the user input)
	//
	uniformity_total_bases -= stats_info->wgs_cov_stats->total_Ns_bases; 

	// now walk through the hash table, find the mode 
	//
	uint32_t count_at_mode=0;
	uint64_t total_area_under_histogram=0;
	for (iter = kh_begin(cov_freq_dist); iter != kh_end(cov_freq_dist); ++iter) {
		if (kh_exist(cov_freq_dist, iter)) {

			total_area_under_histogram += kh_value(cov_freq_dist, iter);

			if (kh_key(cov_freq_dist, iter) == 0) {
				// remove Ns bases and continue
				//
				if (kh_value(cov_freq_dist, iter) > stats_info->wgs_cov_stats->total_Ns_bases) {
					kh_value(cov_freq_dist, iter) -= stats_info->wgs_cov_stats->total_Ns_bases;
				} else {
					fprintf(stderr, "The Ns region %"PRIu32" is larger than calculated one %"PRIu32"\n", stats_info->wgs_cov_stats->total_Ns_bases, kh_value(cov_freq_dist, iter));
				}

				continue;
			}

			if (kh_value(cov_freq_dist, iter) > count_at_mode) {
				stats_info->wgs_cov_stats->mode = kh_key(cov_freq_dist, iter);
				count_at_mode = kh_value(cov_freq_dist, iter);
			}
		}
	}

	// remove Ns
	//
	total_area_under_histogram -= stats_info->wgs_cov_stats->total_Ns_bases;

	// calculate the peak area under histogram
	//
	uint64_t peak_area_under_hist = dynamicCalculateAreaUnderHistogram(stats_info->wgs_cov_stats->mode, cov_freq_dist, user_inputs);
	if (autosome == 1) {
		// Here we need to adjust the total_area_under_histogram as it didn't add Ns for X and Y, but we subtract them anyway
		// Therefore, we need to add them back to make up the loss
		//
		total_area_under_histogram += stats_info->wgs_cov_stats->total_Ns_bases_on_chrX;
		total_area_under_histogram += stats_info->wgs_cov_stats->total_Ns_bases_on_chrY;

		if (primary_chromosomes_only) {
			printf("\nUniformity for primary autosome ONLY, without alt, decoys etc.\n");
			stats_info->wgs_cov_stats->uniformity_metric_primary_autosome_only = (double) peak_area_under_hist / (double) total_area_under_histogram;
		} else {
			printf("\nUniformity for autosome ONLY (including alt, decoy etc.)\n");
			stats_info->wgs_cov_stats->uniformity_metric_autosome_only = (double) peak_area_under_hist / (double) total_area_under_histogram;
		}
	} else {
		// for all chromosome
		//
		if (primary_chromosomes_only) {
			printf("\nUniformity for All Primary Chromosomes, including X and Y chromosomes. But without alt, decoys etc.\n");
			stats_info->wgs_cov_stats->uniformity_metric_all_primary = (double) peak_area_under_hist / (double) total_area_under_histogram;
		} else {
			printf("\nUniformity for All (including X, Y chromosomes, alt and decoys etc.\n");
			stats_info->wgs_cov_stats->uniformity_metric_all = (double) peak_area_under_hist / (double) total_area_under_histogram;
		}
	}

	printf("peak__area_under_histogram\t%"PRIu64"\n", peak_area_under_hist);
	printf("total_area_under_histogram\t%"PRIu64"\n", total_area_under_histogram);
	//printf("uniformity_total_bases is %"PRIu64"\n", uniformity_total_bases);
	
	// clean-up
	//
	//cleanKhashInt(cov_freq_dist);
}

uint64_t dynamicCalculateAreaUnderHistogram(uint32_t peak, khash_t(m32) *cov_freq_dist, User_Input *user_inputs) {

	uint64_t peak_area_under_histogram=0;
	uint8_t counter=0;

	khiter_t iter;
	iter = kh_get(m32, cov_freq_dist, peak);
	if (iter == kh_end(cov_freq_dist)) {
		fprintf(stderr, "ERROR: Peak hash value shouldn't be empty\n");
		exit(EXIT_FAILURE);
	}
	peak_area_under_histogram = kh_value(cov_freq_dist, iter);
	counter = 1;

	uint32_t left  = peak - 1;
	uint32_t right = peak + 1;

	while(1) {
		// the left side could go down to 0. if left side coverage is 0, we should skip it
		//
		uint32_t left_val = 0;
		//if (left > 0) {
			iter = kh_get(m32, cov_freq_dist, left);
			if (iter != kh_end(cov_freq_dist)) {
				left_val = kh_value(cov_freq_dist, iter);
			} else {
				//fprintf(stderr, "Peak left side has no more values\n");
				left_val = 0;
			}
		//}

		uint32_t right_val=0;
		iter = kh_get(m32, cov_freq_dist, right);
		if (iter != kh_end(cov_freq_dist)) {
			right_val = kh_value(cov_freq_dist, iter);
		} else {
			//fprintf(stderr, "Peak right side has no more values\n");
			right_val = 0;
		}

		if (left_val > right_val) {
			peak_area_under_histogram += left_val;
			counter++;
			left--;
		} else if (right_val > left_val) {
			peak_area_under_histogram += right_val;
			counter++;
			right++;
		} else {
			// they are equal
			//
			if ((user_inputs->size_of_peak_area - counter) >= 2) {
				peak_area_under_histogram += right_val * 2;
				left--;
				right++;
				counter += 2;
			} else {
				peak_area_under_histogram += right_val;
				right++;
				counter++;
			}
		}

		if (counter == user_inputs->size_of_peak_area)
			break;
	}

	return peak_area_under_histogram;
}

void set_peak_size_around_mode(Stats_Info *stats_info, User_Input *user_inputs) {
	uint64_t total_genome_non_Ns_bases = stats_info->wgs_cov_stats->total_genome_bases - stats_info->wgs_cov_stats->total_Ns_bases;
	double average_coverage = (double) stats_info->wgs_cov_stats->total_genome_coverage/ (double) total_genome_non_Ns_bases;

	if (!user_inputs->user_set_peak_size_on) {
		if (average_coverage <= 20) {
			user_inputs->size_of_peak_area = 3;
		} else if (average_coverage > 20 && average_coverage <= 25) {
			user_inputs->size_of_peak_area = 4;
		} else if (average_coverage > 25 && average_coverage <= 30) {
			user_inputs->size_of_peak_area = 5;
		} else if (average_coverage > 30 && average_coverage <= 35) {
			user_inputs->size_of_peak_area = 6;
		} else if (average_coverage > 35 && average_coverage <= 55) {
			user_inputs->size_of_peak_area = 7;
        } else if (average_coverage > 55 && average_coverage <= 70) {
            user_inputs->size_of_peak_area = 8;
        } else if (average_coverage > 70 && average_coverage <= 85) {
            user_inputs->size_of_peak_area = 9;
		} else if (average_coverage > 85 && average_coverage <= 100) {
			user_inputs->size_of_peak_area = 10;
		} else if (average_coverage > 100 && average_coverage <= 110) {
			user_inputs->size_of_peak_area = 11;
		} else if (average_coverage > 110 && average_coverage <= 120) {
			user_inputs->size_of_peak_area = 12;
		} else {
            user_inputs->size_of_peak_area = 13;
        }

        fprintf(stderr, "\tThe number of points selected around peak (eg, Mode) area based on average coverage is %d\n", user_inputs->size_of_peak_area);
	} else {
        fprintf(stderr, "\tThe number of points selected around peak (eg, Mode) area based on the user input is %d\n", user_inputs->size_of_peak_area);
    }
}

void outputFreqDistribution(User_Input *user_inputs, khash_t(m32) *cov_freq_dist) {
	// open WGS coverage summary report file handle
	//
	FILE *out_fp = fopen(user_inputs->wgs_cov_report, "a");
	fprintf(out_fp, "\n#Smoothed_Coverage_Frequency_Distribution_for_Whole_Genome\n");
	fprintf(out_fp, "==");
	khiter_t iter;
	for (iter=kh_begin(cov_freq_dist); iter!=kh_end(cov_freq_dist); iter++) {
		if (kh_exist(cov_freq_dist, iter))
			fprintf(out_fp, "%"PRIu32",", kh_value(cov_freq_dist, iter));
	}
	fprintf(out_fp, "\n");
	fclose(out_fp);
}

void printLowCoverageGeneStructure(Low_Coverage_Genes *low_cov_genes) {
	uint32_t i;
	printf("Total Number of Gene Symbol is %"PRIu32"\n", low_cov_genes->total_size);

	for (i=0; i<low_cov_genes->total_size; i++) {
		printf("Gene: %s\tRefSeq: %s\n", low_cov_genes->gene_coverage[i].gene_symbol, low_cov_genes->gene_coverage[i].transcript_name);
	}
}

// To view/print the content of string array before OR after sorting
// Note:  Any arrays will decay to a pointer to the first element when passing to a function.
// Therefore, we will have to pass the size info into the function to make it work!
//
void print_string_array(char** strings_in, size_t length_in) {
	size_t i;
	
	for (i=0; i<length_in; i++) {
		printf("%s\t", strings_in[i]);
	}

	printf("\n");
}
