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
 *      Author:			Peiming (Peter) Huang
 *      Company:		Baylor College of Medicine
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>		// for file access() and getopt()

#include "utils.h"

// define (initialize) global variables declared at the terms.h file
bool N_FILE_PROVIDED = false;
bool TARGET_FILE_PROVIDED = false;

// Before open any file for processing, it is always a good idea to make sure that file exists
//
bool checkFile(char * fName) {
    if(access(fName, F_OK|R_OK) == -1) {
        fprintf(stderr, "No such file as %s;  File not found.\n", fName);
		exit(1);
    }
	return true;
}

// print out the help information to show the general usage of the package
//
void usage() {
	printf("Version %s\n\n", VERSION );
	printf("Usage:\t-o <output directory & file base name\n");
	printf("\t-t <target file>\n");
	printf("\t-b <base quality: to filter out any bases with baseQ less than -b>\n");
	printf("\t-i <BAM/CRAM FILE alignment file (multiple files are not allowed!) >\n");
	printf("\t-m <minimun mapping quality score: to filter out any reads with mapQ less than -m>\n");
	printf("\t-n <the file that contains regions of Ns in the reference genome in bed format>\n");
	printf("\t-T <the number of threads>\n");
	printf("\t[-d] Remove Duplicates and DO NOT use them for statistics\n");
	printf("\t[-s] Remove Supplementary alignments and DO NOT use them for statistics\n");
	printf("\t[-w] Write whole genome coverage\n");
	printf("\t[-h] Print this usage message\n");
}

// some version of sequence file contain chromosome with id in 'chr1' format, as some software can't handle 'chr'.
// we need to remove them before processing them
char * removeChr(char * c) {
	// first we need to see if the string contains 'chr'
	char *find = strstr(c, "chr");

	if (find != NULL) {
		size_t lenC = strlen(c);
		if (lenC == 4) {
			memmove(&c[0], &c[3],(4-3+1));
		} else {
			memmove(&c[0],&c[3],(5-3+1));
		}
	}

    return c;
}

// Calculate the percentage of input numbers
double getPercentage(int num, int dom)
{
    double pc = (double)num/(double)dom;
    pc*=10000.0;pc+=0.5; int ipc = (int)pc; pc = (double)ipc/100;
    return pc;
}

// This is used to check if a string (ie char *) is a number
bool isNumber(const char * inStr) {
    while( *inStr != '\0') {
        if (!isdigit(inStr[0])) {
            return false;
        }
        inStr++;
    }
    return true;
}

// Get command line arguments in and check the sanity of user inputs 
//
void processUserOptions(User_Input *user_inputs, int argc, char *argv[]) {
	int arg;

	//When getopt returns -1, no more options available
	while ((arg = getopt(argc, argv, "b:di:m:n:o:p:st:T:wh")) != -1) {
		printf("current options %c and %s\n", arg, optarg);
		switch(arg) {
			case 'b':
				if (!isNumber(optarg)) {
					fprintf (stderr, "Entered base quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
                user_inputs->min_base_quality = atoi(optarg);
                break;
            case 'd': user_inputs->remove_duplicate = true; break;
            case 'h': usage(); exit(1);
            case 'i':
				user_inputs->bam_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));

                strcpy(user_inputs->bam_file, optarg);
                break;
            case 'm':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered map quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
				user_inputs->min_map_quality = atoi(optarg);
                break;
			case 'n':
				N_FILE_PROVIDED = true;
				user_inputs->n_file = malloc(strlen(optarg)+1 * sizeof(char));
                strcpy(user_inputs->n_file, optarg);
            case 'o':
                user_inputs->out_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->out_file, optarg);
                break;
            case 'p': user_inputs->percentage = atof(optarg); break;
			case 's': user_inputs->remove_supplementary_alignments = true; break;
            case 't':
				TARGET_FILE_PROVIDED = true;
				user_inputs->target_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->target_file, optarg);
                break;
			case 'T':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered number of threads %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
                user_inputs->num_of_threads = atoi(optarg);
            case 'w': user_inputs->wgs_coverage = true; break;
            case '?':
                if (optopt == 'b' || optopt == 'i' || optopt == 'm' || optopt == 'o' || optopt == 't')
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                usage();
                exit(1);
            default: fprintf(stderr, "Non-option argument %c\n", optopt); usage(); exit(1);
        }
    }

	// check the mandatory arguments (will turn this on for the final test/run)
    if (user_inputs->bam_file == NULL) {
        printf("-i\toption is mandatory!\n");
        exit(1);
    }
	
	if (user_inputs->out_file == NULL) {
		printf("-o\toption is mandatory!\n");
		exit(1);
	}

	// Need to check out that all files user provided exist before proceeding
    if (user_inputs->bam_file) checkFile(user_inputs->bam_file);
    if (N_FILE_PROVIDED) checkFile(user_inputs->n_file);
    if (TARGET_FILE_PROVIDED) checkFile(user_inputs->target_file);

	// for cov.fasta file name
	user_inputs->cov_file = (char *) malloc((strlen(user_inputs->bam_file)+15) * sizeof(char));
	strcpy(user_inputs->cov_file, user_inputs->bam_file);
    strcat(user_inputs->cov_file, ".cov.fasta");

	// for wig.fasta file name
	user_inputs->wig_file = (char *) malloc((strlen(user_inputs->bam_file)+15) * sizeof(char));
    strcpy(user_inputs->wig_file, user_inputs->bam_file);
    strcat(user_inputs->wig_file, ".wig.fasta");

	// for whole genome (wgs) file name
	if (user_inputs->wgs_coverage) {
		user_inputs->wgs_file = malloc((strlen(user_inputs->bam_file)+15) * sizeof(char));
		strcpy(user_inputs->wgs_file, user_inputs->bam_file);
		strcat(user_inputs->wgs_file, ".wgs.fasta");
	}
}

User_Input * userInputInit() {
	User_Input * user_inputs = calloc(1, sizeof(User_Input));
	if (user_inputs == NULL) {
		fprintf(stderr, "Memory allocation failed in line %d!\n", __LINE__);
		exit(EXIT_FAILURE);
	}

	user_inputs->min_map_quality  = -1;
	user_inputs->min_base_quality = -1;
	user_inputs->num_of_threads   = 4;
	user_inputs->percentage = 1.0;
	user_inputs->wgs_coverage = false;
	user_inputs->remove_duplicate = false;
	user_inputs->remove_supplementary_alignments = false;
	
	return user_inputs;
}

void userInputDestroy(User_Input *user_inputs) {
	if (user_inputs->target_file != NULL)
		free(user_inputs->target_file);

	if (user_inputs->bam_file != NULL)
		free(user_inputs->bam_file);

	if (user_inputs->cov_file != NULL) 
		free(user_inputs->cov_file);

	if (user_inputs->out_file != NULL)
		free(user_inputs->out_file);

	if (user_inputs->n_file != NULL)
		free(user_inputs->n_file);

	if (user_inputs->wig_file != NULL)
		free(user_inputs->wig_file);

	if (user_inputs->wgs_file != NULL)
		free(user_inputs->wgs_file);
}

void cleanKhashInt(khash_t(m32) *hash_to_clean) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k)
		if (kh_exist(hash_to_clean, k))
			kh_del(m32, hash_to_clean, k);

	kh_destroy(m32, hash_to_clean);
	return;
}

void cleanKhashStr(khash_t(str) *hash_to_clean) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
		if (kh_exist(hash_to_clean, k)) {
			// clean the inner khash if the key exist
			free((char *) kh_key(hash_to_clean, k));
			cleanKhashInt(kh_value(hash_to_clean, k));
			//free(kh_value(hash_to_clean, k));		// already cleaned by cleanKhashInt()
		}
	}

	kh_destroy(str, hash_to_clean);
	return;
}

short getChromIndexFromID(bam_hdr_t *header, char *chrom_id) {
	int i;
	for(i=0; i<header->n_targets; i++) {
		if (strcmp(header->target_name[i], chrom_id) == 0) {
			return i;
		}
	}

	fprintf(stderr, "Can't locate the chromosome name %s at the function get_chrom_index_from_id\n", chrom_id);
	exit(1);
}

// Note: this function should never be used to update any information regarding the Chromosome_Tracking variable
// It is only used for initialization and dynamically allocate memories!
//
void chromosomeTrackingInit(Chromosome_Tracking *chrom_tracking, char *chrom_id, uint32_t chrom_len, int index) {
	// number 25 is used as human has a total of 25 chromosomes including MT
	// if we are going to handle all other non-human contigs or decoy ones, we will have to expand the tracking list
	//
	if ( (index % 25) == 0) {

		uint16_t **c_coverage;
		c_coverage = index == 0 ? calloc(25, sizeof(uint16_t*)) : realloc(chrom_tracking->coverage, index * 2);
		if (!c_coverage) {
			fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->coverage");
			exit(1);
		}
		chrom_tracking->coverage = c_coverage;

		char **c_ids;
		c_ids = index == 0 ? calloc(25, sizeof(char)) : realloc(chrom_tracking->chromosome_ids, index * 2);
		if (!c_ids) {
			fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_ids");
	        exit(1);
		}
		chrom_tracking->chromosome_ids = c_ids;

		uint32_t *c_length;
		c_length = index == 0 ? calloc(25, sizeof(uint32_t)) : realloc(chrom_tracking->chromosome_lengths, index * 2);
		if (!c_length) {
			fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_lengths");
			exit(1);
		}
		chrom_tracking->chromosome_lengths = c_length;

		uint8_t *c_status;
		c_status = index == 0 ? calloc(25, sizeof(uint8_t)) : realloc(chrom_tracking->chromosome_status, index * 2);
		if (!c_status) {
			fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_status");
			exit(1);
		}
		chrom_tracking->chromosome_status = c_status;
	}

	//chrom_tracking->chromosome_ids[index] = calloc(12, sizeof(char));
	//strcpy(chrom_tracking->chromosome_ids[index], chrom_id);
	chrom_tracking->chromosome_ids[index] = strdup(chrom_id);
	if (chrom_tracking->chromosome_ids[index] == NULL) {
		printf("Allocation failed for chrom %s\n", chrom_id);
	} else {
		printf("Allocated chromosome id is %s\n", chrom_tracking->chromosome_ids[index]);
	}

	// As the 0 position will be empty as the position will be 1-based
	// So I used it to store the index information for quick access
	// Also, I need to add 1 to the size to align with the 1-based position
	//
    chrom_tracking->coverage[index] = calloc(chrom_len + 1, sizeof(uint16_t));
	if (!chrom_tracking->coverage[index]) {
		printf("allocation of memory failed for chromosome tracking array number %d\n", index);
		exit(1);
	}

	// need to initialization all of them to 0 (it is supposed to be done by the calloc(), 
	//uint32_t i;
	//for (i=0; i<chrom_len+5; i++) {
	//	printf("before initialization the values are %"PRIu32"\t%d\n", i, chrom_tracking->coverage[index][i]);
		//chrom_tracking->coverage[index][i] = 0;
	//}
	//chrom_tracking->coverage[index][0] = index;

	chrom_tracking->chromosome_lengths[index] = chrom_len;
	chrom_tracking->chromosome_status[index] = 1;
	chrom_tracking->number_tracked++;
}

void chromosomeTrackingDestroy(Chromosome_Tracking *chrom_tracking) {
	int idx = chrom_tracking->number_tracked;
	int i = 0;
	for (i=0; i<idx; i++) {
		if (chrom_tracking->chromosome_ids[i])
			free(chrom_tracking->chromosome_ids[i]);

		if (chrom_tracking->coverage[i])
			free(chrom_tracking->coverage[i]);
	}

	free(chrom_tracking->chromosome_ids);
	free(chrom_tracking->chromosome_lengths);
	free(chrom_tracking->chromosome_status);
	free(chrom_tracking->coverage);

}

uint8_t locateChromosomeIndex(char *chrom_id, Chromosome_Tracking *chrom_tracking) {
	int i;
    for (i = 0; i < chrom_tracking->number_tracked; i++) {
        if (strcmp(chrom_id, chrom_tracking->chromosome_ids[i]) == 0) {
			return i;
        }
    }

	fprintf(stderr, "Something is wrong because the chromosome %s couldn't be found\n", chrom_id);
	return 0;
}

void statsInfoInit(Stats_Info *stats_info) {
	stats_info->target_cov_histogram = kh_init(m16);
    stats_info->genome_cov_histogram = kh_init(m16);
    stats_info->targeted_base_with_N_coverage = kh_init(m16);
    stats_info->genome_base_with_N_coverage   = kh_init(m16);
    stats_info->target_coverage_for_median = kh_init(m16);
    stats_info->genome_coverage_for_median = kh_init(m16);
    memset(stats_info->target_coverage, 0, sizeof(stats_info->target_coverage));
    memset(stats_info->five_prime, 0, sizeof(stats_info->five_prime));
    memset(stats_info->three_prime, 0, sizeof(stats_info->three_prime));

	stats_info->total_genome_bases = 0;
	stats_info->total_buffer_bases = 0;
	stats_info->total_targeted_bases = 0;
	stats_info->total_aligned_bases = 0;
	stats_info->total_target_coverage = 0;
	stats_info->total_genome_coverage = 0;

	stats_info->total_reads_paired = 0;
	stats_info->total_reads_aligned = 0;
	stats_info->total_reads_produced = 0;
	stats_info->total_duplicate_reads = 0;
	stats_info->total_supplementary_reads = 0;
	stats_info->total_paired_reads_with_mapped_mates = 0;

	stats_info->on_target_read_hit_count = 0;
	stats_info->off_target_read_hit_count = 0;
	stats_info->in_buffer_read_hit_count = 0;
	stats_info->hit_target_count = 0;
	stats_info->hit_buffer_only_count = 0;

	stats_info->total_targets = 0;
	stats_info->read_length = 0;
	stats_info->max_coverage = 0;
	stats_info->base_with_max_coverage = 0;
	stats_info->median_genome_coverage = 0;
	stats_info->median_target_coverage = 0;
}

void statsInfoDestroy(Stats_Info *stats_info) {
	kh_destroy(m16, stats_info->target_cov_histogram);
	kh_destroy(m16, stats_info->genome_cov_histogram);
	kh_destroy(m16, stats_info->targeted_base_with_N_coverage);
	kh_destroy(m16, stats_info->genome_base_with_N_coverage);
	kh_destroy(m16, stats_info->target_coverage_for_median);
	kh_destroy(m16, stats_info->genome_coverage_for_median);

	free(stats_info->target_coverage);
	free(stats_info->five_prime);
	free(stats_info->three_prime);

}

void zeroAllNsRegions(char *chrom_id, khash_t(str) *Ns_buffer_hash, Chromosome_Tracking *chrom_tracking) {
    // First, we need to find the index that is used to track current chromosome chrom_id
    uint8_t idx = locateChromosomeIndex(chrom_id, chrom_tracking);

    khiter_t outer_iter, inner_iter;
    for (outer_iter=kh_begin(Ns_buffer_hash); outer_iter!=kh_end(Ns_buffer_hash); outer_iter++) {
        if (kh_exist(Ns_buffer_hash, outer_iter)) {
            if (strcmp(chrom_id, kh_key(Ns_buffer_hash, outer_iter)) == 0) {
                // found the chromosome id key in the hash map
                for (inner_iter=kh_begin(kh_value(Ns_buffer_hash, outer_iter));
                        inner_iter!=kh_end(kh_value(Ns_buffer_hash, outer_iter));
                            inner_iter++) {
                    if (kh_exist(kh_value(Ns_buffer_hash, outer_iter), inner_iter)) {
                        uint32_t pos = kh_key(kh_value(Ns_buffer_hash, outer_iter), inner_iter);

                        // set the coverage variable at the pos in chrom_tracking variable to 0
                        chrom_tracking->coverage[idx][pos] = 0;
                    }
                }
            }
        }
    }
}

void addValueToKhashBucket(khash_t(m16) *hash_in, uint16_t pos_key, uint16_t val) {
    int ret;
    khiter_t k_iter = kh_put(m16, hash_in, pos_key, &ret);
    if (ret == 1) {
        kh_value(hash_in, k_iter) = 0;
    } else if (ret == -1) {
        fprintf(stderr, "can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        exit(1);
    }

    kh_value(hash_in, k_iter) += val;

    return;
}

void addValueToKhashBucket32(khash_t(m32) *hash_in, uint32_t pos_key, uint16_t val) {
    int ret;
    khiter_t k_iter = kh_put(m32, hash_in, pos_key, &ret);
    if (ret == 1) {
        kh_value(hash_in, k_iter) = 0;
    } else if (ret == -1) {
        fprintf(stderr, "can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        exit(1);
    }

    kh_value(hash_in, k_iter) += val;

    return;
}

uint16_t getValueFromKhash(khash_t(m16) *hash16, khash_t(m32) *hash32, uint16_t pos_key) {
	int ret;
    khiter_t k_iter;
	if (hash16 != NULL) k_iter = kh_put(m16, hash16, pos_key, &ret);
	if (hash32 != NULL) k_iter = kh_put(m32, hash32, pos_key, &ret);

    if (ret == 1) {
        if (hash16 != NULL) {
			kh_value(hash16, k_iter) = 0;
			return kh_value(hash16, k_iter);
		}

		if (hash32 != NULL) {
            kh_value(hash32, k_iter) = 0;
            return kh_value(hash32, k_iter);
        }
        if (hash32 != NULL) kh_value(hash32, k_iter) = 0;
    } else if (ret == -1) {
        fprintf(stderr, "can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        exit(1);
    }

	return -1;
}

float calculatePercentage(uint32_t num, uint32_t dom) {
	float val = (double)num/(double)dom;
	int i_val = val * 10000.0 + 0.5;
	return (float)i_val/100.0;
}
