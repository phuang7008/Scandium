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
int READLENGTH = 0;
int MIN_MAP_SCORE  = -1;
int MIN_BASE_SCORE = -1;
int NUM_OF_THREADS = 8;

uint32_t TOTAL_READS_PAIRED   = 0;
uint32_t TOTAL_ALIGNED_BASES  = 0;
uint32_t TOTAL_READS_ALIGNED  = 0;
uint32_t TOTAL_READS_PRODUCED = 0;
uint32_t TOTAL_PAIRED_READS_WITH_MAPPED_MATES = 0;
uint32_t TOTAL_DUPLICATE_READS = 0;
uint32_t TOTAL_SUPPLEMENTARY_ALIGNMENTS = 0;
bool REMOVE_DUPLICATES = false;
bool REMOVE_SUPPLEMENTARY_ALIGNMENTS = false;

float _PERCENTAGE = 1.0;    // only take this proportion of reads into consideration (for scale back experiments), 1 means all

char CURRENT_CHROMOSOME_ID[12];

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
	printf("\t-i <BAM/CRAM FILE alignment file (multiple files are not allowed!) >\n");
	printf("\t-m <minimun mapping quality score: to filter out any reads with mapQ less than -m>\n");
	printf("\t-b <base quality: to filter out any bases with baseQ less than -b>\n");
	printf("\t[-d] Remove Duplicates and DO NOT use them for statistics\n");
	printf("\t[-s] Remove Supplementary alignments and DO NOT use them for statistics\n");
	printf("\t[-w] Write whole genome coverage\n");
	printf("\t[-h] Print this usage message\n");
}

// some version of sequence file contain chromosome with id in 'chr1' format, as some software can't handle 'chr'.
// we need to remove them before processing them
char * removechr(char * c) {
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
bool isnum(const char * inStr) {
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
	int wgs_flag=0, bam_flag=-1, out_flag=-1, target_flag=-1;
	int arg;

	//When getopt returns -1, no more options available
	while ((arg = getopt(argc, argv, "b:di:m:n:o:p:st:wh")) != -1) {
		printf("current options %c and %s\n", arg, optarg);
		switch(arg) {
			case 'b':
				if (!isnum(optarg)) {
					fprintf (stderr, "Entered base quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
                MIN_BASE_SCORE = atoi(optarg);
                break;
            case 'd': REMOVE_DUPLICATES = true; break;
            case 'h': usage(); exit(1);
            case 'i':
                bam_flag = 1;
				user_inputs->bam_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));

                strcpy(user_inputs->bam_file, optarg);
                break;
            case 'm':
                if (!isnum(optarg)) {
                    fprintf (stderr, "Entered map quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
				MIN_MAP_SCORE = atoi(optarg);
                break;
			case 'n':
				if (!isnum(optarg)) {
                    fprintf (stderr, "Entered number of threads %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
				NUM_OF_THREADS = atoi(optarg);
            case 'o':
                out_flag = 1;
                user_inputs->out_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->out_file, optarg);
                break;
            case 'p': _PERCENTAGE=atof(optarg); break;
			case 's': REMOVE_SUPPLEMENTARY_ALIGNMENTS = true; break;
            case 't':
                target_flag = 1;
				user_inputs->target_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->target_file, optarg);
                break;
            case 'w': wgs_flag = 1; break;
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
    /*if (bam_flag == -1 || out_flag == -1 || target_flag == -1) {
        printf("-i\t-o\t-t options are mandatory!\n");
        return 1;
    }*/

	// Need to check out that all files user provided exist before proceeding
    if (target_flag == 1) checkFile(user_inputs->target_file);
    //if (bam_flag == 1) checkFile(user_inputs->bam_file);
    //if (out_flag == 1) checkFile(user_inputs->out_file);

	//char *cov_file = (char *) malloc((strlen(out_file)+15) * sizeof(char));
    //strcpy(cov_file, strcat(out_file, ".cov.fasta"));
    //bed_read(target_file);
}

User_Input * userInputInit() {
	User_Input * user_inputs = calloc(1, sizeof(User_Input));
	if (user_inputs == NULL) {
		fprintf(stderr, "Memory allocation failed in line %d!\n", __LINE__);
		exit(EXIT_FAILURE);
	}
	
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

	if (user_inputs != NULL)
		free(user_inputs);
}

void clean_khash_int(khash_t(m32) *hash_to_clean) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k)
		if (kh_exist(hash_to_clean, k))
			kh_del(m32, hash_to_clean, k);

	kh_destroy(m32, hash_to_clean);
	return;
}

void clean_khash_str(khash_t(str) *hash_to_clean) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k < kh_end(hash_to_clean); ++k) {
		if (kh_exist(hash_to_clean, k)) {
			// clean the inner khash if the key exist
			clean_khash_int(kh_value(hash_to_clean, k));
			free((char *) kh_key(hash_to_clean, k));
		}
	}

	kh_destroy(str, hash_to_clean);
	return;
}

/*
short get_chromosome_length_from_id(char *chrom_id) {
	if (isnum(chrom_id)) {
		return atoi(chrom_id) - 1;
	} else if (strcmp(chrom_id, "X") == 0) {
		return 22;
	} else if (strcmp(chrom_id, "Y") == 0) {
		return 23;
	} else if (strcmp(chrom_id, "M") == 0) {
		return 24;
	} else {
		fprintf(stderr, "chromosome id %s is not part of human genome\n", chrom_id);
		return -1;
	}
}*/

short get_chrom_index_from_id(bam_hdr_t *header, char *chrom_id) {
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
void chromosome_tracking_init(Chromosome_Tracking *chrom_tracking, char *chrom_id, uint32_t chrom_len, int index) {
	// number 25 is used as human has a total of 25 chromosomes including MT
	// if we are going to handle all other non-human contigs or decoy ones, we will have to expand the tracking list
	//
	if ( (index % 25) == 0) {

		unsigned short **c_coverage;
		c_coverage = index == 0 ? calloc(25, sizeof(unsigned short*)) : realloc(chrom_tracking->coverage, index * 2);
		if (!c_coverage) {
			fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->coverage");
			return;
		}
		chrom_tracking->coverage = c_coverage;

		char **c_ids;
		c_ids = index == 0 ? calloc(25, sizeof(char)) : realloc(chrom_tracking->chromosome_ids, index * 2);
		if (!c_ids) {
			fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_ids");
	        return;
		}
		chrom_tracking->chromosome_ids = c_ids;

		uint32_t *c_length;
		c_length = index == 0 ? calloc(25, sizeof(uint32_t)) : realloc(chrom_tracking->chromosome_lengths, index * 2);
		if (!c_length) {
			fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_lengths");
			return;
		}
		chrom_tracking->chromosome_lengths = c_length;

		uint8_t *c_status;
		c_status = index == 0 ? calloc(25, sizeof(uint8_t)) : realloc(chrom_tracking->chromosome_status, index * 2);
		if (!c_status) {
			fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_status");
			return;
		}
		chrom_tracking->chromosome_status = c_status;
	}

	chrom_tracking->chromosome_ids[index] = calloc(12, sizeof(char));
	strcpy(chrom_tracking->chromosome_ids[index], chrom_id);

	// As the 0 position will be empty as the position will be 1-based
	// So I used it to store the index information for quick access
	//
    chrom_tracking->coverage[index] = calloc(chrom_len + 5, sizeof(unsigned short));
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
	chrom_tracking->coverage[index][0] = index;

	chrom_tracking->chromosome_lengths[index] = chrom_len;
	chrom_tracking->chromosome_status[index] = 1;
	chrom_tracking->number_tracked++;
}

void chromosome_tracking_destroy(Chromosome_Tracking *chrom_tracking) {
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
