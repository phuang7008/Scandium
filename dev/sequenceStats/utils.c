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
int MIN_MAP_SCORE  = -1;
int MIN_BASE_SCORE = -1;
int NUM_OF_THREADS = 8;

long TOTAL_READS_PAIRED   = 0;
long TOTAL_ALIGNED_BASES  = 0;
long TOTAL_READS_ALIGNED  = 0;
long TOTAL_READS_PRODUCED = 0;
long TOTAL_PAIRED_READS_WITH_MAPPED_MATES = 0;

long DUPLICATE_READS = 0;
bool REMOVE_DUPLICATES = false;

float _PERCENTAGE = 1.0;    // only take this proportion of reads into consideration (for scale back experiments), 1 means all

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
	printf("\t-i <BAM FILE alignment file (multiple files are not allowed!) >\n");
	printf("\t-m <minimun mapping quality score: to filter out any reads with mapQ less than -m>\n");
	printf("\t-b <base quality: to filter out any bases with baseQ less than -b>\n");
	printf("\t[-d] Remove Duplicates and DO NOT use them for statistics\n");
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
	while ((arg = getopt(argc, argv, "b:di:m:n:o:p:t:wh")) != -1) {
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

void clean_khash(khash_t(32) *hash_to_clean) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k)
		if (kh_exist(hash_to_clean, k))
			kh_del(32, hash_to_clean, k);
}
