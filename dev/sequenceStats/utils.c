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
#include "annotation.h"

// define (initialize) global variables declared at the terms.h file
//
bool N_FILE_PROVIDED = false;
bool HGMD_PROVIDED   = false;
bool TARGET_FILE_PROVIDED  = false;
bool USER_DEFINED_DATABASE = false;

// Before open any file for processing, it is always a good idea to make sure that file exists
//
bool checkFile(char * fName) {
    if(access(fName, F_OK|R_OK) == -1) {
        fprintf(stderr, "No such file as %s;  File not found.\n", fName);
		exit(EXIT_FAILURE);
    }
	return true;
}

uint64_t check_file_size(const char *filename) {
	struct stat st;
	if (stat(filename, &st) == 0)
		return st.st_size;

	fprintf(stderr, "Something is wrong when check the filesize for file: %s\n", filename);
	exit(EXIT_FAILURE);
}

void recordHGMD(Databases *dbs, User_Input *user_inputs, khash_t(khStrInt) *hgmd_genes, khash_t(khStrInt) *hgmd_transcripts) {
	// query from the MySQL database
	//
	char *sql = calloc(100, sizeof(char));
	sprintf(sql, "SELECT gene_symbol, transcript FROM %s", dbs->db_hgmd);

	if (mysql_query(dbs->con,sql))
		finish_with_error(dbs->con);

	dbs->mysql_results = mysql_store_result(dbs->con);
	if (dbs->mysql_results == NULL)
		finish_with_error(dbs->con);

	MYSQL_ROW row;
	while ((row = mysql_fetch_row(dbs->mysql_results))) {
		int absent;

		// check to see if the current gene exists
		//
		khiter_t iter = kh_put(khStrInt, hgmd_genes, row[0], &absent);
		if (absent) {
			kh_key(hgmd_genes, iter) = strdup(row[0]);
		}
		kh_value(hgmd_genes, iter) = 1;

		// check to see if current transcripts exists
		//
		int i=0;
		char *savePtr = row[1];
		char *tokPtr;
		while ((tokPtr = strtok_r(savePtr, ".", &savePtr))) {
			if (i==0) {
				iter = kh_put(khStrInt, hgmd_transcripts, tokPtr, &absent);
				if (absent) {
					kh_key(hgmd_transcripts, iter) = strdup(tokPtr);
				}
				kh_value(hgmd_transcripts, iter) = 1;
			}
			i++;
		}
	}

	if (dbs->mysql_results) { 
		mysql_free_result(dbs->mysql_results);
		dbs->mysql_results = NULL;
	}
	free(sql);
}

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

void stringArrayDestroy(stringArray *arrayIn) {
    uint16_t i;
    for(i=0; i<arrayIn->size; i++) {
        if (arrayIn->theArray[i]) {
            free(arrayIn->theArray[i]);
        }
    }

    if (arrayIn->theArray != NULL)
        free(arrayIn->theArray);
}

// I found this on the internet
//
char* stristr( const char* str1, const char* str2 ) {
    const char* p1 = str1 ;
    const char* p2 = str2 ;
    const char* r = *p2 == 0 ? str1 : 0 ;

    while ( *p1 != 0 && *p2 != 0 ) {
        if ( tolower( (unsigned char)*p1 ) == tolower( (unsigned char)*p2 ) ) {
            if ( r == 0 ) {
                r = p1 ;
            }

            p2++ ;
        } else {
            p2 = str2 ;
            if ( r != 0 ) {
                p1 = r + 1 ;
            }

            if ( tolower( (unsigned char)*p1 ) == tolower( (unsigned char)*p2 ) ) {
                r = p1 ;
                p2++ ;
            } else {
                r = 0 ;
            }
        }

        p1++;
    }

    return *p2 == 0 ? (char*)r : 0 ;
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
			uint32_t start, end;
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
							fprintf(stderr, "Memory re-allocation failed at processLowCovRegionFromKhash()!\n");
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
				fprintf(stderr, "Memory re-allocation failed at the processLowCovRegionFromKhash()!\n");
				exit(EXIT_FAILURE);
			}
		}

		*output[0] = '\0';	// set to the null terminator, so we could use strcat() all the way

		for (i=0; i<size_L; i++) {
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

uint32_t processLowCovRegionFromStrArray(stringArray *merged_low_cov_regions, char **output) {
	uint32_t low_cov_region_size=0, stringLength=0;
	uint32_t i;

	for (i=0; i<merged_low_cov_regions->size; i++) {
		stringLength += strlen(merged_low_cov_regions->theArray[i]) + 1;
	}

	*output = realloc(*output, (stringLength + 50) * sizeof(char));
	if (*output == NULL) {
		fprintf(stderr, "Memory re-allocation failed at processLowCovRegionFromStrArray()\n");
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
		uint32_t start, end, k=0;

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
		kh_value(gene_transcripts, iter) = calloc(1, sizeof(stringArray));
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
				fprintf(stderr, "Memory re-allocation failed at addToGeneTranscriptKhashTable()!\n");
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

void annotationWrapperDestroy(Annotation_Wrapper *annotation_wrapper) {
	uint16_t i;

	// Free the inner most part first
	//
	for (i=0; i<annotation_wrapper->real_size; i++) {
		if (annotation_wrapper->annotations[i].gene != NULL)
			free(annotation_wrapper->annotations[i].gene);

		if (annotation_wrapper->annotations[i].Synonymous != NULL)
			free(annotation_wrapper->annotations[i].Synonymous);

		//free(annotation_wrapper->annotations[i].prev_genes);

		if (annotation_wrapper->annotations[i].exon_info != NULL)
			free(annotation_wrapper->annotations[i].exon_info);
	}

	// free the next level
	//
	if (annotation_wrapper->annotations != NULL)
		free(annotation_wrapper->annotations);

	if (annotation_wrapper != NULL)
		free(annotation_wrapper);
}

// print out the help information to show the general usage of the package
//
void usage() {
	printf("Version %s\n\n", VERSION_ );
	printf("Usage:  coverage -i bam/cram -o output_dir [options ...]\n");
	printf("Note: this is a multi-threading program. Each thread need 3gb of memory. So please allocate them accordingly!\n");
	printf("Note: for example: 3 threads would use 8-9gb of memory, while 4 threads would need 12 gb of memory, etc.\n\n");
	printf("Mandatory:\n");
	printf("\t-i <BAM/CRAM alignment file (multiple files are not allowed!). It Is Mandatory >\n");
	printf("\t-o <output directory. It Is Mandatory>\n\n");

	printf("The Followings Are Optional:\n");
	printf("\t-b <minimal base quality: to filter out any bases with baseQ less than b. Default 0>\n");
	printf("\t-f <the file that contains user defined database (for annotation only)>\n");
	printf("\t-g <the percentage used for gVCF blocking: Default 10 for 1000%%>\n");
	printf("\t-m <minimal mapping quality score: to filter out any reads with mapQ less than m. Default 0>\n");
	printf("\t-n <the file that contains regions of Ns in the reference genome in bed format>\n");
	printf("\t-p <the percentage (fraction) of reads used for this analysis. Default 1.0 (ie, 100%%)>\n");
	printf("\t-t <target file. If this is specified, all of the output file names related to this will contain .Capture_>\n");

	printf("\t-B <the Buffer size immediate adjacent to a target region. Default: 100>\n");
	printf("\t-D <the version of human genome database (either hg19 [or hg37], or hg38). Default:hg19>\n");
	printf("\t-H <the high coverage cutoff value. Any coverages larger than it will be outputted. Default=10000>\n");
	printf("\t-L <the low coverage cutoff value. Any coverages smaller than it will be outputted. Default=20>\n");
	printf("\t-T <the number of threads (Note: when used with HPC's msub, make sure number of processors:ppn matches to number of threads). Default 3>\n");

	printf("The Followings are for the block region output (used for Uniformity analysis)\n");
	printf("\t-l <the lower bound for the block region output. Default: 1>\n");
	printf("\t-u <the upper bound for the block region output. Default: 150>\n\n");

	printf("The Followings Are Flags\n");
	printf("\t[-a] Turn off the annotation information for genes, exons and transcript. Default: ON\n");
	printf("\t[-d] Remove Duplicates is OFF! Specify this flag only when you want to keep Duplicates reads. Default: ON\n");
	printf("\t[-s] Remove Supplementary alignments and DO NOT use them for statistics. Default: off\n");
	printf("\t[-w] Write whole genome coverage related reports (all of the output file names related to this will have .WGS_ in them). This flag doesn't produce the WGS Coverage.fasta file, use -W for that. Default: off\n");

	printf("\t[-G] Write/Dump the WIG formatted file. Default: off\n");
	printf("\t[-M] Use HGMD annotation. Default: off\n");
	printf("\t[-W] Write/Dump the WGS Coverage.fasta file (both -w and -W needed). Default: off\n");
	printf("\t[-h] Print this help/usage message\n");
}

// some version of sequence file contain chromosome with id in 'chr1' format, as some software can't handle 'chr'.
// we need to remove them before processing them
char * removeChr(char * c) {
	// first we need to see if the string contains 'chr'
	char *find = strstr(c, "chr");

	if (find != NULL) {
		size_t lenC = strlen(c);
		memmove(&c[0], &c[3],(lenC-3+1));
	}

    return c;
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
	//
	while ((arg = getopt(argc, argv, "ab:B:c:dD:f:g:GH:i:L:l:m:Mn:o:p:st:T:u:wWy:h")) != -1) {
		//printf("User options for %c is %s\n", arg, optarg);
		switch(arg) {
			case 'a':
				user_inputs->annotation_on = false; break;
			case 'b':
				if (!isNumber(optarg)) {
					fprintf (stderr, "Entered base quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                user_inputs->min_base_quality = atoi(optarg);
                break;
			case 'B': user_inputs->target_buffer_size = atoi(optarg); break;
            case 'd': user_inputs->remove_duplicate = false; break;
            case 'D': 
				strcpy(user_inputs->database_version, optarg); 

				// change all to lower case
				int i;
				for(i = 0; user_inputs->database_version[i]; i++){
					user_inputs->database_version[i] = tolower(user_inputs->database_version[i]);
				}

				if (strcmp(user_inputs->database_version, "hg37") == 0)
					strcpy(user_inputs->database_version, "hg19");

				break;
			case 'f':
				USER_DEFINED_DATABASE = true;
				user_inputs->user_defined_database_file = (char*) malloc((strlen(optarg)+1) * sizeof(char));
				strcpy(user_inputs->user_defined_database_file, optarg);
				break;
			case 'g': user_inputs->gVCF_percentage = atoi(optarg); break;
            case 'G': user_inputs->Write_WIG = true; break;
            case 'h': usage(); exit(EXIT_FAILURE);
			case 'H':
				if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered High coverage cutoff value %s is not a number\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                user_inputs->high_coverage_to_report = atoi(optarg);
                break;
            case 'i':
				user_inputs->bam_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));

                strcpy(user_inputs->bam_file, optarg);
                break;
			case 'L':
				if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered Lower coverage cutoff value %s is not a number\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                user_inputs->low_coverage_to_report = atoi(optarg);
                break;
			case 'l':
				if (!isNumber(optarg)) {
					fprintf (stderr, "Entered lower_bound value %s is not a number\n", optarg);
					usage();
					exit(EXIT_FAILURE);
				}
				user_inputs->lower_bound = atoi(optarg);
				break;
            case 'm':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered map quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
				user_inputs->min_map_quality = atoi(optarg);
                break;
			case 'M':
				HGMD_PROVIDED = true;
				break;
			case 'n':
				N_FILE_PROVIDED = true;
				user_inputs->n_file = malloc(strlen(optarg)+1 * sizeof(char));
                strcpy(user_inputs->n_file, optarg);
				break;
            case 'o':
                user_inputs->output_dir = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->output_dir, optarg);
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
                    exit(EXIT_FAILURE);
                }
                user_inputs->num_of_threads = atoi(optarg);
				break;
			case 'u':
				if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered upper_bound value %s is not a number\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                user_inputs->upper_bound = atoi(optarg);
				break;
            case 'w': user_inputs->wgs_coverage = true; break;
            case 'W': user_inputs->Write_WGS = true; break;
            case '?':
					  if (optopt == 'b' || optopt == 'B' || optopt == 'D' || optopt == 'g'
							  || optopt == 'H' || optopt == 'i' || optopt == 'L' || optopt == 'l' || optopt == 'm'
							  || optopt == 'n' || optopt == 'o' || optopt == 'p' || optopt == 't'
							  || optopt == 'T' || optopt == 'u')
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                usage();
                exit(EXIT_FAILURE);
            default: fprintf(stderr, "Non-option argument %c\n", optopt); usage(); exit(EXIT_FAILURE);
        }
    }

	// don't proceed if the user doesn't specify either -t or -w or both
	//
	if (!user_inputs->wgs_coverage && !TARGET_FILE_PROVIDED) {
		printf("\nYou specify neigher -t (for Capture Project) nor -w (for WGS analysis)\n");
		printf("Please specify either -t or -w or both before proceed. Thanks!\n\n");
		usage();
		exit(EXIT_FAILURE);
	}

	// check the mandatory arguments (will turn this on for the final test/run)
    if (user_inputs->bam_file == NULL) {
        printf("\n-i\toption is mandatory!\n\n");
		usage();
        exit(EXIT_FAILURE);
    }

	// check database version
	if ( (strcmp(user_inputs->database_version, "hg19") != 0) && (strcmp(user_inputs->database_version, "hg37") != 0)
			&& strcmp(user_inputs->database_version, "hg38") != 0) {
		printf("\n-D\toption is not correct! It should be either hg19 or hg37 or hg38! All in lower case, please! Thanks!\n\n");
		usage();
		exit(EXIT_FAILURE);
	}
	
	if (user_inputs->output_dir == NULL) {
		printf("\n-o\toption is mandatory!\n\n");
		usage();
		exit(EXIT_FAILURE);
	} else {
		// check to see if the directory exist!
		DIR* dir = opendir(user_inputs->output_dir);
		if (dir) {
			/* Directory exists */
			closedir(dir);
		} else if (ENOENT == errno) {
			printf("\nThe output directory doesn't exist! Please double check the output directory and try again. Thanks!!\n\n");
			usage();
			exit(EXIT_FAILURE);
		} else {
			/* opendir() failed for some other reason, such as permission */
			printf("\nCan't open the output directory! Please check to see if the permission is set correctly. Thanks!\n\n");
			usage();
			exit(EXIT_FAILURE);
		}
	}

	if (user_inputs->upper_bound < user_inputs->lower_bound) {
		printf("\nThe value for -u should be larger than the value for -l option \n\n");
		usage();
		exit(EXIT_FAILURE);
	}

	// Need to check out that all files user provided exist before proceeding
    if (user_inputs->bam_file) checkFile(user_inputs->bam_file);
    if (N_FILE_PROVIDED) checkFile(user_inputs->n_file);
    if (TARGET_FILE_PROVIDED) checkFile(user_inputs->target_file);

	// need to get the basename from BAM/CRAM filename
	//char *tmp_basename = basename(strdup(user_inputs->bam_file));
	char *tmp_basename = basename(user_inputs->bam_file);
	if (!tmp_basename || strlen(tmp_basename) == 0) {
		printf("\nSomething went wrong for extracting the basename from the input BAM/CRAM file\n");
		exit(EXIT_FAILURE);
	}

	char string_to_add[350];

	// For all Capture related output files
	//
	if (TARGET_FILE_PROVIDED) {
		// output average coverage for all target regions (capture)
		sprintf(string_to_add, ".Capture_AllSites_REPORT.txt");
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_all_site_file, string_to_add);
		writeHeaderLine(user_inputs->capture_all_site_file, 1);

		// For capture coverage summary report
		sprintf(string_to_add, ".Capture_Coverage_Summary_Report.csv");
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_cov_report, string_to_add);

		// for cov.fasta file name
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_cov_file, ".Capture_cov.fasta");

		// for target regions have no coverage at all
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->missed_targets_file, ".Capture_missed_targets.txt");
		writeHeaderLine(user_inputs->missed_targets_file, 2);

		// for off target good hit wig.fasta file name
		if (user_inputs->Write_WIG)
			createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_wig_file, ".Capture_off_target_good_hits.wig.fasta");

		// output low coverage regions for target (capture)
		sprintf(string_to_add, ".Capture_below%dx_REPORT.txt", user_inputs->low_coverage_to_report);
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_low_cov_file, string_to_add);
		writeHeaderLine(user_inputs->capture_low_cov_file, 1);

		// output too high coverage regions for target (capture)
		//
		sprintf(string_to_add, ".Capture_above%dx_REPORT.txt", user_inputs->high_coverage_to_report);
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_high_cov_file, string_to_add);
		writeHeaderLine(user_inputs->capture_high_cov_file, 1);

		// output range block file for Uniformity Analysis
		//
		sprintf(string_to_add, ".Capture_between%dx_%dx_REPORT.txt", user_inputs->lower_bound, user_inputs->upper_bound);
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_range_file, string_to_add);

		// for low coverage gene/exon/transcript reports
		//
		if (user_inputs->annotation_on) {
			sprintf(string_to_add, ".Capture_below%dx_Gene_pct.txt", user_inputs->low_coverage_to_report);
			createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_gene_pct_file, string_to_add);
			writeHeaderLine(user_inputs->low_cov_gene_pct_file, 5);

			sprintf(string_to_add, ".Capture_below%dx_Exon_pct.txt", user_inputs->low_coverage_to_report);
			createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_exon_pct_file, string_to_add);
			writeHeaderLine(user_inputs->low_cov_exon_pct_file, 4);

			sprintf(string_to_add, ".Capture_below%dx_Transcript_pct.txt", user_inputs->low_coverage_to_report);
			createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_transcript_file, string_to_add);
			writeHeaderLine(user_inputs->low_cov_transcript_file, 3);
		}
	}

	// For whole Genome report
	//
	if (user_inputs->wgs_coverage) {
		// output WGS coverage summary report
		sprintf(string_to_add, ".WGS_Coverage_Summary_Report.csv");
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_cov_report, string_to_add);

		// output low coverage regions for WGS
		sprintf(string_to_add, ".WGS_below%dx_REPORT.txt", user_inputs->low_coverage_to_report);
    	createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_low_cov_file, string_to_add);
		writeHeaderLine(user_inputs->wgs_low_cov_file, 1);

	    // output too high coverage regions for the whole genome
		//
		sprintf(string_to_add, ".WGS_above%dx_REPORT.txt", user_inputs->high_coverage_to_report);
	    createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_high_cov_file, string_to_add);
		writeHeaderLine(user_inputs->wgs_high_cov_file, 1);

		// output the range block file for Uniformity Analysis
		//
		sprintf(string_to_add, ".WGS_between%dx_%dx_REPORT.txt", user_inputs->lower_bound, user_inputs->upper_bound);
	    createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_range_file, string_to_add);

		// for whole genome (wgs) file name
    	if (user_inputs->Write_WGS) {
        	createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_cov_file, ".WGS_cov.fasta");
    	    //printf("Create wgs file name %s\n", user_inputs->wgs_file);
	    }
	}

	// KEEP the following Please!
	// free the memory
	// However, here is the quote from the basename() official site
	// Both dirname() and basename() return pointers to null-terminated strings. (Do not pass these pointers to free(3))
	//if (tmp_basename) {
	//	free(tmp_basename);
	//	tmp_basename=NULL;
	//}

	// string_to_add is declared at the stack, so no need to free it!
}

//Here I need to pass in file_in name string as reference, otherwise, it will be by value and will get segmentation fault
void createFileName(char *output_dir, char *base_name, char **file_in, char *string_to_append) {
	*file_in = calloc(strlen(output_dir)+strlen(base_name)+50,  sizeof(char));
	strcpy(*file_in, output_dir);
	strcat(*file_in, "/");
	strcat(*file_in, base_name);
	strcat(*file_in, string_to_append);

	// need to write the version number to every output file
	FILE *out_fp = fopen(*file_in, "w");
	fprintf(out_fp, "%s\n", VERSION_);
	fclose(out_fp);
}

// need to write the header line for some of the output files
void writeHeaderLine(char *file_in, uint8_t type) {
	FILE *out_fp = fopen(file_in, "a");

	if (type == 1) {
		// for the coverage annotation report (for example: below20x, above10000x coverage reports)
		fprintf(out_fp, "##This file will be produced when the user specifies the target file (For Capture only) and need detailed annotations\n");
		fprintf(out_fp, "##%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Chrom", "Start", "End", "Length", "Coverage", "Gene_Symbol", "Synonymon", "Prev_Gene_Symbol", "RefSeq", "CCDS", "VEGA", "miRNA", "Others (SNP, Pseudo-Gene etc.)");
	} else if (type == 2) {
		// for capture missed target file
		fprintf(out_fp, "##This file will be produced when the user specifies the target file (For Capture only)\n");
		fprintf(out_fp, "##It contains the target regions that do not get covered\n");
		fprintf(out_fp, "##%s\t%s\t%s\n", "Chrom", "Start", "End");
	} else if (type == 3) {
		// for gene percentage coverage annotation reports
		fprintf(out_fp, "##This file will be produced when the user specifies the target file (For Capture only)\n");
		fprintf(out_fp, "##It contains the percentage of coverage information related to a Gene/RefSeq pair that intersects with a specific set of target regions from the input target bed file\n");
		fprintf(out_fp, "##%s\t%s\t%s\t%s\t%s\t%s\n", "Chrom", "Gene_Symbol", "RefSeq", "Length", "Exon_Count", "Percentage_Of_Coverage");
	} else if (type == 4) {
		// for exon percentage coverage annotation reports
		fprintf(out_fp, "##This file will be produced when the user specifies the target file (For Capture only)\n");
		fprintf(out_fp, "##It contains the percentage of coverage information related to a group of Gene/RefSeq/Exon that intersects with a specific target region from the input target bed file\n");
		fprintf(out_fp, "##%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Chrom", "Gene_Symbol", "RefSeq", "Exon_ID", "Start", "End", "Percentage_Of_Coverage", "Regions_With_Low_Coverage");
	} else if (type == 5) {
		// for gene/transcript percentage coverage reports
		fprintf(out_fp, "##This file will be produced when the user specifies the target file (For Capture only)\n");
		fprintf(out_fp, "##It contains the percentage of coverage information related to all RefSeq transcripts under a Gene that intersects with a specific set of target regions from the input target bed file\n");
		fprintf(out_fp, "##%s\t%s\t%s\n", "Gene_Symbol", "RefSeq_List(Percentage_Of_Coverage)", "Average_Percentage_Of_Coverage");
	}

	fclose(out_fp);
}

void outputUserInputOptions(User_Input *user_inputs) {
	fprintf(stderr, "The following are the options you have choosen:\n");
	fprintf(stderr, "\tInput bam/cram file: %s\n", user_inputs->bam_file);
	fprintf(stderr, "\tOutput directory is: %s\n", user_inputs->output_dir);
	fprintf(stderr, "\tThe version of official gene annotation is: %s\n", user_inputs->database_version);

	if (user_inputs->target_file)
		fprintf(stderr, "\tThe capture/target bed file is: %s\n", user_inputs->target_file);

	if (user_inputs->n_file)
		fprintf(stderr, "\tThe file that contains all Ns regions is: %s\n", user_inputs->n_file);

	if (USER_DEFINED_DATABASE) {
		fprintf(stderr, "\tUser provided database is %s.\n", user_inputs->user_defined_database_file);
		fprintf(stderr, "\t\tAll annotations will be based on this file!\n");
	}

	fprintf(stderr, "\tThe minimum mapping quality is: %d\n", user_inputs->min_map_quality);
	fprintf(stderr, "\tThe minimum base quality is: %d\n", user_inputs->min_base_quality);
	fprintf(stderr, "\tThe number of thread used is: %d\n", user_inputs->num_of_threads);
	fprintf(stderr, "\tThe percentage of reads used for analysis is: %.1f%%\n", user_inputs->percentage*100);
	fprintf(stderr, "\tThe coverage number used for low coverage report is:  < %d\n", user_inputs->low_coverage_to_report);
	fprintf(stderr, "\tThe coverage number used for high coverage report is: > %d\n", user_inputs->high_coverage_to_report);
	fprintf(stderr, "\tThe percentage used for gVCF block grouping is %d\n", user_inputs->gVCF_percentage);
	fprintf(stderr, "\tThe buffer size around a target region is %d\n", user_inputs->target_buffer_size);

	fprintf(stderr, "\tThe range block file will be produced\n");
	fprintf(stderr, "\t\tThe range is between %d and %d inclusive! \n", user_inputs->lower_bound, user_inputs->upper_bound);

	if (user_inputs->annotation_on) {
		fprintf(stderr, "\tThe detailed gene annotation is ON\n");
	} else {
		fprintf(stderr, "\tThe detailed gene annotation is OFF\n");
	}

	if (user_inputs->Write_WGS) {
		fprintf(stderr, "\tThe whole genome coverage dump is ON (this will create a WGS cov.fasta file\n");
	} else {
		fprintf(stderr, "\tThe whole genome coverage dump is OFF\n");
	}

	if (user_inputs->Write_WIG) {
		fprintf(stderr, "\tThe WIG file creation is ON\n");
	} else {
		fprintf(stderr, "\tThe WIG file creation is OFF\n");
	}

	if (user_inputs->wgs_coverage) {
		fprintf(stderr, "\tThe whole genome analysis is ON\n");
	} else {
		fprintf(stderr, "\tThe whole genome analysis is OFF\n");
	}

	if (user_inputs->remove_duplicate) {
		fprintf(stderr, "\tRemove duplicate reads is ON\n");
	} else {
		fprintf(stderr, "\tRemove duplicate reads is OFF\n");
	}

	if (user_inputs->remove_supplementary_alignments) {
		fprintf(stderr, "\tRemove supplementaty alignments is ON\n");
	} else {
		fprintf(stderr, "\tRemove supplementaty alignments is OFF\n");
	}

	fprintf(stderr, "User Input Options ===> DONE!\n\n");
	//printf("The  is: %d\n", user_inputs->);

}

User_Input * userInputInit() {
	User_Input * user_inputs = calloc(1, sizeof(User_Input));
	if (!user_inputs) {
		fprintf(stderr, "Memory allocation failed in line %d!\n", __LINE__);
		exit(EXIT_FAILURE);
	}

	user_inputs->min_map_quality  = 0;
	user_inputs->min_base_quality = 0;
	user_inputs->low_coverage_to_report = 20;
	user_inputs->high_coverage_to_report = 10000;
	user_inputs->lower_bound = 1;
	user_inputs->upper_bound = 150;
	user_inputs->target_buffer_size = 100;
	user_inputs->gVCF_percentage = 10;
	user_inputs->num_of_threads   = 3;
	user_inputs->percentage = 1.0;
	user_inputs->annotation_on = true;
	user_inputs->wgs_coverage = false;
	user_inputs->Write_WIG = false;
	user_inputs->Write_WGS = false;
	user_inputs->remove_duplicate = true;
	user_inputs->remove_supplementary_alignments = false;

	user_inputs->database_version = calloc(10, sizeof(char));
	strcpy(user_inputs->database_version, "hg37");
	
	return user_inputs;
}

void userInputDestroy(User_Input *user_inputs) {

	if (user_inputs->database_version)
		free(user_inputs->database_version);

	if (user_inputs->bam_file)
		free(user_inputs->bam_file);

	if (user_inputs->output_dir)
		free(user_inputs->output_dir);

	if (user_inputs->n_file)
		free(user_inputs->n_file);

	// whole genome output files clean-up
	//
	if (user_inputs->wgs_wig_file)
		free(user_inputs->wgs_wig_file);

	if (user_inputs->wgs_cov_file)
		free(user_inputs->wgs_cov_file);

	if (user_inputs->wgs_cov_report)
		free(user_inputs->wgs_cov_report);

	if (user_inputs->wgs_low_cov_file)
		free(user_inputs->wgs_low_cov_file);

	if (user_inputs->wgs_high_cov_file)
		free(user_inputs->wgs_high_cov_file);

	if (user_inputs->wgs_range_file)
		free(user_inputs->wgs_range_file);

	// Capture (target) output files clean-up
	//
	if (user_inputs->target_file)
		free(user_inputs->target_file);

	if (user_inputs->capture_cov_file) 
		free(user_inputs->capture_cov_file);

	if (user_inputs->capture_cov_report)
		free(user_inputs->capture_cov_report);

	if (user_inputs->capture_all_site_file)
		free(user_inputs->capture_all_site_file);

	if (user_inputs->capture_low_cov_file)
		free(user_inputs->capture_low_cov_file);

	if (user_inputs->capture_high_cov_file)
		free(user_inputs->capture_high_cov_file);

	if (user_inputs->capture_range_file)
		free(user_inputs->capture_range_file);

	if (user_inputs->missed_targets_file)
		free(user_inputs->missed_targets_file);

	if (user_inputs->low_cov_gene_pct_file)
		free(user_inputs->low_cov_gene_pct_file);

	if (user_inputs->low_cov_exon_pct_file)
		free(user_inputs->low_cov_exon_pct_file);

	if (user_inputs->low_cov_transcript_file)
		free(user_inputs->low_cov_transcript_file);

	// For User-Defined Database output files clean-up
	//
	if (user_inputs->user_defined_database_file)
		free(user_inputs->user_defined_database_file);

	if (user_inputs)
		free(user_inputs);
}

void fetchTotalGenomeBases(bam_hdr_t *header, Stats_Info *stats_info) {
	int i;
	for ( i = 0; i < header->n_targets; i++)
		stats_info->cov_stats->total_genome_bases += header->target_len[i];
}

void cleanKhashInt(khash_t(m32) *hash_to_clean) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k)
		if (kh_exist(hash_to_clean, k))
			kh_del(m32, hash_to_clean, k);

	if (hash_to_clean) kh_destroy(m32, hash_to_clean);
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

void cleanKhashStrInt(khash_t(khStrInt) *hash_to_clean) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
		if (kh_exist(hash_to_clean, k)) {
			// clean key if the key exist
			//
			if (kh_key(hash_to_clean, k)) free((char *) kh_key(hash_to_clean, k));
		}
	}

	if (hash_to_clean) kh_destroy(khStrInt, hash_to_clean);
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

// Note: this function should never be used to update any information regarding the Chromosome_Tracking variable
// It is only used for initialization and dynamically allocate memories!
//
Chromosome_Tracking * chromosomeTrackingInit(bam_hdr_t *header) {
	// number 25 is used as human has a total of 25 chromosomes including MT
	// However, it seems that everything is considered. 
	// So we need to fetch the total number of targets from the bam/cram header
	// if we are going to handle all other non-human contigs or decoy ones, we will have to expand the tracking list
	//
	uint32_t i=0;
	Chromosome_Tracking *chrom_tracking = calloc(1, sizeof(Chromosome_Tracking));
	if (!chrom_tracking) {
		fprintf(stderr, "Memory allocation for Chromosome Tracking variable failed\n");
		exit(EXIT_FAILURE);
	}

	chrom_tracking->coverage = calloc(header->n_targets, sizeof(uint32_t*));
	if (!chrom_tracking->coverage) {
		fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->coverage");
		exit(EXIT_FAILURE);
	}

	//since the thread finishes unevenly, some chromosome might be processed before its predecessor,
	//hence we need to initialize them here
	for(i=0; i<header->n_targets; i++)
		chrom_tracking->coverage[i] = NULL;

	chrom_tracking->chromosome_ids = calloc(header->n_targets, sizeof(char*));
	if (!chrom_tracking->chromosome_ids) {
		fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_ids");
	       exit(EXIT_FAILURE);
	}
	for(i=0; i<header->n_targets; i++)
		// need to increase the size if we need to analysis any other chromosomes (such as decoy)
		// but at this moment, we will just set it to NULL
		//
        chrom_tracking->chromosome_ids[i] = NULL;

	chrom_tracking->chromosome_lengths = calloc(header->n_targets, sizeof(uint32_t));
	if (!chrom_tracking->chromosome_lengths) {
		fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_lengths");
		exit(EXIT_FAILURE);
	}
	for(i=0; i<header->n_targets; i++)
        chrom_tracking->chromosome_lengths[i] = 0;

	chrom_tracking->chromosome_status = calloc(header->n_targets, sizeof(uint8_t));
	if (!chrom_tracking->chromosome_status) {
		fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_status");
		exit(EXIT_FAILURE);
	}
	for(i=0; i<header->n_targets; i++)
        chrom_tracking->chromosome_status[i] = 0;

	chrom_tracking->number_tracked = header->n_targets;
	chrom_tracking->more_to_read = true;

	return chrom_tracking;
}

void chromosomeTrackingUpdate(Chromosome_Tracking *chrom_tracking, char *chrom_id, uint32_t chrom_len, int index) {
	uint8_t id_len = strlen(chrom_id);
	chrom_tracking->chromosome_ids[index] = calloc(id_len+1, sizeof(char));
	strcpy(chrom_tracking->chromosome_ids[index], chrom_id);
	if (chrom_tracking->chromosome_ids[index] == NULL) {
		printf("Allocation failed for chrom %s\n", chrom_id);
		exit(EXIT_FAILURE);
	}

	// As the 0 position will be empty as the position will be 1-based
	// So I used it to store the index information for quick access
	// Also, I need to add 1 to the size to align with the 1-based position
	//
	chrom_tracking->coverage[index] = calloc(chrom_len + 1, sizeof(uint32_t));
	if (!chrom_tracking->coverage[index]) {
		fprintf(stderr, "Memory allocation failed for chromosome_tracking->coverage");
		exit(EXIT_FAILURE);
	}

	if (chrom_tracking->chromosome_status[index] == 0) 
		chrom_tracking->chromosome_status[index] = 1;
	chrom_tracking->chromosome_lengths[index] = chrom_len;
}

void chromosomeTrackingDestroy(Chromosome_Tracking *chrom_tracking) {
	uint32_t i = 0;
	for (i=0; i<chrom_tracking->number_tracked; i++) {
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

// This function is used to dynamically allocate the memory and then copy everything in
//
void dynamicStringAllocation(char *str_in, char **storage_str) {
	if (!str_in) { printf("String is null\n"); }

	if (*storage_str) {
		if (strlen(str_in) > strlen(*storage_str)) {
			*storage_str = realloc(*storage_str, (strlen(str_in) + 2)*sizeof(char));
			if (*storage_str == NULL) {
				fprintf(stderr, "Dynamic Memory allocation failed\n");
				exit(EXIT_FAILURE);
			}
		}
	} else {
		if (str_in == NULL || strlen(str_in) == 0) {
			*storage_str = calloc(5, sizeof(char));
		} else {
			*storage_str = calloc(strlen(str_in) + 2, sizeof(char));
		}
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

int32_t locateChromosomeIndexForChromTracking(char *chrom_id, Chromosome_Tracking *chrom_tracking) {
    int32_t i=0;
    for (i = 0; i < chrom_tracking->number_tracked; i++) {
        if (chrom_tracking->chromosome_ids[i]) {
            if (strcmp(chrom_id, chrom_tracking->chromosome_ids[i]) == 0) {
                return i;
            }
        }
    }

    //fprintf(stderr, "Something is wrong because the chromosome %s couldn't be found\n", chrom_id);
    return -1;
}

void statsInfoInit(Stats_Info *stats_info) {
	if (!stats_info) {
		fprintf(stderr, "Memory allocation failed for Stats_Info\n");
		exit(EXIT_FAILURE);
	}

	stats_info->target_cov_histogram = kh_init(m32);
    stats_info->genome_cov_histogram = kh_init(m32);

    stats_info->targeted_base_with_N_coverage = kh_init(m32);
    stats_info->genome_base_with_N_coverage   = kh_init(m32);

    stats_info->target_coverage_for_median = kh_init(m32);
    stats_info->genome_coverage_for_median = kh_init(m32);

	stats_info->cov_stats = calloc(1, sizeof(Coverage_Stats));
	coverageStatsInit(stats_info->cov_stats);

	// initializing all numbers to 0
	int i = 0;
	for (i=0; i<PRIMER_SIZE; i++) {
		stats_info->five_prime[i] = 0;
		stats_info->three_prime[i] = 0;
	}

	for (i=0; i<101; i++)
		stats_info->target_coverage[i] = 0;
}

void coverageStatsInit(Coverage_Stats * cov_stats) {
	if (!cov_stats) {
		fprintf(stderr, "Memory allocation failed for Coverage_Stats\n");
		exit(EXIT_FAILURE);
	}

	cov_stats->total_genome_bases = 0;
	cov_stats->total_buffer_bases = 0;
	cov_stats->total_targeted_bases = 0;
	cov_stats->total_aligned_bases = 0;
	cov_stats->total_target_coverage = 0;
	cov_stats->total_genome_coverage = 0;

	cov_stats->total_reads_paired = 0;
	cov_stats->total_reads_aligned = 0;
	cov_stats->total_reads_produced = 0;
	cov_stats->total_duplicate_reads = 0;
	cov_stats->total_supplementary_reads = 0;
	cov_stats->total_paired_reads_with_mapped_mates = 0;

	cov_stats->on_target_read_hit_count = 0;
	cov_stats->off_target_read_hit_count = 0;
	cov_stats->in_buffer_read_hit_count = 0;
	cov_stats->hit_target_count = 0;
	cov_stats->hit_target_buffer_only_count = 0;
	cov_stats->non_target_good_hits = 0;

	cov_stats->total_targets = 0;
	cov_stats->read_length = 0;
	cov_stats->max_coverage = 0;
	cov_stats->base_with_max_coverage = 0;
	cov_stats->median_genome_coverage = 0;
	cov_stats->median_target_coverage = 0;
}

void statsInfoDestroy(Stats_Info *stats_info) {
	kh_destroy(m32, stats_info->target_cov_histogram);
	kh_destroy(m32, stats_info->genome_cov_histogram);
	kh_destroy(m32, stats_info->targeted_base_with_N_coverage);
	kh_destroy(m32, stats_info->genome_base_with_N_coverage);
	kh_destroy(m32, stats_info->target_coverage_for_median);
	kh_destroy(m32, stats_info->genome_coverage_for_median);

	if (stats_info->cov_stats) { 
		free(stats_info->cov_stats);
		stats_info->cov_stats = NULL;
	}
	if (stats_info) { free(stats_info); stats_info=NULL; }
}

void zeroAllNsRegions(char *chrom_id, Bed_Info *Ns_info, Chromosome_Tracking *chrom_tracking, Target_Buffer_Status *target_buffer_status) {
	// First, we need to find the index that is used to track current chromosome chrom_id
	//
	uint32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
   	uint32_t i=0,j=0, chrom_len=0;

	// here we need to find out the length of the current chromsome to prevent out of bound segmentation error
	// get the index for the target_buffer_status
	//
    for(i=0; i<target_buffer_status[0].num_of_chromosomes; i++) {
        if (strcmp(target_buffer_status[i].chrom_id, chrom_id) == 0) {
            chrom_len = target_buffer_status[i].size;

            break;
        }
    }

	for (i=0; i<Ns_info->size; i++) {
		if (strcmp(Ns_info->coords[i].chrom_id, chrom_id) == 0) {
			//printf("%s\t%"PRIu32"\t%"PRIu32"\n", Ns_info->coords[i].chr, Ns_info->coords[i].start, Ns_info->coords[i].end);
			for (j=Ns_info->coords[i].start; j<=Ns_info->coords[i].end; j++) {
				if (j>=chrom_len) continue;

				//printf("value of j is %d\n", j);
				chrom_tracking->coverage[idx][j] = 0;
			}
		}
	}
	printf("Finished for zero all N zeros\n");
}
void addValueToKhashBucket16(khash_t(m16) *hash_in, uint16_t pos_key, uint16_t val) {
    int ret;
    khiter_t k_iter = kh_put(m16, hash_in, pos_key, &ret);
    if (ret == 1) {
        kh_value(hash_in, k_iter) = 0;
    } else if (ret == -1) {
        fprintf(stderr, "can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
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
        fprintf(stderr, "can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        exit(EXIT_FAILURE);
    }

    kh_value(hash_in, k_iter) += val;
	//if (pos_key == 100) printf("added for combined value is %"PRIu32"\n", kh_value(hash_in, k_iter));

    return;
}

uint32_t getValueFromKhash32(khash_t(m32) *hash32, uint32_t pos_key) {
	int ret;
    khiter_t k_iter;
	if (hash32 != NULL) {
		k_iter = kh_put(m32, hash32, pos_key, &ret);

		if (ret == -1) {
        	fprintf(stderr, "can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
        	exit(EXIT_FAILURE);
    	}

		// this is needed as if the bucket is never used, the value will be whatever left there, and the value is undefined!
    	if (ret == 1)
            kh_value(hash32, k_iter) = 0;

        return kh_value(hash32, k_iter);
    }

	return 0;
}

uint16_t getValueFromKhash(khash_t(m16) *hash16, uint32_t pos_key) {
    int ret;
    khiter_t k_iter;
    if (hash16 != NULL) {
		k_iter = kh_put(m16, hash16, pos_key, &ret);

		if (ret == -1) {
			fprintf(stderr, "can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
			exit(EXIT_FAILURE);
		}

		//if (ret == 1)
        //    kh_value(hash16, k_iter) = 0;

        return kh_value(hash16, k_iter);
    }

    return 0;
}

float calculatePercentage(uint32_t num, uint32_t dom) {
	float val = (double)num/(double)dom;
	int i_val = val * 10000.0 + 0.5;
	return (float)i_val/100.0;
}

void combineCoverageStats(Stats_Info *stats_info, Coverage_Stats *cov_stats) {
	stats_info->cov_stats->total_reads_produced  += cov_stats->total_reads_produced;
	stats_info->cov_stats->total_reads_aligned   += cov_stats->total_reads_aligned;
	stats_info->cov_stats->total_reads_paired    += cov_stats->total_reads_paired;
	stats_info->cov_stats->total_aligned_bases   += cov_stats->total_aligned_bases;
	stats_info->cov_stats->total_duplicate_reads += cov_stats->total_duplicate_reads;
	stats_info->cov_stats->total_supplementary_reads += cov_stats->total_supplementary_reads;
	stats_info->cov_stats->on_target_read_hit_count  += cov_stats->on_target_read_hit_count;
	stats_info->cov_stats->in_buffer_read_hit_count  += cov_stats->in_buffer_read_hit_count;
	stats_info->cov_stats->off_target_read_hit_count += cov_stats->off_target_read_hit_count;

	stats_info->cov_stats->total_paired_reads_with_mapped_mates += cov_stats->total_paired_reads_with_mapped_mates;

	if (stats_info->cov_stats->read_length == 0) 
		stats_info->cov_stats->read_length = cov_stats->read_length;

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
			fprintf(stderr, "Memory re-allocation failed at lowCoverageGeneHashBucketKeyInit()!\n");
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
		fprintf(stderr, "gc1 or gc2 shouldn't be NULL\n");
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
			gc2->low_cov_regions = calloc(1, sizeof(stringArray));
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
				fprintf(stderr, "Memory re-allocation failed at copyGeneCoverageLowCovRegions()!\n");
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
void getOneLowCovRegions(khash_t(khStrInt) *low_cov_regions_hash, stringArray *mergedArray) {
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

void mergeLowCovRegions(khash_t(khStrInt) *low_cov_regions_hash, stringArray *mergedArray, uint32_t size_in,  uint32_t cds_t_start, uint32_t cds_t_end) {
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
					kh_value(starts, k_iter)++;
				}

				if (i==1) {
					uint32_t end = (uint32_t) strtol(tokPtr, NULL, 10);
					if (end > cds_t_end)
						end = cds_t_end;

					k_iter = kh_put(m32, ends, end, &ret);
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
	uint32_t p_start, num_of_items=0;
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

// the following comparison is used to compare the uint32_t array
//
int compare(const void * val1, const void * val2) {
	uint32_t tmp_val1 = *((uint32_t*) val1);
	uint32_t tmp_val2 = *((uint32_t*) val2);

	if (tmp_val1 == tmp_val2) return 0;
	else if (tmp_val1 < tmp_val2) return -1;
	else return 1;
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
