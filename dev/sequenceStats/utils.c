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
#include <dirent.h>		// for checking output directory
#include <libgen.h>		// for function basename()
#include "terms.h"
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

uint64_t check_file_size(const char *filename) {
	struct stat st;
	if (stat(filename, &st) == 0)
		return st.st_size;

	fprintf(stderr, "Something is wrong when check the filesize for file: %s\n", filename);
	exit(1);
}

// Split a string to a word array
// Right now we will use static way:
// The first dimension is for the number of strings to be stored
// (the max is 96-100, but I will use 300 instead)
// The second dimension is used to store the string like OTTHUMT00000367763_exon_16 or NM_198384_exon_16 etc.
// (the max is 26-30, but we will use 50 instead)
// In the future, we might want to use a way that is more dynamic!
//
void splitStringToArray(char *stringPtr, char (*arrayPtr)[50]) {
	//int sizeP1 = 300, sizeP2 = 50;
	int sizeP1 = 300;
	char *tokPtr;
	arrayPtr = calloc(sizeP1, sizeof(char*));
	size_t i, j, k, l;

	tokPtr = strtok(stringPtr, ";");
	//arrayPtr[0] = calloc(sizeP2, sizeof(char));
	strcpy(arrayPtr[0], tokPtr);

	for(i = 1; (tokPtr = strtok(NULL, " ")) != NULL; i++) {
		//arrayPtr[i] = calloc(sizeP2, sizeof(char));
		strcpy(arrayPtr[i], tokPtr);
	}

	for(j = 0; j <= i; j++) {
		for(k = j + 1; k <= i; k++) {
			if(strcmp(arrayPtr[j], arrayPtr[k]) == 0) {
				for(l = k; l < i; l++)
					strcpy(arrayPtr[l], arrayPtr[l + 1]);

				k = j;
				i--;
			}
		}
	}

	//for(l = 0; l <= i; l++)
	//	printf("%s ", wordTable[l]);
}

// It will remove duplicated words from a string array
// The first dimension is the array size
// The second dimension is the size of the corresponding string (30 is OK, but we will use 50)
//
void removeSameWordsFromArray (char **arrayPtr, uint16_t arraySize, char **str_in_and_out) {
	// Before we start, let's sort the string array first!
	//
	qsort(arrayPtr, arraySize, sizeof(char *), compare3);

	size_t i, j, k;

	for (j=0; j<arraySize; j++) {
		for (k=j+1; k<arraySize; k++) {
			if (strcmp(arrayPtr[j], arrayPtr[k]) == 0) {
				for (i=k; i<arraySize-1; i++)
					strcpy(arrayPtr[i], arrayPtr[i+1]);
				k=j;
				arraySize--;
			}
		}
	}

	// find the total size needed for string 
	//
	uint32_t str_len_needed=0;

	for (i=0; i<=arraySize; i++) {
		str_len_needed += strlen(arrayPtr[i]);
	}

	// reallocate the memory 
	//
	char *tmp = realloc(*str_in_and_out, str_len_needed);
	if (!tmp) {
		fprintf(stderr, "Memory re-allocation for string failed in checkExonRegion\n");
		exit(1);
	}
	*str_in_and_out = tmp;

	for(i=0; i<arraySize; i++) {
		strcat(*str_in_and_out, arrayPtr[i]);
	}
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
	printf("\t-g <the percentage used for gVCF blocking: Default 5 for 500%%>\n");
	printf("\t-m <minimal mapping quality score: to filter out any reads with mapQ less than m. Default 0>\n");
	printf("\t-n <the file that contains regions of Ns in the reference genome in bed format>\n");
	printf("\t-p <the percentage (fraction) of reads used for this analysis. Default 1.0 (ie, 100%%)>\n");
	printf("\t-t <target file. If this is specified, all of the output file names related to this will contain .Capture_>\n");
	printf("\t-y <type of annotation: 1 for dynamic or 2 for static. (Default 1: dynamic)>\n\n");
	printf("\t-c <database catetory if you choose static annotation type: 1:VCRome+PKv2, 2:eMerge or 3: right_10K. (Default 1: VCRome+PKv2)>\n\n");

	printf("\t-B <the Buffer size immediate adjacent to a target region. Default: 100>\n");
	printf("\t-D <the version of human genome database (either hg19 [or hg37], or hg38). Default:hg19>\n");
	printf("\t-H <the high coverage cutoff value. Any coverages larger than it will be outputted. Default=10000>\n");
	printf("\t-L <the low coverage cutoff value. Any coverages smaller than it will be outputted. Default=20>\n");
	printf("\t-T <the number of threads (Note: when used with HPC's msub, make sure number of processors:ppn matches to number of threads). Default 3>\n");
	printf("\t-U <the upper bound cutoff value when reporting the range information (it's associated with -H to form an interval. Hence, -H is the lower bound. It is very useful for plotting the coverage graphs if needed). Default -1 (not set) >\n\n");

	printf("The Followings Are Flags\n");
	printf("\t[-a] write the annotation information for genes, exons and transcript. Default off\n");
	printf("\t[-d] Remove Duplicates! Specify this flag only when you want to use Duplicates reads. Default: ON (not specified)\n");
	printf("\t[-s] Remove Supplementary alignments and DO NOT use them for statistics. Default off\n");
	printf("\t[-w] Write whole genome coverage related reports (all of the output file names related to this will have .WGS_ in them). This flag doesn't produce the WGS Coverage.fasta file, use -W for that. Default off\n");

	printf("\t[-G] Write/Dump the WIG formatted file. Default off\n");
	printf("\t[-W] Write/Dump the WGS Coverage.fasta file (both -w and -W needed). Default off\n");
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
	//
	while ((arg = getopt(argc, argv, "ab:B:c:dD:g:GH:i:L:m:n:o:p:st:T:U:wWy:h")) != -1) {
		//printf("User options for %c is %s\n", arg, optarg);
		switch(arg) {
			case 'a':
				user_inputs->annotation_on = true; break;
			case 'b':
				if (!isNumber(optarg)) {
					fprintf (stderr, "Entered base quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
                user_inputs->min_base_quality = atoi(optarg);
                break;
			case 'B': user_inputs->target_buffer_size = atoi(optarg); break;
			case 'c': user_inputs->database_category = atoi(optarg); break;
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
			case 'g': user_inputs->gVCF_percentage = atoi(optarg); break;
            case 'G': user_inputs->Write_WIG = true; break;
            case 'h': usage(); exit(1);
			case 'H':
				if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered map quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
                user_inputs->high_coverage_to_report = atoi(optarg);
                break;
            case 'i':
				user_inputs->bam_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));

                strcpy(user_inputs->bam_file, optarg);
                break;
			case 'L':
				if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered map quality filter score %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
                user_inputs->low_coverage_to_report = atoi(optarg);
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
                    exit(1);
                }
                user_inputs->num_of_threads = atoi(optarg);
				break;
			case 'U':
				if (!isNumber(optarg)) {
                    fprintf (stderr, "Entered number of threads %s is not a number\n", optarg);
                    usage();
                    exit(1);
                }
                user_inputs->upper_bound_to_report = atoi(optarg);
				break;
            case 'w': user_inputs->wgs_coverage = true; break;
            case 'W': user_inputs->Write_WGS = true; break;
			case 'y': user_inputs->annotation_type = atoi(optarg); break;
            case '?':
					  if (optopt == 'b' || optopt == 'B' || optopt == 'c' || optopt == 'D' || optopt == 'g' 
							  || optopt == 'H' || optopt == 'i' || optopt == 'L' || optopt == 'm' 
							  || optopt == 'n' || optopt == 'o' || optopt == 'p' || optopt == 't'
							  || optopt == 'T' || optopt == 'U' || optopt == 'y')
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

	// don't proceed if the user doesn't specify either -t or -w or both
	//
	if (!user_inputs->wgs_coverage && !TARGET_FILE_PROVIDED) {
		printf("\nYou specify neigher -t (for Capture Project) nor -w (for WGS analysis)\n");
		printf("Please specify either -t or -w or both before proceed. Thanks!\n\n");
		usage();
		exit(1);
	}

	// check the mandatory arguments (will turn this on for the final test/run)
    if (user_inputs->bam_file == NULL) {
        printf("\n-i\toption is mandatory!\n\n");
		usage();
        exit(1);
    }

	// check database version
	if ( (strcmp(user_inputs->database_version, "hg19") != 0) && (strcmp(user_inputs->database_version, "hg37") != 0)
			&& strcmp(user_inputs->database_version, "hg38") != 0) {
		printf("\n-D\toption is not correct! It should be either hg19 or hg37 or hg38! All in lower case, please! Thanks!\n\n");
		usage();
		exit(1);
	}
	
	if (user_inputs->output_dir == NULL) {
		printf("\n-o\toption is mandatory!\n\n");
		usage();
		exit(1);
	} else {
		// check to see if the directory exist!
		DIR* dir = opendir(user_inputs->output_dir);
		if (dir) {
			/* Directory exists */
			closedir(dir);
		} else if (ENOENT == errno) {
			printf("\nThe output directory doesn't exist! Please double check the output directory and try again. Thanks!!\n\n");
			usage();
			exit(1);
		} else {
			/* opendir() failed for some other reason, such as permission */
			printf("\nCan't open the output directory! Please check to see if the permission is set correctly. Thanks!\n\n");
			usage();
			exit(1);
		}
	}

	if ((user_inputs->upper_bound_to_report > 0) && 
		(user_inputs->upper_bound_to_report < user_inputs->high_coverage_to_report)) {
		printf("\n-U option should be larger than -H option (default -H option is 10000)\n\n");
		usage();
		exit(1);
	}

	// Need to check out that all files user provided exist before proceeding
    if (user_inputs->bam_file) checkFile(user_inputs->bam_file);
    if (N_FILE_PROVIDED) checkFile(user_inputs->n_file);
    if (TARGET_FILE_PROVIDED) checkFile(user_inputs->target_file);

	// need to get the basename from BAM/CRAM filename
	char *tmp_basename = basename(strdup(user_inputs->bam_file));
	if (!tmp_basename || strlen(tmp_basename) == 0) {
		printf("\nSomething went wrong for extracting the basename from the input BAM/CRAM file\n");
		exit(1);
	}

	char string_to_add[350];

	// For all Capture related output files
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
		if (user_inputs->upper_bound_to_report == -1) {
			sprintf(string_to_add, ".Capture_above%dx_REPORT.txt", user_inputs->high_coverage_to_report);
		} else {
			sprintf(string_to_add, ".Capture_between%dx_%dx_REPORT.txt", user_inputs->high_coverage_to_report, user_inputs->upper_bound_to_report);
		}
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_high_cov_file, string_to_add);
		writeHeaderLine(user_inputs->capture_high_cov_file, 1);

		// for low coverage gene/exon/transcript reports
		//
		if (user_inputs->annotation_on) {
			sprintf(string_to_add, ".Capture_below%dx_Gene_pct.txt", user_inputs->low_coverage_to_report);
			createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_gene_pct_file, string_to_add);
			writeHeaderLine(user_inputs->low_cov_gene_pct_file, 3);

			sprintf(string_to_add, ".Capture_below%dx_Exon_pct.txt", user_inputs->low_coverage_to_report);
			createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_exon_pct_file, string_to_add);
			writeHeaderLine(user_inputs->low_cov_exon_pct_file, 4);

			sprintf(string_to_add, ".Capture_below%dx_Transcript_pct.txt", user_inputs->low_coverage_to_report);
			createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_transcript_file, string_to_add);
			writeHeaderLine(user_inputs->low_cov_transcript_file, 5);
		}
	}

	// For whole Genome report
	if (user_inputs->wgs_coverage) {
		// output WGS coverage summary report
		sprintf(string_to_add, ".WGS_Coverage_Summary_Report.csv");
		createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_cov_report, string_to_add);

		// output low coverage regions for WGS
		sprintf(string_to_add, ".WGS_below%dx_REPORT.txt", user_inputs->low_coverage_to_report);
    	createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_low_cov_file, string_to_add);
		writeHeaderLine(user_inputs->wgs_low_cov_file, 1);

	    // output too high coverage regions for target (capture)
		if (user_inputs->upper_bound_to_report == -1) {
			sprintf(string_to_add, ".WGS_above%dx_REPORT.txt", user_inputs->high_coverage_to_report);
		} else {
			sprintf(string_to_add, ".WGS_between%dx_%dx_REPORT.txt", user_inputs->high_coverage_to_report, user_inputs->upper_bound_to_report);
		}
	    createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_high_cov_file, string_to_add);
		writeHeaderLine(user_inputs->wgs_high_cov_file, 1);

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
		fprintf(out_fp, "##%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Chrom", "Start", "End", "Length", "Coverage", "Gene_Symbol", "Prev_Gene_Symbol", "Synonymon", "RefSeq", "CCDS", "VEGA", "miRNA");
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

	fprintf(stderr, "\tThe minimum mapping quality is: %d\n", user_inputs->min_map_quality);
	fprintf(stderr, "\tThe minimum base quality is: %d\n", user_inputs->min_base_quality);
	fprintf(stderr, "\tThe number of thread used is: %d\n", user_inputs->num_of_threads);
	fprintf(stderr, "\tThe percentage of reads used for analysis is: %.1f%%\n", user_inputs->percentage*100);
	fprintf(stderr, "\tThe coverage number used for low coverage report is:  < %d\n", user_inputs->low_coverage_to_report);
	fprintf(stderr, "\tThe coverage number used for high coverage report is: > %d\n", user_inputs->high_coverage_to_report);
	fprintf(stderr, "\tThe percentage used for gVCF block grouping is %d\n", user_inputs->gVCF_percentage);
	fprintf(stderr, "\tThe buffer size around a target region is %d\n", user_inputs->target_buffer_size);

	if (user_inputs->upper_bound_to_report > 1) {
		fprintf(stderr, "\tThe range block file will be produced\n");
		fprintf(stderr, "\t\tThe range is between %d and %d inclusive! \n", user_inputs->high_coverage_to_report, user_inputs->upper_bound_to_report);
	}

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
	user_inputs->upper_bound_to_report = -1;
	user_inputs->target_buffer_size = 100;                                                                    
	user_inputs->gVCF_percentage = 5;
	user_inputs->num_of_threads   = 3;
	user_inputs->percentage = 1.0;
	user_inputs->annotation_on = false;
	user_inputs->wgs_coverage = false;
	user_inputs->Write_WIG = false;
	user_inputs->Write_WGS = false;
	user_inputs->remove_duplicate = true;
	user_inputs->remove_supplementary_alignments = false;

	user_inputs->annotation_type = 1;
	user_inputs->database_category = 1;
	user_inputs->database_version = calloc(10, sizeof(char));
	strcpy(user_inputs->database_version, "hg37");
	
	return user_inputs;
}

void userInputDestroy(User_Input *user_inputs) {

	if (user_inputs->database_version)
		free(user_inputs->database_version);

	if (user_inputs->target_file)
		free(user_inputs->target_file);

	if (user_inputs->bam_file)
		free(user_inputs->bam_file);

	if (user_inputs->capture_cov_file) 
		free(user_inputs->capture_cov_file);

	if (user_inputs->missed_targets_file)
		free(user_inputs->missed_targets_file);

	if (user_inputs->output_dir)
		free(user_inputs->output_dir);

	if (user_inputs->n_file)
		free(user_inputs->n_file);

	if (user_inputs->wgs_wig_file)
		free(user_inputs->wgs_wig_file);

	if (user_inputs->wgs_cov_file)
		free(user_inputs->wgs_cov_file);

	if (user_inputs->wgs_cov_report)
		free(user_inputs->wgs_cov_report);

	if (user_inputs->capture_cov_report)
		free(user_inputs->capture_cov_report);

	if (user_inputs->capture_all_site_file)
		free(user_inputs->capture_all_site_file);

	if (user_inputs->capture_low_cov_file)
		free(user_inputs->capture_low_cov_file);

	if (user_inputs->capture_high_cov_file)
		free(user_inputs->capture_high_cov_file);

	if (user_inputs->low_cov_gene_pct_file)
		free(user_inputs->low_cov_gene_pct_file);

	if (user_inputs->low_cov_exon_pct_file)
		free(user_inputs->low_cov_exon_pct_file);

	if (user_inputs->low_cov_transcript_file)
		free(user_inputs->low_cov_transcript_file);

	if (user_inputs->wgs_low_cov_file)
        free(user_inputs->wgs_low_cov_file);

	if (user_inputs->wgs_high_cov_file)
        free(user_inputs->wgs_high_cov_file);

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

	//printf("before clean hash int\n");
	if (hash_to_clean) kh_destroy(m32, hash_to_clean);
	//printf("after clean hash int\n");
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

uint32_t getChromIndexFromID(bam_hdr_t *header, char *chrom_id) {
	uint32_t i=0;
	for(i=0; i<header->n_targets; i++) {
		if (strcmp(removeChr(header->target_name[i]), chrom_id) == 0) {
			return i;
		}
	}

	fprintf(stderr, "Can't locate the chromosome name %s at the function get_chrom_index_from_id\n", chrom_id);
	exit(1);
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
		exit(1);
	}

	chrom_tracking->coverage = calloc(header->n_targets, sizeof(uint32_t*));
	if (!chrom_tracking->coverage) {
		fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->coverage");
		exit(1);
	}

	//since the thread finishes unevenly, some chromosome might be processed before its predecessor,
	//hence we need to initialize them here
	for(i=0; i<header->n_targets; i++)
		chrom_tracking->coverage[i] = NULL;

	chrom_tracking->chromosome_ids = calloc(header->n_targets, sizeof(char*));
	if (!chrom_tracking->chromosome_ids) {
		fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_ids");
	       exit(1);
	}
	for(i=0; i<header->n_targets; i++)
		// need to increase the size if we need to analysis any other chromosomes (such as decoy)
		// but at this moment, we will just set it to NULL
		//
        chrom_tracking->chromosome_ids[i] = NULL;

	chrom_tracking->chromosome_lengths = calloc(header->n_targets, sizeof(uint32_t));
	if (!chrom_tracking->chromosome_lengths) {
		fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_lengths");
		exit(1);
	}
	for(i=0; i<header->n_targets; i++)
        chrom_tracking->chromosome_lengths[i] = 0;

	chrom_tracking->chromosome_status = calloc(header->n_targets, sizeof(uint8_t));
	if (!chrom_tracking->chromosome_status) {
		fprintf(stderr, "Memory allocation for %s failed\n", "chrom_tracking->chromosome_status");
		exit(1);
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
		exit(1);
	}

	// As the 0 position will be empty as the position will be 1-based
	// So I used it to store the index information for quick access
	// Also, I need to add 1 to the size to align with the 1-based position
	//
	chrom_tracking->coverage[index] = calloc(chrom_len + 1, sizeof(uint32_t));
	if (!chrom_tracking->coverage[index]) {
		fprintf(stderr, "Memory allocation failed for chromosome_tracking->coverage");
		exit(1);
	}

	/*
	// need to initialization all of them to 0 (it is supposed to be done by the calloc(), 
	uint32_t i=0;
	for (i=0; i<chrom_len+1; i++) {
	//	printf("before initialization the values are %"PRIu32"\t%d\n", i, chrom_tracking->coverage[index][i]);
		if (i >= chrom_len) printf("I am at %"PRIu32"\n", i);
		chrom_tracking->coverage[index][i] = 0;
	} */
	//printf("After the initializatioin to 0\n");

	if (chrom_tracking->chromosome_status[index] == 0) 
		chrom_tracking->chromosome_status[index] = 1;
	chrom_tracking->chromosome_lengths[index] = chrom_len+1;
	//chrom_tracking->number_tracked++;
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
	char *tmp;
	if (!str_in) { printf("String is null\n"); }

	if (strlen(str_in) == 0) { 
		//printf("String is empty\n"); 
		strcpy(str_in, "."); 
	}

	if (*storage_str) {
		if (strlen(str_in) > strlen(*storage_str)) {
			tmp = realloc(*storage_str, strlen(str_in) + 2);
			if (!tmp) {
				fprintf(stderr, "Dynamic Memory allocation failed\n");
				exit(1);
			}
			*storage_str = tmp;
		}
	} else {
		*storage_str = calloc(strlen(str_in) + 2, sizeof(char));
	}

	strcpy(*storage_str, str_in);

	//if (tmp) { free(tmp); tmp=NULL; }
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

Stats_Info * statsInfoInit() {
	Stats_Info *stats_info = malloc(sizeof(Stats_Info));
	if (!stats_info) {
		fprintf(stderr, "Memory allocation failed for Stats_Info\n");
		exit(1);
	}

	stats_info->target_cov_histogram = kh_init(m32);
    stats_info->genome_cov_histogram = kh_init(m32);

    stats_info->targeted_base_with_N_coverage = kh_init(m32);
    stats_info->genome_base_with_N_coverage   = kh_init(m32);

    stats_info->target_coverage_for_median = kh_init(m32);
    stats_info->genome_coverage_for_median = kh_init(m32);

	stats_info->cov_stats = coverageStatsInit();

	// initializing all numbers to 0
	int i = 0;
	for (i=0; i<PRIMER_SIZE; i++) {
		stats_info->five_prime[i] = 0;
		stats_info->three_prime[i] = 0;
	}

	for (i=0; i<101; i++)
		stats_info->target_coverage[i] = 0;

	return stats_info;
}

Coverage_Stats * coverageStatsInit() {
	Coverage_Stats *cov_stats = malloc(sizeof(Coverage_Stats));
	if (!cov_stats) {
		fprintf(stderr, "Memory allocation failed for Coverage_Stats\n");
		exit(1);
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

	return cov_stats;
}

void statsInfoDestroy(Stats_Info *stats_info) {
	kh_destroy(m32, stats_info->target_cov_histogram);
	kh_destroy(m32, stats_info->genome_cov_histogram);
	kh_destroy(m32, stats_info->targeted_base_with_N_coverage);
	kh_destroy(m32, stats_info->genome_base_with_N_coverage);
	kh_destroy(m32, stats_info->target_coverage_for_median);
	kh_destroy(m32, stats_info->genome_coverage_for_median);

	//free(stats_info->five_prime);
	//free(stats_info->three_prime);

	free(stats_info->cov_stats);
	if (stats_info) { free(stats_info); stats_info=NULL; }
}

/*
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
*/

void zeroAllNsRegions(char *chrom_id, Bed_Info *Ns_info, Chromosome_Tracking *chrom_tracking, Target_Buffer_Status *target_buffer_status) {
	// First, we need to find the index that is used to track current chromosome chrom_id
	uint32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
   	uint32_t i=0,j=0, chrom_len=0;

	// here we need to find out the length of the current chromsome to prevent out of bound segmentation error
	// get the index for the target_buffer_status
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
        exit(1);
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
        exit(1);
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
        	exit(1);
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
			exit(1);
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

void getDatabaseName(User_Input *user_inputs, char* db_name_in_out) {
	if (user_inputs->annotation_type == 1) {
		// dynamic gene annotation
		//
		if (strcmp(user_inputs->database_version, "hg38") == 0) {
			strcpy(db_name_in_out, "Gene_RefSeq_CDS38");
		} else {
			strcpy(db_name_in_out, "Gene_RefSeq_CDS37");
		}
	} else {
		// static gene annotation
		//
		if (strcmp(user_inputs->database_version, "hg38") == 0) {
			if (user_inputs->database_category == 1) {
				strcpy(db_name_in_out, "Gene_RefSeq_CDS38");
			} else if (user_inputs->database_category == 2) {
				strcpy(db_name_in_out, "eMerge_CDS38");
			} else {
				strcpy(db_name_in_out, "Righ10K_CDS38");
			}
		} else {
			if (user_inputs->database_category == 1) {
				strcpy(db_name_in_out, "Gene_RefSeq_CDS37");
			} else if (user_inputs->database_category == 2) {
				strcpy(db_name_in_out, "eMerge_CDS37");
			} else {
				strcpy(db_name_in_out, "Righ10K_CDS37");
			}
		}
	}
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

void printLowCoverageGeneStructure(Low_Coverage_Genes *low_cov_genes) {
	uint32_t i;

	printf("Total Number of Gene Symbol is %"PRIu32"\n", low_cov_genes->total_size);

	for (i=0; i<low_cov_genes->total_size; i++) {
		printf("Gene: %s\tRefSeq: %s\n", low_cov_genes->gene_coverage[i].gene_symbol, low_cov_genes->gene_coverage[i].gene_name);
	}
}

// this compare is to compare the refseq_name inside Gene_Coverage variable
//
int compare(const void *gene_coverage1, const void *gene_coverage2) {
	Gene_Coverage *gc1 = (Gene_Coverage *) gene_coverage1;
	Gene_Coverage *gc2 = (Gene_Coverage *) gene_coverage2;

	int compare_result = strcmp(gc1->gene_name, gc2->gene_name);
	if (compare_result == 0) {
		return gc1->cds_target_start - gc2->cds_target_start;
	} else {
		return compare_result;
	}
}

// the following comparison is used to compare the gene_symbol in Transcript_Coverage_Percentage
//
int compare2(const void *transcript_cov_pct1, const void *transcript_cov_pct2) {
	Transcript_Coverage_Percentage *tcp1 = (Transcript_Coverage_Percentage *) transcript_cov_pct1;
	Transcript_Coverage_Percentage *tcp2 = (Transcript_Coverage_Percentage *) transcript_cov_pct2;

	int compare_result = strcmp(tcp1->gene_symbol, tcp2->gene_symbol);
    if (compare_result == 0) {
        return strcmp(tcp1->gene_name, tcp2->gene_name);
    } else {
        return compare_result;
    }
}

// The following comparison is used to compare the refseq_name array variable
//
int compare3(const void *refseq_array1, const void *refseq_array2) {
	char** rs_array1 = (char**) refseq_array1;
	char** rs_array2 = (char**) refseq_array2;

	return strcmp(*rs_array1, *rs_array2);
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
