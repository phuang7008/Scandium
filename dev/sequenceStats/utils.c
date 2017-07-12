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
	printf("\t-p <the percentage (fraction) of reads used for this analysis>\n");
	printf("\t-T <the number of threads>\n");
	printf("\t[-d] Remove Duplicates and DO NOT use them for statistics\n");
	printf("\t[-G] Write/Dump the WIG formatted file\n");
	printf("\t[-s] Remove Supplementary alignments and DO NOT use them for statistics\n");
	printf("\t[-w] Write whole genome coverage report, this flag doesn't produce the WGS Coverage.fasta file, use -W for that\n");
	printf("\t[-W] Write/Dump the WGS Coverage.fasta file\n");
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
	while ((arg = getopt(argc, argv, "b:dGi:m:n:o:p:st:T:wWh")) != -1) {
		printf("User options for %c is %s\n", arg, optarg);
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
            case 'G': user_inputs->Write_WIG = true; break;
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
				break;
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
				break;
            case 'w': user_inputs->wgs_coverage = true; break;
            case 'W': user_inputs->Write_WGS = true; break;
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
	createFileName(user_inputs->bam_file, &user_inputs->cov_file, ".cov.fasta");

	// for target regions have no coverage at all
	createFileName(user_inputs->bam_file, &user_inputs->missed_targets_file, ".missedTargets.txt");

	// for off target good hit wig.fasta file name
	if (user_inputs->Write_WIG)
		createFileName(user_inputs->bam_file, &user_inputs->wig_file, ".wig.fasta");

	// for whole genome (wgs) file name
	if (user_inputs->wgs_coverage && user_inputs->Write_WGS) {
		createFileName(user_inputs->bam_file, &user_inputs->wgs_file, ".wgs.fasta");
		//printf("Create wgs file name %s\n", user_inputs->wgs_file);
	}
}

//Here I need to pass in file_in name string as reference, otherwise, it will be by value and will get segmentation fault
void createFileName(char *base_name, char **file_in, char *string_to_append) {
	*file_in = calloc(strlen(base_name)+30,  sizeof(char));
	strcpy(*file_in, base_name);
	strcat(*file_in, string_to_append);
}

User_Input * userInputInit() {
	User_Input * user_inputs = calloc(1, sizeof(User_Input));
	if (!user_inputs) {
		fprintf(stderr, "Memory allocation failed in line %d!\n", __LINE__);
		exit(EXIT_FAILURE);
	}

	user_inputs->min_map_quality  = -1;
	user_inputs->min_base_quality = -1;
	user_inputs->num_of_threads   = 4;
	user_inputs->percentage = 1.0;
	user_inputs->wgs_coverage = false;
	user_inputs->Write_WIG = false;
	user_inputs->Write_WGS = false;
	user_inputs->remove_duplicate = false;
	user_inputs->remove_supplementary_alignments = false;
	
	return user_inputs;
}

void userInputDestroy(User_Input *user_inputs) {
	if (user_inputs->target_file)
		free(user_inputs->target_file);

	if (user_inputs->bam_file)
		free(user_inputs->bam_file);

	if (user_inputs->cov_file) 
		free(user_inputs->cov_file);

	if (user_inputs->missed_targets_file)
		free(user_inputs->missed_targets_file);

	if (user_inputs->out_file)
		free(user_inputs->out_file);

	if (user_inputs->n_file)
		free(user_inputs->n_file);

	if (user_inputs->wig_file)
		free(user_inputs->wig_file);

	if (user_inputs->wgs_file)
		free(user_inputs->wgs_file);

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

void cleanKhashStr(khash_t(str) *hash_to_clean) {
	khint_t k;
	for (k = kh_begin(hash_to_clean); k != kh_end(hash_to_clean); ++k) {
		if (kh_exist(hash_to_clean, k)) {
			// clean key if the key exist
			if (kh_key(hash_to_clean, k)) free((char *) kh_key(hash_to_clean, k));

			// clean Temp_Coverage_Array
			free(kh_value(hash_to_clean, k)->cov_array);
			free(kh_value(hash_to_clean, k));
			//kh_del(str, hash_to_clean, k);
		}
	}
	//printf("before clean hash string\n");

	if (hash_to_clean) kh_destroy(str, hash_to_clean);
	//printf("after clean hash string\n");
}

uint32_t getChromIndexFromID(bam_hdr_t *header, char *chrom_id) {
	uint32_t i=0;
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
		//need to increase the size if we need to analysis any other chromosomes (such as decoy)
        //chrom_tracking->chromosome_ids[i] = calloc(25, sizeof(char));
        chrom_tracking->chromosome_ids[i] = NULL;

	//chrom_tracking->chromosome_lengths = calloc(25, sizeof(uint32_t));
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
	chrom_tracking->chromosome_ids[index] = calloc(50, sizeof(char));
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

uint32_t locateChromosomeIndex(char *chrom_id, Chromosome_Tracking *chrom_tracking) {
	uint32_t i=0;
    for (i = 0; i < chrom_tracking->number_tracked; i++) {
		if (chrom_tracking->chromosome_ids[i]) {
			if (strcmp(chrom_id, chrom_tracking->chromosome_ids[i]) == 0) {
				return i;
			}
        }
    }

	fprintf(stderr, "Something is wrong because the chromosome %s couldn't be found\n", chrom_id);
	return 0;
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

void zeroAllNsRegions(char *chrom_id, Bed_Info *Ns_info, Chromosome_Tracking *chrom_tracking) {
	// First, we need to find the index that is used to track current chromosome chrom_id
	uint32_t idx = locateChromosomeIndex(chrom_id, chrom_tracking);
   	uint32_t i=0,j=0;

	for (i=0; i<Ns_info->size; i++) {
		if (strcmp(Ns_info->coords[i].chrom_id, chrom_id) == 0) {
			//printf("%s\t%"PRIu32"\t%"PRIu32"\n", Ns_info->coords[i].chr, Ns_info->coords[i].start, Ns_info->coords[i].end);
			for (j=Ns_info->coords[i].start; j<=Ns_info->coords[i].end; j++) {
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
