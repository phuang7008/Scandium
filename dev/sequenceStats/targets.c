/*
 * =====================================================================================
 *
 *		Filename:		targets.c
 *
 *		Description:	For the base coverage calculation
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

#include "targets.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>		// for file access()
#include <time.h>

// It will open the bed-formatted file for the first time. and count the number of lines (items) within it.
//
uint32_t getLineCount(char *bed_file) {
    FILE *fp = fopen(bed_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("target file %s open failed!", bed_file);
        exit(1);
    }

	// Extract characters from file and store in character c
    char c;         // To store a character read from file
    long count = 0; // To store total line number in the target file

    for (c = getc(fp); c != EOF; c = getc(fp)) {
        if (c == '\n')      // Increment count if this character is newline
            count = count + 1;
    }

    fclose(fp);     // Close the file

    return count;
}

// It will open the bed-formatted file and then record all the start and stop positions along with chromosome ids
// for chromosome X and Y are characters, I will use array of chars (ie array of string) to store chromosome ids
//
void loadBedFiles(char * bed_file, Bed_Coords * coords) {
	FILE *fp = fopen(bed_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("target file %s open failed!", bed_file);
        exit(1);
    }

	int num = 0;		// to keep track the tokens from string split using strtok()
	int count = 0;		// index for each target line (item)
	ssize_t read;
	size_t len = 0;
	char *p_token, *line;	// here p_ means a point. so p_token is a point to a token

    while((read = getline(&line, &len, fp)) != -1) {
		//printf("%s\n", line);
		if (line[0] == '\0') continue;		// skip if it is a blank line
		if (line[0] == '#')  continue;		// skip if it is a comment line

		p_token = strtok(line, " \t");
		while (p_token != NULL) {
			if (num == 0) 
				strcpy(coords[count].chr,  p_token);
			else if (num == 1)
				coords[count].start = atoi(p_token);
			else if (num == 2) {
				coords[count].end = atoi(p_token);
				num = 0;
				p_token = NULL;
				break;		// stop 'while' loop as we don't need anything after stop position
			}

			p_token = strtok(NULL, " \t");	// In subsequent calls, strtok expects a null pointer
			num++;
		}
		count++;
	}

	fclose(fp);
	if (line != NULL)
		free(line);
	if (p_token != NULL)
		free(p_token);
}

void processBedFiles(char *bed_file, Bed_Info *bed_info, khash_t(str) *bed_buffer_hash, Stats_Info *stats_info, bam_hdr_t *header, short type) {
	// First, let's get the total number of lines(items or count) within the target file
	bed_info->size = getLineCount(bed_file);

	// Now initialize the storage for arrays that are used to store target coordinates
    // Do need to remember to free these allocated memories at the end of the program!!!
    bed_info->coords = calloc(bed_info->size, sizeof(Bed_Coords));
	
    // load target file again and store the target information (starts, stops and chromosome ids)
    loadBedFiles(bed_file, bed_info->coords);

    // Now we are going to generate target-buffer lookup table for all the loaded targets
    // we will store targets and buffers information based on chromosome ID
    // we will have 22 + X + Y = 24 chromosomes, Here X=23 and Y will be 24
	// Right now, we only need to do it for target file. 
	if (type == 1)
		generateBedBufferLookupTable(bed_info, bed_buffer_hash, stats_info, header, type);
}

void outputForDebugging(Bed_Info *bed_info, khash_t(str) *bed_buffer_hash) {
	// Check to see if target and buffer is set correctly
    khiter_t outer_iter, inner_iter;
    
    for (outer_iter=kh_begin(bed_buffer_hash); outer_iter!=kh_end(bed_buffer_hash); outer_iter++) {
        if (kh_exist(bed_buffer_hash, outer_iter)) {
            // find the key
            printf("Current key is %s\n", kh_key(bed_buffer_hash, outer_iter));
            if (strcmp("1", kh_key(bed_buffer_hash, outer_iter)) == 0) {
                for (inner_iter=kh_begin(kh_value(bed_buffer_hash, outer_iter));
                            inner_iter!=kh_end(kh_value(bed_buffer_hash, outer_iter));
                                inner_iter++) {
					if (kh_exist(kh_value(bed_buffer_hash, outer_iter), inner_iter)) {
						uint32_t pos = kh_key(kh_value(bed_buffer_hash, outer_iter), inner_iter);
						//if ( (pos >= 1219000) && (pos<=1220290) )  {        // for region: 1220086 1220186
							printf("%"PRIu32"\t%d %s\n",pos, kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter), kh_key(bed_buffer_hash, outer_iter));
						//}
                    }
                }
				printf("Finished writing...\n");
            }
        }
    }

	// Output stored coordinates for verification
    uint32_t i=0;
    for (i=0; i<bed_info->size; i++) {
        printf("%s\t%"PRIu32"\t%"PRIu32"\n", bed_info->coords[i].chr, bed_info->coords[i].start, bed_info->coords[i].end);
    }
}

void generateBedBufferLookupTable(Bed_Info * bed_info, khash_t(str) *bed_buffer_hash, Stats_Info *stats_info, bam_hdr_t *header, short type) {
	int ret=0, absent=0;
	uint32_t i=0, j=0, chrom_len=0;
	khiter_t outer_iter, inner_iter;
	char cur_chrom_id[50];
	strcpy(cur_chrom_id, "something");

	for (i = 0; i < bed_info->size; i++) {
		//printf("Processing id %d\n", i);

		// check to see if we need to initialize the hash for current chromosome
		if (strcmp(bed_info->coords[i].chr, cur_chrom_id) != 0) {
			khash_t(m32) *tmp_hash = kh_init(m32);
			//outer_iter = kh_put(str, bed_buffer_hash, bed_info->coords[i].chr, &absent);
			//if (absent)		// insert the key if there is no key
			//	kh_key(bed_buffer_hash, outer_iter) = strdup(bed_info->coords[i].chr);
			outer_iter = kh_put(str, bed_buffer_hash, strdup(bed_info->coords[i].chr), &absent);
			kh_value(bed_buffer_hash, outer_iter) = tmp_hash;

			strcpy(cur_chrom_id, bed_info->coords[i].chr);
			//printf("hash is added for chromosome: %s\n", coords[i].chr);

			// get the current chromosome id's length from header file
    		for (j=0; j<header->n_targets; j++) {
		        if (strcmp(cur_chrom_id, header->target_name[j]) == 0) {
        		    chrom_len = header->target_len[j];
					//printf("Get Chromosome %s length %"PRIu32"\n", header->target_name[j], chrom_len);
					break;
        		}
    		}
		}

		// now get the hash table for current chromosome id
		outer_iter = kh_put(str, bed_buffer_hash, bed_info->coords[i].chr, &ret);
		if (ret)
			fprintf(stderr, "Something went wrong, has key for %s is not added", bed_info->coords[i].chr);

		// type == 1 for target, type == 2 for Ns regions in the reference sequence
		uint8_t c_type = type == 1 ? 1 : 3;

		// for positions on targets or Ns
		for (j=bed_info->coords[i].start; j<=bed_info->coords[i].end; j++) {
			if (j > chrom_len) continue;

			inner_iter = kh_put(m32, kh_value(bed_buffer_hash, outer_iter), j, &ret);	// add the key
			kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) = c_type;		// add the value

			if (type == 1)
				stats_info->cov_stats->total_targeted_bases += 1;

			if (type == 2)
				stats_info->cov_stats->total_Ns_bases += 1;
		}

		if (type == 1) {
			// for positions on the buffer at the left side
			for (j=bed_info->coords[i].start-BUFFER; j < bed_info->coords[i].start; j++) {
				if (j < 0) continue;

				inner_iter = kh_put(m32, kh_value(bed_buffer_hash, outer_iter), j, &ret);
				if (ret == 1) {
					kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) = 2;
                    //stats_info->cov_stats->total_buffer_bases += 1;
				} else if (kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) == 2) {
					continue;
				} else if (kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) != 1) {
					kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) = 2;
					//stats_info->cov_stats->total_buffer_bases += 1;
				}
			}

			// for positions on the buffer at the right side
			for (j=bed_info->coords[i].end+1; j <= bed_info->coords[i].end+BUFFER; j++ ) {
				if (j >= chrom_len) continue;

				inner_iter = kh_put(m32, kh_value(bed_buffer_hash, outer_iter), j, &ret);
				if (ret == 1) {
					kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) = 2;
                    //stats_info->cov_stats->total_buffer_bases += 1;
				} else if (kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) == 2) {
					continue;
				} else if (kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) != 1) {
					kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) = 2;
					//stats_info->cov_stats->total_buffer_bases += 1;
				}
			}
		}
	}

	// now need to find out the number of buffer bases
	if (type == 1) {
		for(outer_iter = kh_begin(bed_buffer_hash); outer_iter!=kh_end(bed_buffer_hash); outer_iter++) {
			if (kh_exist(bed_buffer_hash, outer_iter)) {
				for (inner_iter = kh_begin(kh_value(bed_buffer_hash, outer_iter)); 
						inner_iter != kh_end(kh_value(bed_buffer_hash, outer_iter));
							inner_iter++) {
					if (kh_exist(kh_value(bed_buffer_hash, outer_iter), inner_iter)) {
						if (kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) == 2)
							stats_info->cov_stats->total_buffer_bases += 1;
					}
				}
			}
		}
	}
}

void cleanBedInfo(Bed_Info *bed_info) {
	if (bed_info) {
		if (bed_info->coords)
			free(bed_info->coords);
		free(bed_info);
	}
}
