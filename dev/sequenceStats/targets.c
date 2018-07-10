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
#include "utils.h"

// It will open the bed-formatted file for the first time. and count the number of lines (items) within it.
//
uint32_t getLineCount(char *bed_file) {
    FILE *fp = fopen(bed_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("target file %s open failed!", bed_file);
        exit(EXIT_FAILURE);
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
        exit(EXIT_FAILURE);
    }

	uint32_t count = 0;		// index for each target line (item)
	ssize_t read;
	size_t len = 0;
	char *p_token=NULL, *line=NULL;	// here p_ means a point. so p_token is a point to a token
	char *savePtr;

    while((read = getline(&line, &len, fp)) != -1) {
		//printf("%s\n", line);
		if (line[0] == '\0') continue;		// skip if it is a blank line
		if (line[0] == '#')  continue;		// skip if it is a comment line

		savePtr = line;

		// to keep track the tokens from string split using strtok_r()
		//
		uint32_t i = 0;

		// loop through the rest of items, but here we are only interested in start and end position
		//
		while ((p_token = strtok_r(savePtr, "\t", &savePtr))) {
			// get the first item, which is chromosome id
			//
			if (i == 0)
				strcpy(coords[count].chrom_id,  p_token);

			if (i == 1)
				coords[count].start = (uint32_t) strtol(p_token, NULL, 10);
			
			if (i == 2)
				coords[count].end = (uint32_t) strtol(p_token, NULL, 10);

			i++;
		}
		count++;
	}

	if (line != NULL) free(line);

	fclose(fp);
}

void processBedFiles(User_Input *user_inputs, Bed_Info *bed_info, Stats_Info *stats_info, Target_Buffer_Status *target_buffer_status,  bam_hdr_t *header, short type) {
	// First, let's get the total number of lines(items or count) within the target file
	//
	if (type == 1) {
		bed_info->size = getLineCount(user_inputs->target_file);
	} else {
		bed_info->size = getLineCount(user_inputs->n_file);
	}

	// Now initialize the storage for arrays that are used to store target coordinates
    // Do need to remember to free these allocated memories at the end of the program!!!
	//
    bed_info->coords = calloc(bed_info->size, sizeof(Bed_Coords));
	
	// load target file or Ns bed file again and store the information (starts, stops and chromosome ids)
	//
	if (type == 1) {
		loadBedFiles(user_inputs->target_file, bed_info->coords);
	} else {
		loadBedFiles(user_inputs->n_file, bed_info->coords);
	}

    // Now we are going to generate target-buffer lookup table for all the loaded targets
    // we will store targets and buffers information based on chromosome ID
	//
	generateBedBufferStats(bed_info, stats_info, target_buffer_status, header, user_inputs, type);
}

void outputForDebugging(Bed_Info *bed_info) {
	// Output stored coordinates for verification
    uint32_t i=0;
    for (i=0; i<bed_info->size; i++) {
        printf("%s\t%"PRIu32"\t%"PRIu32"\n", bed_info->coords[i].chrom_id, bed_info->coords[i].start, bed_info->coords[i].end);
    }
}

// Here type 1 refers to target bed, while teyp 2 means Ns region bed
// For values stored in the status_array:
// 1: target		2: buffer		3: Ns		4: 1+3 (target+Ns overlaps)		5: 2+3 (buffer+Ns overlaps)
//
void generateBedBufferStats(Bed_Info * bed_info, Stats_Info *stats_info, Target_Buffer_Status *target_buffer_status, bam_hdr_t *header, User_Input *user_inputs, short type) {
	uint32_t i=0, j=0, k=0, chrom_len=0;
	int idx = -1;
	char cur_chrom_id[50];
	strcpy(cur_chrom_id, "something");

	for (i = 0; i < bed_info->size; i++) {

		//if (i >= 100986)
		//	printf("alt in\n");

        if (strcmp(bed_info->coords[i].chrom_id, cur_chrom_id) != 0) {
			strcpy(cur_chrom_id, bed_info->coords[i].chrom_id);

			// get the index for the target_buffer_status
			for(k=0; k<header->n_targets; k++) {
				if (strcmp(target_buffer_status[k].chrom_id, cur_chrom_id) == 0) {
					idx = k;
					chrom_len = target_buffer_status[k].size;
					target_buffer_status[k].index = k;
					break;
				}
			}

        }

		if (idx == -1) return;

		uint8_t c_type = type == 1 ? 1 : 3;

		// for positions on targets or Ns
		for (j=bed_info->coords[i].start; j<=bed_info->coords[i].end; j++) {
			if (j >= chrom_len) continue;

			if (type == 1) {
			   	if (target_buffer_status[idx].status_array[j] == 3) { 
					target_buffer_status[idx].status_array[j] = 4;
				} else {
					target_buffer_status[idx].status_array[j] = c_type;
				}
			}

			// The old Java version include the end position in the bed file. 
			// So I will include it here as well. 
			// In the future, need to make sure if the bed file format include end position or not!
			//
			//if (type == 2 && j < bed_info->coords[i].end) {
			if (type == 2 && j <= bed_info->coords[i].end) {
				if (j >= chrom_len) continue;

				stats_info->cov_stats->total_Ns_bases += 1;

				if (target_buffer_status[idx].status_array[j] == 1) {
					target_buffer_status[idx].status_array[j] = 4;
				} else {
					target_buffer_status[idx].status_array[j] = c_type;
				}
			}
		}

		if (type == 1) {
			// for the buffer positions at the left side
			//
			for (j=bed_info->coords[i].start-user_inputs->target_buffer_size; j < bed_info->coords[i].start; j++) {
				if (j < 0) continue;
				if (j >= chrom_len) continue;

				if (target_buffer_status[idx].status_array[j] == 0) {
					target_buffer_status[idx].status_array[j] = 2;
				} else if (target_buffer_status[idx].status_array[j] == 3) {
					target_buffer_status[idx].status_array[j] = 5;
				} else {
					// target_buffer_status[idx].status_array[j] == 1
					continue;
				}
			}

			// for the buffer positions at the right side
			//
			for (j=bed_info->coords[i].end+1; j <= bed_info->coords[i].end+user_inputs->target_buffer_size; j++ ) {
				if (j >= chrom_len) continue;
				
				if (target_buffer_status[idx].status_array[j] == 0) {
					target_buffer_status[idx].status_array[j] = 2;
				} else if (target_buffer_status[idx].status_array[j] == 3) {
					target_buffer_status[idx].status_array[j] = 5;
				} else {
					// target_buffer_status[idx].status_array[j] == 1
					continue;
				}
			}
		}
	}

	for (i=0; i<header->n_targets; i++) {
		if (target_buffer_status[i].index == -1)
			continue;

		for (j=0; j<target_buffer_status[i].size; j++) {
		   	// update the target/buffer stats here
			//
			if (type == 1) {
	        	if (target_buffer_status[i].status_array[j] == 1 || target_buffer_status[i].status_array[j] == 4)
		        	stats_info->cov_stats->total_targeted_bases += 1;
	
				if (target_buffer_status[i].status_array[j] == 2 || target_buffer_status[i].status_array[j] == 5)
		        	stats_info->cov_stats->total_buffer_bases += 1;
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
