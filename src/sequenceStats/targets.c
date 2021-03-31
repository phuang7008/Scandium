/*
 * =====================================================================================
 *
 *        Filename:        targets.c
 *
 *        Description:    The implementation of the target.h file
 *
 *      Version:        1.0
 *      Created:        02/06/2017 04:45:04 PM
 *      Revision:        none
 *      Compiler:        gcc
 *
 *      Author:            Peiming (Peter) Huang (phuang@bcm.edu)
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */

#include "targets.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>        // for file access()
#include <time.h>
#include "utils.h"

// It will open the bed-formatted file for the first time. and count the number of lines (items) within it.
//
uint32_t getLineCount(char *bed_file) {
    FILE *fp = fopen(bed_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("ERROR: target file \n%s\n open failed!", bed_file);
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
uint32_t loadBedFiles(User_Input *user_inputs, char *bed_file, Bed_Coords * coords, khash_t(khStrInt)* wanted_chromosome_hash) {

    FILE *fp = fopen(bed_file, "r");;

    if (fp == NULL) {       // ensure the target file open correctly!
        printf("ERROR: target or N region bed file \n%s\n open failed!", bed_file);
        exit(EXIT_FAILURE);
    }

    // setup a variable to store chromosomes that have been seen
    // this is used to check if the bed file is sorted!
    // if the file is not sorted, exit and give an error message!
    //
    khash_t(khStrInt) *seen_chromosome_hash = kh_init(khStrInt);
    char *prev_chr_id = calloc(150, sizeof(char));
    strcpy(prev_chr_id, "nothing99999");
    uint32_t prev_start=0;
    bool sorted=true;

    uint32_t count=0;        // index for each target line (item)
    uint32_t total_size=0;    // used for target bed input file merge and uniq checking
    ssize_t read;
    size_t len = 0;
    char *p_token=NULL, *line=NULL;    // here p_ means a point. so p_token is a point to a token
    char *savePtr;

    while((read = getline(&line, &len, fp)) != -1) {
        //printf("%s\n", line);
        if (line[0] == '\0') continue;        // skip if it is a blank line
        if (line[0] == '#')  continue;        // skip if it is a comment line

        savePtr = line;

        // to keep track the tokens from string split using strtok_r()
        //
        uint32_t i = 0;

        // loop through the rest of items, but here we are only interested in start and end position
        //
        while ((p_token = strtok_r(savePtr, "\t", &savePtr))) {
            // get the first item, which is chromosome id
            //
            if (i == 0) {
                // skip if the chromosome is not going to be processed
                //
                if (wanted_chromosome_hash != NULL) {
                    khiter_t iter_p = kh_get(khStrInt, wanted_chromosome_hash, p_token);
                    if (iter_p == kh_end(wanted_chromosome_hash))
                        break;
                }

                checkChromosomeID(user_inputs, p_token);       // check chromosome id for versioning
                strcpy(coords[count].chrom_id,  p_token);
            }

            if (i == 1)
                coords[count].start = (uint32_t) strtol(p_token, NULL, 10);
            
            if (i == 2)
                coords[count].end = (uint32_t) strtol(p_token, NULL, 10);

            i++;
        }

        if (coords[count].chrom_id == NULL)
            continue;

        //if (coords[count].chrom_id && strcmp(coords[count].chrom_id, '\0') == 0)
        if (coords[count].chrom_id && strlen(coords[count].chrom_id) == 0)
            continue;

        // checking if the bed file is sorted
        //
        if (strcmp(coords[count].chrom_id, prev_chr_id) == 0) {
            if (prev_start >= coords[count].start)
                sorted = false;
        } else {
            // its a new chromosome
            // need to check if we have seen the chromosome id before
            //
            int absent;
            khiter_t iter = kh_put(khStrInt, seen_chromosome_hash, coords[count].chrom_id, &absent);
            if (absent) {
                // haven't seen it, just add to the seen_chromosome_hash
                //
                kh_key(seen_chromosome_hash, iter) = strdup(coords[count].chrom_id);
                kh_value(seen_chromosome_hash, iter) = 1;

                strcpy(prev_chr_id, coords[count].chrom_id);
            } else {
                sorted = false;
            }
        }

        if (!sorted) {
            fprintf(stderr, "ERROR: The input bed file \n%s\n", bed_file);
            fprintf(stderr, "\tis not sorted between chrom %s and %s!\n", prev_chr_id, coords[count].chrom_id);
            fprintf(stderr, "\tis not sorted between %"PRIu32" and %"PRIu32"!\n", prev_start, coords[count].start);
            fprintf(stderr, "Please make sure your input bed file is sorted, merged and unique!\n");
            exit(EXIT_FAILURE);
        }

        // reset the prev_start position and sum the total bed file size
        //
        prev_start = coords[count].start;
        total_size += coords[count].end - coords[count].start;

        count++;
    }

    if (line != NULL) free(line);
    if (prev_chr_id != NULL) free(prev_chr_id);
    if (seen_chromosome_hash != NULL) 
        cleanKhashStrInt(seen_chromosome_hash);

    fclose(fp);

    return total_size;
}

void processBedFiles(User_Input *user_inputs, Bed_Info *bed_info, Stats_Info *stats_info, Target_Buffer_Status *target_buffer_status, khash_t(khStrInt)* wanted_chromosome_hash, char* bedfile_name, int number_of_chromosomes, short target_file_index, short type) {
    // First, let's get the total number of lines(items or count) within the target file
    //
    bed_info->size = getLineCount(bedfile_name);

    // Now initialize the storage for arrays that are used to store target coordinates
    // Do need to remember to free these allocated memories at the end of the program!!!
    //
    bed_info->coords = calloc(bed_info->size, sizeof(Bed_Coords));
    
    // load target file or Ns bed file again and store the information (starts, stops and chromosome ids)
    //
    uint32_t total_size = loadBedFiles(user_inputs, bedfile_name, bed_info->coords, wanted_chromosome_hash);
    //printf("total size is %"PRIu32"\n", total_size);

    // Now we are going to generate target-buffer lookup table for all the loaded targets
    // we will store targets and buffers information based on chromosome ID
    //
    generateBedBufferStats(bed_info, stats_info, target_buffer_status, user_inputs, wanted_chromosome_hash, number_of_chromosomes, target_file_index, type);

    // Here we need to check if the bed file is merged and uniqued by comparing two different ways of addition of bases
    //
    if (type == 1) {
        //printf("Total target is %"PRIu32"\n", stats_info->cov_stats->total_targeted_bases);
        if (total_size != stats_info->capture_cov_stats[target_file_index]->total_targeted_bases) {
            printf("\nERROR: The target bed file \n%s\n needs to be bedtools sorted, merged and uniqued!\n", bedfile_name);
            printf("\ttotal size: %"PRIu32"\n", total_size);
            printf("\ttotal target bases: %"PRIu32"\n", stats_info->capture_cov_stats[target_file_index]->total_targeted_bases);
            //printf("\tThe chromosome ids in the capture file MUST appear in the chromosome input file (--chr_list option)!\n");
            exit(EXIT_FAILURE);
        }
    } else {
        //printf("Total Ns region is %"PRIu32"\n", stats_info->cov_stats->total_Ns_bases);
        if (total_size != stats_info->wgs_cov_stats->total_Ns_bases) {
            printf("\nERROR: The Ns-region bed file \n%s\n needs to be bedtools sorted, merged and uniqued!\n", bedfile_name);
            printf("\ttotal size: %"PRIu32"\n", total_size);
            printf("\ttotal N bases: %"PRIu32"\n", stats_info->wgs_cov_stats->total_Ns_bases);
            //printf("\tThe chromosome ids listed in the Ns-region file MUST appear in the chromosome input file (--chr_list option)!\n");
            exit(EXIT_FAILURE);
        }
    }
}

void outputForDebugging(Bed_Info *bed_info) {
    // Output stored coordinates for verification
    uint32_t i=0;
    for (i=0; i<bed_info->size; i++) {
        printf("%s\t%"PRIu32"\t%"PRIu32"\n", bed_info->coords[i].chrom_id, bed_info->coords[i].start, bed_info->coords[i].end);
    }
}

// Here type 1 refers to target bed, while type 2 means Ns region bed
// For values stored in the status_array:
// 1: target        2: buffer        3: Ns        4: 1+3 (target+Ns overlaps)        5: 2+3 (buffer+Ns overlaps)
//
void generateBedBufferStats(Bed_Info * bed_info, Stats_Info *stats_info, Target_Buffer_Status *target_buffer_status, User_Input *user_inputs, khash_t(khStrInt)* wanted_chromosome_hash, int number_of_chromosomes, short target_file_index, short type) {
    uint32_t i=0, j=0, k=0, chrom_len=0;
    int idx = -1;
    char cur_chrom_id[50];
    strcpy(cur_chrom_id, "something");

    for (i = 0; i < bed_info->size; i++) {

        // skip if the chromosome is not going to be processed
        //
        if (wanted_chromosome_hash != NULL) {
            khiter_t iter_p = kh_get(khStrInt, wanted_chromosome_hash, bed_info->coords[i].chrom_id);
            if (iter_p == kh_end(wanted_chromosome_hash))
                continue;
        }

        if (strcmp(bed_info->coords[i].chrom_id, cur_chrom_id) != 0) {
            strcpy(cur_chrom_id, bed_info->coords[i].chrom_id);

            // get the index for the target_buffer_status
            //
            for(k=0; k<(uint32_t)number_of_chromosomes; k++) {
                if (strcmp(target_buffer_status[k].chrom_id, cur_chrom_id) == 0) {
                    idx = k;
                    chrom_len = target_buffer_status[k].size;
                    target_buffer_status[k].index = k;
                    break;
                }
            }

        }

        if (idx == -1) return;

        // for positions on targets or Ns
        // since the chromosome position array is 0 based, we will just use the 0-based bed file as is
        //
        for (j=bed_info->coords[i].start; j<bed_info->coords[i].end; j++) {
            if (j >= chrom_len) continue;

            if (type == 1) {
                setTargetBufferStatus(target_buffer_status, idx, j, target_file_index, 1);
            }

            if (type == 2 && j < bed_info->coords[i].end) {
                if (j >= chrom_len) continue;

                stats_info->wgs_cov_stats->total_Ns_bases += 1;

                if ((strcmp(bed_info->coords[i].chrom_id, "chrX") == 0) || (strcmp(bed_info->coords[i].chrom_id, "X") == 0))
                    stats_info->wgs_cov_stats->total_Ns_bases_on_chrX += 1;

                if ((strcmp(bed_info->coords[i].chrom_id, "chrY") == 0) || (strcmp(bed_info->coords[i].chrom_id, "Y") == 0))
                    stats_info->wgs_cov_stats->total_Ns_bases_on_chrY += 1;
            }
        }

        // If it is for target bed file, we need to handle the buffer positions on both sides of a target
        //
        if (type == 1) {
            // for the buffer positions at the left-hand side
            //
            processBufferRegions(bed_info->coords[i].start-user_inputs->target_buffer_size, bed_info->coords[i].start,
                    idx, target_buffer_status, target_file_index);

            // for the buffer positions at the right hand side
            //
            processBufferRegions(bed_info->coords[i].end+1, bed_info->coords[i].end+user_inputs->target_buffer_size+1,
                    idx, target_buffer_status, target_file_index);
        }
    }

    if (type == 1) {
        uint8_t cur_target_buffer_bit = getTargetBufferBit(target_file_index);

        for (i=0; i<(uint32_t)number_of_chromosomes; i++) {
            if (target_buffer_status[i].index == -1)
                continue;

            for (j=0; j<target_buffer_status[i].size; j++) {
                // update the target/buffer stats here
                //
                if (target_buffer_status[i].target_status_array[j] & cur_target_buffer_bit) {
                    stats_info->capture_cov_stats[target_file_index]->total_targeted_bases += 1;
                    if (target_buffer_status[i].buffer_status_array[j] & cur_target_buffer_bit)
                        target_buffer_status[i].buffer_status_array[j] ^= cur_target_buffer_bit;
                }

                if (target_buffer_status[i].buffer_status_array[j] & cur_target_buffer_bit)
                    stats_info->capture_cov_stats[target_file_index]->total_buffer_bases += 1;
            }
        }
    }
}

uint8_t getTargetBufferBit(uint8_t target_file_index) {
    if (target_file_index == 0) { return TRT_BFR_1; }
    else if (target_file_index == 1){ return TRT_BFR_2; }
    else if (target_file_index == 2) { return TRT_BFR_3; }
    else if (target_file_index == 3) { return TRT_BFR_4; }
    else if (target_file_index == 4) { return TRT_BFR_5; }
    else if (target_file_index == 5) { return TRT_BFR_6; }
    else if (target_file_index == 6) { return TRT_BFR_7; }
    else if (target_file_index == 7) { return TRT_BFR_8; }
    else { printf("ERROR: Something is wrong as target file index > 8\n"); exit(EXIT_FAILURE); }
}

void setTargetBufferStatus(Target_Buffer_Status *target_buffer_status, int chrom_idx, uint32_t pos_idx, uint8_t target_file_index, uint8_t type) {
    uint8_t target_buffer_bit = getTargetBufferBit(target_file_index);

    if (type == 1) {
        if ((target_buffer_status[chrom_idx].target_status_array[pos_idx] & target_buffer_bit) == 0) 
            target_buffer_status[chrom_idx].target_status_array[pos_idx] |= target_buffer_bit;
    } else {
        if ((target_buffer_status[chrom_idx].buffer_status_array[pos_idx] & target_buffer_bit) == 0)
            target_buffer_status[chrom_idx].buffer_status_array[pos_idx] |= target_buffer_bit;
    }
}

void processBufferRegions(uint32_t start, uint32_t end, int chrom_idx, Target_Buffer_Status *target_buffer_status, short target_file_index) {
    uint32_t j;
    for (j=start; j < end; j++) {
        //if (j < 0) continue;      // Removed! As it is always False!
        if (j >= target_buffer_status[chrom_idx].size) continue;

        setTargetBufferStatus(target_buffer_status, chrom_idx, j, target_file_index, 2);
        
    }
}
