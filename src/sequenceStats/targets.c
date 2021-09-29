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

void processBedFiles(User_Input *user_inputs, Bed_Info *bed_info, khash_t(khStrInt)* wanted_chromosome_hash, char* bedfile_name) {
    // First, let's get the total number of lines(items or count) within the target file
    //
    bed_info->size = getLineCount(bedfile_name);

    // Now initialize the storage for arrays that are used to store target coordinates
    // Do need to remember to free these allocated memories at the end of the program!!!
    //
    bed_info->coords = calloc(bed_info->size, sizeof(Bed_Coords));
    
    // load target file or Ns bed file again and store the information (starts, stops and chromosome ids)
    //
    loadBedFiles(user_inputs->database_version, bedfile_name, bed_info->coords, wanted_chromosome_hash, user_inputs->database_version);
    //printf("total size is %"PRIu32"\n", total_size);
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
void generateBedBufferStats(Bed_Info * bed_info, Stats_Info *stats_info, Target_Buffer_Status *target_buffer_status, int32_t target_buffer_index, User_Input *user_inputs, char* chrom_id, short target_file_index, short type) {
    if (target_buffer_index == -1) return;

    uint32_t i=0, j=0;
    uint32_t chrom_len = target_buffer_status[target_buffer_index].size;  // used to check the boundary
    target_buffer_status[target_buffer_index].index = target_buffer_index;
    uint32_t prev_end = 0;

    for (i = 0; i < bed_info->size; i++) {

        if (strcmp(bed_info->coords[i].chrom_id, chrom_id) != 0)
            continue;
        
        // check if bed file is sorted and merged
        //
        if ((prev_end != 0) && (prev_end > bed_info->coords[i].start)) {
            char error_message[250];
            sprintf(error_message, "The file %s is not sorted or merged.\n", user_inputs->target_files[target_file_index]);
            exitWithFailure(error_message);
        } else {
            prev_end = bed_info->coords[i].end;
        }

        // for positions on targets or Ns
        // since the chromosome position array is 0 based, we will just use the 0-based bed file as is
        //
        for (j=bed_info->coords[i].start; j<bed_info->coords[i].end; j++) {
            if (j >= chrom_len) continue;

            if (type == 1) {
                setTargetBufferStatus(target_buffer_status, target_buffer_index, j, target_file_index, 1);
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
                    target_buffer_index, target_buffer_status, target_file_index);

            // for the buffer positions at the right hand side
            //
            processBufferRegions(bed_info->coords[i].end+1, bed_info->coords[i].end+user_inputs->target_buffer_size+1,
                    target_buffer_index, target_buffer_status, target_file_index);
        }
    }

    if (type == 1) {
        uint8_t cur_target_buffer_bit = getTargetBufferBit(target_file_index);

        for (j=0; j<target_buffer_status[target_buffer_index].size; j++) {
            // update the target/buffer stats here
            //
            if (target_buffer_status[target_buffer_index].target_status_array[j] & cur_target_buffer_bit) {
                // the following addition need to be in CRITICAL step 
                //
                stats_info->capture_cov_stats[target_file_index]->total_targeted_bases += 1;
                if (target_buffer_status[target_buffer_index].buffer_status_array[j] & cur_target_buffer_bit)
                    target_buffer_status[target_buffer_index].buffer_status_array[j] ^= cur_target_buffer_bit;
            }

            if (target_buffer_status[target_buffer_index].buffer_status_array[j] & cur_target_buffer_bit)
                // the following addition need to be in CRITICAL step
                //
                stats_info->capture_cov_stats[target_file_index]->total_buffer_bases += 1;
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

void setTargetBufferStatus(Target_Buffer_Status *target_buffer_status, int target_buffer_index, uint32_t pos_idx, uint8_t target_file_index, uint8_t type) {
    uint8_t target_buffer_bit = getTargetBufferBit(target_file_index);

    if (type == 1) {
        if ((target_buffer_status[target_buffer_index].target_status_array[pos_idx] & target_buffer_bit) == 0) 
            target_buffer_status[target_buffer_index].target_status_array[pos_idx] |= target_buffer_bit;
    } else {
        if ((target_buffer_status[target_buffer_index].buffer_status_array[pos_idx] & target_buffer_bit) == 0)
            target_buffer_status[target_buffer_index].buffer_status_array[pos_idx] |= target_buffer_bit;
    }
}

void processBufferRegions(uint32_t start, uint32_t end, int target_buffer_index, Target_Buffer_Status *target_buffer_status, short target_file_index) {
    uint32_t j;
    for (j=start; j < end; j++) {
        //if (j < 0) continue;      // Removed! As it is always False!
        if (j >= target_buffer_status[target_buffer_index].size) continue;

        setTargetBufferStatus(target_buffer_status, target_buffer_index, j, target_file_index, 2);
        
    }
}
