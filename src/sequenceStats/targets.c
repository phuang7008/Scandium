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
    uint32_t total_size = loadBedFiles(user_inputs->database_version, bedfile_name, bed_info->coords, wanted_chromosome_hash);
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
