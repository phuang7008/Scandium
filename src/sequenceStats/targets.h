/*
 * =====================================================================================
 *
 *      Filename:		targets.h
 *
 *		Description:	For the Capture/Target related functionalities
 *
 *      Version:		1.0
 *      Created:		02/06/2017 04:45:04 PM
 *      Revision:		none
 *      Compiler:		gcc
 *
 *      Author:			Peiming (Peter) Huang (phuang@bcm.edu)
 *      Company:		Baylor College of Medicine
 *
 * =====================================================================================
 */

#ifndef TARGETS_H
#define TARGETS_H

#include <stdbool.h>
#include <stdio.h>
#include "htslib/sam.h"
#include "terms.h"
#include "user_inputs.h"

/**
 * generate target and buffer lookup table for quick access. The khash.h file is used for hash table setup
 * @param bed_info: the bed information that is stored for the future usage
 * @param stat_info, statistical information for the reads/bases
 * @param user_inputs, contains all the user inputs information, including target_buffer_size
 * @param type: either target bed (type 1) or Ns regions in the reference sequences (type 2)
 */
void generateBedBufferStats(Bed_Info * bed_info, Stats_Info *stats_info, Target_Buffer_Status *target_status, int32_t target_buffer_index, User_Input *user_inputs, char* chrom_id, short target_file_index, short type);

/**
 * process bed-formatted file and populate the coordinates and lookup hash table
 * @param user_inputs: contains all the user inputs including the target or Ns bed file names
 * @param bed_info: the storage of bed coordinates and the size of the bed file
 * @param stats_info: a variable that contains various statistical information
 * @param target_buffer_status: a variable to store target/buffer info
 * @param type: either target bed (type 1) or Ns regions in the reference sequences (type 2)
 */
void processBedFiles(User_Input *user_inputs, Bed_Info *bed_info, khash_t(khStrInt)* wanted_chromosome_hash, char* bedfile_name);

uint8_t getTargetBufferBit(uint8_t target_file_index);

void setTargetBufferStatus(Target_Buffer_Status *target_buffer_status, int chrom_idx, uint32_t pos_idx, uint8_t target_file_index, uint8_t type);

void processBufferRegions(uint32_t start, uint32_t end, int chrom_idx, Target_Buffer_Status *target_buffer_status, short target_file_index);

/**
 * just to output some information for debugging
 */
void outputForDebugging(Bed_Info *bed_info);

/**
 * It is used to store both target and buffer information
 * @param chromosome_id: the current chromosome to load
 * @param size: the size of the chromosome
 * @param target_buffer_regions: the 2D char array that store the flags either targets('T') or their surrounding buffers ('B')
 * @param target_coords: the struct array that stores the coordinates of targets
 * @return
 */
void getTargetAndBufferPositions(char * chromosome_id, int size, char *target_buffer_regions, Bed_Coords *coords);

void processTargetAnnotationFile(User_Input *user_inputs);

#endif //TARGETS_H
