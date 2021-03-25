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
 * open the input file (bed-formatted) and count how many items within the bed (such as target) files
 * @param bed_file: file in bed-formatted such as (target file or the Ns regions of reference sequences in bed format)
 * @return count	=>	Total number of targets or Ns regions
 */
uint32_t getLineCount(char *bed_file);

/**
 * open the input bed-formated file and then load the coordinate regions into memory and leave them there
 * Note: the input file should be in .bed format
 * @param user_inputs: a variable contains all the user inputs including the target or Ns bed file names
 * @param bed_file: input name for the bed file
 * @param coords: the Bed_Coords structure to store coordinates of each bed section
 * @return the total number of bases covered in the target/N-regions bed file
 */
uint32_t loadBedFiles(User_Input *user_inputs, char *bed_file, Bed_Coords * coords, khash_t(khStrInt)* wanted_chromosome_hash);

/**
 * generate target and buffer lookup table for quick access. The khash.h file is used for hash table setup
 * @param bed_info: the bed information that is stored for the future usage
 * @param stat_info, statistical information for the reads/bases
 * @param header, the bam/sam/cram header pointer that hold the length info of each chromosome
 * @param user_inputs, contains all the user inputs information, including target_buffer_size
 * @param type: either target bed (type 1) or Ns regions in the reference sequences (type 2)
 */
void generateBedBufferStats(Bed_Info * bed_info, Stats_Info *stats_info, Target_Buffer_Status *target_status, User_Input *user_inputs, khash_t(khStrInt)* wanted_chromosome_hash, short target_file_index, short type);

/**
 * process bed-formatted file and populate the coordinates and lookup hash table
 * @param user_inputs: contains all the user inputs including the target or Ns bed file names
 * @param bed_info: the storage of bed coordinates and the size of the bed file
 * @param stats_info: a variable that contains various statistical information
 * @param target_buffer_status: a variable to store target/buffer info
 * @param type: either target bed (type 1) or Ns regions in the reference sequences (type 2)
 */
void processBedFiles(User_Input *user_inputs, Bed_Info *bed_info, Stats_Info *stats_info, Target_Buffer_Status *target_buffer_status, khash_t(khStrInt)* wanted_chromosome_hash, char* bedfile_name, short target_file_index, short type);

uint8_t getTargetBufferBit(uint8_t target_file_index);

void setTargetBufferStatus(Target_Buffer_Status *target_buffer_status, int chrom_idx, uint32_t pos_idx, uint8_t target_file_index, uint8_t type);

void processBufferRegions(uint32_t start, uint32_t end, int chrom_idx, Target_Buffer_Status *target_buffer_status, short target_file_index);

/**
 * just to output some information for debugging
 */
void outputForDebugging(Bed_Info *bed_info);

/**
 * this function is used to free all the memories that are allocated by the program
 * @param bed_info: the declared Bed_Info variable
 */
void cleanBedInfo(Bed_Info *bed_info);

/**
 * It is used to store both target and buffer information
 * @param chromosome_id: the current chromosome to load
 * @param size: the size of the chromosome
 * @param target_buffer_regions: the 2D char array that store the flags either targets('T') or their surrounding buffers ('B')
 * @param target_coords: the struct array that stores the coordinates of targets
 * @return
 */
void getTargetAndBufferPositions(char * chromosome_id, int size, char *target_buffer_regions, Bed_Coords *coords);

void TargetBufferStatusInit(Target_Buffer_Status *target_buffer_status, bam_hdr_t *header);

void TargetBufferStatusInit2(Target_Buffer_Status *target_buffer_status, khash_t(khStrInt)* wanted_chromosome_hash);

void TargetBufferStatusDestroy(Target_Buffer_Status *target_buffer_status);

#endif //TARGETS_H
