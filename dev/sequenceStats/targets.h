/*
 * =====================================================================================
 *
 *      Filename:		targets.h
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

#ifndef TARGETS_H
#define TARGETS_H

#include <stdbool.h>
#include "htslib/sam.h"
#include "terms.h"

/**
 * open the input target file and count how many items within the target file, finally it will return the count
 * Note: the target file should be in bed format
 * @param targetFile
 * @return count	=>	Total number of targets
 */
long getTargetCount(char *target_file);

/**
 * open the input target file and then load the target regions into memory and leave them there
 * Note: the target file should be in .bed format
 * @param targetFile
 * @param targetStarts	=> integer array to store all the target starts
 * @param targetStopa	=> integer array to store all the target stops
 * @param targetChrs	=> char array to store all the target chromosome ids (chars as well). so it is an array of string
 */
void loadTargets(char * target_file, Target_Coords * targets);

/**
 * generate target and buffer lookup table for quick access. The khash.h file is used for hash table setup
 * @param target_buffer_lookup_table, which is a type of khash_t(32)
 */
void generateTargetBufferLookupTable(Target_Coords * target_coords, long target_line_count, khash_t(32) *target_buffer_hash[]);

/**
 * Gets the target regions from the target file, and writes over the coverage fasta files, as well as determines many of the coverage statistics.
 * @param chromo  Current chromosome
 * @param COVERAGE The array which contains the coverage of every base in the genome
 * @param covFasta The filewriter for the target-specific coverage
 * @param missTraget  Write a "wig" format file (good for ucsc) which shows you where all off-target regions with high coverage are
 * @param wgCoverage A filewriter for the whole genome coverage... if null, this file won't be written
 */
void getTargetsAndWriteCoverage(char * chromo,  short * COVERAGE, char * covFasta, char * missTraget, char * wgCoverage);

/**
 * 
 * @param chromosome_id: the current chromosome to load
 * @param size: the size of the chromosome
 * @param target_buffer_regions: the 2D char array that store the flags either targets('T') or their surrounding buffers ('B')
 * @param target_coords: the struct array that stores the coordinates of targets
 * @return
 */
void getTargetAndBufferPositions(char * chromosome_id, int size, char *target_buffer_regions, Target_Coords *target_coords);


#endif //TARGETS_H
