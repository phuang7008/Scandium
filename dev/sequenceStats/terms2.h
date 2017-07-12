/*
 *
 * ===================================================================================
 *
 *		Filename:		terms.h
 *
 *		Description:	define constants to be used in Capture Statistics
 *
 *		Version:		1.0
 *		Created:		01/30/2017 04:45:04 PM
 *		Revision:		none 
 *		Compiler:		gcc
 *
 *		Author:			Peiming (Peter) Huang
 *		Company:		Baylor College of Medicine
 *
 *		=====================================================================================
 */

#ifndef TERMS_H
#define TERMS_H

#include <inttypes.h>
#include <stdbool.h>	// for bool definition
#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>     // for getopt()
#include <math.h>
#include <zlib.h>

// Users defined header files
#include "htslib/sam.h"

// The followings are defined as macro/constants. The program should never try to change their values
#define VERSION "SeqStats v1.0 2017-02-05"

#define PRIMER_SIZE			1000	//upstream or downstream of a target
#define BUFFER				100		//buffer region around a target

#define TARGET_LINE_SIZE	150		//the number of chars to be read for each line within a target file

#define DYMMY				"dummy"	//
#define WRITE_EXOME			0		//write whole genome coverage statistics

// We need to declared the followings as glabal since the program will change these values!!!
// The naming convention for this type of data is CAPTICAL_WORD1_WORD2_WORD3...
extern bool N_FILE_PROVIDED;		// this is the file that contains regions of all Ns in the reference
extern bool TARGET_FILE_PROVIDED;

/**
 * define a structure that holds the file strings from user inputs
 */
typedef struct {
	//file info
	char * target_file;
	char * bam_file;
	char * cov_file;
	char * n_file;		// provide the regions with Ns in the reference genome in bed format
	char * out_file;

	//misc
	uint8_t min_map_quality;
	uint8_t min_base_quality;
	short num_of_threads;
	float percentage;				// percentage (fraction) of total bam reads will be used for analysis
	bool remove_duplicate;
	bool remove_supplementary_alignments;
	bool wgs_coverage;
} User_Input;

/**
 * define a structure that holds the target coordinates
 */
typedef struct {
    char chr[15];    // some chromosome ID would be quite long, not sure if they will be used though as we are only interested in human ones
    uint32_t start;
    uint32_t end;
} Bed_Coords;

/**
 * define a structure that holds the size and coordinates info of a bed file
 */
typedef struct {
	Bed_Coords *coords;
	uint32_t size;
} Bed_Info;

/**
 * define a structure to hold a chunk of read buff to be processed by each thread
 */
typedef struct {
	bam1_t **chunk_of_reads;
	uint32_t size;
} Read_Buffer;

/**
 * define a strcuture that hold the chromosome ids status in array 
 * (that is all of the member variables are array, except number_tracked)
 */
typedef struct {
	size_t number_tracked;			// the number of chromosomes we are tracking so far!
	unsigned short **coverage;		// the coverage count info for each base on each chromosome will be stored here!

    char **chromosome_ids;
    uint32_t *chromosome_lengths;
    uint8_t  *chromosome_status;	// 0 pending, 1 working, 2 finish processing, 3.done writing!
} Chromosome_Tracking;

#include "htslib/khash.h"

// Instantiate a hash map containing integer keys
// m32 means the key is 32 bit, while the value is a char type
//KHASH_MAP_INIT_INT(32, char)

// m32 means the key is 32 bit integer, while the value is of unsigned short type
KHASH_MAP_INIT_INT(m32, unsigned short)
//KHASH_MAP_INIT_INT(m32)
KHASH_MAP_INIT_INT(m16, unsigned short)
KHASH_MAP_INIT_INT(m8, unsigned short)

/**
 * define a khash like structure that has string as key and khash_t(m32) as values
 */
KHASH_MAP_INIT_STR(str, khash_t(m32)*)
//KHASH_SET_INIT_STR(str)

/**
 * define a structure to store various information, such as coverage histogram etc
 */
typedef struct {
    khash_t(m16) *target_cov_histogram;				//target coverage histogram
    khash_t(m16) *genome_cov_histogram;				//coverage histogram for the whole Genome

    khash_t(m16) *targeted_base_with_N_coverage;	// here N stands for 1, 5, 10, 15, 20, 30, 40, 50, 60, 100
    khash_t(m16) *genome_base_with_N_coverage;		// here N stands for 1, 5, 10, 15, 20, 30, 40, 50, 60, 100

	khash_t(m16) *target_coverage_for_median;		//Used for calculating the median coverage.
	khash_t(m16) *genome_coverage_for_median;		//Used for calculating the median coverage for Whole Genome.

	khash_t(m32) *target_coverage;					//stores data about the coverage across a target
	khash_t(m32) *five_prime;						//stores data about coverage upstream of the target 
	khash_t(m32) *three_prime;						//stores data about coverage downstream of the target 

	//base stats
	uint32_t total_genome_bases;		//total number of bases in Genome
	uint32_t total_buffer_bases;		//total number of bases in the buffer region
	uint32_t total_targeted_bases;		//total number of bases targeted
	uint32_t total_Ns_bases;			//total number of bases that are N (unknown)
	uint32_t total_aligned_bases;		//total number of aligned bases
	uint32_t total_target_coverage;		//total number of read bases aligned to the target
	uint32_t total_genome_coverage;		//total number of read bases aligned to the Genome

	//read stats
	uint32_t total_reads_paired;			//total number of reads with mate pairs (if any)
	uint32_t total_reads_aligned;			//total reads aligned to a target region
	uint32_t total_reads_produced;			//total reads contains in the bam
	uint32_t total_duplicate_reads;			//total number of duplicate reads
	uint32_t total_supplementary_reads;		//total number of reads with supplementary flag set
	uint32_t total_paired_reads_with_mapped_mates; //total number of aligned reads which have mapped mates

	//read stats on target/buffer
    uint32_t on_target_read_hit_count;		//total number of reads which align to a target region
    uint32_t off_target_read_hit_count;		//total number of reads which do not align to a target region
    uint32_t in_buffer_read_hit_count;		//total number of reads which align to the buffer region
	uint32_t hit_target_count;				//total targets with at least 1 read aligned to them
	uint32_t hit_buffer_only_count;			//total targets with no hits, except in buffer

	//misc
	uint32_t total_targets;					//total taregted regions
	uint16_t read_length;					//the sequenced READ Length
	uint16_t max_coverage;
	uint16_t base_with_max_coverage;
} Stats_Info;

#endif //TERMS_H
