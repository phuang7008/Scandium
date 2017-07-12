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
extern int READLENGTH;
extern int MIN_MAP_SCORE;
extern int MIN_BASE_SCORE;
extern int NUM_OF_THREADS;

extern uint32_t TOTAL_READS_PAIRED;
extern uint32_t TOTAL_ALIGNED_BASES;
extern uint32_t TOTAL_READS_ALIGNED;
extern uint32_t TOTAL_READS_PRODUCED;
extern uint32_t TOTAL_PAIRED_READS_WITH_MAPPED_MATES;
extern uint32_t TOTAL_DUPLICATE_READS;
extern uint32_t TOTAL_SUPPLEMENTARY_ALIGNMENTS;
extern bool REMOVE_DUPLICATES;
extern bool REMOVE_SUPPLEMENTARY_ALIGNMENTS;

extern float _PERCENTAGE;	// percentage (fraction) of total bam reads will be used for analysis

extern char CURRENT_CHROMOSOME_ID[12];

/**
 * define a structure that holds the file strings from user inputs
 */
typedef struct {
	char * target_file;
	char * bam_file;
	char * cov_file;
	char * out_file;
} User_Input;

/**
 * define a structure that holds the target coordinates
 */
typedef struct {
    char chr[12];    // some chromosome ID would be quite long, not sure if they will be used though
    uint32_t start;
    uint32_t end;
} Target_Coords;

/**
 * define a structure to hold a chunk of read buff to be processed by each thread
 */
typedef struct {
	bam1_t **chunk_of_reads;
	uint32_t size;
} Read_Buffer;

/**
 * define a strcuture that hold the working chromosome ids and complete chromosome ids
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
// 32 means the key is 32 bit, while the value is a char type
//KHASH_MAP_INIT_INT(32, char)

// 32 means the key is 32 bit, while the value is of unsigned shortt type
KHASH_MAP_INIT_INT(m32, unsigned short)

/**
 * define a structure to hold coverage hash table, which uses position as key and coverage counts as values!
 */
typedef struct {
	khash_t(m32) *cov_hash;
	char *chrom_id;
} Coverage_Hash;

/**
 * define a khash like structure that has string as key and khash_t(m32) as values
 */
//KHASH_SET_INIT_STR(str)
KHASH_MAP_INIT_STR(str, khash_t(m32)*)

#endif //TERMS_H
