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
#include <unistd.h>     // for getopt() and usleep()
#include <math.h>
#include <zlib.h>
#include <my_global.h>
#include <mysql.h>

// Users defined header files
#include "htslib/sam.h"

// The followings are defined as macro/constants. The program should never try to change their values
#define VERSION_ "SeqStats v1.0 2017-02-05"

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
	char * bam_file;
	char * capture_cov_file;
	char * n_file;					// provide the regions with Ns in the reference genome in bed format
	char * out_file;
	char * target_file;
	char * wgs_wig_file;				// output the off target good hit in wig formatted
	char * wgs_cov_file;				// output the whole genome coverage information
	char * missed_targets_file;		// for target regions that have no coverage at all
	char * capture_low_cov_file;	// for target regions with lower coverage and their detailed annotation
	char * capture_high_cov_file;	// for target regions with high overage without detailed annotation
	char * capture_all_site_file;	// for all target regions with average coverage and the detailed annotation
	char * low_cov_gene_pct_file;	// for percentage of a gene with low coverage bases
	char * low_cov_exon_pct_file;	// for percentage of an exon with low coverage bases
	char * low_cov_transcript_file;	// for low cov percentage of every gene in capture file with all different transcripts
	char * wgs_low_cov_file;		// for the regions within the whole genome that have lower coverage with detailed annotation
	char * wgs_high_cov_file;		// for the regions within the whole genome that have high coverage without detailed annotation

	//misc
	int8_t min_map_quality;
	int8_t min_base_quality;
	uint8_t low_coverage_to_report;		// default 20, users are allowed to change it as an input option
	uint16_t high_coverage_to_report;	// default 10000, to report regions with higher coverage as users specified
	int16_t upper_bound_to_report;		// default -1. It is used to set the coverage upper bound to report. it is used with high_coverage_to_report for a range to report
	short num_of_threads;
	float percentage;					// percentage (fraction) of total bam reads will be used for analysis
	bool remove_duplicate;
	bool remove_supplementary_alignments;
	bool annotation_on;
	bool wgs_coverage;
	bool Write_WIG;
	bool Write_WGS;
} User_Input;

/**
 * define a structure that holds the target coordinates
 */
typedef struct {
    char chrom_id[15];    // some chromosome ID would be quite long
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
 * define a strcuture for quick lookup of target information
 */
typedef struct {
	char chrom_id[50];				// which chromosome it is tracking
	uint32_t size;					// size of current chromosome
	int32_t index;					// -1 means this chromosome is not important/processed, and we should skip it!
	uint8_t *status_array;			// the status for each chrom position, 1 for target, 2 for buffer and 3 for Ns
} Target_Buffer_Status;

/**
 * define a structure to hold a chunk of read buff to be processed by each thread
 */
typedef struct {
	bam1_t **chunk_of_reads;
	uint32_t size;
} Read_Buffer;

/**
 * data structure to store temp coverage results
 */
typedef struct {
	uint32_t size;
	uint32_t * cov_array;
} Temp_Coverage_Array;

/**
 * define a strcuture that hold the chromosome ids status in array 
 * (that is all of the member variables are array, except number_tracked)
 */
typedef struct {
	uint32_t number_tracked;		// the number of chromosomes we are tracking so far!
	uint32_t **coverage;			// the coverage count info for each base on each chromosome will be stored here!

    char **chromosome_ids;
    uint32_t *chromosome_lengths;
    uint8_t  *chromosome_status;	// 0 pending, 1 working, 2 finish coverage calculation, 3. done summary-report/annotation processing, 4.done writing!
	bool more_to_read;
} Chromosome_Tracking;

/**
 * define a structure that holds the exonic, inter-genic and intronic regions
 * single '*' for 1-D pointer, double '**' for 2-D array, while three '***' for 2-D array of strings
 *
 * chromosome_ids array list:					'1', '2', '3', '4', '5', ...... 'X', 'Y', 'MT'
 * size_r: num_of_regions on each chrom:        55,  96,  183, 66,  9,   ...... 32,  15,  2
 * region starts/endds on one chrom:			0:  31938		32954
 *												1:  65938		78932
 *												2:  101343		123908
 *												.
 *												54: 1498903		2809823
 * gene_name for each region on one chrom:		0: 'ABL2'
 *												1: 'BCKL1'	'DSSL3'
 *												.
 *
 */
typedef struct {
	uint16_t chrom_list_size;	// the size of first dimention 
	char **chromosome_ids;		
	uint32_t *size_r;			// this is the total number of regions on each chromosome
	uint32_t **starts;			// each region has its start and end position
	uint32_t **ends;
	uint32_t prev_search_loc_index;		// this is to avoid search starts from beginning
	uint32_t prev_search_chrom_index;	// need to make sure we skip search for the previous same chromosome

	char ***gene;
	char ***Synonymous;
	char ***prev_genes;
	char ***exon_info;
} Regions_Skip_MySQL;

/**
 * define a structure for Gene RefSeq Exons for the calculation of gene/exon Percentage Coverage Reports
 */
typedef struct {
	char *gene_symbol;
	char *refseq_name;
	//char *chrom_id;
	uint32_t refseq_start;		// for the RefSeq transcript
	uint32_t refseq_end;
    uint32_t exon_target_start;	// the intersect between the official exon coordinates from RefSeq DB and the target regions from PKv2
    uint32_t exon_target_end;
	uint16_t exon_id;
	uint16_t exon_count;
    uint16_t num_of_low_cov_bases;
	char *low_cov_regions;
} Gene_Coverage;

typedef struct {
	Gene_Coverage *gene_coverage;
	uint32_t num_of_refseq;
} Low_Coverage_Genes;

/**
 * define a data structure to store the information related to transcript percentage coverage
 */
typedef struct {
	char *gene_symbol;
	char *refseq_name;
	float refseq_cov_percentage;
} Transcript_Coverage_Percentage;

typedef struct {
	Transcript_Coverage_Percentage *transcript_cov_pct;
	uint32_t num_of_refseq;
} Transcript_Coverage;

#include "htslib/khash.h"

// Instantiate a hash map containing integer keys
// m32 means the key is 32 bit, while the value is a char type
//KHASH_MAP_INIT_INT(32, char)

// m32 means the key is 32 bit integer, while the value is of unsigned short type (ie uint16_t)
KHASH_MAP_INIT_INT(m32, uint32_t)
//KHASH_MAP_INIT_INT(m32)
KHASH_MAP_INIT_INT(m16, uint16_t)
KHASH_MAP_INIT_INT(m8, uint16_t)

/**
 * define a khash like structure that has string as key and khash_t(m32) as values
 */
KHASH_MAP_INIT_STR(str, Temp_Coverage_Array*)
//KHASH_MAP_INIT_STR(str, khash_t(m32)*)
//KHASH_SET_INIT_STR(str)

/**
 * define a coverage statistics structure
 */
typedef struct {
	//base stats
	uint64_t total_genome_bases;		//total number of bases in Genome
	uint32_t total_buffer_bases;		//total number of bases in the buffer region
	uint32_t total_targeted_bases;		//total number of bases targeted
	uint32_t total_Ns_bases;			//total number of bases that are N (unknown)
	uint64_t total_aligned_bases;		//total number of aligned bases
	uint64_t total_target_coverage;		//total number of read bases aligned to the target
	uint64_t total_genome_coverage;		//total number of read bases aligned to the Genome

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
	uint32_t hit_target_buffer_only_count;	//total targets with no hits, except in buffer
	uint32_t non_target_good_hits;			//regions that have high coverage but are not in the target

	//misc
	uint32_t total_targets;					//total taregted regions
	uint16_t read_length;					//the sequenced READ Length
	uint32_t max_coverage;
	uint32_t base_with_max_coverage;
	uint16_t median_genome_coverage;
	uint16_t median_target_coverage;
} Coverage_Stats;

/**
 * define a structure to store various information, such as coverage histogram etc
 */
typedef struct {
    khash_t(m32) *target_cov_histogram;             //target coverage histogram
    khash_t(m32) *genome_cov_histogram;             //coverage histogram for the whole Genome

    khash_t(m32) *targeted_base_with_N_coverage;    // here N stands for 1, 5, 10, 15, 20, 30, 40, 50, 60, 100
    khash_t(m32) *genome_base_with_N_coverage;      // here N stands for 1, 5, 10, 15, 20, 30, 40, 50, 60, 100

    khash_t(m32) *target_coverage_for_median;       //Used for calculating the median coverage.
    khash_t(m32) *genome_coverage_for_median;       //Used for calculating the median coverage for Whole Genome.

    uint32_t target_coverage[101];                  //stores data about the coverage across a target
    uint32_t five_prime[PRIMER_SIZE];              //stores data about coverage upstream of the target 
    uint32_t three_prime[PRIMER_SIZE];             //stores data about coverage downstream of the target 

    Coverage_Stats *cov_stats;
} Stats_Info;
#endif //TERMS_H
