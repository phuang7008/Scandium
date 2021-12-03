/*
 *
 * ===================================================================================
 *
 *        Filename:        terms.h
 *
 *        Description:    define constants/struct to be used in Sequencing Statistics
 *
 *        Version:        1.0
 *        Created:        01/30/2017 04:45:04 PM
 *        Revision:        none 
 *        Compiler:        gcc
 *
 *        Author:            Peiming (Peter) Huang (phuang@bcm.edu)
 *        Company:        Baylor College of Medicine
 *
 *        =====================================================================================
 */

#ifndef TERMS_H
#define TERMS_H

#include <inttypes.h>   // for PRIu32 and PRIu64 
#include <stdbool.h>    // for bool definition
#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>     // for getopt() and usleep()
#include <math.h>
#include <zlib.h>
#include <mysql.h>
#include <sys/stat.h>
#include <stdio.h>      // for file read and write

// Users defined header files
#include "htslib/sam.h"

#include "data_structure.h"
#include "coverage_tracking.h"

// The followings are defined as macro/constants. The program should never try to change their values
#define VERSION_ "##Scandium v3.0.2"

#define PRIMER_SIZE         1000    // upstream or downstream of a target

#define TARGET_LINE_SIZE    150     // the number of chars to be read for each line within a target file

/* @abstract the bit flag for the 1st capture input files's target/buffer*/
#define TRT_BFR_1    1
/* @abstract the bit flag for the 2nd capture input file's target/buffer*/
#define TRT_BFR_2    2
/* @abstract the bit flag for the 3rd capture input file's target/buffer*/
#define TRT_BFR_3    4
/* @abstract the bit flag for the 4th capture input file's target/buffer*/
#define TRT_BFR_4    8
/* @abstract the bit flag for the 5th capture input file's target/buffer*/
#define TRT_BFR_5    16
/* @abstract the bit flag for the 6th capture input file's target/buffer*/
#define TRT_BFR_6    32
/* @abstract the bit flag for the 7th capture input file's target/buffer*/
#define TRT_BFR_7    64
/* @abstract the bit flag for the 8th capture input file's target/buffer*/
#define TRT_BFR_8    128

// We need to declared the followings as glabal since the program will change these values!!!
// The naming convention for this type of data is CAPTICAL_WORD1_WORD2_WORD3...
extern bool N_FILE_PROVIDED;            // this is the file that contains regions of all Ns in the reference
extern bool TARGET_FILE_PROVIDED;
extern bool HGMD_PROVIDED;
extern bool USER_DEFINED_DATABASE;
extern int  khStrStr;
//extern int  khStrInt;
extern int  khStrFloat;
extern int  khStrStrArray;
extern int  khStrLCG;
extern int  khStrGTP;

/**
 * define a structure that holds the file strings from user inputs
 */
typedef struct {
    // General file info
    char * bam_file;
    char * n_file;                  // provide the regions with Ns in the reference genome in bed format
    char * output_dir;              // output directory (mandatory)
    char * reference_file;          // reference file name for cram input file
    char * chromosome_bed_file;     // a file contains chromosome ids and regions need to be processed in bed format
    char ** target_files;
    char ** user_defined_annotation_files;  // users can provide multiple annotation files
    char ** annotation_file_basenames;      // the basename w/o extension of user defined annotation files
    char ** target_file_basenames;          // the basename w/o extension of target bed files
    char * target_annotation_list_file;     // file contains multiple lines, 
                                            // each line has a target file and its corresponding annotation file separated by a tab
    char * database_version;            // either hg19 (hg37) or hg38
    char * user_name;                   // user name for the database
    char * passwd;                      // password for the user

    // For whole genome (WGS) related outputs
    char * wgs_cov_file;            // output the whole genome coverage count information
    char * wgs_cov_report;          // output the whole genome coverage summary report    
    char * wgs_low_cov_file;        // for the regions within the whole genome that have lower coverage with detailed annotation
    char * wgs_high_cov_file;       // for the regions within the whole genome that have high coverage without detailed annotation
    char * wgs_uniformity_file;     // for whole genome uniformity data file

    // For Capture related outputs
    //char * missed_targets_file;       // for target regions that have no coverage at all
    char ** capture_wig_files;          // output the off target good hit in wig formatted
    char ** capture_cov_files;          // output the capture target regions coverage count information
    char ** capture_cov_bedfiles;       // output the capture target regions coverage count information
    char ** capture_cov_reports;        // output the capture target regions coverage summary report
    char ** capture_low_cov_files;      // for target regions with lower coverage and their detailed annotation
    char ** capture_high_cov_files;     // for target regions with high overage without detailed annotation
    char ** capture_all_site_files;     // for all target regions with average coverage and the detailed annotation
    char ** low_cov_gene_pct_files;     // for percentage of a gene with low coverage bases
    char ** low_cov_exon_pct_files;     // for percentage of an exon with low coverage bases
    char ** low_cov_transcript_files;   // for low cov percentage of every gene in capture file with all different transcripts

    //misc
    int8_t min_map_quality;
    int8_t min_base_quality;
    uint16_t num_of_target_files;
    uint16_t num_of_annotation_files;
    uint16_t low_coverage_to_report;    // default 20, users are allowed to change it as an input option
    uint32_t high_coverage_to_report;   // default 10000, to report regions with higher coverage as users specified
    uint16_t lower_bound;               // default 1. Used with -u option (upper_bound) for the uniformity output
    uint16_t upper_bound;               // default 150. Used with -l option (lower_bound) for the uniformity output
    uint16_t gVCF_percentage;           // default 5 for 500%. For gVCF formula: BLOCKAVE_Xp, where X=gVCF_percentage
    uint16_t target_buffer_size;        // default 100. For regions immediate adjacent to any target regions
    unsigned short num_of_threads;
    float percentage;                   // percentage (fraction) of total bam reads will be used for analysis
    uint8_t size_of_peak_area;          // the points around peak area to pick for uniformity calculation
    bool user_set_peak_size_on;
    bool remove_duplicate;
    bool remove_supplementary_alignments;
    bool wgs_annotation_on;
    bool above_10000_on;                // output bases/regions with coverage > 10000
    bool excluding_overlapping_bases;   // to exclude overlapping reads/bases for double counting
    bool wgs_coverage;
    bool Write_Capture_cov_fasta;
    bool Write_WGS_cov_fasta;           // need two different flags, one for Capture and one for WGS
    bool Write_WIG;
    //bool primary_chromosomes_only;    // do we need all chromosomes including decoy, alt etc or primary only

    // developer testing options
    bool non_MC_tag_ON;                 // use non_MC_tag approach for overlap base removal even though the bam file has MC_tags
} User_Input;


/**
 * define a structure to hold the database name infomation
 */
typedef struct {
    char* db_annotation;    // Database contains annotation column (annotations are partitioned)
    char* db_coords;        // it contains cds_target_start and cds_start columns, for gene/transcript/exon percentage calculation
    char* db_introns;       // the intronic regions
    char* db_hgmd;          // For official HGMD gene/transcripts
    MYSQL *con;
    MYSQL_RES *mysql_results;
} Databases;

typedef struct {
    char *chrom_id;
    uint32_t num_of_cds;
} User_Defined_Database;

typedef struct {
    uint32_t num_of_chroms;             // the number of chromosomes in the user specified annotation file
    uint64_t num_of_lines;
    User_Defined_Database *ud_database_per_chrom;
} User_Defined_Database_Wrapper;

/*
 * store the raw user-defined-annotation with the structure like the following
 * num_of_chromosomes                25 + alt + hla + decoy etc.
 * chrom_id     (array)                "1"        "2"         "3"    ...    "7"        "Y"        "alt" ...
 * annotation_size (array)             0         5            17            12         22         35  ...
 * annotations (array of arrays)                              "7    117292897    117292985    CFTR|ENST00000600166|cds_1|gene
 *                                                            "7    117304742    117304914    CFTR|ENST00000600166|cds_2|gene
 *                                                            "7    117305513    117305618    CFTR|ENST00000600166|cds_3|gene
 *                                                            "7    117355812    117355913    CFTR|ENST00000600166|cds_4|gene
 *                                                            .....
 */
typedef struct {
    uint32_t *annotation_size;
    uint32_t num_of_chromosomes;        // the number of chromosomes in the user specified annotation file
    char **chrom_id;
    char *** annotations;
    int32_t number_of_unique_genes;
    int32_t number_of_unique_transcripts;
} Raw_User_Defined_Database;

/**
 * define a structure that holds the exonic, inter-genic and intronic regions
 * single '*' for 1-D pointer, double '**' for 2-D array, while three '***' for 2-D array of strings
 *
 * chromosome_ids array list:                    '1', '2', '3', '4', '5', ...... 'X', 'Y', 'MT'
 * size_r: num_of_regions on each chrom:         55,  96,  183, 66,  9,   ...... 32,  15,  2
 * region starts/ends on one chrom:               0:  31938          32954
 *                                                1:  65938          78932
 *                                                2:  101343         123908
 *                                                .
 *                                                54: 1498903        2809823
 * transcript_name for a region per chrom:        0: 'ABL2'
 *                                                1: 'BCKL1'         'DSSL3'
 *                                                .
 *
 */
typedef struct {
    uint16_t chrom_list_size;   // the size of first dimension 
    char **chromosome_ids;        
    uint32_t *size_r;           // this is the total number of regions on each chromosome
    uint32_t **starts;          // each region has its start and end position
    uint32_t **ends;
    uint32_t prev_search_loc_index;     // this is to avoid search starts from beginning
    uint32_t prev_search_chrom_index;   // need to make sure we skip search for the previous same chromosome

    char ***gene;
    char ***Synonymous;
    char ***prev_genes;
    char ***exon_info;
} Regions_Skip_MySQL;

/** define annotation structure that is used to store multiple hits from search results
 *  This is used to handle the case where there might be multiple genes within a low coverage region
 *  And we need to include all of them in our annotation
 */
typedef struct {
    char *gene;
    char *Synonymous;
    //char *prev_genes;        // not used!
    char *exon_info;
} Annotation;

/* define an annotation array structure to hold the annotation array
 * n_regions, iter)
 */
typedef struct {
    Annotation *annotations;
    uint16_t allocated_size;    // the allocated size through memory management
    uint16_t real_size;         // the number of exon info currently stored in this array
} Annotation_Wrapper;

/**
 * define a structure for Gene RefSeq CDS Exons for the calculation of gene/transcript/cds Percentage Coverage Reports
 * For exon_target_start and exon_target_end, they are different things for different variables
 * For the refseq_cds_genes, they just refer to the cds exon coordinates
 * For the low_cov_genes, they are the intersect between the official cds exon coordinates from 
 *        the RefSeq DB and the target regions from user inputs
 */
typedef struct {
    char *gene_symbol;
    char *transcript_name;
    bool targeted;
    uint32_t cds_start;             // for the RefSeq transcript coding region (whole) only
    uint32_t cds_end;               // for the RefSeq transcript coding region (whole) only
    uint32_t cds_target_start;      // for a single CDS exon only    
    uint32_t cds_target_end;        // for a single CDS exon only
    int16_t  exon_id;               // for SNP, it is -1. So it should be signed int
    uint16_t exon_count;
    uint32_t cds_length;
    uint16_t num_of_low_cov_bases;
    StringArray * low_cov_regions;    // it will be an array of string structure
} Gene_Coverage;

typedef struct {
    Gene_Coverage *gene_coverage;
    uint32_t capacity;
    uint32_t total_size;
} Low_Coverage_Genes;

typedef struct {
    char* transcript_name;
    float percentage;
    bool  HGMD;
} Transcript_Percentage;

typedef struct {
    Transcript_Percentage *transcript_percentage;
    uint16_t capacity;
    uint16_t size;
    bool has_HGMD;
} Gene_Transcript_Percentage;

#include "htslib/khash.h"

/**
 * define a khash like structure that has string as key and various structures as values
 * Note: name part of init must be unique for the key, value types.
 * In our case, 33/32/31 (defined in main.c) are arbitrary symbolic names for hashtables
 * that contains string keys and Low_Coverage_Genes* and int values.
 */
KHASH_MAP_INIT_STR(khStrStr, char*)

KHASH_MAP_INIT_STR(khStrStrArray, StringArray*)

//KHASH_MAP_INIT_STR(khStrInt, uint32_t)

KHASH_MAP_INIT_STR(str, Temp_Coverage_Array*)

KHASH_MAP_INIT_STR(khStrLCG, Low_Coverage_Genes*)

KHASH_MAP_INIT_STR(khStrGTP, Gene_Transcript_Percentage*)

//KHASH_MAP_INIT_STR(khStrRSM, Regions_Skip_MySQL*)

/**
 * define a coverage statistics structure for general read stats
 */
typedef struct {
    //read stats
    uint64_t total_reads_produced;          // total reads contains in the bam
    uint64_t total_reads_aligned;           // total reads aligned to a genomic region
    uint32_t total_duplicate_reads;         // total number of duplicate reads
    uint64_t total_reads_paired;            // total number of reads with mate pairs (if any)
    uint64_t total_reads_proper_paired;     // total number of reads with mate pairs (if any)
    uint32_t total_chimeric_reads;          // total number of duplicate reads
    uint32_t total_supplementary_reads;     // total number of reads with supplementary flag set
    uint64_t total_paired_reads_with_mapped_mates; // total number of aligned reads which have mapped mates

    uint16_t read_length;                   // the sequenced READ Length, it is taken from => read_buff_in->chunk_of_reads[i]->core.l_qseq
} Read_Coverage_Stats;

typedef struct {
    //base stats
    uint64_t total_genome_bases;            // total number of bases in Genome
    uint32_t total_Ns_bases;                // total number of bases that are N (unknown bases)
    uint32_t total_Ns_bases_on_chrX;        // total number of bases that are N (unknown bases) on X chromosome
    uint32_t total_Ns_bases_on_chrY;        // total number of bases that are N (unknown bases) on Y chromosome
    uint64_t total_mapped_bases;            // total number of mapped bases
    uint64_t total_uniquely_aligned_bases;  // aka. Reads Usable - where "Usable" is uniquely aligned, non-duplicate, on-target reads
    uint64_t total_genome_coverage;         // total number of read bases aligned to the Genome (used to calculate average coverage)
    uint64_t base_quality_20;               // total number of aligned bases with quality >= 20
    uint64_t base_quality_30;               // total number of aligned bases with quality >=30
    uint32_t total_overlapped_bases;        // total number of overlapped bases from pair-end reads

    //misc
    uint32_t wgs_max_coverage;
    uint32_t base_with_wgs_max_coverage;
    uint16_t median_genome_coverage;
    uint32_t mode;
    double uniformity_metric_all;                   // it contains all the alt, decoy etc.
    double uniformity_metric_all_primary;           // it doesn't contain alt, decoy chromosomes.
    double uniformity_metric_autosome_only;         // it contains all the alt, decoy etc.
    double uniformity_metric_primary_autosome_only; // it doesn't contain alt, decoy chromosomes.

    uint32_t genome_cov_histogram[1001];            // coverage histogram array
    khash_t(m32) *genome_base_with_N_coverage;      // here N stands for 1, 5, 10, 15, 20, 30, 40, 50, 60, 100
    khash_t(m32) *genome_coverage_for_median;       // Used for calculating the median coverage for Whole Genome.
} WGS_Coverage_Stats;

/**
 * define a coverage statistics structure for Capture
 */
typedef struct {
    //read stats on target/buffer
    uint64_t total_target_coverage;         // total number of read bases aligned to the target (used to calculate average coverage)
    uint32_t total_targeted_bases;          // total number of bases targeted
    uint32_t total_buffer_bases;            // total number of bases in the buffer region
    uint32_t on_target_read_hit_count;      // total number of reads which align to a target region
    uint32_t off_target_read_hit_count;     // total number of reads which do not align to a target region
    uint32_t in_buffer_read_hit_count;      // total number of reads which align to the buffer region
    uint32_t hit_target_count;              // total targets with at least 1 read aligned to them
    uint32_t hit_target_buffer_only_count;  // total targets with no hits, except in buffer
    uint32_t non_target_good_hits;          // regions that have high coverage but are not in the target

    //misc
    uint32_t total_targets;                 // total taregted regions
    uint32_t target_max_coverage;
    uint32_t base_with_target_max_coverage;
    uint16_t median_target_coverage;
    uint32_t target_coverage[101];          // stores data about the coverage across a target
    uint32_t five_prime[PRIMER_SIZE];       // stores data about coverage upstream of the target 
    uint32_t three_prime[PRIMER_SIZE];      // stores data about coverage downstream of the target

    uint32_t target_cov_histogram[1001];        // coverage histogram array
    khash_t(m32) *target_base_with_N_coverage;  // here N stands for 1, 5, 10, 15, 20, 30, 40, 50, 60, 100
    khash_t(m32) *target_coverage_for_median;   // Used for calculating the median coverage for Whole Genome.

} Capture_Coverage_Stats;

/**
 * define a structure to store various information, such as coverage histogram etc
 */
typedef struct {
    Read_Coverage_Stats     *read_cov_stats;
    WGS_Coverage_Stats      *wgs_cov_stats;
    Capture_Coverage_Stats  **capture_cov_stats;
} Stats_Info;

#endif //TERMS_H
