/*
 * =====================================================================================
 *
 *       Filename:  reports.h
 *
 *    Description:  The header file for the reporting of sequencing statistics
 *
 *        Version:  1.0
 *        Created:  02/24/2017 03:47:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */

#ifndef REPORTS_H
#define REPORTS_H

#include "htslib/sam.h"
#include "terms.h"
#include "annotation.h"

/**
 * To write to the coverage fasta files, as well as determines many of the coverage statistics. 
 * It will also process target info if the target file is available
 * @param chrom_id:  Current chromosome id
 * @param target_info: the target information stored in a struct array
 * @param chrom_tracking: The array which contains the coverage of every base in the genome
 * @param user_inputs: variable that contains all the user input info
 * @param intronic_regions: the official intronic regions fetched from MySQL to be checked against
 * @param exon_regions: the official exonic regions fetched from MySQL to be checked against
 */
void writeCoverage(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL **exon_regions);

/**
 * To compile base related statistics
 * @param stats_info: this variable contains all the base and read related statistics (here we will store them)
 * @param cov_val: the coverage count 
 * @param target: 1 means we need to add target related statistics
 * @param wgs: 1 means we need to add whole genome statistics
 */
void addBaseStats(Stats_Info *stats_info, uint32_t cov_val, uint8_t target, uint8_t wgs);

/**
 * Writes all the statistical information to an output file.
 * @param stats_info, the statistical information to be outputted
 * @param user_inputs: the flag info to indicated if users have specify to output both WGS or Capture(targets) only
 */
void writeReport(Stats_Info *stats_info, User_Input *user_inputs);

/**
 * Write the general information out
 * @param fp: opened file handle for writing
 * @param stats_info: a pointer to the stats storage
 * @param average_coverage: the average coverage of the current sequencing run
 * @param user_inputs: contains all the user_inputs including target_buffer_size
 * @param type: type of output, 1 for the whole genome, 2 for the target only
 */
void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage, User_Input *user_inputs, uint8_t type);

void generateCaptureStats(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, int32_t chrom_idx);

void produceReportsOnThresholds(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL **exon_regions);

/**
 * Write off target wig file for off target statistics
 * This method is destructive to the data structure, no further work can be done after this method has ran.  
 * Works out whether reads are on or off target and how far off target they are
 * The original script uses 500 away from target as off target, here I am using 100 Buffer size as the off target
 * @param chrom_tracking where it contains all the coverage information
 * @param chrom_id: the chrom_id needs to be taken care of
 * @param target_bed_info: contains the coordinates of all the targets
 * @param user_input: the base file name needed to produce the name of the off-target wig file name
 * @param stats_info: the function will increment a non_traget_good_hits member
 */
void produceOffTargetWigFile(Chromosome_Tracking *chrom_tracking, char *chrom_id, Bed_Info *target_bed_info, User_Input *user_inputs, Stats_Info *stats_info);

/**
 * Write all Capture Regions out with detailed annotation including the average coverage
 * @param begin: the start position of the region to be checked
 * @param length: the length of the regions to be inspected
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param chrom_id: the chromosome id to be handed
 * @maram fh_all_sites: the opened file handle for all capture sites annotation report file
 * @param intronic_regions: the official intronic regions fetched from MySQL to be checked against            
 * @param exon_regions: the official exonic regions fetched from MySQL to be checked against 
 */
void produceCaptureAllSitesReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char * chrom_id, FILE **fh_all_sites, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

void writeAnnotations(char *chrom_id, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

/**
 * This is a wrapper function to help generate the coverage range info using writeCoverageRanges()
 * This will not generate annotation information for the speed reason
 * @param chrom_id: the chromosome id to be handed
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param user_inputs: contains all the user_inputs options
 */
void coverageRangeInfoForGraphing(char *chrom_id, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs);

/**
 * This is used to generate average coverage information for a range of position based on different binning strategies
 * @param begin: the start position of the region to check
 * @param length: the length of the regions to be inspected
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param chrom_idx: the chromosome index 
 * @param user_inputs: contains all the user_inputs options
 * @param fh_uniformity: an opened uniformity data report file handle
 */
void writeCoverageRanges(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, uint16_t chrom_idx, User_Input *user_inputs, FILE *fh_uniformity);

/**
 * output those regions with lower than or higher than user specified coverage values
 * @param begin: the start position of the region to check
 * @param length: the length of the regions to be inspected
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param chrom_id: the chromosome id to be handled
 * @param user_inputs: contains all the user_inputs options 
 * @maram fh_low: the opened file handle for lower coverage report file
 * @param fh_high: the opend file handle for higher coverage report file
 * @param intronic_regions: the official intronic regions fetched from MySQL to be checked against            
 * @param exon_regions: the official exonic regions fetched from MySQL to be checked against
 * @param type: mode of processing => 1 for WGS, 2 for Capture
 * @return the end position of the region with lower or higher base coverage
 */
uint32_t writeLow_HighCoverageReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char *chrom_id, User_Input *user_inputs, FILE *fh_low, FILE *fh_high, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions, uint8_t type);

void produceAnnotation(uint8_t type, int32_t chrom_idr, char *chrom_id, uint32_t start, uint32_t end, User_Input *user_inputs, FILE *out_fh, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

bool fetchAnnotation(char *chrom_id, uint32_t start, uint32_t end, FILE *out_fh, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

/*
 * this is the real function that is used to fetch the detailed gene annotations
 * @param start: the start position to be annotationed
 * @param end: the end position to be annotationed
 * @param chrom_id: the chromosome id to be handled
 * @param intronic_regions: the official intronic regions fetched from MySQL to be checked against            
 * @param exon_regions: the official exonic regions fetched from MySQL to be checked against
 * @param type: for faster search or slow search (not working for faster search)
 * @return the detailed annotations
 */
char * getRegionAnnotation(uint32_t start, uint32_t end, char *chrom_id, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions, uint8_t type);

/*
 * find out the detailed gene/transcript coverage percentage
 * @param chrom_id: the chromosome id to be handled
 * @param target_info: detailed targets in bed format
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param user_inputs: contains all the user_inputs options
 * @param low_cov_gene_hash: a khash_t(khStrLCG) variable that contains all the low coverage information
 */
void calculateGenePercentageCoverage(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, khash_t(khStrLCG) *low_cov_gene_hash);

/*
 * write the gene/transcript coverage percentage to a file
 * @param chrom_id: the chromosome id to be handled
 * @param target_info: detailed targets in bed format
 * @param user_inputs: contains all the user_inputs options
 * @param transcript_hash: a khash_t(khStrLCG) variable that contains all transcript information
 * @param gene_transcripts: a khash_t(khStrStrArray) variable that contains all genes and their transcripts
 * @param hgmd_genes: a lookup hash table that contains all HGMD genes
 * @param hgmd_transcripts: a lookup hash table that contains all HGMD transcripts
 */
void storeGenePercentageCoverage(char *chrom_id, User_Input *user_inputs, khash_t(khStrLCG) *transcript_hash, khash_t(khStrStrArray) *gene_transcripts, khash_t(khStrInt) *hgmd_transcripts, khash_t(khStrGTP) **gene_transcript_percentage_hash, uint8_t file_index);

/* it is used to store percentage for current gene's transcripts for later usage
 * @param hgmd_transcripts_hash: a lookup hash table that contains all HGMD transcripts
 * @param transcript_name: the transcript name to be handled
 * @param percentage: the percentage of coverage for current transcript
 */
void storeHGMD_TranscriptPercentageInfo(khash_t(khStrGTP) *gene_transcript_percentage_hash, char *gene_symbol, char* transcript_name, float percentage, bool HGMD_on);

void outputGeneCoverage(khash_t(khStrGTP) *gene_transcript_percentage_hash, User_Input *user_inputs, uint8_t file_index);

#endif // REPORTS_H
