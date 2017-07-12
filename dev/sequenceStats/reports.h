/*
 * =====================================================================================
 *
 *       Filename:  reports.h
 *
 *    Description:  The header file for sequence stats analysis
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

//#include <my_global.h>
//#include <mysql.h>
//#include <stdbool.h>
#include "htslib/sam.h"
#include "terms.h"
#include "annotation.h"

/**
 * To write to the coverage fasta files, as well as determines many of the coverage statistics. 
 * It will also process target info if the target file is available
 * @param chrom_id:  Current chromosome id
 * @param Ns_buffer_hash: it contains all the regions with Ns 
 * @param target_info: the target information stored in a struct array
 * @param chrom_tracking: The array which contains the coverage of every base in the genome
 * @param user_inputs: variable that contains all the user input info
 * @param stats_info: variable used to store all the statistical information regarding bases and reads
 * @param cov_fp: the opened file handler for cov.fasta file for writing
 * @param wig_fp: the opened file handler for wig.fasta file in "wig" format file (good for ucsc) which shows you where all off-target regions with high coverage are
 * @param wgs_fp: the opened file handler for the file that contains whole genome coverage
 */
void writeCoverage(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

/**
 * To compile base related statistics
 * @param stats_info: this variable contains all the base and read related statistics (here we will store them)
 * @param cov_val: the coverage count 
 * @param target: 1 means we need to add target related statistics
 * @param wgs: 1 means we need to add whole genome statistics
 * @return boolean
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
 * @param fp: opened file handle
 * @param type: type of output, 1 for the whole genome, 2 for the target only
 */
void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage, uint8_t type);

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
 * @param begin: the start position of the region to check
 * @param length: the length of the regions to be inspected
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param chrom_idx: the chromosome index
 * @maram fh_all_sites: the opened file handle for all capture sites annotation report file
 */
void produceCaptureAllSitesReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char * chrom_id, User_Input *user_inputs, FILE *fh_all_sites, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

void writeAnnotations(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

/**
 * output those regions with lower than or higher than user specified coverage values
 * @param begin: the start position of the region to check
 * @param length: the length of the regions to be inspected
 * @param chrom_tracking: contains the coverage information for the current chromosome
 * @param chrom_idx: the chromosome index
 * @maram fh_low: the opened file handle for lower coverage report file
 * @param fh_high: the opend file handle for higher coverage report file
 * @return the end position of the region with lower or higher base coverage
 */
uint32_t writeLow_HighCoverageReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char *chrom_id, User_Input *user_inputs, FILE *fh_low, FILE *fh_high, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

char * getRegionAnnotation(uint32_t start, uint32_t end, char *chrom_id, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions);

#endif // REPORTS_H
