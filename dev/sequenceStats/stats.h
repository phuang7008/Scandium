/*
 * =====================================================================================
 *
 *       Filename:  stats.h
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

#ifndef STATS_H
#define STATS_H

#include <stdbool.h>
#include "htslib/sam.h"
#include "terms.h"

/**
 * This is used to initialize an array of read buffers of type bam1_t (pointers)
 * @param read_buff_in: the array of read_buff_in to be initialized
 * @param size: the size of the read_buff_in array
 */
void readBufferInit(Read_Buffer *read_buff_in);

/**
 * This function is used to delete all the memory allocated for read_buff_in of type bam1_t
 * @param read_buff_in: the array of read_buff_in to be destroyed
 * @param size: the size of the read_buff_in array
 */
void readBufferDestroy(Read_Buffer *read_buff_in);

/**
 * This function is used to read one chunk (10,000 -- 1,000,000) of aligned reads from bam into the memory 
 * @param sffh: samFile file hander: this is actually refer to the opened bam file handler
 * @param header: bam file header information
 * @param more_to_read: a boolean flag that is used to signal is bam file reading is finished (shared among all threads)
 * @param read_buff_in: the array of Read_Buffer to be used to store chunk of alignment reads to be processed by individual thread
 */
uint32_t readBam(samFile *sffh, bam_hdr_t *header, bool more_to_read,  Read_Buffer *read_buff_in);

/**
 * This function is used to process one chunk of aligned reads from bam file and process it using current thread
 * @param user_inputs: variable that contains all the user input info
 * @param stats_info: variable used to store all the statistical information regarding bases and reads
 * @param coverage_hash: the hash table that is used to store temp calculation results
 * @param header: bam file header information
 * @param read_buff_in: the array of Read_Buffer that will be processed by current function
 */
void processBamChunk(User_Input *user_inputs, Stats_Info *stats_info, khash_t(str) *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in);
//void process_chunk_of_bam(int thread_id, Chromosome_Tracking *chrom_tracking, Coverage_Hash *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in);

/**
 * This function is used to process individual aligned read and put the results into a hash table where the key is the chromosome location
 * while the value is the count
 * @param user_inputs: variable that contains all the user input info
 * @param stats_info: variable used to store all the statistical information regarding bases and reads
 * @param coverage_hash: the hash table used to store the coverage count
 * @param chrom_id
 * @param rec: the individual alignment record to be processed
 */
void processRecord(User_Input *user_inputs, Stats_Info *stats_info, khash_t(str) *coverage_hash, char *chrom_id, bam1_t *rec);
//void processRecord(Coverage_Hash *coverage_hash, int hash_index, bam1_t *rec);

/**
 * This function is used to combine each thread coverage results and put them onto a big array with length of each chromosome
 * This function should be critical as no more than one thread should write to the same array at the same time
 * @param chrom_tracking: the variable that is used to tracking the status of all chromosomes
 * @param coverage_hash: the coverage hash table whose contents will write into the chrom_coverage
 */
void combineThreadResults(Chromosome_Tracking *chrom_tracking, khash_t(str) *coverage_hash, bam_hdr_t *header);

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
void writeCoverage(char *chrom_id, khash_t(str) *Ns_buffer_hash, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info);

/**
 * To compile base related statistics
 * @param stats_info: this variable contains all the base and read related statistics (here we will store them)
 * @param cov_val: the coverage count 
 * @param target: 1 means we need to add target related statistics
 * @param wgs: 1 means we need to add whole genome statistics
 * @return boolean
 */
bool addBaseStats(Stats_Info *stats_info, uint16_t cov_val, uint8_t target, uint8_t wgs);

/**
 * Writes all the statistical information to an output file.
 * @param fname output file name
 */
void writeReport(Stats_Info *stats_info, char *output_file);

/**
 * Write the general information out
 * @param fp: opened file handle
 * @param type: type of output, 1 for the whole genome, 2 for the target only
 */
void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, float average_coverage, uint8_t type);

#endif // STATS_H
