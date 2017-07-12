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
void read_buff_init(Read_Buffer *read_buff_in);

/**
 * This function is used to delete all the memory allocated for read_buff_in of type bam1_t
 * @param read_buff_in: the array of read_buff_in to be destroyed
 * @param size: the size of the read_buff_in array
 */
void read_buff_destroy(Read_Buffer *read_buff_in);

/**
 * This function is used to read one chunk (10,000 -- 1,000,000) of aligned reads from bam into the memory 
 * @param sffh: samFile file hander: this is actually refer to the opened bam file handler
 * @param header: bam file header information
 * @param more_to_read: a boolean flag that is used to signal is bam file reading is finished (shared among all threads)
 * @param read_buff_in: the array of Read_Buffer to be used to store chunk of alignment reads to be processed by individual thread
 */
uint32_t read_bam(samFile *sffh, bam_hdr_t *header, bool more_to_read,  Read_Buffer *read_buff_in);

/**
 * This function is used to process one chunk of aligned reads from bam file and process it using current thread
 * @param header: bam file header information
 * @param read_buff_in: the array of Read_Buffer that will be processed by current function
 * @param coverage_hash: the hash table that is used to store temp calculation results
 */
void process_chunk_of_bam(int thread_id, Chromosome_Tracking *chrom_tracking, khash_t(str) *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in);
//void process_chunk_of_bam(int thread_id, Chromosome_Tracking *chrom_tracking, Coverage_Hash *coverage_hash, bam_hdr_t *header, Read_Buffer *read_buff_in);

/**
 * This function is used to process individual aligned read and put the results into a hash table where the key is the chromosome location
 * while the value is the count
 * @param coverage_hash: the hash table used to store the coverage count
 * @rec: the individual alignment record to be processed
 */
void processRecord(khash_t(str) *coverage_hash, char *chrom_id, bam1_t *rec);
//void processRecord(Coverage_Hash *coverage_hash, int hash_index, bam1_t *rec);

/**
 * This function is used to combine each thread coverage results and put them onto a big array with length of each chromosome
 * This function should be critical as no more than one thread should write to the same array at the same time
 * @param one_chromosome_array: an array for one chromosome (the array type is short int)
 * @param coverage_hash: the coverage hash table whose contents will write into the one_chromosome_array
 */
void combine_thread_results(short *one_chromosome_array, khash_t(32) * coverage_hash);

#endif // STATS_H
