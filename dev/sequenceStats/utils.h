/*
 * =====================================================================================
 *
 *      Filename:       utilss.h
 *
 *      Description:    For the general utility functionalities
 *
 *      Version:        1.0
 *      Created:        02/06/2017 04:45:04 PM
 *      Revision:       none
 *      Compiler:       gcc
 *
 *      Author:         Peiming (Peter) Huang
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>	// for bool return type
#include <ctype.h>		// for isdigit()
#include "terms.h"

/**
 * makes sure a file exists
 * @param fName (file name to check)
 * @return bool => return true if file exist, otherwise, exit(1)
 */
bool checkFile(char * fName);

/**
 * produces the usage information for end users
 */
void usage();

/**
 * removes the "chr" prefix from cName (chromosome name) passed in
 * @param cName
 * @return char * without 'chr' prefix
 */
char * removeChr(char * cName);

/**
 * converts fractions into percentages with 2 decimal positions
 * @param num
 * @param dom
 * @return
 */
double getPercentage(int num, int dom);

/**
 * check if the input string is a number
 * @param inStr
 * @return bool
 */
bool isNum(const char *inStr);

/**
 * process user input options and ensure they are correct!
 * @param argv[]: an array of strings that are used to decide user input options
 */
void processUserOptions(User_Input *user_inputs, int argc, char *argv[]);

/**
 * this is a helper function that is used to create a file name from the existing base file name
 * @param base_name: the file name to be created based on the base_name
 * @param file_in: the file_name to be created into
 * @param string_to_append: the string to be attached to the base name
 */
void createFileName(char *base_name, char **file_in, char *string_to_append);

/**
 * this method is used to initialize the user input structure
 * @return returns an instance of User_Input upon the memory allocation success. On failure, the entire program will exit!
 */
User_Input * userInputInit();

/**
 * This method is used to clean up the memory after the User_Input variable is done!
 * @param User_Input: a user defined struct
 */
void userInputDestroy(User_Input *user_inputs);

/**
 * calculate the total number of genome bases from the bam/cram header info
 * @param header: bam/cram/sam header info
 * @param stats_info: to store the final number of total genome bases
 */
void fetchTotalGenomeBases(bam_hdr_t *header, Stats_Info *stats_info);

/**
 * for Coverage_Stats variable initialization
 * @return an instance of Coverage_Stats upon successful memory allocation
 */
Coverage_Stats * coverageStatsInit();

/**
 * This function is used to clean the khash_t (int key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 */
void cleanKhashInt(khash_t(m32) *hash_to_clean);

/**
 * This function is used to clean the khash_t (string key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 */
void cleanKhashStr(khash_t(str) *hash_to_clean, uint8_t type);

/**
 * This function is used to dynamically allocate string and grow it accordingly
 * @param storage_str, the string to be dynamically allocated
 * @param str_in, the string to be copied
 * @return allocated string pointer
 */
void dynamicStringAllocation(char *str_in, char **storage_str);

/**
 * Find the corresponding chromosome index from input chromosome id
 * @param header: the bam/cram/sam header
 * @param chrom_id
 * @return the index of the corresponding chromosome id
 */
uint32_t getChromIndexFromID(bam_hdr_t *header, char *chrom_id);
//uint32_t getChromIndexFromID(uint16_t chr_count, char *chrom_id);

/**
 * initialize the Chromosome_Tracking variable
 * @return an instance of Chromosome_Tracking upon successful
 */
Chromosome_Tracking * chromosomeTrackingInit();

/**
 * This function is used to update all members for the Chromosome_Tracking variable
 * @param chrom_tracking: a Chromosome_Tracking used for tracking
 * @param chrom_id: the current chromosome id
 * @param chrom_len: the length of current chromosome
 * @param index: the members in Chromosome_Tracking variable are stored in arrays, using index will help locate the chromosome info
 */
void chromosomeTrackingUpdate(Chromosome_Tracking *chrom_tracking, char *chrom_id, uint32_t chrom_len, int index);

/**
 * This function is used to update members for the Chromosome_Tracking variable
 * @param chrom_tracking: a Chromosome_Tracking used for tracking
 * @param chrom_id: the current chromosome id
 * @param chrom_len: the length of current chromosome
 * @param index: the members in Chromosome_Tracking variable are stored in arrays, using index will help locate the chromosome info
 * @param status: the status of current chromosomd id
 */
//void chromosome_tracking_update(Chromosome_Tracking *chrom_tracking, char *chrom_id, uint32_t chrom_len, int index, int status);

/**
 * To clean up the allocated memory for chromosome tracking
 * @param chrom_tracking: the tracking variable to be cleaned
 */
void chromosomeTrackingDestroy(Chromosome_Tracking * chrom_tracking);

/**
 * To locate the index of a chromosome id in an array give the chromosome id information
 * @param chrom_id
 * @param chrom_tracking: the Chromosome_Tracking variable to track the status of chromosome processed
 */
int32_t locateChromosomeIndexForChromTracking(char *chrom_id, Chromosome_Tracking *chrom_tracking);

int32_t locateChromosomeIndexForRegionSkipMySQL(char *chrom_id, Regions_Skip_MySQL *regions_in);

/**
 * Initialize the member of the Stats_Info variable
 * @return an instance of Stats_Info upon successful
 */
Stats_Info * statsInfoInit();

/**
 * to destroy everything allocated for stats_info
 * @param stats_info
 */
void statsInfoDestroy(Stats_Info *stats_info);

/**
 * it will set the coverage for all of the Ns regions in the genome to zero
 * @param chrom_id
 * @param Ns_buffer_hash
 * @param chrom_tracking
 */
//void zeroAllNsRegions(char *chrom_id, khash_t(str) *Ns_buffer_hash, Chromosome_Tracking *chrom_tracking);
void zeroAllNsRegions(char *chrom_id, Bed_Info *Ns_info, Chromosome_Tracking *chrom_tracking);

/**
 * To add value into a hash table by the key
 * @param hash_in: the hash table to be modified
 * @param pos_key: the key for a specific position
 * @param val: value to be added
 */
void addValueToKhashBucket16(khash_t(m16) *hash_in, uint16_t pos_key, uint16_t val);
void addValueToKhashBucket32(khash_t(m32) *hash_in, uint32_t pos_key, uint32_t val);

/**
 * Get value from the hash table by the key
 * @param hash_in
 * @param key
 * @return the value pointed by the key
 */
uint32_t getValueFromKhash32(khash_t(m32) *hash32, uint32_t pos_key);
uint16_t getValueFromKhash16(khash_t(m16) *hash16, uint32_t pos_key);

/** converts fractions into percentages with 2 decimal positions
 * @param num: the numeric value 
 * @param dom: the donominator value
 * @return a float value with 2 decimal points
 */
float calculatePercentage(uint32_t num, uint32_t dom);

/** 
 * combine all the coverage stats from individual thread to stats_info
 * @param stats_info
 * @param cov_stats
 */
void combineCoverageStats(Stats_Info *stats_info, Coverage_Stats *cov_stats);

void printLowCoverageGeneStructure(Low_Coverage_Genes *low_cov_genes);

//int compare(const Gene_Coverage *gene_coverage1, const Gene_Coverage *gene_coverage2);
int compare(const void *gene_coverage1, const void *gene_coverage2);

int compare2(const void *transcript_cov_pct1, const void *transcript_cov_pct2);

#endif //UTILS_H
