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
 *      Author:         Peiming (Peter) Huang (phuang@bcm.edu)
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>	// for bool return type
#include <ctype.h>		// for isdigit()
#include <errno.h>
#include "terms.h"

/* 
 * Split a string to a word hash
 * Note: arrayPtr should be inialized outside the function before passing it to the call!
 * We use a dynamic way for memory allocation!

 * The following is an example what two original strings look like and we need to combine them togetner!
 *    NM_022834_exon_0;NM_199121_exon_0    CCDS27_exon_0;CCDS28_exon_0 OTTHUMT00000008291_exon_0;OTTHUMT00000008294_exon_0 .
 *    NM_001242361_exon_29;NR_031580_exon_0    CCDS55599_exon_29;CCDS559_exon_29   OTTHUMT00000023045_exon_29  hsa-mir-761

 * Here the first delimiter is "\t", while the second delimiter is ";"

 * The StringArray will store all the information about the detailed exon annotation 
 * The following illustrates the StringArray structure
 * StringArray ->
                   RefSeq              CCDS              VEGA/Gencode          miRNA
                     0                   1                    2                  3   
               char **theArray    char **theArray        char **theArray     char **theArray
                  -> theArray[0] -> theArray[1] -> theArray[2] ... -> theArray[size]

 * @param stringPtr, a string pointer. It's content will be splitted into words
 * @param arrayPtr, a StringArray pointer. It is used to store string array information
 * @param index the array index of khashArrayPtr
*/ 
void splitStringToKhash(char *stringPtr, khash_t(khStrInt) **khashArrayPtr, uint8_t index);

/*
 * it is used to calculation the size of low coverage regions from a StrInt Hash table
 * In addition, it will generate the low coverage regions in sorted order in string format
 * @param low_cov_regions: a StrInt hash table contain the start and end position of low coverage region
 * @param output: output string array
 * @return total low coverage region size
 */
uint32_t processLowCovRegionFromKhash(khash_t(khStrInt) *low_cov_regions, char **output);

/* this is used to process low coverage regions for output
 * @param low_cov_regions: a string array that contains all low coverage regions
 * @param output: combined all low coverage regions
 * @return number of low coverage regions
 */
uint32_t processLowCovRegionFromStrArray(StringArray *low_cov_regions, char **output);

/* 
 * It is used to clean the kash_t (char* as key, but string array as value) hash table
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories including keys and char* array
 */
void cleanKhashStrStrArray(khash_t(khStrStrArray) * hash_to_clean);

/**
 * This function is used to clean the khash_t (string key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 * @param type: if it is set to 1, it will clean Temp_Coverage_Array struct variable as well
 */
void cleanKhashStr(khash_t(str) *hash_to_clean, uint8_t type);

/**
 * This is used to free the memory used by the hash_to_clean
 * @param hash_to_clean: a khash_t variable to be cleaned
 */
void cleanKhashStrStr(khash_t(khStrStr) * hash_to_clean);

/**
 *
 */
void cleanGeneTranscriptPercentage(khash_t(khStrGTP) *gene_transcript_percentage_hash);

/**
 * This function is used to dynamically allocate string and copy value accordingly
 * @param str_in, the string to be copied
 * @param storage_str, the string to be dynamically allocated/expanded
 */
void dynamicStringAllocation(char *str_in, char **storage_str);

/**
 * This function is used to dynamically allocate string and grow it accordingly
 * @param str_in, the string to be added
 * @param storage_str, the string to be dynamically allocated/expanded
 */
void dynamicStringExpansion(char *str_in, char **storage_str);

int32_t locateChromosomeIndexForRegionSkipMySQL(char *chrom_id, Regions_Skip_MySQL *regions_in);

int32_t findChromsomeIndex(Chromosome_Tracking *chrom_tracking, bam_hdr_t *header, int32_t index);

int32_t findTargetBufferIndex(Target_Buffer_Status *target_buffer_status, int32_t number_of_chromosomes, char* chrom_id);

char* createTmpFileName(char* base_file_name, char* string_to_add);

/**
 * Initialize the member of the Stats_Info variable
 * @param stats_info: an instance of Stats_Info to store the whole genome coverage stat info
 */
void statsInfoInit(Stats_Info *stats_info, User_Input *user_inputs);

/**
 * for XXX_Coverage_Stats variable initialization
 * @param xxx_cov_stats: an instance of XXX_Coverage_Stats to store the coverage statistics for current chrom
 */
void WGSCoverageStatsInit(WGS_Coverage_Stats * wgs_cov_stats);
void readCoverageStatsInit(Read_Coverage_Stats * read_cov_stats);
void captureCoverageStatsInit(Capture_Coverage_Stats * capture_cov_stats);

/**
 * to destroy everything allocated for stats_info
 * @param stats_info
 */
void statsInfoDestroy(Stats_Info *stats_info, User_Input *user_inputs);

/**
 * To add value into a hash table by the key
 * @param hash_in: the hash table to be modified
 * @param pos_key: the key for a specific position
 * @param val: value to be added
 */
void addValueToKhashBucket16(khash_t(m16) *hash_in, uint16_t pos_key, uint16_t val);
void addValueToKhashBucket32(khash_t(m32) *hash_in, uint32_t pos_key, uint32_t val);
void addValueToKhashBucketStrStr(khash_t(khStrStr) *hash_in, char *key, char * val);

/**
 * Get value from the hash table by the key
 * @param hash_in
 * @param key
 * @return the value pointed by the key
 */
uint32_t getValueFromKhash32(khash_t(m32) *hash32, uint32_t pos_key);
uint16_t getValueFromKhash16(khash_t(m16) *hash16, uint32_t pos_key);
char * getValueFromKhashStrStr(khash_t(khStrStr) *hash_in, char* key);

/** converts fractions into percentages with 2 decimal positions
 * @param num: the numeric value 
 * @param dom: the donominator value
 * @return a float value with 2 decimal points
 */
float calculatePercentage64(uint64_t num, uint64_t dom);
float calculatePercentage32(uint32_t num, uint32_t dom);
float calculatePercentage32_64(uint32_t num, uint64_t dom);

/** 
 * combine all the coverage stats from each individual thread and store them in the stats_info
 * @param stats_info: a storage place for all summarized sequencing stats
 * @param cov_stats: detailed coverage stats
 */
void combineCoverageStats(Stats_Info *stats_info, Stats_Info *tmp_stats_info, User_Input *user_inputs);

void copyReadCoverageStats(Read_Coverage_Stats *read_cov_stats, Read_Coverage_Stats *tmp_read_cov_stats);
void copyWGSCoverageStats(WGS_Coverage_Stats *wgs_cov_stats, WGS_Coverage_Stats *tmp_wgs_cov_stats);
void copyCaptureCoverageStats(Capture_Coverage_Stats **capture_cov_stats, Capture_Coverage_Stats **tmp_capture_cov_stats, User_Input *user_inputs);

/* obtain the hash key for the integer value passed in                                                                
 * the keys are defined every 1000 position for quick lookup
 * For example: 3241645 (pos) -> 3241000 (key)
 * @param position_in, the position for the key                                                               
 * @return hash key                                                                                           
 */
uint32_t getHashKey(uint32_t position_in);

/*
 * To initialize each khash_t(khStrLCG) bucket
 * @param low_cov_gene_hash: a khash_t hash that is used to store low coverage gene info
 * @param key_in: the hash key used for a specific bucket. If the key doesn't exist, initialize it accordingly
 */
void lowCoverageGeneHashBucketKeyInit(khash_t(khStrLCG) *low_cov_gene_hash, char *key_in);

/*
 * To initialize the Gene_Coverage structure
 * @param gc: a instance of Gene_coverage to be initialized
 */
void geneCoverageInit(Gene_Coverage *gc);

/*
 * Add gene transcript info to the corresponding hash table
 * @param gene_symbol: gene symbol name
 * @param transcript_name: transcript name
 * @param gene_transcripts: a hash table used to store transcripts for each everything gene encountered
 * @param seen_transcripts: if the transcript has been seen before, don't add it again
 */
void addToGeneTranscriptKhashTable(char *gene_symbol, char *transcript_name, khash_t(khStrStrArray) *gene_transcripts, khash_t(khStrInt) *seen_transcript);

/*
 * Used to copy the low coverage region (a char*) to the Gene_Coverage
 * @param gc1: the source of low coverage region
 * @param gc2: the destination of low coverage region
 * @param copy_gene: also copy all othere gene related information such as exon start and end etc
 */ 
void copyGeneCoverageLowCovRegions(Gene_Coverage* gc1, Gene_Coverage* gc2, bool copy_gene);

/*
 * to merge to low coverage regions if they overlap
 * @param low_cov_regions_hash: the low coverage hash table that contains the overlapped regions
 * @param mergedArray: a temp array to store merged low coverage regions
 * @param size_in: need to know the size of regions which is a string array
 * @param cds_t_start: the targeted cds start position as we are not interested in anything before that
 * @param cds_t_end: the targeted cds end position as we are not interested in anything beyond that
 */
void mergeLowCovRegions(khash_t(khStrInt) *low_cov_regions_hash, StringArray *mergedArray, uint32_t size_in, uint32_t cds_t_start, uint32_t cds_t_end);

/* it is used for reporting purpose. If there is only one low coverage regions, we don't have to do the merge.
 * We can just print its content accordingly
 * @param low_cov_regions_hash: a hash table that contains the low coverage region
 * @param mergedArray: just copy everything into thie mergedArray for outputting
 */
void getOneLowCovRegions(khash_t(khStrInt) *low_cov_regions_hash, StringArray *mergedArray);

/*
 * help function that is used to print the low coverage gene info for debugging
 * @param low_cov_genes: a struct used to store all the information related to low coverage genes
 */
void printLowCoverageGeneStructure(Low_Coverage_Genes *low_cov_genes);

/*
 * calculate uniformity metrics
 * @param stats_info: used to store the uniformity metrics
 * @param user_inputs: it store the file name -> open for writing
 * @param primary_chromosome_hash: lookup hash to store primary chromosomes only
 * @param autosome: it is used to indicate if it is only going to handle non-sex chromosomes
 * @param primary_chromosomes_only: it is used to indicate if it is only going to handle primary chromosome without alt or decoys
 */
void calculateUniformityMetrics(Stats_Info *stats_info, User_Input *user_inputs, khash_t(khStrInt) *wanted_chromosome_hash, khash_t(m32) *cov_freq_dist, bool autosome, bool primary_chromosomes_only);

/*
 * Instead of taking points around Mode evenly, it will always try to pick the higher points around Mode
 * @param peak: the Mode or Mean or Median. But here we are going to use Mode
 * @param cov_freq_dist: it is hash table that is used to store the sequencing coverage frequency distribution for quick lookup
 * @param user_inputs: the peak size (default 7) could be defined by the end user and store at user_inputs
 */
uint64_t dynamicCalculateAreaUnderHistogram(uint32_t peak, khash_t(m32) *cov_freq_dist, User_Input *user_inputs);

void outputFreqDistribution(User_Input *user_inputs, khash_t(m32) *cov_freq_dist);

void set_peak_size_around_mode(Stats_Info *stats_info, User_Input *user_inputs);

/**
 * It is used to print a string array before (OR after sorting) for viewing and comparison.
 * @param strings_in: the string array to be printed!
 * @param length_in: specifies the size of input string array
 */
void print_string_array(char** strings_in, size_t length_in);

#endif //UTILS_H
