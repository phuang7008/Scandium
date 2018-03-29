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
 * @return bool => return true if file exist, otherwise, exit(EXIT_FAILURE)
 */
bool checkFile(char * fName);

/**
 * check the size of input bam/cram file
 */
uint64_t check_file_size(const char *filename);

/* Split a string to a word array
 * Note: arrayPtr should be inialized outside the function before passing it to the call!
 * We use a dynamic way for memory allocation!

 * The following is an example what two original strings look like and we need to combine them togetner!
 *    NM_022834_exon_0;NM_199121_exon_0    CCDS27_exon_0;CCDS28_exon_0 OTTHUMT00000008291_exon_0;OTTHUMT00000008294_exon_0 .
 *    NM_001242361_exon_29;NR_031580_exon_0    CCDS55599_exon_29;CCDS559_exon_29   OTTHUMT00000023045_exon_29  hsa-mir-761

 * Here the first delimiter is "\t", while the second delimiter is ";"

 * The stringArray will store all the information about the detailed exon annotation 
 * The following illustrates the stringArray structure
 * stringArray ->
                   RefSeq              CCDS              VEGA/Gencode          miRNA
                     0                   1                    2                  3   
               char **theArray    char **theArray        char **theArray     char **theArray
                  -> theArray[0] -> theArray[1] -> theArray[2] ... -> theArray[size]

 * @param stringPtr, a string pointer. It's content will be splitted into words
 * @param arrayPtr, a stringArray pointer. It is used to store string array information
*/ 
void splitStringToArray(char *stringPtr, stringArray *arrayPtr);

/**
 * to remove duplated words from a give string
 * @param array1, a stringArrat pointer that contains an string array with duplicates to be removed!
 * @param array2, a stringArray pointer that contains an string array that its string content need to be removed from array1
 */
void removeDuplicatesFromStringArray(stringArray *array1, stringArray *array2);

void stringArrayDestroy(stringArray *arrayIn);

/*
 * find if a string contains a substr (case insensative)
 * @param string1, the string to be checked
 * @param substr2, the substring to search for
 */
char* stristr( const char* string1, const char* substr2 );

void annotationWrapperDestroy(Annotation_Wrapper *annotation_wrapper);

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
 * Output the User Input Options to the end user so he/she can double-check if all the options are correct!
 * user_inputs, the object contains all the user input options
 */
void outputUserInputOptions(User_Input *user_inputs);

/**
 * this is a helper function that is used to create a file name from the existing base file name
 * @param output_dir: the directory name to store all the output reports/files
 * @param base_name: the file name to be created based on the base_name
 * @param file_in: the file_name to be created into the output_directory
 * @param string_to_append: the string to be attached to the base name
 */
void createFileName(char *output_dir, char *base_name, char **file_in, char *string_to_append);

void writeHeaderLine(char *file_in, uint8_t type);

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
 * @param cov_stats: an instance of Coverage_Stats to store the coverage statistics for current chrom
 */
void coverageStatsInit(Coverage_Stats * cov_stats);

/**
 * This function is used to clean the khash_t (int key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 */
void cleanKhashInt(khash_t(m32) *hash_to_clean);

void cleanKhashStrInt(khash_t(khStrInt) *hash_to_clean);

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
 * @param stats_info: an instance of Stats_Info to store the whole genome coverage stat info
 */
void statsInfoInit(Stats_Info *stats_info);

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
void zeroAllNsRegions(char *chrom_id, Bed_Info *Ns_info, Chromosome_Tracking *chrom_tracking, Target_Buffer_Status *target_buffer_status);

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

void getDatabaseName(User_Input *user_inputs, char* db_name_in_out);

/** 
 * combine all the coverage stats from individual thread to stats_info
 * @param stats_info
 * @param cov_stats
 */
void combineCoverageStats(Stats_Info *stats_info, Coverage_Stats *cov_stats);

void printLowCoverageGeneStructure(Low_Coverage_Genes *low_cov_genes);

/**
 * It is used to sort the input Gene_Coverage variable based on gene_symbol
 * @param gene_coverage1, the Gene_Coverage input variable 1
 * @param gene_coverage2, the Gene_Coverage input variable 2
 */
int compare(const void *gene_coverage1, const void *gene_coverage2);

/**
 * It is used to sort the input Transcript_Cov_PCT variable based on refseq_name
 * @param transcript_cov_pct1, the Transcript_Cov_PCT input variable 1
 * @param transcript_cov_pct2, the Transcript_Cov_PCT input variable 2
 */
int compare2(const void *transcript_cov_pct1, const void *transcript_cov_pct2);

/**
 * The following comparison is used to compare the refseq_name array variable (string array)
 * @param refseq_array1, the string array input 1
 * @param refseq_array2, the string array input 2
 */
int compare3(const void *refseq_array1, const void *refseq_array2);

/**
 * It is used to print a string array before (OR after sorting) for viewing and comparison.
 * strings_in: the string array to be printed!
 */
void print_string_array(char** strings_in, size_t length_in);

#endif //UTILS_H
