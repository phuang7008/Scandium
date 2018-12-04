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

/*
 * It is used to fetch all the HGMD info from the database
 * @param dbs: a Databases variable used to store all info related to the opened MySQL database
 * @param hgmd_genes: a variable that is used to store all the HGMD genes
 * @param hgmd_transcripts: a variable that is used to store all the HGMD transcripts
 */
void recordHGMD(Databases *dbs, User_Input *user_inputs, khash_t(khStrInt) *hgmd_genes, khash_t(khStrInt) *hgmd_transcripts);

/* 
 * Split a string to a word hash
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
void splitStringToKhash(char *stringPtr, khash_t(khStrInt) **khashArrayPtr, uint8_t index);

void stringArrayDestroy(stringArray *arrayIn);

/*
 * find if a string contains a substr (case insensative)
 * @param string1, the string to be checked
 * @param substr2, the substring to search for
 */
char* stristr( const char* string1, const char* substr2 );

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
uint32_t processLowCovRegionFromStrArray(stringArray *low_cov_regions, char **output);

/* a helper function that is used to clean-up the annotation wrapper
 * @param annotation_wrapper: the variable to be cleaned
 */
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
 * check if the input string is a number
 * @param inStr: a variable need to be checked
 * @return true or false
 */
bool isNum(const char *inStr);

/**
 * process user input options and ensure they are correct!
 * @param argv[]: an array of strings that are used to decide user input options
 */
void processUserOptions(User_Input *user_inputs, int argc, char *argv[]);

/**
 * Output the User Input Options to the end user so he/she can double-check if all the options are correct!
 * @param user_inputs: a variable contains all the user input options
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

/*
 * for output only. So that users will understand the meaning of each output column
 * @param file_in: the output file handle
 * @param type: different type of files will have different header
 */
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
 * calculate the total number of genome bases and load chromosome info from the bam/cram header info
 * @param header: bam/cram/sam header info
 * @param stats_info: to store the final number of total genome bases
 */
void loadGenomeInfoFromBamHeader(khash_t(khStrInt) *wanted_chromosome_hash, bam_hdr_t *header, Stats_Info *stats_info, User_Input *user_inputs);

/**
 * for Coverage_Stats variable initialization
 * @param cov_stats: an instance of Coverage_Stats to store the coverage statistics for current chrom
 */
void coverageStatsInit(Coverage_Stats * cov_stats);

/**
 * This function is used to clean the khash_t (uint32_t key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 */
void cleanKhashInt(khash_t(m32) *hash_to_clean);

/*
 * It is used to clean the kash_t (char* as key, with uint32_t as value) hash table
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories including keys
 */
void cleanKhashStrInt(khash_t(khStrInt) *hash_to_clean);

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
 *
 */
void cleanGeneTranscriptPercentage(khash_t(khStrGTP) *gene_transcript_percentage_hash);

/**
 * This function is used to dynamically allocate string and grow it accordingly
 * @param str_in, the string to be added
 * @param storage_str, the string to be dynamically allocated/expanded
 */
void dynamicStringAllocation(char *str_in, char **storage_str);

/**
 * initialize the Chromosome_Tracking variable, this approach will process all chromosomes
 * @param header: a bam header that contains all the chromosome information
 * @return an instance of Chromosome_Tracking upon successful
 */
void chromosomeTrackingInit1(uint32_t num_of_chroms, Chromosome_Tracking *chrom_tracking, khash_t(khStrInt) *wanted_chromosome_hash, bam_hdr_t *header);

/**
 * Initialize the chromosome_tracking variable using user specified region file
 * This approach will process those chromosomes specified by user only
 * @param wanted_chromosome_hash: a kh_hash table stores chromosomes to be processed
 * @return an instance of Chromosome_Tracking upon successful 
 */
uint32_t chromosomeTrackingInit2(khash_t(khStrInt) *wanted_chromosome_hash, Chromosome_Tracking *chrom_tracking, bam_hdr_t *header);

/**
 * This function is used to update all members for the Chromosome_Tracking variable
 * @param chrom_tracking: a Chromosome_Tracking used for tracking
 * @param chrom_id: the current chromosome id
 * @param chrom_len: the length of current chromosome
 * @param index: the members in Chromosome_Tracking variable are stored in arrays, using index will help locate the chromosome info
 */
void chromosomeTrackingUpdate(Chromosome_Tracking *chrom_tracking, char *chrom_id, uint32_t chrom_len, int index);

/**
 * To clean up the allocated memory for chromosome tracking
 * @param chrom_tracking: the tracking variable to be cleaned
 */
void chromosomeTrackingDestroy(Chromosome_Tracking * chrom_tracking);

/**
 * To locate the index of a chromosome id in an array give the chromosome id information
 * @param chrom_id: current chromosome id to be handled
 * @param chrom_tracking: the Chromosome_Tracking variable to track the status of chromosome processed
 * @return a index at the tracking array
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
 * @param chrom_id: current chromosome id to be handled
 * @param Ns_info: the detailed Ns info
 * @param chrom_tracking: a storage used to track each chromosome in details
 * @param target_buffer_status: it is used to tell which regions are targets and which regions are buffer. 
 *		Sometimes, Ns regions will overlap with the buffer regions
 */
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
float calculatePercentage64(uint64_t num, uint64_t dom);
float calculatePercentage32(uint32_t num, uint32_t dom);

/** 
 * combine all the coverage stats from each individual thread and store them in the stats_info
 * @param stats_info: a storage place for all summarized sequencing stats
 * @param cov_stats: detailed coverage stats
 */
void combineCoverageStats(Stats_Info *stats_info, Coverage_Stats *cov_stats);

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
void mergeLowCovRegions(khash_t(khStrInt) *low_cov_regions_hash, stringArray *mergedArray, uint32_t size_in, uint32_t cds_t_start, uint32_t cds_t_end);

/* it is used for reporting purpose. If there is only one low coverage regions, we don't have to do the merge.
 * We can just print its content accordingly
 * @param low_cov_regions_hash: a hash table that contains the low coverage region
 * @param mergedArray: just copy everything into thie mergedArray for outputting
 */
void getOneLowCovRegions(khash_t(khStrInt) *low_cov_regions_hash, stringArray *mergedArray);

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

/*
 * it is a help function that is used to load all chromosomes need to be processed
 * @param primary_chromosome_hash: a hash table that is used to store the primary chromosomes for quick lookup
 * @param user_inputs: different version of human genomes define chromosome id differently. For other genomes, users will have to make the adjustment
 */
void loadWantedChromosomes(khash_t(khStrInt) *primary_chromosome_hash, User_Input *user_inputs, Stats_Info *stats_info);

/*
 * the following comparison is used to compare the int array
 * @param val1: int array 1 used for comparison
 * @param val2: int array 2 used for comparison
 */
int compare(const void * val1, const void * val2);

/**
 * It is used to print a string array before (OR after sorting) for viewing and comparison.
 * @param strings_in: the string array to be printed!
 * @param length_in: specifies the size of input string array
 */
void print_string_array(char** strings_in, size_t length_in);

//int8_t strcasecmp(char const *a, char const *b);

/* need to check for the naming convention between the chromosome bed file and the input bam/cram files
 * @param header: the bam/cram header to store the bam/cram header info
 * @param wanted_chromosome_hash: a kh_hash table to store the chromosome ids that need to be processed
 */
void checkNamingConvention(bam_hdr_t *header, khash_t(khStrInt)* wanted_chromosome_hash);

#endif //UTILS_H
