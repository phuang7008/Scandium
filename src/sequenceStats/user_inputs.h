/*
 * =====================================================================================
 *
 *      Filename:       user_inputs.h
 *
 *      Description:    For the functionalities of user inputs
 *
 *      Version:        1.0
 *      Created:        03/16/2020 04:45:04 PM
 *      Revision:       none
 *      Compiler:       gcc
 *
 *      Author:         Peiming (Peter) Huang (phuang@bcm.edu)
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef USER_INPUTS_H
#define USER_INPUTS_H

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
void recordHGMD(Databases *dbs, khash_t(khStrInt) *hgmd_genes, khash_t(khStrInt) *hgmd_transcripts);

/* a helper function that is used to clean-up the annotation wrapper
 * @param annotation_wrapper: the variable to be cleaned
 */
void annotationWrapperDestroy(Annotation_Wrapper *annotation_wrapper);

/**
 * produces the usage information for end users
 */
void usage();

/**
 * check if the input string is a number
 * @param inStr: a variable need to be checked
 * @return true or false
 */
bool isNumber(const char *inStr);

/**
 * check if the input string is a float number
 * @param str, the string to be coverted
 * @param dest, the value to be copied to
 * @return true or false
 */
bool isFloat(const char *str, float *dest);

/**
 * Handle an input with multiple user defined annotation files
 * @param optarg: the annotation files at the command line separated by comma ','
 * @param user_inputs: a variable to store all of the user provided inputs
 */
void formAnnotationFileArray(char* optarg, User_Input *user_inputs);

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
 * This function is used to check if the input bam file named as .cram; or the input cram file named as .bam
 * @param user_inputs: it store the file name information
 * @param std: a samFile open file hand pointer
 */ 
void checkFileExtension(User_Input *user_inputs, samFile* sfd);

/**
 * This function is used to obtain the input filename extension
 * @param filename: the input filename
 * @return char*: the file extension
 */
const char* getFileExtension(const char* filename);

/**
 * Get the base name of the file including the extension, but without the full path
 * @param filepath: a full file path
 * @return char*: the basename from the file path
 */
const char* baseFilename(char const *filepath);

/**
 * Obtain the base file name without the extension from a filepath
 * @param user_inputs: variable that stores the user inputs related information
 */
void getBaseFilenameWithoutExtension(User_Input *user_inputs);

/*
 * The method is added for cram file input as cram file processing needs reference sequence
 * @param fn_ref: the file path to the reference sequences
 * @return the file path to the reference sequence's fai file
 */
char* getReferenceFaiPath(const char *fn_ref);

void setupMySQLDB(Databases **dbs, User_Input *user_inputs);

/*
 * since hg37 and hg38 used different chromosome nameing convention, 
 * we need to ensure that the input files follow this convention
 * @param user_inputs: variable that stores the user inputs
 * @param chrom_id: chromosome id to be checked
 */
void checkChromosomeID(User_Input *user_inputs, char* chrom_id);

/*
 * it is a help function that is used to load all chromosomes need to be processed
 * @param primary_chromosome_hash: a hash table that is used to store the primary chromosomes for quick lookup
 * @param user_inputs: different version of human genomes define chromosome id differently. For other genomes, users will have to make the adjustment
 */
void loadWantedChromosomes(khash_t(khStrInt) *primary_chromosome_hash, User_Input *user_inputs, Stats_Info *stats_info);

//int8_t strcasecmp(char const *a, char const *b);

/* need to check for the naming convention between the chromosome bed file and the input bam/cram files
 * @param header: the bam/cram header to store the bam/cram header info
 * @param wanted_chromosome_hash: a kh_hash table to store the chromosome ids that need to be processed
 */
void checkNamingConvention(bam_hdr_t *header, khash_t(khStrInt)* wanted_chromosome_hash);

#endif //USER_INPUTS_H
