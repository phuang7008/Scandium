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
char * removechr(char * cName);

/**
 * removes the "chr" prefix from cName (chromosome name) passed in
 * @param cName
 * @return char * without 'chr' prefix
 */
char * removechr(char * cName);

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
bool isnum(const char *inStr);

/**
 * process user input options and ensure they are correct!
 * @param argv[]: an array of strings that are used to decide user input options
 */
void processUserOptions(User_Input *user_inputs, int argc, char *argv[]);

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
 * This function is used to clean the khash_t (int key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 */
void clean_khash_int(khash_t(m32) *hash_to_clean);

/**
 * This function is used to clean the khash_t (string key) hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the allocated memories
 */
void clean_khash_str(khash_t(str) *hash_to_clean);

/**
 * Find the corresponding chromosome index from input chromosome id
 * @param header: the bam/cram/sam header
 * @param chrom_id
 * @return the index of the corresponding chromosome id
 */
short get_chrom_index_from_id(bam_hdr_t *header, char *chrom_id);

/**
 * This function is used to initialize all members for the Chromosome_Tracking variable
 * @param chrom_tracking: a Chromosome_Tracking used for tracking
 * @param chrom_id: the current chromosome id
 * @param chrom_len: the length of current chromosome
 * @param index: the members in Chromosome_Tracking variable are stored in arrays, using index will help locate the chromosome info
 */
void chromosome_tracking_init(Chromosome_Tracking *chrom_tracking, char *chrom_id, uint32_t chrom_len, int index);

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
void chromosome_tracking_destroy(Chromosome_Tracking * chrom_tracking);

#endif //UTILS_H
