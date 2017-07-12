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
 * This function is used to clean the khash_t hash table used by the users
 * @param hash_to_clean: loop through the hash table to clean all the memories by the hash
 */
void clean_khash(khash_t(32) *hash_to_clean);

#endif //UTILS_H
