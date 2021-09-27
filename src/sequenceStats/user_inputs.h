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

#include <ctype.h>		// for isdigit()
#include <errno.h>
#include <getopt.h>
#include <stdbool.h>	// for bool return type
#include "terms.h"
#include "utility.h"

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
 * Handle an input with multiple user defined annotation files
 * @param optarg: the annotation files at the command line separated by comma ','
 * @param user_inputs: a variable to store all of the user provided inputs
 * @param type: the type of file array: 1 for target file, 2 for annotation file
 */
//void formTargetAnnotationFileArray(char* optarg, User_Input *user_inputs, uint8_t type);
void formTargetAnnotationFileArray(khash_t(khStrStr) *capture_files, khash_t(khStrStr) *annotation_files, User_Input *user_inputs);

void readTargetAnnotationFilesIn(User_Input *user_inputs, char* file_in, int type);

void checkInputCaptureAndAnnotationFiles(User_Input *user_inputs);

void checkRepeatedCaptureFiles(User_Input *user_inputs);

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

void setupOutputReportFiles(User_Input *user_inputs);

/*
 * for output only. So that users will understand the meaning of each output column
 * @param file_in: the output file handle
 * @param user_inputs: it contains all of the options from user's command line inputs
 * @param annotation_file_index: the file index in annotation file array
 * @param type: different type of files will have different header
 */
void writeHeaderLine(char *file_in, User_Input *user_inputs, int annotation_file_index, uint8_t type);

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
 * Obtain the base file name without the extension from a filepath
 * @param user_inputs: variable that stores the user inputs related information
 * @param type: the type of files used => 1 for target files, 2 for annotation files
 */
void getBaseFilenameWithoutExtension(User_Input *user_inputs, uint8_t type);

/*
 * This method will free all the memories allocated for storing file names (array)
 * @param f_size: the file array size
 * @param file_array: the input file array
 */
void cleanCreatedFileArray(uint8_t f_size, char **file_array);

void setupMySQLDB(Databases **dbs, User_Input *user_inputs);


#endif //USER_INPUTS_H
