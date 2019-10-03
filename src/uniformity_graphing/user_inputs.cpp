/*
 * =====================================================================================
 *
 *       Filename:  user_inputs.cpp
 *
 *    Description:  the detailed inplementation of user input options
 *
 *        Version:  1.0
 *        Created:  01/05/2018 02:32:03 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */

#include <iostream>
#include <unistd.h>			// for file access() and getopt()
#include <cstdlib>			// for exit()
#include "user_inputs.hpp"

// constructor
User_Inputs::User_Inputs() {
	input_file  = "";
	output_file = "";
	chrom_id    = "";
	reference_version = "hg37";

	upper_bound = 150;
	buffer_size = 20;

	smoothing_type = 1;
	percentage = 2;
}

void User_Inputs::usage() {
	std::cout << "USAGE: draw_graph [options]" << std::endl;
	std::cout << "\t-i: the input file name [Mandatory]" << std::endl;
	std::cout << "\t-o: the output file name [Mandatory]" << std::endl;
	std::cout << "\t-c: the chromosome id [Mandatory]" << std::endl;
	std::cout << std::endl;
	std::cout << "\t-b: the buffer size for smoothing noises (used for type 2 approach)! (Default: 20)" << std::endl;
	std::cout << "\t-p: the percentage used for the (type 1) gVCF BLOCKAVE_'p'p formula (should be >=2)! (Default: 2 as 200%)" << std::endl;
	std::cout << "\t-r: the reference version (for human)! (Default: hg37)" << std::endl;
	std::cout << "\t-t: the type of method used for smoothing (type 1: gVCF, type 2: lower-upper range)! (Default: 1)" << std::endl;
	std::cout << "\t-u: the upper bound used to generate the range file (Default: 150)" << std::endl;
}

// Get command line arguments in and check the sanity of user inputs 
//
void User_Inputs::processUserOptions(int argc, char* argv[]) {
	int arg;

	if (argc == 1) {
		//fprintf (stderr, "No options provided\n");
		usage();
		exit(1);
	}

	//When getopt returns -1, no more options available
	//
	 while ((arg = getopt(argc, argv, "b:c:i:o:p:r:t:u:h")) != -1) {
		//printf("User options for %c is %s\n", arg, optarg);
		switch(arg) {
			case 'b':
				if (!isNumber(optarg)) {
					fprintf (stderr, "Entered buffer size %s is not a number\n", optarg);
					usage();
					exit(1);
				}
				set_buffer_size(std::atol(optarg));
				break;
			case 'c':
				set_chrom_id(optarg);
				break;
			case 'i':
				input_file=optarg;
				//printf("%s\n", optarg);
				break;
			case 'o':
				set_output_file(optarg);
				break;
			case 'p':
				percentage=atoi(optarg);
				break;
			case 'r':
				reference_version=optarg;
				break;
			case 't':
				smoothing_type=std::atol(optarg);
				break;
			case 'u':
				set_upper_bound(std::atol(optarg));
				break;
			case '?':
				if (optopt == 'b' || optopt == 'c' || optopt == 'i' || optopt == 'o' 
						|| optopt == 'p' || optopt == 'r' || optopt == 't' || optopt == 'u')
					std::cerr << "Option " << optopt << " requires an argument." << std::endl;
				else if (isprint (optopt))
					std::cerr << "Unknown option " <<  optopt << std::endl;
				else
					std::cerr << "Unknown option character " << optopt << std::endl;
				usage();
				exit(1);
			default: 
				std::cerr << "Non-option argument " << optopt << std::endl;
				usage(); 
				exit(1);
		}
	}

	// check if the Mandatory options are available
	//
	check_inputs(input_file, "The input file is not available!");
	check_inputs(output_file, "Please specify the output file name!");
	check_inputs(chrom_id, "Please specify the chromosome ID to be drawn!");

	if (reference_version == "hg19")
		reference_version = "hg37";
}

void User_Inputs::check_inputs(std::string input_to_check, std::string error_message) {
	if (input_to_check.length() == 0) {
		std::cerr << error_message << std::endl;
		usage();
		exit(1);
	}
}

// This is used to check if a string (ie char *) is a number                                                  
//                                                                                                            
bool User_Inputs::isNumber(const char * inStr) {                                                                           
    while( *inStr != '\0') {                                                                                  
        if (!isdigit(inStr[0])) {                                                                             
            return false;                                                                                     
        }                                                                                                     
        inStr++;                                                                                              
    }                                                                                                         
    return true;                                                                                              
}
