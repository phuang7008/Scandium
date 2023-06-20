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
	chromosome_length = 0;
	chromosome_centromere_file = "";
}

void User_Inputs::usage() {
	std::cout << "USAGE: draw_graph [options]" << std::endl;
	std::cout << "\t-i: the input file name [Mandatory]" << std::endl;
	std::cout << "\t-o: the output file name [Mandatory]" << std::endl;
	std::cout << "\t-c: the chromosome id [Mandatory]" << std::endl;
	std::cout << std::endl;
	std::cout << "\t-r: the reference version (for human)! (Default: hg37)" << std::endl;
	std::cout << "\t-m: the chromosome centromere region file (for human hg38 and above)!" << std::endl;
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
	 while ((arg = getopt(argc, argv, "c:i:o:m:r:h")) != -1) {
		//printf("User options for %c is %s\n", arg, optarg);
		switch(arg) {
			case 'c':
				chrom_id=optarg;
				break;
			case 'i':
				input_file=optarg;
				break;
			case 'm':
				chromosome_centromere_file=optarg;
				break;
			case 'o':
				//set_output_file(optarg);
				output_file=optarg;
				break;
			case 'r':
				reference_version=optarg;
				break;
			case '?':
				if (optopt == 'c' || optopt == 'i' || optopt == 'm' || optopt == 'o' || optopt == 'r')
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
	check_inputs(input_file, "The input file is not available!\n");
	check_inputs(output_file, "Please specify the output file name!\n");
	check_inputs(chrom_id, "Please specify the chromosome ID to be drawn!\n");

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
