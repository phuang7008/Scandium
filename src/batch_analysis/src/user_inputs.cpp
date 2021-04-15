/*
 * =====================================================================================
 *
 *       Filename:  user_inputs.cpp
 *
 *    Description:  the detailed inplementation of user input options
 *
 *        Version:  1.0
 *        Created:  08/23/2017 01:32:03 PM
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
	mode = 1;
	annotation_source = 1;
	reference_version = "";
	output_dir = "";
	gene_list_file = "";
	designed_targets_bed_file = "";
	user_defined_annotation_file = "";
	low_cov_annotation_file  = "";
	output_low_cov_gene_file = "";
	output_low_cov_exon_file = "";
	output_exon_bed_file  = "";
	output_low_cov_transcript_file = "";
}

void User_Inputs::usage() {
	cout << "USAGE: analysis [options]" << endl;
	cout << "\t-d: database version (either hg19 or hg38) [Mandatory]" << endl;
	cout << "\t-i: input file that contains the low coverage regions in bed format [Mandatory]" << endl;
	cout << "\t-o: output directory [Mandatory]" << endl;
	cout << "\t-t: designed targets in bed file format [Mandatory]" << endl;
	cout << endl;
	cout << "\t-f: user-defined annotation database. If no -f option, no annotations will be provided!" << endl;
	cout << "\t-g: a file contains list of genes of interest. If no -g option, all genes with low coverage will be processed!" << endl;
}

// Get command line arguments in and check the sanity of user inputs 
//
void User_Inputs::process_user_inputs(int argc, char* argv[]) {
	int arg;

	//When getopt returns -1, no more options available
	while ((arg = getopt(argc, argv, "a:d:f:g:hi:m:o:t:u:")) != -1) {
		//printf("User options for %c is %s\n", arg, optarg);
		switch(arg) {
			case 'a':
				annotation_source=atoi(optarg);
				break;
			case 'd':
				reference_version=optarg;
				break;
			case 'f':
				user_defined_annotation_file=optarg;
				break;
			case 'g':
				gene_list_file=optarg;
				break;
			case 'i':
				low_cov_annotation_file=optarg;
			case 'm':
				mode=atoi(optarg);
				break;
			case 'o':
				output_dir=optarg;
				break;
			case 't':
				designed_targets_bed_file=optarg;
				break;
			case 'u':
				design_uncovered_exon_file=optarg;
				break;
			case '?':
				if (optopt == 'd' || optopt == 'f' || optopt == 'i' || optopt == 'm' || optopt == 'g' || optopt == 'o' || optopt == 'u')
					cerr << "Option " << optopt << " requires an argument." << endl;
				else if (isprint (optopt))
					cerr << "Unknown option " <<  optopt << endl;
				else
					cerr << "Unknown option character " << optopt << endl;
				usage();
				exit(1);
			default: 
				cerr << "Non-option argument " << optopt << endl; 
				usage(); 
				exit(1);
		}
	}

	// check if the Mandatory options are available
	//
	check_inputs(low_cov_annotation_file, "Low Coverage Regions with Annotatioin File is missing!");
	check_inputs(output_dir, "Output Directory is not specified!");
	check_inputs(reference_version, "Please specify the Database Version you are using!");
	check_inputs(designed_targets_bed_file, "Please specify the Designed Targets in bed format!");

	// now we need to generate the other file names based on the input file name
	// but first get the base file name first without the path information
	//
	// string base_filename = input_file.substr(input_file.find_last_of("\\/") + 1);
	string base_filename = low_cov_annotation_file;
	const size_t last_slash_idx = base_filename.find_last_of("\\/");
	if (string::npos != last_slash_idx) {
		base_filename.erase(0, last_slash_idx + 1);
	}

	output_low_cov_gene_file = output_dir + "/" + base_filename + "_" + "Gene_File.txt";
	output_low_cov_exon_file = output_dir + "/" + base_filename + "_" + "Exon_File.txt";
	output_low_cov_transcript_file = output_dir + "/" + base_filename + "_" + "Transcript_File.txt";
	output_exon_bed_file = output_dir + "/" + "Official_Exons.bed";
	//output_exon_bed_file = output_dir + "/" + "Official_Merged_Exons.bed";

	// delete these files if they exist, as the mode for writing is append
	// Therefore, we don't want to have files append forever
	//
	remove_exising_output_file(output_low_cov_gene_file.c_str());
	remove_exising_output_file(output_low_cov_exon_file.c_str());
	remove_exising_output_file(output_low_cov_transcript_file.c_str());
}

void User_Inputs::check_inputs(string input_to_check, string error_message) {
	if (input_to_check.length() == 0) {
		cerr << error_message << endl;
		usage();
		exit(1);
	}
}

void User_Inputs::remove_exising_output_file(const char* file_name) {
	if ( access(file_name, F_OK) != -1) {
		remove(file_name);
	}
}

void User_Inputs::output_user_options() {
	cerr << endl;
	cerr << "The following are the options you have choosen:" << endl;
	cerr << "\tOutput directory: " << output_dir << endl;
	cerr << "\tDatabase version used: " << reference_version << endl;
	cerr << "\tDesigned targeted bed file: " << designed_targets_bed_file << endl;
	cerr << "\tLow coverage annotation file: " << low_cov_annotation_file << endl;
	cerr << "\tUser-Defined annotation file: " << user_defined_annotation_file << endl;
	
	if (get_design_uncovered_exon_file().length() > 0) {
		cerr << "\tDesigned un-covered exon file: " << design_uncovered_exon_file << endl;
	}

	if (gene_list_file.length() > 0) {
		cerr << "\tGene list file: " << gene_list_file << endl;
	}

	if (annotation_source == 1) {
		cerr << "\tAnnotation Source: RefSeq" << endl;
	} else if (annotation_source == 2) {
		cerr << "\tAnnotation Source: CCDS" << endl;
	} else if (annotation_source == 3) {
		cerr << "\tAnnotation Source: VEGA" << endl;
	} else if (annotation_source == 4) {
		cerr << "\tAnnotation Source: Gencode" << endl;
	} else {
		cerr << "\tAnnotation Source: All (RefSeq + CCDS + VEGA(hg19) or Gencode(hg38) " << endl;
	}

	cerr << endl;
	//cerr << "" << << endl;
}
