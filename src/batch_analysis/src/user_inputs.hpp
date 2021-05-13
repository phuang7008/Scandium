/*
 * =====================================================================================
 *
 *       Filename:  user_inputs.hpp
 *
 *    Description:  get user inputs and process user options
 *
 *        Version:  1.0
 *        Created:  08/23/2017 11:36:54 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Mecidine
 *
 * =====================================================================================
 */
#ifndef USER_INPUTS_HPP
#define USER_INPUTS_HPP

#include <string>

using namespace std;

class User_Inputs {
	public:
		User_Inputs();
		~User_Inputs() {}

		// access functions
		//
		string get_reference_version() const { return reference_version; }
		string get_output_dir() const { return output_dir; }
		string get_gene_list_file() const { return gene_list_file; }
		string get_exon_bed_file() const { return output_exon_bed_file; }
		string get_low_cov_annotation_file() const { return low_cov_annotation_file; }
		string get_designed_targets_bed_file() const { return designed_targets_bed_file; }
		string get_design_uncovered_exon_file() const { return design_uncovered_exon_file; }
		string get_user_defined_annotation_file() const { return user_defined_annotation_file; }

		// output file list
		//
		string get_output_low_cov_gene_file() const { return output_low_cov_gene_file; }
		string get_output_low_cov_exon_file() const { return output_low_cov_exon_file; }
		string get_output_low_cov_transcript_file() const { return output_low_cov_transcript_file; }

		void usage();
		void process_user_inputs(int argc, char* argv[]);
		void check_inputs(string input_to_check, string error_message);
		void remove_exising_output_file(const char* file_name);
		void output_user_options();

	private:
		string reference_version;
		string output_dir;
		string gene_list_file;					// list of genes that users intested in
		string design_uncovered_exon_file;		// the text file that contains exons which are not covered by capture design
		string designed_targets_bed_file;		// the text file that contains regions which are covered by capture design
		string low_cov_annotation_file;			// the text file that contains the low cov regions with detailed annotations
		string user_defined_annotation_file;	// the text file that contains the detailed user-defined annotations
		string output_exon_bed_file;			// exon bed file used to intersect the low coverage regions
		string output_low_cov_gene_file;
		string output_low_cov_exon_file;
		string output_low_cov_transcript_file;
};

#endif
