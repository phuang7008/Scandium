/*
 * =====================================================================================
 *
 *       Filename:  user_inputs.hpp
 *
 *    Description:  fetch user inputs
 *
 *        Version:  1.0
 *        Created:  01/05/2018 01:57:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef USER_INPUTS_HPP
#define USER_INPUTS_HPP

#include <string>

class User_Inputs {

	public:
		User_Inputs();
		~User_Inputs() {}

		// access function
		//
		std::string get_input_file()  { return input_file; }
		std::string get_output_file() { return output_file; }
		std::string get_chrom_id()    { return chrom_id; }
		std::string get_reference_version() { return reference_version; }
		std::uint64_t get_chromosome_length() { return chromosome_length;}
		std::string get_chromosome_centromere_file() { return chromosome_centromere_file; }

		// setters
		//
		void set_input_file(std::string in_file) { input_file = in_file; }
		void set_output_file(std::string o_file) { output_file = o_file; }
		void set_chrom_id(std::string cid) { chrom_id = cid; }
		void set_reference_version (std::string ref) { reference_version = ref; }
		void set_chromosome_length (uint64_t chr_length) { chromosome_length = chr_length; }
		void set_chromosome_centromere_file (std::string chr_centromere_file) { chromosome_centromere_file = chr_centromere_file; }

		void usage();
		bool isNumber(const char * inStr);
		void processUserOptions(int argc, char *argv[]);
		void check_inputs(std::string input_to_check, std::string error_message);

	private:
		std::string input_file;
		std::string output_file;
		std::string chrom_id;			// the chromosome id to be drawn
		std::string reference_version;
		std::string chromosome_centromere_file;
		uint64_t chromosome_length;
};

#endif
