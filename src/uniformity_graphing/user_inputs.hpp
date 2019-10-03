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
 *         Author:  Dr. Fritz Mehner (mn), mehner@fh-swf.de
 *        Company:  FH Südwestfalen, Iserlohn
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
		uint16_t get_upper_bound()    { return upper_bound; }
		uint16_t get_buffer_size()    { return buffer_size; }
		uint8_t  get_smoothing_type() { return smoothing_type; }
		uint8_t  get_percentage()     { return percentage; }

		// setters
		//
		void set_input_file(std::string in_file) { input_file = in_file; }
		void set_output_file(std::string o_file) { output_file = o_file; }
		void set_chrom_id(std::string cid) { chrom_id = cid; }
		void set_reference_version (std::string ref) { reference_version = ref; }
		void set_upper_bound(uint16_t upper)  { upper_bound = upper; }
		void set_buffer_size(uint16_t buffer) { buffer_size = buffer; }
		void set_smoothing_type(uint8_t type) { smoothing_type = type; }
		void set_percentage(uint8_t pct) { percentage = pct; }

		void usage();
		bool isNumber(const char * inStr);
		void processUserOptions(int argc, char *argv[]);
		void check_inputs(std::string input_to_check, std::string error_message);

	private:
		std::string input_file;
		std::string output_file;
		std::string chrom_id;			// the chromosome id to be drawn
		std::string reference_version;

		uint16_t upper_bound;		// upper bound used to generate the range file
		uint16_t buffer_size;		// the smoothing size around the upper_bound for type 2
		uint8_t smoothing_type;		// type 1: gVCF; type 2: range based
		uint8_t percentage;			// the percentage for the type 1 gVCF formula BLOCKAVE_'p'p
};

#endif
