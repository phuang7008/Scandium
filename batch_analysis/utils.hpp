/*
 * =====================================================================================
 *
 *       Filename:  utils.hpp
 *
 *    Description:  general utility class to process user inputs, check errors etc.
 *
 *        Version:  1.0
 *        Created:  08/23/2017 11:32:11 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */
#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <sstream>		// used for string split
#include <vector>
#include <iterator>		// used for string split
#include <unordered_map>
#include <sys/stat.h>   // used for checking file stats call stat()

using namespace std;

#include "hts_data.hpp"
#include "user_inputs.hpp"

class Utils {
	public:
		Utils();
		~Utils() {}

		// the following is a helper function that will help me with the debugging
		// convert char* to std::string object
		//
		std::string& charToString(const char* str);

		// for string split
		//
		template<typename T>
		void split(const string &s, char delim, T result) {
			stringstream ss;
			ss.str(s);
			string item;
			while (getline(ss, item, delim)) {
				*(result++) = item;
			}
		}

        // the following is also a help function
        //
        bool isFileEmpty(User_Inputs *user_inputs);

        // the following is used to check if a file path is a file (true) or directory (false)
        //
        bool isFile(User_Inputs *user_inputs);

		std::vector<string> split(const string &s, char delim) {
			vector<string> elems;
			split(s, delim, back_inserter(elems));
			return elems;
		}

		// To join together strings within a vector separated by delim
		//
		string join( const std::vector<std::string>& s, const char* const delim);

		// Read in the detailed information for the low coverage regions, with annotations
		//
		void read_low_cov_annotation_file(User_Inputs *user_inputs, HTS_Data *hts_data);

		// Read in Designed Targets
		//
		void read_target_bed_file(User_Inputs *user_inputs, HTS_Data *hts_data);

		// Read in the genes of interest
		//
		void store_gene_list(User_Inputs *user_inputs, HTS_Data *hts_data);

		// Output the detailed genes, transcripts and exons information to their corresponding files
		//
		void output_detailed_annotations(string chrom_id, HTS_Data *hts_data, User_Inputs *user_inputs);

	//private:

};

#endif
