/*
 * =====================================================================================
 *
 *       Filename:  utils.hpp
 *
 *    Description:  This is used for graphing in C++
 *
 *        Version:  1.0
 *        Created:  01/04/2018 10:39:16 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef UTILS_HPP
#define UTILS_HPP

#include <inttypes.h>   // for PRIu32 and PRIu64 
#include <unistd.h>     // for file access() and getopt()
#include <stdbool.h>    // for bool definition

#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>

#include "user_inputs.hpp"
#include "coverage.hpp"

// Data structures
//
class Utils {
	public:
		// constructor and destruction!
		//
		Utils();
		~Utils() { delete centromere_regions; }

		// accessors
		//
		uint32_t get_init_vector_size() { return init_vector_size; }
		uint32_t get_init_block_size()  { return init_block_size;  }

		// setters
		//
		void set_init_vector_size(uint32_t v_size) { init_vector_size = v_size; }
		void set_init_block_size(uint32_t b_size)  { init_block_size  = b_size; }

		// methods
		//
		void processInputFile(User_Inputs *user_inputs, std::vector<Coverage> *cov_raw);
		void read_in_centromere_regions(User_Inputs *user_inputs);
		void check_for_centromere_regions(unsigned int start, unsigned int stop, unsigned int length, uint32_t coverage, unsigned int max_start, unsigned int min_end, std::vector<Coverage> *cov_raw, bool overalp);

		// the following is a helper function that will help me with the debugging
		// convert char* to std::string object
		//
		std::string& charToString(const char* str);

	private:
		uint32_t init_vector_size;
		uint32_t init_block_size;
		std::map<unsigned int, unsigned int > *centromere_regions;
};

#endif
