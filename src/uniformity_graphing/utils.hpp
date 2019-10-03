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

#include "user_inputs.hpp"
#include "coverage.hpp"

// Data structures
//
class Utils {
	public:
		// constructor and destruction!
		//
		Utils();
		~Utils() {}

		// accessors
		//
		uint32_t get_init_vector_size() { return init_vector_size; }
		uint32_t get_init_block_size()  { return init_block_size;  }
		std::map<std::string, uint32_t> get_chrom_length_map() { return chrom_length_map; }

		// setters
		//
		void set_init_vector_size(uint32_t v_size) { init_vector_size = v_size; }
		void set_init_block_size(uint32_t b_size)  { init_block_size  = b_size; }

		// there are 2 types of processInputFile() funciton
		// one for type 1: gVCF
		// another one for type 2: lower_bound - upper_bound range type
		// type 3 is to process input file without smoothing
		//
		void processInputFile1(User_Inputs *user_inputs, std::vector<Coverage> *cov_raw, std::vector<Coverage> *cov_sm);
		void processInputFile2(User_Inputs *user_inputs, std::vector<Coverage> *cov_raw, std::vector<Coverage> *cov_sm);
		void processInputFile3(User_Inputs *user_inputs, std::vector<Coverage> *cov_raw, std::vector<Coverage> *cov_sm);
		void combinedBlockInfo(User_Inputs *user_inputs, std::vector<Coverage> *sm, std::vector<Coverage> *block);
		void setupChromLengthMap(User_Inputs *user_inputs);

		// the following is a helper function that will help me with the debugging
		// convert char* to std::string object
		//
		std::string& charToString(const char* str);

	private:
		uint32_t init_vector_size;
		uint32_t init_block_size;
		std::map<std::string, uint32_t> chrom_length_map;
};

#endif
