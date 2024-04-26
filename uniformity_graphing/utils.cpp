/*
 * =====================================================================================
 *
 *       Filename:  utils.cpp
 *
 *    Description:  The detailed implementation for utils.hpp
 *
 *        Version:  1.0
 *        Created:  01/04/2018 10:39:16 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang. phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77036
 *
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <limits>   // define numeric_limits
#include <cstdlib>
#include <sstream>

#include "utils.hpp"

Utils::Utils() {
	init_vector_size=50000; 
	init_block_size=10;
}

void Utils::processInputFile(User_Inputs *user_inputs, std::vector<Coverage> *cov_raw) {
	// variables for raw coverage data
	//
	unsigned int start, stop, length, new_size;
	uint32_t coverage;		// need to compare with negative value...
	std::string chrom_id;
	unsigned int chr_length = 0;

	// open file for reading
	//
	std::ifstream infile(user_inputs->get_input_file());

	// need to check for the comment lines and skip them
	//
	std::string line_c;
	while (std::getline(infile, line_c)) {
		if (line_c[0] == '#') {
			continue;
		} else {
			// process the first data line here
			// Insert the string into a stream
			//
			std::stringstream ss(line_c);
			std::string buf;				// Have a buffer string
			int i = 0;
			bool to_skip=false;

			while (ss >> buf) {
				if (i == 0) chrom_id = buf;

				// only process regions that the user interested in
				//
				if (chrom_id != user_inputs->get_chrom_id()) {
					to_skip = true;
					break;
				}

				if (i == 1) start = (unsigned int) std::stol(buf);
				if (i == 2) stop  = (unsigned int) std::stol(buf);
				if (i == 3) {
					length = (unsigned int) std::stol(buf);
					chr_length += (unsigned int) std::stol(buf);
				}
				if (i == 4) coverage = (uint32_t) std::stol(buf);

                if (stop > chr_length) chr_length = stop;

				i++;
			}

			if (to_skip) continue;

			// need to exclude the centromere regions only if users specify the centromere file
			//
			if (user_inputs->get_chromosome_centromere_file().length() > 0) {
				bool overlap = false;
				unsigned int max_start=0, min_end=0;	// the overlap region

				// Here keyword 'auto' stands for "automatic type inference" based on the initialization statement 
				//
				for (auto it=centromere_regions->begin(); it!=centromere_regions->end(); ++it) {
					/*
					 * centromere_begin (it->first) ---------- centromere end (it->second)
					 *             current_region_start (start) +++++++++++ current_region_end (stop)
					 *
					 *					centromere_begin (it->first) ---------- centromere end (it->second)
					 * current_region_start (start) ++++++++++++ current_region_end (stop)
					 *
					 * No overlap! So continue everything as is!
					 */
					if (it->first > stop || start > it->second) {
						continue;
					} else {
						overlap = true;
						(it->first >= start) ?  max_start = it->first : max_start = start;                                      
						(it->second >= stop) ?  min_end = stop : min_end = it->second;
					}
				}

				check_for_centromere_regions(start, stop, length, coverage, max_start, min_end, cov_raw, overlap);
			} else {
				Coverage cov_storage;	// position Coverage array

				// don't need to exclude the centromere regions
				//
				//cov_storage.set_chrom_id(chrom_id);
				cov_storage.set_start(start);
				cov_storage.set_stop(stop);
				cov_storage.set_length(length);
				cov_storage.set_cov(coverage);

				cov_raw->push_back(cov_storage);
			}
		}

        // need to dynamically increase the capacity of the vector
		//
        if (cov_raw->capacity() <= cov_raw->size() + 5) {
            new_size = cov_raw->capacity() + get_init_vector_size();
            cov_raw->resize(new_size);
        }
	}

	infile.close();
	user_inputs->set_chromosome_length(chr_length);
}

void Utils::check_for_centromere_regions(unsigned int start, unsigned int stop, unsigned int length, uint32_t coverage, unsigned int max_start, unsigned int min_end, std::vector<Coverage> *cov_raw, bool overlap) {

	Coverage cov_storage;   // position Coverage array

	if (!overlap) {
		cov_storage.set_start(start);
		cov_storage.set_stop(stop);
		cov_storage.set_length(length);
		cov_storage.set_cov(coverage);
		cov_raw->push_back(cov_storage);

	} else {
		/* 
		 * they overlap.
		 */

		/*
		 * the max_start and min_end will divide the region into three parts
		 *
		 *		start -------- max_start ------------- min_end ------ stop
		 *		                         <- overlap ->
		 */
		if (start < max_start) {
			cov_storage.set_start(start);
			cov_storage.set_stop(max_start - 1);
			cov_storage.set_length(max_start - 1 - start);
			cov_storage.set_cov(coverage);
			cov_raw->push_back(cov_storage);
		}

		// for overlaped region
		//
		cov_storage.set_start(max_start);
		cov_storage.set_stop(min_end);
		cov_storage.set_length(min_end - max_start);
		cov_storage.set_cov(0);
		cov_raw->push_back(cov_storage);

		if (stop > min_end) {
			cov_storage.set_start(min_end + 1);
			cov_storage.set_stop(stop);
			cov_storage.set_length(min_end + 1 - stop);
			cov_storage.set_cov(coverage);
			cov_raw->push_back(cov_storage);
		}
	}
}

// used for debug only as char* need to be converted into std::string
//
std::string& Utils::charToString(const char* str) {
     return *(new std::string(str));
}

void Utils::read_in_centromere_regions(User_Inputs *user_inputs) {
	centromere_regions = new std::map<unsigned int, unsigned int>;

	// open file for reading
	//
	std::ifstream infile(user_inputs->get_chromosome_centromere_file());	

	std::string line_c;
	while (std::getline(infile, line_c)) {
		// need to check for the comment lines and skip them
		//
		if (line_c[0] == '#') {
			continue;
		} else {
			std::stringstream ss(line_c);
			std::string buf, chrom_id;
			unsigned int start, stop;
			int i = 0;
			
			while (ss >> buf) {
				if (i == 0) chrom_id = buf;
				if (i == 1) start = std::stol(buf);
				if (i == 2) stop  = std::stol(buf);
				i++;
			}

			// skip those regions that are not part of the same chromosome users interested in
			//
			if (chrom_id != user_inputs->get_chrom_id()) continue;

			// now insert start and stop information to the map<unsigned int, unsigned int>
			//
			centromere_regions->insert(std::make_pair(start, stop));
		}
	}
}
