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

#include "utils.hpp"

Utils::Utils() {
	init_vector_size=50000; 
	init_block_size=10;
}

void Utils::processInputFile1(User_Inputs *user_inputs, std::vector<Coverage> *cov_raw, std::vector<Coverage> *cov_sm) {
	// variables for raw coverage data
	//
	uint32_t start, stop, length, new_size;
	int32_t coverage;		// need to compare with negative value...
	std::string chrom_id;
	Coverage cov_storage;	// position Coverage array

	// check lowest_value_allowed, anything below 70 won't be smoothed!
	//
	int8_t lowest_value_allowed;
	lowest_value_allowed = user_inputs->get_upper_bound() / user_inputs->get_percentage();
	if (lowest_value_allowed < 70) lowest_value_allowed = 70;

	// check to see if we need to group nearby bases and coverages
	// 1: no, 2: yes
	//
	uint32_t flag=1;
	std::vector<Coverage> *block_to_be_smoothed = new std::vector<Coverage>;
	block_to_be_smoothed->reserve(10);
	int32_t lower_bound = user_inputs->get_upper_bound();
	int32_t upper_bound = user_inputs->get_upper_bound();

	// open file for reading
	//
	std::ifstream infile(user_inputs->get_input_file());

	// need to check for the comment lines and skip them
	//
	/*string line_c;
	while (getline(infile, line_c)) {
		if (line_c[0] == '#') {
			continue;
		} else {
			// process the first data line here
			//

			break;
		}
	}*/


	while (infile >> chrom_id >> start >> stop >> length >> coverage) {
		// process current input chromosome only
		//
		if (chrom_id != user_inputs->get_chrom_id())
			continue;

        cov_storage.set_chrom_id(chrom_id);
        cov_storage.set_start(start);
        cov_storage.set_stop(stop);
        cov_storage.set_length(length);
        cov_storage.set_cov(coverage);

        // need to dynamically increase the capacity of the vector
		//
        if (cov_raw->capacity() == cov_raw->size() + 2) {
            new_size = cov_raw->capacity() + get_init_vector_size();
            cov_raw->resize(new_size);
        }
		
        cov_raw->push_back(cov_storage);

		if (coverage > lowest_value_allowed) {
			// update lower_bound or upper_bound
			//
			if (coverage > upper_bound) 
				upper_bound = coverage;

			if (coverage < lower_bound)
				lower_bound = coverage;

			// check the gVCF formula constrain
			//
			if (lower_bound*user_inputs->get_percentage() >= upper_bound) {
				if (flag == 1) flag = 2;

				// check if we need to dynamically increase the block size
				//
				if (block_to_be_smoothed->capacity() == block_to_be_smoothed->size() + 2) {
					new_size = block_to_be_smoothed->capacity() + get_init_block_size();
					block_to_be_smoothed->resize(new_size);
				}
				block_to_be_smoothed->push_back(cov_storage);
			} else {
				// doesn't satify the BLOCKAVE_'p'p constrain
				// break and output block info if block exist
				//
				if (flag == 2) {
					// output previous block
					//
					combinedBlockInfo(user_inputs, cov_sm, block_to_be_smoothed);

					// clean-up
					//
					block_to_be_smoothed->clear();
					flag = 1;
				}

				// check if we need to increase the size for smoothed coverage vector                                     
				//
				if (cov_sm->capacity() == cov_sm->size() + 2) {                                                                   
					uint32_t new_size = cov_sm->capacity() + get_init_vector_size();                                          
					cov_sm->resize(new_size);                                                                                 
				}
				cov_sm->push_back(cov_storage);
				//cov_storage.print_coverage();
				lower_bound = user_inputs->get_upper_bound();
				upper_bound = user_inputs->get_upper_bound();
			}
		} else {
			// first check if any previous block exist
			//
			if (flag == 2) {
				// output previous block
				//
				combinedBlockInfo(user_inputs, cov_sm, block_to_be_smoothed);

				// clean-up
				//
				block_to_be_smoothed->clear();
				flag = 1;
			}

			// check if we need to increase the size for smoothed coverage vector
			//
			if (cov_sm->capacity() == cov_sm->size() + 2) {                                                                    
				uint32_t new_size = cov_sm->capacity() + get_init_vector_size();                                       
				cov_sm->resize(new_size);
			}
			cov_sm->push_back(cov_storage);
			//cov_storage.print_coverage();
			lower_bound = user_inputs->get_upper_bound();
			upper_bound = user_inputs->get_upper_bound();
		}
	}

	infile.close();
	delete block_to_be_smoothed;
}

void Utils::processInputFile2(User_Inputs *user_inputs, std::vector<Coverage> *cov_raw, std::vector<Coverage> *cov_sm) {
	// variables for raw coverage data
	//
	uint32_t start, stop, length, new_size;
	int32_t coverage;	// need to compare with negative value...
	std::string chrom_id;
	Coverage cov_storage;

	// when used to check to see if we need to group nearby bases and coverages 
	// 1: no, 2: yes
	// when used to check the first char each read-in line
	// -1: '#', 1: others
	//
	int8_t flag=-1;
	std::vector<Coverage> *block_to_be_smoothed = new std::vector<Coverage>;
	block_to_be_smoothed->reserve(10);
	int32_t lower_bound = user_inputs->get_upper_bound() - user_inputs->get_buffer_size();
	int32_t upper_bound = user_inputs->get_upper_bound() + user_inputs->get_buffer_size();

	// open file for reading
	//
	std::ifstream infile(user_inputs->get_input_file());
	if (!infile.is_open()) {
		std::cerr << "open file " << user_inputs->get_input_file() << " failed" << std::endl;
		exit(EXIT_FAILURE);
	}

	// as file start with comments, we need to skip them accordingly
	//
	char first;

	while (flag == -1) {
		// check the first character
		//
		first = infile.get();
		if (first == '#') {
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		} else {
			flag = 1;
		}
	}

	chrom_id = first;
	cov_storage.set_chrom_id(chrom_id);
	infile >> start >> stop >> length >> coverage;
	cov_storage.set_start(start);
	cov_storage.set_stop(stop);
	cov_storage.set_length(length);
	cov_storage.set_cov(coverage);

	while (infile >> chrom_id >> start >> stop >> length >> coverage) {
		// process current chrom_id only
		//
		if (chrom_id != user_inputs->get_chrom_id()) 
			continue;

		cov_storage.set_chrom_id(chrom_id);
		cov_storage.set_start(start);
		cov_storage.set_stop(stop);
		cov_storage.set_length(length);
		cov_storage.set_cov(coverage);

		// need to dynamically increase the capacity of the vector
		//
		if (cov_raw->capacity() == cov_raw->size() + 2) {
			new_size = cov_raw->capacity() + get_init_vector_size();
			cov_raw->resize(new_size);
		}

		cov_raw->push_back(cov_storage);

		// now we need to do smoothing
		//
		if ( (lower_bound <= coverage) && (coverage <= upper_bound) ) {
			if (flag == 1) flag = 2;

			// check if we need to dynamically increase the block size
			//
			if (block_to_be_smoothed->capacity() == block_to_be_smoothed->size() + 2) {
				new_size = block_to_be_smoothed->capacity() + get_init_block_size();
				block_to_be_smoothed->resize(new_size);
			}
			block_to_be_smoothed->push_back(cov_storage);
		} else {
			// reset the flag
			//
			if (flag == 2) {
				// output the previous block
				//
				combinedBlockInfo(user_inputs, cov_sm, block_to_be_smoothed);
				// clean-up
				//
				block_to_be_smoothed->clear();
				flag = 1;
			}

			// check if we need to increase the size for smoothed coverage vector
			//
			if (cov_sm->capacity() == cov_sm->size() + 2) {                                                           
				uint32_t new_size = cov_sm->capacity() + get_init_vector_size();                                  
				cov_sm->resize(new_size);
			}
			cov_sm->push_back(cov_storage);

			//cov_storage.print_coverage();
		}
	}

	infile.close();
	delete block_to_be_smoothed;
}

void Utils::processInputFile3(User_Inputs *user_inputs, std::vector<Coverage> *cov_raw, std::vector<Coverage> *cov_sm) {
	
}

void Utils::combinedBlockInfo(User_Inputs *user_inputs, std::vector<Coverage> *sm, std::vector<Coverage> *block) {
	uint32_t i=0, total_length=0, total_coverage=0;
	Coverage coverage;
	coverage.set_start(block->at(0).get_start());
	coverage.set_stop(block->at(block->size()-1).get_stop());
	coverage.set_chrom_id(block->at(0).get_chrom_id());

	for (i=0; i<block->size(); i++) {
		// the items in the block should be ordered by their start position,
		// so we will be OK to loop through them one by one
		//
		total_length   += block->at(i).get_length();
		total_coverage += block->at(i).get_cov() * block->at(i).get_length();
	}

	uint32_t average = uint32_t ((double)total_coverage / (double)total_length + 0.5);
	coverage.set_cov(average);
	coverage.set_length(total_length);

	// check if we need to increase the size for smoothed coverage vector
	//
	if (sm->capacity() == sm->size() + 2) {
		uint32_t new_size = sm->capacity() + get_init_vector_size();
		sm->resize(new_size);
	}

	sm->push_back(coverage);
	//coverage.print_coverage();
}

// used for debug only as char* need to be converted into std::string
//
std::string& Utils::charToString(const char* str) {
     return *(new std::string(str));
}

void Utils::setupChromLengthMap(User_Inputs *user_inputs) {
	// This is not initialization, I can NOT use { {"1", 248956422}, {"2", 242193529}, ...} notation
	//
	if (user_inputs->get_reference_version() == "hg38") {
		chrom_length_map["chr1"] = 248956422;
	   	chrom_length_map["chr2"] = 242193529;
	   	chrom_length_map["chr3"] = 198295559;
	   	chrom_length_map["chr4"] = 190214555;
		chrom_length_map["chr5"] = 181538259;
	   	chrom_length_map["chr6"] = 170805979;
		chrom_length_map["chr7"] = 159345973;
		chrom_length_map["chr8"] = 145138636;
		chrom_length_map["chr9"] = 138394717;
	   	chrom_length_map["chr10"] = 133797422; 
		chrom_length_map["chr11"] = 135086622;
		chrom_length_map["chr12"] = 133275309;
		chrom_length_map["chr13"] = 114364328;
		chrom_length_map["chr14"] = 107043718;
		chrom_length_map["chr15"] = 101991189;
		chrom_length_map["chr16"] = 90338345;
		chrom_length_map["chr17"] = 83257441; 
		chrom_length_map["chr18"] = 80373285;
		chrom_length_map["chr19"] = 58617616;
		chrom_length_map["chr20"] = 64444167;
		chrom_length_map["chr21"] = 46709983;
		chrom_length_map["chr22"] = 50818468;
		chrom_length_map["chrX"] = 156040895;
		chrom_length_map["chrY"] = 57227415;
		chrom_length_map["chrM"] = 16569;
	} else {
		chrom_length_map["1"] = 249250621;
		chrom_length_map["2"] = 243199373;
		chrom_length_map["3"] = 198022430;
		chrom_length_map["4"] = 191154276;
		chrom_length_map["5"] = 180915260;
		chrom_length_map["6"] = 171115067;
		chrom_length_map["7"] = 159138663;
		chrom_length_map["8"] = 146364022;
		chrom_length_map["9"] = 141213431;
		chrom_length_map["10"] = 135534747;
		chrom_length_map["11"] = 135006516;
		chrom_length_map["12"] = 133851895;
		chrom_length_map["13"] = 115169878;
		chrom_length_map["14"] = 107349540;
		chrom_length_map["15"] = 102531392;
		chrom_length_map["16"] = 90354753;
		chrom_length_map["17"] = 81195210;
		chrom_length_map["18"] = 78077248;
		chrom_length_map["19"] = 59128983;
		chrom_length_map["20"] = 63025520;
		chrom_length_map["21"] = 48129895;
		chrom_length_map["22"] = 51304566;
		chrom_length_map["X"] = 155270560;
		chrom_length_map["Y"] = 59373566;
		chrom_length_map["MT"] = 16569;
	}
}
