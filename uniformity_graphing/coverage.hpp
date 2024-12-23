/*
 * =====================================================================================
 *
 *       Filename:  coverage.hpp
 *
 *    Description:  the structure/class for sequencing coverage
 *
 *        Version:  1.0
 *        Created:  01/05/2018 03:30:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef COVERAGE_HPP
#define COVERAGE_HPP

#include <string>
#include <iostream>

class Coverage {
	public:
		Coverage() { start=0; stop=0; length=0; cov=0; }
		~Coverage() {}

		// accessors
		//
		uint64_t get_start()  { return start; }
		uint64_t get_stop()   { return stop; }
		uint64_t get_length() { return length; }
		uint32_t get_cov()    { return cov; }
		//std::string get_chrom_id() { return chrom_id; }

		// setters
		//
		void set_start(uint64_t st) { start = st; }
		void set_stop(uint64_t sp)  { stop = sp; }
		void set_length(uint64_t len) { length = len; }
		void set_cov(uint32_t cv) { cov = cv; }
		//void set_chrom_id(std::string cid) { chrom_id = cid; }

		void print_coverage() {
			std::cout << get_start() << "\t" << get_stop() << "\t";
			std::cout << get_length() << "\t" << get_cov() << std::endl;
		}

	private:
		uint64_t start;
		uint64_t stop;
		uint64_t length;
		uint32_t cov;
		//std::string chrom_id;
};

#endif
