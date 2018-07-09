/*
 * =====================================================================================
 *
 *       Filename:  exon.hpp
 *
 *    Description:  define the structure to hold the exon information
 *
 *        Version:  1.0
 *        Created:  08/30/2017 01:05:21 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef EXON_HPP
#define EXON_HPP

#include <string>

using namespace std;

// define a struct for holding the exon data from query result
//
struct Exon {
	private:
		uint32_t exon_start;
		uint32_t exon_end;
		uint16_t exon_id;
		//uint16_t exon_count;	we don't need this one as it is the size() of the exon vector
		uint32_t num_of_low_qual_bases;
		vector<string> *low_cov_regions;
		//enum Covered_By_Design {yes=1, no=0};

	public:
		Exon() { 
			exon_start=0; exon_end=0; exon_id=0; num_of_low_qual_bases=0; 
			low_cov_regions = new vector<string>;
			low_cov_regions->reserve(1);
		}
		~Exon() {}

		// Access functions
		//
		uint32_t get_exon_start() const { return exon_start; }
		uint32_t get_exon_end()   const { return exon_end; }
		uint16_t get_exon_id()    const { return exon_id;  }
		//uint16_t get_exon_count() const { return exon_count; }
		uint32_t get_num_of_low_qual_bases() const { return num_of_low_qual_bases; }
		vector<string> * get_low_cov_regions() const { return low_cov_regions; }

		// setters
		//
		void set_exon_start(uint32_t start) { exon_start = start; }
		void set_exon_end(uint32_t end)     { exon_end = end; }
		void set_exon_id (uint16_t e_id)    { exon_id  = e_id; }
		//void set_exon_count(uint16_t e_cnt) { exon_count = e_cnt; }
		void set_num_of_low_qual_bases(uint32_t num_of_bases) { num_of_low_qual_bases = num_of_bases; }
		void add_low_cov_regions(string low_cov_region)  { low_cov_regions->push_back(low_cov_region); }
};

#endif
