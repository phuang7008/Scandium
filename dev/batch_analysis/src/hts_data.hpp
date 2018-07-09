/*
 * =====================================================================================
 *
 *       Filename:  hts_data.hpp
 *
 *    Description:  To declare all the major data and data structure used by the program
 *
 *        Version:  1.0
 *        Created:  08/25/2017 09:46:53 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef HTS_DATA_HPP
#define HTS_DATA_HPP

#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

#include "exon.hpp"

class HTS_Data {
	public:
		HTS_Data() : initial_vector_size(30) { 
			gene_list_to_check = new unordered_map<string, string>;
			gene_transcript_exon_map = new unordered_map<string, unordered_map<string, vector<Exon::Exon> > >;
			target_annotation_map = new unordered_map<string, unordered_map<string, unordered_map<string, vector<Exon::Exon> > > >;
			low_cov_exon_map = new unordered_map<string, unordered_map<string, unordered_map<string, vector<Exon::Exon> > > >;
			chromosome_list = new unordered_map<string, int>;
		}

		~HTS_Data() { 
			delete gene_list_to_check; 
			delete low_cov_exon_map; 
			delete target_annotation_map; 
			delete gene_transcript_exon_map;
			delete chromosome_list;
		}

		// access functions
		//
		uint32_t get_initial_vector_size () const { return initial_vector_size; }
		unordered_map<string, int> * get_chromosome_list() const { return chromosome_list; }

		unordered_map<string, string> * get_gene_list_to_check() const { return gene_list_to_check; }

		unordered_map<string, unordered_map<string, unordered_map<string, vector<Exon::Exon> > > > * get_low_cov_exon_map() const { return low_cov_exon_map; }
		unordered_map<string, unordered_map<string, unordered_map<string, vector<Exon::Exon> > > > * get_target_annotation_map() const { return target_annotation_map; }
		unordered_map<string, unordered_map<string, vector<Exon::Exon> > > * get_gene_transcript_exon_map() const { return gene_transcript_exon_map; }

		//unordered_map<string, int> *get_exon_count_per_transcript_map() const { return exon_count_per_transcript_map; }

		void add_to_chromosome_list(string chrom_id) { (*chromosome_list)[chrom_id] = 1; }

		void reset_gene_transcript_exon_map() { 
			delete gene_transcript_exon_map; 
			gene_transcript_exon_map = NULL;
			gene_transcript_exon_map = new unordered_map<string, unordered_map<string, vector<Exon::Exon> > >; 
		}

	private:
		const uint32_t initial_vector_size;
		unordered_map<string, int> *chromosome_list;
		unordered_map<string, string> * gene_list_to_check;

		// Store MySQL query results for Gene Exon Annotation information (official gene exon annotation)
		// NOTE: the first key is gene_symbol, while the second key is transcript name
		// The following is the schema:
		// gene_transcript_exon_map -----> gene_symbol key -----> transcript_name1 key -----> vector<Exon::Exon> -----> exon1
		//                                                                                                       -----> exon2
		//                                                                                                       -----> ...
		//                                                                                                       -----> exonn
		//                                                 -----> transcript_name2 key
		//                                                 -----> transcript_name3 key
		//
		unordered_map<string, unordered_map<string, vector<Exon::Exon> > > * gene_transcript_exon_map;

		// Store designed targets. The key info is
		// target_map -----> chrom_id key -----> gene_symbol key -----> transcript_name1 key -----> vector<Exon::Exon> -----> exon1
		//																											   -----> exon2
		//																											   -----> ...
		//
		unordered_map<string, unordered_map<string, unordered_map<string, vector<Exon::Exon> > > > * target_annotation_map;

		// Store the information from input file regarding low coverage regions
		// NOTE: the first key is chrom_id, second key is gene_symbol, while the third key is the transcript name 
		// The following is the schema:
		// low_cov_exon_map -----> chrom_id key -----> gene_symbol1 key (same as the schema above)
		//                                      -----> gene_symbol2 key (same as the schema above)
		//                                      -----> ...
		//
		unordered_map<string, unordered_map<string, unordered_map<string, vector<Exon::Exon> > > > * low_cov_exon_map;

};

#endif
