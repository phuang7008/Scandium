/*
 * =====================================================================================
 *
 *       Filename:  utils.cpp
 *
 *    Description:  the detailed implementation of the header file
 *
 *        Version:  1.0
 *        Created:  08/23/2017 04:42:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030 
 *
 * =====================================================================================
 */
#include <fstream>
#include <iostream>
#include <iomanip>		// for std::setprecision
#include "utils.hpp"

Utils::Utils() { }

bool Utils::isFileEmpty(User_Inputs *user_inputs) {
    ifstream infile(user_inputs->get_low_cov_annotation_file());
    if (infile.peek() == std::ifstream::traits_type::eof()) {
        infile.close();
        return true;
    } else {
        infile.close();
        return false;
    }
}

bool Utils::isFile(User_Inputs *user_inputs) {
    struct stat strct;
    if ( stat(user_inputs->get_low_cov_annotation_file().c_str(), &strct) == 0 ) {
        if (strct.st_mode & S_IFDIR ) {
            // it is a directory
            //
            return false;
        } else if( strct.st_mode & S_IFREG ) {
            return true;
        } else {
            return false;
        }
    } else {
        cout << "check stats failed for file: " << user_inputs->get_low_cov_annotation_file()<< endl;
        return false;
    }
}

// Need to read in the target exons with low covorage regions here
//
void Utils::read_low_cov_annotation_file(User_Inputs *user_inputs, HTS_Data *hts_data) {
    // check if user provides the low coverage annotation file
    //
    if (user_inputs->get_low_cov_annotation_file() == "") return;

    // check if the input file is a directory
    //
    if (!isFile(user_inputs)) {
        cout << endl << "ERROR:" << endl;
        cout << "    the input low coverage annotation path is not a file: " << endl;
        cout << "    " << user_inputs->get_low_cov_annotation_file() << endl;
        cout << "    Please check the file and try again. Thanks!" << endl;
        cout << endl;
        return;
    }

    // check the file size
    //
    if (isFileEmpty(user_inputs)) {
        cout << "NOTE:" << endl;
        cout << "    The low coverage annotation file is empty:" << endl;
        cout << "    " << user_inputs->get_low_cov_annotation_file() << endl;
        cout << "    Return directly" << endl;
        cout << endl;
        return;
    }

	ifstream infile(user_inputs->get_low_cov_annotation_file());
	string chrom_id, exon_annotation;
	uint32_t exon_target_start, exon_target_end, low_cov_start, low_cov_end, size_low_cov;

	// Here is the example
	// chr11   2170753 2170770 chr11   2170689 2170770 TH|NM_199293|cds_1|gene 17
	//
	while (infile >> chrom_id >> low_cov_start >> low_cov_end >> chrom_id 
				  >> exon_target_start >> exon_target_end >> exon_annotation >> size_low_cov) {

		// need to skip if there is no annotations
		// chr8    27578719        27578724        .       -1      -1              .       0
		//
		if (exon_annotation == "." ) continue;

		// if the end of low cov region is > exon end, we need to cut the end of low cov region to be the same as exon_target end
		// chrY    1206432 1206599 chrY    1206432 1206445 CRLF2|NM_001012288-2_cds_1      13
		//
		if (low_cov_end > exon_target_end) 
			low_cov_end = exon_target_end;

		// if the start of low cov region < exon start, adjust the start of low cov region to be the same as exon target start
		// chrY	1288758	1288888	chrY	1288814	1288888	CSF2RA|NM_006140-2_cds_5	74
		//
		if (low_cov_start < exon_target_start)
			low_cov_start = exon_target_start;

		// if the chrom_id key is not found in the low_cov_exon_map iterator, create an empty container
		//
		if (hts_data->get_low_cov_exon_map()->find(chrom_id) == hts_data->get_low_cov_exon_map()->end()) {
			hts_data->get_low_cov_exon_map()->insert(std::make_pair(chrom_id, unordered_map<string, unordered_map<string, vector<Exon> > >()));
		}

		// Now process the exon annotation first. Here is an example:
		// TH|NM_199293|cds_1|gene 17
		// 
		vector<string> annotations = split(exon_annotation, '|');

		// Here we need to check to see if the current gene_symbol key is present. If not, we need to create an empty container
		//
		if (hts_data->get_low_cov_exon_map()->at(chrom_id).find(annotations[0]) == hts_data->get_low_cov_exon_map()->at(chrom_id).end()) {
			hts_data->get_low_cov_exon_map()->at(chrom_id).insert(std::make_pair(annotations[0], unordered_map<string, vector<Exon> >()));
			//cout << "Processing gene symbol " << annotations[0] << endl;
		}

		//cerr << annotations[0] << "\t" << annotations[2] << endl;
		vector<string> exon_id_info = split(annotations[2], '_');

		uint16_t exon_id = atoi(exon_id_info.back().c_str());
		exon_id_info.pop_back();		// pop the 'exon_id' out
		exon_id_info.pop_back();		// pop the 'exon' out

		string transcript_name(annotations[1]);

		// check to see if transcript_name key exists. If not, create an empty container
		//
		if (hts_data->get_low_cov_exon_map()->at(chrom_id).at(annotations[0]).find(transcript_name) == hts_data->get_low_cov_exon_map()->at(chrom_id).at(annotations[0]).end()) {
			hts_data->get_low_cov_exon_map()->at(chrom_id).at(annotations[0]).insert(std::make_pair(transcript_name, vector<Exon>()));

			// set the initial size
			// This is needed because according to the following:
			// When calling push_back() all iterators and references to vector's elements are invalidated 
			// if the new size() > capacity(). To avoid this, call reserve first.
			//
			hts_data->get_low_cov_exon_map()->at(chrom_id).at(annotations[0]).at(transcript_name).reserve(9);
		}

		// need to check if we have see this exon before to append or to add
		//
		bool found = false;
		char *lc_region = static_cast<char*> (calloc(50, sizeof(char)));

		// If the exon exist, we append the new low coverage regions to the exising one
		//
		for (auto en= hts_data->get_low_cov_exon_map()->at(chrom_id).at(annotations[0]).at(transcript_name).begin(); 
				  en!=hts_data->get_low_cov_exon_map()->at(chrom_id).at(annotations[0]).at(transcript_name).end(); ++en) {
			if (en->get_exon_start() == exon_target_start && en->get_exon_end() == exon_target_end) {
				// for debugging
				//
				//cerr << en->get_exon_start() << endl;
				//for (auto in=en->get_low_cov_regions()->begin(); in!=en->get_low_cov_regions()->end(); ++in) {
				//	cerr << *in << endl;
				//}
			
				// Found! Now need to update the low coverage information
				// 
				int32_t low_cov_count = en->get_num_of_low_qual_bases();
				low_cov_count += low_cov_end - low_cov_start;
				en->set_num_of_low_qual_bases(low_cov_count);

				sprintf(lc_region, "%d-%d", low_cov_start, low_cov_end);
				string lc_string(lc_region);
				en->add_low_cov_regions(lc_string);

				// for debugging
				//
				//for (auto in=en->get_low_cov_regions()->begin(); in!=en->get_low_cov_regions()->end(); ++in) {
				//	cerr << *in << endl;
				//}

				found = true;
			}
		}

		if (!found) {
			// The low_cov_exon is the low coverage regions
			//
			Exon low_cov_exon;
			low_cov_exon.set_exon_start(exon_target_start);
			low_cov_exon.set_exon_end(exon_target_end);
			low_cov_exon.set_num_of_low_qual_bases(low_cov_end - low_cov_start);
			low_cov_exon.set_exon_id(exon_id);

			sprintf(lc_region, "%d-%d", low_cov_start, low_cov_end);
			string lc_string(lc_region);
			low_cov_exon.add_low_cov_regions(lc_string);

			// insert the element into the right place
			//
			hts_data->get_low_cov_exon_map()->at(chrom_id).at(annotations[0]).at(transcript_name).push_back(low_cov_exon);
		}
		free(lc_region);
	}

	infile.close();
}

// if user specifies a list of genes interested in, we need to read them in here
// Note: gene names should be listed in the file, as one gene per line!
//
void Utils::store_gene_list(User_Inputs *user_inputs, HTS_Data *hts_data) {
	// return if user didn't specify the gene list
	//
	if (user_inputs->get_gene_list_file().length() == 0)
		return;

	ifstream infile(user_inputs->get_gene_list_file());
	string cur_gene;

	while (infile >> cur_gene) {
		// Insert the key-value here
		//
		(*hts_data->get_gene_list_to_check())[cur_gene] = cur_gene;
		//hts_data->get_gene_list_to_check()->insert({cur_gene, cur_gene});
	}
	infile.close();
}

// combines vector of strings together, separated by the delim char
//
string Utils::join( const std::vector<std::string>& s, const char* const delim) {
	switch (s.size())
    {
        case 0:
            return "";
        case 1:
            return s[0];
        default:
            std::ostringstream os; 
            std::copy(s.begin(), s.end()-1, std::ostream_iterator<std::string>(os, delim));
            os << *s.rbegin();
            return os.str();
    }
}

// used for debug only as char* need to be converted into std::string
//
std::string& Utils::charToString(const char* str) { 
	return *(new std::string(str)); 
}

/*
 * The following function is used to output the detailed annotation for gene/transcripts/exons respectively
 */
void Utils::output_detailed_annotations(string chrom_id, HTS_Data *hts_data, User_Inputs *user_inputs) {
	// open output files for writing
	//
	ofstream gene_pct_fs(user_inputs->get_output_low_cov_gene_file(), std::ios_base::app);
	ofstream exon_pct_fs(user_inputs->get_output_low_cov_exon_file(), std::ios_base::app);
	ofstream transcript_pct_fs(user_inputs->get_output_low_cov_transcript_file(), std::ios_base::app);

	for (auto git=hts_data->get_gene_transcript_exon_map()->begin(); 
				 git!=hts_data->get_gene_transcript_exon_map()->end(); ++git) {

		string gene_symbol=git->first;

		// Next check if gene list is present. If it is, process it accordingly. Otherwise, just skip it!
		// 
		//cerr << "The size of gene list is " << hts_data->get_gene_list_to_check()->size() << endl;

		if ((user_inputs->get_gene_list_file().length() > 0 ) || (hts_data->get_gene_list_to_check()->size() > 0)) {
			if (hts_data->get_gene_list_to_check()->find(gene_symbol) == hts_data->get_gene_list_to_check()->end())
				continue;
		}

		// output the gene related information: gene_symbol first
		//
		gene_pct_fs << gene_symbol << "\t";
		float total_transcript_cov_rate=0.0;
		uint16_t num_of_transcripts=0;

		for (auto it=hts_data->get_gene_transcript_exon_map()->at(gene_symbol).begin(); 
				 it!=hts_data->get_gene_transcript_exon_map()->at(gene_symbol).end(); ++it) {

			string transcript_name=it->first;

			// NOTE: some gene doesn't have a gene symbol or transcript name, therefore, we will skip them
			//
			if (transcript_name.compare(".") == 0) continue;

			vector<Exon> tmp_exons = it->second;	
			uint32_t total_exon_size=0;
			uint32_t total_exon_low_qual_bases=0;

			for (auto eit=tmp_exons.begin(); eit!=tmp_exons.end(); ++eit) {
				// handle current exon
				//
				uint32_t exon_length  = eit->get_exon_end() - eit->get_exon_start();
				float exon_cov_rate = ((float) exon_length - eit->get_num_of_low_qual_bases())*100 / (float) exon_length;

				exon_pct_fs << chrom_id << "\t" << gene_symbol << "\t" << transcript_name << "\texon_" << eit->get_exon_id() << "\t";
				exon_pct_fs << eit->get_exon_start() << "\t" << eit->get_exon_end() << "\t" << std::fixed << std::setprecision(3) << exon_cov_rate << "\t";

				// for low coverage region output
				//
				vector<string> *low_cov_regions = eit->get_low_cov_regions();

				if (low_cov_regions->size() == 0) {
					exon_pct_fs << ".";
				} else {
					bool first_print=true;
					for (auto lt=low_cov_regions->begin(); lt!=low_cov_regions->end(); ++lt) {
						if (!first_print) {
							exon_pct_fs << ";";
						} else {
							first_print=false;
						}

						exon_pct_fs << *lt;
					}
				}
				exon_pct_fs << endl;

				total_exon_size += exon_length;
				total_exon_low_qual_bases += eit->get_num_of_low_qual_bases();
			}

			// handle current transcript
			//
			++num_of_transcripts;
			float transcript_cov_rate = ((float) total_exon_size - total_exon_low_qual_bases)*100 / (float) total_exon_size;
			total_transcript_cov_rate += transcript_cov_rate;
			if (num_of_transcripts > 1)
				gene_pct_fs << ";";

			gene_pct_fs << transcript_name << "(" ;
			gene_pct_fs << std::fixed << std::setprecision(3) <<transcript_cov_rate << ")";
			
			// now output the transcript related information
			//
			transcript_pct_fs << chrom_id << "\t" << gene_symbol << "\t" << transcript_name << "\t" << total_exon_size << "\t" ;
			//transcript_pct_fs << (*hts_data->get_exon_count_per_transcript_map())[transcript_name] << "\t";
			transcript_pct_fs << tmp_exons.size() << "\t";
			transcript_pct_fs << std::fixed << std::setprecision(3) << transcript_cov_rate << endl;

			// now erase it
			//
			//it = hts_data->get_gene_transcript_exon_map()->at(gene_symbol).erase(it);

		}

		hts_data->get_gene_transcript_exon_map()->at(gene_symbol).clear();

		float ave_transcript_cov_rate = total_transcript_cov_rate / (float) num_of_transcripts;
		gene_pct_fs << std::fixed << std::setprecision(2) << "\t" << ave_transcript_cov_rate << endl;

		// erase the key and value
		//
		//git = hts_data->get_gene_transcript_exon_map()->erase(git);
	}

	hts_data->get_gene_transcript_exon_map()->clear();
}
