/*
 * =====================================================================================
 *
 *       Filename:  analysis.cpp
 *
 *    Description:  for batch analysis of genes with low coverage across multiple runs
 *
 *        Version:  1.0
 *        Created:  08/23/2017 11:20:38 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Peiming (Peter) Huang, mehner@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */

#include<iostream>
#include<fstream>
#include<string>

#include "hts_data.hpp"
#include "user_inputs.hpp"
#include "utils.hpp"
#include "annotation.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	cerr << endl << "Start Batch Analysis ..." << endl;

	// process user inputs
	User_Inputs user_inputs=User_Inputs();
	user_inputs.process_user_inputs(argc, argv);
	user_inputs.output_user_options();

	Annotation annotation=Annotation(user_inputs.get_db_version());
	if (user_inputs.get_mode() == 1) {
		annotation.fetch_and_dump_exon_info(&user_inputs);
	} else {
		// declare HTS_Data object on the heap as the data wil grow!
		//
		HTS_Data *hts_data = new HTS_Data();

		// process low coverage annotation file and genes to check list
		//
		Utils utils=Utils();
		utils.read_low_cov_annotation_file(&user_inputs, hts_data);
		utils.store_gene_list(&user_inputs, hts_data);

		// get the low_cov_exon_map, which we just read in
		//
		unordered_map <string, unordered_map<string, unordered_map<string, vector<Exon::Exon> > > > *low_cov_exon_map = hts_data->get_low_cov_exon_map();
		//cout << low_cov_exon_map->size() << endl;

		// Here keyword 'auto' stands for "automatic type inference" based on the initialization statement
		// the following is just for debugging
		//
		for (auto it=low_cov_exon_map->begin(); it != low_cov_exon_map->end(); ++it) {
			// iterator like pointers, we need to use pointer dereference to access them!
			//
			//cout << "size for " << it->first << " : " << it->second.size() << endl;
			continue;
		}

		// need to find out the chrom_list, either from the database or from the user_defined database
		//
		if (user_inputs.get_user_defined_annotation_file() == "") {
			annotation.get_chrom_list_from_mysql(hts_data, &user_inputs);
		} else {
			annotation.get_chrom_list_from_user_defined_db(hts_data, &user_inputs);
		}

		// loop through chrom_list, as we are going to process one chromosome at a time!
		//
		unordered_map<string, int>* chrom_list = hts_data->get_chromosome_list();
		for (auto it=chrom_list->begin(); it!=chrom_list->end(); ++it) {
			//cout << "Current chromosome is " << it->first << endl;

			if (user_inputs.get_user_defined_annotation_file() == "") {
				// fetch the official gene exon information from MySQL database and update low coverage regions accordingly
				//
				//fprintf(stderr, "Use MySQL to fetch gene annotations\n");
				annotation.fetch_gene_exon_info_from_mysql(it->first, hts_data, &user_inputs);
			} else {
				// fetch the gene exon information from user-defined database
				//
				//fprintf(stderr, "Use user-defined database for gene annotations\n");
				annotation.fetch_gene_exon_info_from_user_defined_db(it->first, hts_data, &user_inputs, &utils);
			}

			// output the detailed gene, transcript and exon information
			//
			utils.output_detailed_annotations(it->first, hts_data, &user_inputs);

			// now clean the gene_exon_map
			//
			hts_data->reset_gene_transcript_exon_map();

		}

		delete hts_data;
	}

	return 0;
}
