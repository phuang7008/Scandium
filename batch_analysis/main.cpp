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

	Annotation annotation=Annotation(user_inputs.get_reference_version());

	// declare HTS_Data object on the heap as the data wil grow!
	//
	HTS_Data *hts_data = new HTS_Data();

	// process low coverage annotation file and genes to check list
	//
	Utils utils=Utils();
	utils.store_gene_list(&user_inputs, hts_data);

	// get the low_cov_exon_map, which we just read in
	//
	unordered_map <string, unordered_map<string, unordered_map<string, vector<Exon> > > > *low_cov_exon_map = hts_data->get_low_cov_exon_map();
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

	// need to find out the chrom_list, from the user_defined annotations
	//
	if (user_inputs.get_user_defined_annotation_file() == "") {
        cout << "Since there is no user-defined annotation provided, the program will exit!" << std::endl;
	} else {
        cout << "User provided an annotation file, so will proceed accordingly" << std::endl;
	
        utils.read_low_cov_annotation_file(&user_inputs, hts_data);
	    annotation.get_chrom_list_from_user_defined_db(hts_data, &user_inputs);
	    unordered_map<string, int>* chrom_list = hts_data->get_chromosome_list();

	    // loop through chrom_list, as we are going to process one chromosome at a time!
	    //
	    for (auto it=chrom_list->begin(); it!=chrom_list->end(); ++it) {
		    cout << "Current chromosome is " << it->first << endl;

			// fetch the gene exon information from user-defined database
			//
			fprintf(stderr, "Use user-defined database for gene annotations\n");
			annotation.fetch_gene_exon_info_from_user_defined_db(it->first, hts_data, &user_inputs, &utils);

		    // output the detailed gene, transcript and exon information
    		//
	    	utils.output_detailed_annotations(it->first, hts_data, &user_inputs);

		    // now clean the gene_exon_map
		    //
		    hts_data->reset_gene_transcript_exon_map();
        }
	}

	delete hts_data;

	return 0;
}
