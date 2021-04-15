/*
 * =====================================================================================
 *
 *       Filename:  annotation.cpp
 *
 *    Description:  the detailed inplementation of the header file
 *
 *        Version:  1.0
 *        Created:  08/23/2017 04:05:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "annotation.hpp"

Annotation::Annotation(string reference_version_in) {
	reference_version = reference_version_in;
	//make_mysql_connection();
}

/*void Annotation::make_mysql_connection() {
	con = mysql_init(NULL);
	if (con == NULL)
		finish_with_error();

	if (mysql_library_init(0, NULL, NULL)) {
		cerr << "Could not initialize MySQL Library!" << endl;
		exit(1);
	}

	if (mysql_real_connect(con, "sug-esxa-db1", "phuang468", "phuang468", "GeneAnnotations", 0, NULL, 0) == NULL) {
		finish_with_error();
	}
}*/

/*void Annotation::finish_with_error() {
	cerr << mysql_error(con) << endl;
	mysql_close(con);
	exit(1);
}*/

/*void Annotation::process_mysql_query(char* sql) {
	if (mysql_query(con,sql)) 
		finish_with_error();

	result = mysql_store_result(con);

	if (result == NULL)
		finish_with_error();
}*/

// this function is used to fetch the detailed information on each exon region from the Gene_Exon database
// In this database, one exon is one row and gene_symbol and transcript_name are splitted already!
//
/*void Annotation::fetch_gene_exon_info_from_mysql(const string chrom_id, HTS_Data *hts_data, User_Inputs *user_inputs) {
	char *sql = static_cast<char*> (calloc(250, sizeof(char)));		// in C++, you have to force the type conversion

	string string_contains;
    if (user_inputs->get_annotation_source() == 1) {
        string_contains="N%";
    } else if (user_inputs->get_annotation_source() == 2) {
        string_contains="CCDS%";
    } else if (user_inputs->get_annotation_source() == 3) {
        string_contains="OTT%";
    } else if (user_inputs->get_annotation_source() == 4) {
        string_contains="ENST%";
    } else {
        string_contains="%%";
    }

	if (db_version.compare("hg38") == 0) {
		sprintf(sql, "SELECT exon_start, exon_end, exon_id, exon_count, gene_symbol, transcript_name FROM Gene_Exon38 WHERE chrom='%s' and transcript_name like '%s' order by exon_start", chrom_id.c_str(), string_contains.c_str());
	} else if ( (db_version.compare("hg37") == 0) || (db_version.compare("hg19") == 0) ) {
		sprintf(sql, "SELECT exon_start, exon_end, exon_id, exon_count, gene_symbol, transcript_name FROM Gene_Exon37 WHERE chrom='%s' and transcript_name like '%s' order by exon_start", chrom_id.c_str(), string_contains.c_str());
	} else {
		cerr << "The db version " << db_version << " is not available!" << endl;
	}

	//cerr << sql << endl;

	process_mysql_query(sql);
	MYSQL_ROW row;

	while ((row = mysql_fetch_row(result))) {
		Exon exon;
		exon.set_exon_start(atol(row[0]));
		exon.set_exon_end(atol(row[1]));
		exon.set_exon_id(atol(row[2]));
		string gene_symbol(row[4]);
		string transcript_name(row[5]);

		// if the gene symbol is not found in the unordered_map, create an empty container first
		//
		if (hts_data->get_gene_transcript_exon_map()->find(gene_symbol) == hts_data->get_gene_transcript_exon_map()->end()) {
			hts_data->get_gene_transcript_exon_map()->insert(std::make_pair(gene_symbol, unordered_map<string, vector<Exon> >()));
		}

		// if the transcript name is not found in the unordered_map, create an empty container first
		//
		if (hts_data->get_gene_transcript_exon_map()->at(gene_symbol).find(transcript_name) == hts_data->get_gene_transcript_exon_map()->at(gene_symbol).end()) {
			hts_data->get_gene_transcript_exon_map()->at(gene_symbol).insert(std::make_pair(transcript_name, vector<Exon>()));
			hts_data->get_gene_transcript_exon_map()->at(gene_symbol).at(transcript_name).reserve(hts_data->get_initial_vector_size());
		}

		// now we need to update the num_of_low_qual_bases!
		// first get the target_exon corresponding to the current chromosome id,  gene_symbol and transcript_name
		//
		if (hts_data->get_low_cov_exon_map()->find(chrom_id) != hts_data->get_low_cov_exon_map()->end()) {
			if (hts_data->get_low_cov_exon_map()->at(chrom_id).find(gene_symbol) != hts_data->get_low_cov_exon_map()->at(chrom_id).end()) {
				if (hts_data->get_low_cov_exon_map()->at(chrom_id).at(gene_symbol).find(transcript_name) != 
						hts_data->get_low_cov_exon_map()->at(chrom_id).at(gene_symbol).end()) {

					// iterate through exon vector
					//
					for (auto it=hts_data->get_low_cov_exon_map()->at(chrom_id).at(gene_symbol).at(transcript_name).begin(); 
							  it!=hts_data->get_low_cov_exon_map()->at(chrom_id).at(gene_symbol).at(transcript_name).end(); ++it) {
						if (it->get_exon_id() == exon.get_exon_id()) {
							// iterate through low cov regions vector
							//
							for (auto lc=it->get_low_cov_regions()->begin(); lc!=it->get_low_cov_regions()->end(); ++lc) {
								exon.add_low_cov_regions(*lc);
							}

							exon.set_num_of_low_qual_bases(exon.get_num_of_low_qual_bases() + it->get_num_of_low_qual_bases());

							//cerr << gene_symbol << "\t" << transcript_name << "\t" << exon.get_exon_id() << "\t" << exon.get_num_of_low_qual_bases() << endl;
						}
					}
				}
			}
		}

		// insert the exon
		//
		hts_data->get_gene_transcript_exon_map()->at(gene_symbol).at(transcript_name).push_back(exon);
	}

	if (result) {
		mysql_free_result(result);
		result = NULL;
	}

	free(sql);
	sql=NULL;
}*/

void Annotation::fetch_gene_exon_info_from_user_defined_db(const string chrom_in, HTS_Data *hts_data, User_Inputs *user_inputs, Utils *utils) {
	ifstream infile(user_inputs->get_user_defined_annotation_file());
	string chrom_id, exon_annotation;                                                                     
	uint32_t exon_start, exon_end;

	// here is an example:
	// chr6	32610994	32610995	PD_HLA-DRB5_rs112485576|PD_HLA-DRB5_rs112485576|snp_1|snp
	// chr1	7969344	7969404	PARK7|NM_007262|cds_3|gene
	//
	while (infile >> chrom_id >> exon_start >> exon_end >> exon_annotation) {
		if (chrom_id.compare(chrom_in) == 0 && exon_annotation.length() > 0) {
			Exon exon;
			exon.set_exon_start(exon_start);
			exon.set_exon_end(exon_end);

			vector<string> exon_info = utils->split(exon_annotation, '|');
			string gene_symbol(exon_info[0]);

			vector<string> exon_id_info = utils->split(exon_info[2], '_');
			uint16_t exon_id = atoi(exon_id_info.back().c_str());
			exon.set_exon_id(exon_id);

			exon_id_info.pop_back();	 // remove exon_id 
			exon_id_info.pop_back();     // remove the 'cds' field from the exon_info vector
			string transcript_name( exon_info[1] );

			//hts_data->get_exon_count_per_transcript_map()->at(transcript_name) += 1;
			//(*hts_data->get_exon_count_per_transcript_map())[transcript_name]++;

			// if the gene symbol unordered_map is not found, create an empty container first
			//
			if (hts_data->get_gene_transcript_exon_map()->find(gene_symbol) == hts_data->get_gene_transcript_exon_map()->end()) {
				hts_data->get_gene_transcript_exon_map()->insert(std::make_pair(gene_symbol, unordered_map<string, vector<Exon> >()));
			}

			// if the transcript name unordered_map is not found, create an empty container first
			//
			if (hts_data->get_gene_transcript_exon_map()->at(gene_symbol).find(transcript_name) == hts_data->get_gene_transcript_exon_map()->at(gene_symbol).end()) {
				hts_data->get_gene_transcript_exon_map()->at(gene_symbol).insert(std::make_pair(transcript_name, vector<Exon>()));
				hts_data->get_gene_transcript_exon_map()->at(gene_symbol).at(transcript_name).reserve(hts_data->get_initial_vector_size());
			}

			// now we need to update the num_of_low_qual_bases!
			// first get the target_exon corresponding to the current chromosome id,  gene_symbol and transcript_name
			//
			if (hts_data->get_low_cov_exon_map()->find(chrom_id) != hts_data->get_low_cov_exon_map()->end()) {
				if (hts_data->get_low_cov_exon_map()->at(chrom_id).find(gene_symbol) != hts_data->get_low_cov_exon_map()->at(chrom_id).end()) {
					if (hts_data->get_low_cov_exon_map()->at(chrom_id).at(gene_symbol).find(transcript_name) !=
						hts_data->get_low_cov_exon_map()->at(chrom_id).at(gene_symbol).end()) {

						// iterate through exon vector
						//
						for (auto it=hts_data->get_low_cov_exon_map()->at(chrom_id).at(gene_symbol).at(transcript_name).begin();
								  it!=hts_data->get_low_cov_exon_map()->at(chrom_id).at(gene_symbol).at(transcript_name).end(); ++it) {
							if (it->get_exon_id() == exon.get_exon_id()) {
								// iterate through low cov regions vector
								//
								for (auto lc=it->get_low_cov_regions()->begin(); lc!=it->get_low_cov_regions()->end(); ++lc) {
									exon.add_low_cov_regions(*lc);
								}

								exon.set_num_of_low_qual_bases(exon.get_num_of_low_qual_bases() + it->get_num_of_low_qual_bases());

								//cerr << gene_symbol << "\t" << transcript_name << "\t" << exon.get_exon_id() << "\t" << exon.get_num_of_low_qual_bases() << endl;
							}
						}
					}
				}
			}

			// insert the exon
			//
			hts_data->get_gene_transcript_exon_map()->at(gene_symbol).at(transcript_name).push_back(exon);
		}
	}

	infile.close();
}

// This function is used to obtain the exon detailed information from one of the following two sources:
// 1). the MySQL database (here the annotation is combined)
// 2). user defined database
// The output looks like: chr1    13220   14409   DDX11L1|NR_046018_exon_2
// users will have to split the annotation later for further analysis
//
/*void Annotation::fetch_and_dump_exon_info(User_Inputs *user_inputs) { 
	// open a file stream for writing
	//
	ofstream os_file(user_inputs->get_exon_bed_file());
	if (user_inputs->get_user_defined_annotation_file() == "") {
		// now user defined database, so we search the MySQL database
		//
		string database="Gene_Exon";

		char *sql = static_cast<char*> (calloc(250, sizeof(char)));     // in C++, you have to force the type conversion

		string string_contains;
		if (user_inputs->get_annotation_source() == 1) {
			string_contains="N%";
		} else if (user_inputs->get_annotation_source() == 2) {
			string_contains="CCDS%";
		} else if (user_inputs->get_annotation_source() == 3) {
			string_contains="OTT%";
		} else if (user_inputs->get_annotation_source() == 4) {
			string_contains="ENST%";
		} else {
			string_contains="%%";
		}

		if (db_version.compare("hg38") == 0) {
			sprintf(sql, "SELECT chrom, exon_start, exon_end, exon_id, exon_count, gene_symbol, transcript_name FROM %s%d WHERE transcript_name like '%s' order by chrom, exon_start", database.c_str(), 38, string_contains.c_str());
	    } else if ( (db_version.compare("hg37") == 0) || (db_version.compare("hg19") == 0) ) {
		    sprintf(sql, "SELECT chrom, exon_start, exon_end, exon_id, exon_count, gene_symbol, transcript_name FROM %s%d WHERE transcript_name like '%s' order by chrom, exon_start", database.c_str(), 37, string_contains.c_str());
	    } else {
		    cerr << "The db version " << db_version << " is not available!" << endl;
	    }

		//cerr << sql << endl;

	    process_mysql_query(sql);
	    MYSQL_ROW row;
	    Exon exon;

		while ((row = mysql_fetch_row(result))) {
			char *tmp_str = static_cast<char*> (calloc(100, sizeof(char)));
			//sprintf(tmp_str, "%s|%s_exon_%s_%s", row[5], row[6], row[3], row[4]);
			sprintf(tmp_str, "%s|%s_exon_%s", row[5], row[6], row[3]);
			os_file << row[0] << "\t" << row[1] << "\t" << row[2] << "\t" << tmp_str << endl;
			free(tmp_str);
		}

		free(sql);
	} else {
		// user-defined database is provided, read everything in!
		//
		ifstream infile(user_inputs->get_user_defined_annotation_file());
		string chrom_id, exon_annotation;
		uint32_t exon_start, exon_end;
		while (infile >> chrom_id >> exon_start >> exon_end >> exon_annotation) {
			os_file << chrom_id << "\t" << exon_start << "\t" << exon_end << "\t" << exon_annotation << endl;
		}
		infile.close();
	}

	os_file.close();
}*/

/*void Annotation::get_chrom_list_from_mysql(HTS_Data *hts_data, User_Inputs *user_inputs) {
	string database="Gene_Exon";
	char *sql = static_cast<char*> (calloc(250, sizeof(char)));     // in C++, you have to force the type conversion

	string string_contains;
	if (user_inputs->get_annotation_source() == 1) { 
		string_contains="N%";
	} else if (user_inputs->get_annotation_source() == 2) {
		string_contains="CCDS%";
	} else if (user_inputs->get_annotation_source() == 3) {
		string_contains="OTT%";
	} else if (user_inputs->get_annotation_source() == 4) {
		string_contains="ENST%";
	} else {
		string_contains="%%";
	}

	if (db_version.compare("hg38") == 0) {
		sprintf(sql, "SELECT DISTINCT chrom FROM %s%d WHERE transcript_name like '%s'", database.c_str(), 38, string_contains.c_str());
	} else if ( (db_version.compare("hg37") == 0) || (db_version.compare("hg19") == 0) ) {
		sprintf(sql, "SELECT DISTINCT chrom FROM %s%d WHERE transcript_name like '%s'", database.c_str(), 37, string_contains.c_str());
	} else {
		cerr << "The db version " << db_version << " is not available!" << endl;
	}

	process_mysql_query(sql);
	MYSQL_ROW row;

	while ((row = mysql_fetch_row(result))) {
		hts_data->add_to_chromosome_list(row[0]);
	}
}*/

void Annotation::get_chrom_list_from_user_defined_db(HTS_Data *hts_data, User_Inputs *user_inputs) {

	ifstream infile(user_inputs->get_user_defined_annotation_file());
	string chrom_id, chrom_prev, exon_annotation;
	uint32_t exon_start, exon_end;
	chrom_prev = "None";

	// here is an example:
	// 1       1043537 1043732 AGRN|NM_001305275_cds_8
	//
	while (infile >> chrom_id >> exon_start >> exon_end >> exon_annotation) {
		if (chrom_id != chrom_prev) {
			hts_data->add_to_chromosome_list(chrom_id);
			chrom_prev = chrom_id;
		} else {
			continue;
		}
	}

	infile.close();
}
