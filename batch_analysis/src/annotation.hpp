/*
 * =====================================================================================
 *
 *       Filename:  annotation.hpp
 *
 *    Description:  process user-defined annotation and produce the gene annotation in details
 *
 *        Version:  1.0
 *        Created:  08/23/2017 03:51:10 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */
#ifndef ANNOTATION_HPP
#define ANNOTATION_HPP

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

#include "hts_data.hpp"
#include "utils.hpp"

class Annotation {
	public:
		Annotation(string reference_version_in);
		~Annotation() {}

		// access functions
		//MYSQL * get_mysql_connection() const { return con; }

		//void process_mysql_query(char* sql);
		//void finish_with_error();
		//void fetch_gene_exon_info_from_mysql(const string chrom_id, HTS_Data *hts_data, User_Inputs *user_inputs);
		void fetch_gene_exon_info_from_user_defined_db(const string chrom_id, HTS_Data *hts_data, User_Inputs *user_inputs, Utils *utils);
		//void fetch_and_dump_exon_info(User_Inputs *user_inputs);

		//void get_chrom_list_from_mysql(HTS_Data *hts_data, User_Inputs *user_inputs);
		void get_chrom_list_from_user_defined_db(HTS_Data *hts_data, User_Inputs *user_inputs);

	private:
		//void make_mysql_connection();

		//MYSQL *con;
		//MYSQL_RES *result;
		string reference_version;
};

#endif
