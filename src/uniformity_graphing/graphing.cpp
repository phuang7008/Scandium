/*
 * =====================================================================================
 *
 *       Filename:  graphing.cpp
 *
 *    Description:  This is used to test the graphing speed in C++
 *
 *        Version:  1.0
 *        Created:  01/04/2018 10:39:16 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX
 *
 * =====================================================================================
 */
#include <vector>
#include <iostream>
//#include <inttypes.h>   // for PRIu32 and PRIu64
#include <cinttypes>

#include "utils.hpp"
#include "user_inputs.hpp"

using namespace std;

/* The argument is the WGS_uniformity file produced by the Scandium run
 * 1). smoothing
 * 2). graph generation
 * 3). plot to a file
 * For gnuplot graph color: https://www2.uni-hamburg.de/Wiss/FB/15/Sustainability/schneider/gnuplot/colors.htm
 */
int main(int argc, char *argv[]) {
	// get user options
	//
	User_Inputs *user_inputs = new User_Inputs();
	user_inputs->processUserOptions(argc, argv);
	
	Utils *utils = new Utils();

	// declare Coverage_Array with size of 2
	// the first one is for the original data from the WGS_uniformity file
	// the second one is WGS_uniformity data after smoothing
	// Here I declare them as vector pointers as their sizes will grow
	//
	vector<Coverage> *cov_raw = new vector<Coverage>;
	vector<Coverage> *cov_sm  = new vector<Coverage>;
	cov_raw->reserve(utils->get_init_vector_size());
	cov_sm->reserve(utils->get_init_vector_size());

	if (user_inputs->get_smoothing_type() == 1) {
		// gVCF type smoothing
		//
		utils->processInputFile1(user_inputs, cov_raw, cov_sm);
	} else {
		// range bound smoothing
		//
		utils->processInputFile2(user_inputs, cov_raw, cov_sm);
	}

	// graphing
	//
	utils->setupChromLengthMap(user_inputs);
	map<string, uint32_t> len_map = utils->get_chrom_length_map();
	
	uint32_t length;
	//map<string,uint32_t>::iterator it = utils->get_chrom_length_map().find(user_inputs->get_chrom_id());
	map<string, uint32_t>::iterator it = len_map.find(user_inputs->get_chrom_id());
	//if (it == utils->get_chrom_length_map().end()) {
	if (it == len_map.end()) {
		cerr << "Can't find the chromosome length for chromosome id " << user_inputs->get_chrom_id() << endl;
		exit(1);
	} else {
		length = it->second;
	}
	cerr << "Total chromosome length is " << length << endl;
	cerr << "output file name is " << user_inputs->get_output_file().c_str() << endl;

	// One can avoid having to write to a file by sending gnuplot the plot '-' command
    // followed by data points followed by the letter "e"
	//
	FILE *gnuplot_pipe = popen("gnuplot -persist", "w");

	// the following line is added to get rid of the error message 
	// "Could not find/open font when opening font "arial", using internal non-scalable font"
	// export GDFONTPATH=/usr/share/fonts/liberation
	// export GNUPLOT_DEFAULT_GDFONT=LiberationSans-Regular
	//
	fprintf(gnuplot_pipe, "set term png\n");	

	fprintf(gnuplot_pipe, "set terminal png size 2000,700\n");
	fprintf(gnuplot_pipe, "set output '%s.png'\n", user_inputs->get_output_file().c_str());
	fprintf(gnuplot_pipe, "set format x \"%%.0s*10^%%T\"\n");	// to set the format to 15*10^8
	fprintf(gnuplot_pipe, "set xrange [0:%" PRIu32 "]\n", length);
	fprintf(gnuplot_pipe, "set yrange [1:200]\n");
	fprintf(gnuplot_pipe, "set title 'Uniformity Graphing'\n");
	fprintf(gnuplot_pipe, "set xlabel 'Chromsome Position'\n");
	fprintf(gnuplot_pipe, "set ylabel 'Coverage'\n");
	fprintf(gnuplot_pipe, "plot '-'\n");

	uint32_t i, j;

	
	FILE * temp = fopen("data.temp", "w");
	for (i = 0; i < cov_sm->size(); i++) {
		int32_t val=-1;
		for (j = cov_sm->at(i).get_start(); j < cov_sm->at(i).get_stop(); j++) {
			if (cov_sm->at(i).get_chrom_id() != user_inputs->get_chrom_id())
				continue;

			if (val == -1) {
				fprintf(temp, "%" PRIu32" %" PRIu32 "\n", j, cov_sm->at(i).get_cov());
				val = 1;
			} else {
				fprintf(temp, "%" PRIu32" %" PRIu32 "\n", j, 0);
			}
		}
	}
	fclose(temp);
	

	for (i = 0; i < cov_sm->size(); i++) {
		int32_t val=-1;

		for (j = cov_sm->at(i).get_start(); j < cov_sm->at(i).get_stop(); j++) {
			if (cov_sm->at(i).get_chrom_id() != user_inputs->get_chrom_id())
				continue;

			if (val == -1) {
				fprintf(gnuplot_pipe, "%" PRIu32" %" PRIu32 "\n", j, cov_sm->at(i).get_cov());
				val = 1;
			} else {
				fprintf(gnuplot_pipe, "%" PRIu32" %" PRIu32 "\n", j, 0);
			}
		}
	}

	fprintf(gnuplot_pipe, "e\n");    // finally, e
	fflush(gnuplot_pipe);    // flush the pipe-gnuplot to update the plot
	pclose(gnuplot_pipe);

	// clean-up
	//
	delete user_inputs;
	delete utils;
	delete cov_raw;
	delete cov_sm;
	return 0;
}
