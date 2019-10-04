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
 * 1). read data from the input uniformity file
 * 2). graph generation
 * 3). plot to a file
 * For gnuplot graphing: https://www2.uni-hamburg.de/Wiss/FB/15/Sustainability/schneider/gnuplot/colors.htm
 */
int main(int argc, char *argv[]) {
	// get user options
	//
	User_Inputs *user_inputs = new User_Inputs();
	user_inputs->processUserOptions(argc, argv);
	
	Utils *utils = new Utils();

	// declare Coverage_Array
	// Here I declare it as vector pointers as their sizes will grow
	//
	vector<Coverage> *cov_raw = new vector<Coverage>;
	cov_raw->reserve(utils->get_init_vector_size());

	// check if a centromere region file is provided!
	//
	if (user_inputs->get_chromosome_centromere_file().length() > 0)
		utils->read_in_centromere_regions(user_inputs);

	utils->processInputFile(user_inputs, cov_raw);

	// graphing
	//
	uint32_t chr_length = user_inputs->get_chromosome_length();
	cerr << "Total chromosome length is " << chr_length << endl;
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
	fprintf(gnuplot_pipe, "set xrange [0:%" PRIu32 "]\n", chr_length);
	fprintf(gnuplot_pipe, "set yrange [1:200]\n");
	fprintf(gnuplot_pipe, "set title 'Uniformity Graphing'\n");
	fprintf(gnuplot_pipe, "set xlabel 'Chromsome Position'\n");
	fprintf(gnuplot_pipe, "set ylabel 'Coverage'\n");
	fprintf(gnuplot_pipe, "plot '-' with points pointtype 7 ps 0.75 lc rgb \"#32CD32\"\n");
	//fprintf(gnuplot_pipe, "plot '-' with points pointtype 26\n");

	uint32_t i, j;

	
	/*FILE * temp = fopen("data.temp", "w");
	for (i = 0; i < cov_raw->size(); i++) {
		int32_t val=-1;
		for (j = cov_raw->at(i).get_start(); j < cov_raw->at(i).get_stop(); j++) {
			if (val == -1) {
				fprintf(temp, "%" PRIu32" %" PRIu32 "\n", j, cov_raw->at(i).get_cov());
				val = 1;
			} else {
				fprintf(temp, "%" PRIu32" %" PRIu32 "\n", j, 0);
			}
		}
	}
	fclose(temp);
	*/

	for (i = 0; i < cov_raw->size(); i++) {
		int32_t val=-1;

		for (j = cov_raw->at(i).get_start(); j < cov_raw->at(i).get_stop(); j++) {
			if (val == -1) {
				fprintf(gnuplot_pipe, "%" PRIu32" %" PRIu32 "\n", j, cov_raw->at(i).get_cov());
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
	return 0;
}