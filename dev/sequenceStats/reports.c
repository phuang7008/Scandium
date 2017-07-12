/*
 * =====================================================================================
 *
 *       Filename:  reports.c
 *
 *    Description:  the detailed implementation of sequencing statistics
 *
 *        Version:  1.0
 *        Created:  02/22/2017 01:55:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang (phuang@bcm.edu)
 *        Company:  Baylor College of medicine
 *
 * =====================================================================================
 */

#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reports.h"
#include "utils.h"

void writeCoverage(char *chrom_id, Bed_Info *Ns_bed_info, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, MYSQL *con, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions, Regions_Skip_MySQL *all_site_reports) {

    // First, we need to find the index that is used to track current chromosome chrom_id
    int32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	if (idx == -1) return;

    // for the whole genome, we need to use the file that contains regions of all Ns in the reference
    // As they will be not used, so we are going to set the count info in these regions to 0
    //if (N_FILE_PROVIDED)
    //    zeroAllNsRegions(chrom_id, Ns_bed_info, chrom_tracking);

    // write to the file that contains whole genome info
	uint32_t i=0;
    if(user_inputs->wgs_coverage) {

		if (user_inputs->Write_WGS) {
			FILE *wgs_fp = fopen(user_inputs->wgs_cov_file, "a");
			printf("Whole Genome output for chrom id %s is on\n", chrom_id);

			// no need to add newline as the next line will take care of it!
			fprintf(wgs_fp, ">%s 1 %"PRIu32, chrom_id, chrom_tracking->chromosome_lengths[idx]-1);

	        for (i = 0; i < chrom_tracking->chromosome_lengths[idx]; i++) {
				if(i%100==0) fputc('\n', wgs_fp);
				fprintf(wgs_fp, "%d ", chrom_tracking->coverage[idx][i]);

				addBaseStats(stats_info, (uint32_t) chrom_tracking->coverage[idx][i], 0, 1);
			}

			fputc('\n', wgs_fp);
			fclose(wgs_fp);
		} else {
			for (i = 0; i < chrom_tracking->chromosome_lengths[idx]; i++)
				addBaseStats(stats_info, (uint32_t) chrom_tracking->coverage[idx][i], 0, 1);
		}
    }

	// if the target bed file is available, we will need to handle it here and write the results to cov.fasta file
	if (TARGET_FILE_PROVIDED) {
		FILE * cov_fp = fopen(user_inputs->capture_cov_file, "a");
		FILE * missed_target_fp = fopen(user_inputs->missed_targets_file, "a");

		for(i = 0; i < target_info->size; i++) {
            if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
                continue;

            stats_info->cov_stats->total_targets++;
            uint32_t start = target_info->coords[i].start;
            uint32_t end = target_info->coords[i].end;
            int length = end - start + 1;
            bool collect_target_cov = length > 99 ? true : false ;  // TODO: is it the min length of the exon? Check.
			//if (collect_target_cov) 
			//	fprintf(stderr, "%s %"PRIu32" %"PRIu32" %d\n", chrom_id, start, end, length);

            uint32_t j=0;
            if (collect_target_cov) {
                for(j = 0; j < PRIMER_SIZE; j++) {
                    if ( (start - j) < 0 || (end + j) >= chrom_tracking->chromosome_lengths[idx])
                        continue;

					stats_info->five_prime[j]  += chrom_tracking->coverage[idx][start-j];
					stats_info->three_prime[j] += chrom_tracking->coverage[idx][end+j];
                }
			}

			//In the original java code, the following is used for the pc and pc2 definition
            //short pc[101], pc2[101];
			//After discussion with Divya and Qiaoyan, here is what we understand.
			//pc and pc2 is used to setup the 'target_coverage' variable which is defined as an array of 101 in the origain java code
			//here I will use percentage_bin and percentage_count_per_bin
			uint32_t percentage_bin[101], percentage_count_per_bin[101];
			for(j=0; j<101; j++) {
				percentage_bin[j] = 0;
				percentage_count_per_bin[j] = 0;
			}

            bool target_hit = false;
            fprintf(cov_fp, ">%s %"PRIu32" %"PRIu32"\n", chrom_id, start, end);

            bool space_it = false;
            if(end - start > 10000) space_it = true;

            for(j = 0; j < length; j++) {
                if (j+start >= chrom_tracking->chromosome_lengths[idx])
                    continue;

                if (space_it && j%100 == 0) fputc('\n', cov_fp);    // enter a new line after every 100 bases

				uint32_t cov = chrom_tracking->coverage[idx][j+start];
				addBaseStats(stats_info, cov, 1, 0);

				if (!target_hit && (cov > 0))
					target_hit = true;

                // output to the cov.fasta file
                fprintf(cov_fp, "%d ", cov);

                if (collect_target_cov) {
					float num, den;
					memcpy(&num, &j, sizeof(num));
					memcpy(&den, &length, sizeof(den));
                    int percentage_pos = (int)((num*100)/den + 0.5);
                    //int percentage_pos = lround((num/den)*100);
                    percentage_bin[percentage_pos] += cov;
                    percentage_count_per_bin[percentage_pos] += 1;
					//fprintf(stderr, "%"PRIu32" %d %d\n", j, length, percentage_pos);
                }
            }

            // output a newline char to the cov.fasta file 
            fputc('\n', cov_fp);

			if (collect_target_cov) {
				for (j = 0; j < 101; j++) {
					if(percentage_count_per_bin[j] != 0) {
						float num, den;
				        memcpy(&num, &percentage_bin[j], sizeof(num));
					    memcpy(&den, &percentage_count_per_bin[j], sizeof(den));
						int d = (int) (num/den+0.5);
	                    //int d = lround(num/den);
		                percentage_bin[j] = d;	
					}
                }

				for(j = 0; j < 101; j++)
					stats_info->target_coverage[j] += percentage_bin[j];
            }

			if (target_hit) {
				stats_info->cov_stats->hit_target_count += 1;
			} else {
				// need to write to the missed target file
				fprintf(missed_target_fp, "%s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, start, end);
				bool hit = false;
				for (j = start - BUFFER; j < start && !hit; j++) {
					if (j >= chrom_tracking->chromosome_lengths[idx])
						continue;

					if ( chrom_tracking->coverage[idx][j] > 0 )
						hit = true;
				}

				if (hit)
					stats_info->cov_stats->hit_target_buffer_only_count += 1;
			}

			// now need to report those Capture regions with low or too high coverages
	        // NOTE: the bed format is different here, the end position is included!
    	    //
        	FILE *capture_low_x_fp    = fopen(user_inputs->capture_low_cov_file, "a");
	        FILE *capture_high_x_fp   = fopen(user_inputs->capture_high_cov_file, "a");
	        FILE *capture_all_site_fp = fopen(user_inputs->capture_all_site_file, "a");

			// For All Sites Report
			produceCaptureAllSitesReport(start, length, chrom_tracking, chrom_id, user_inputs, capture_all_site_fp, con, all_site_reports);

			// For low coverage and high coverage Report
			writeLow_HighCoverageReport(start, length, chrom_tracking, chrom_id, user_inputs, capture_low_x_fp, capture_high_x_fp, con, inter_genic_regions, intronic_regions, exon_regions);

	        fclose(capture_low_x_fp);
    	    fclose(capture_high_x_fp);
    	    fclose(capture_all_site_fp);
		}
		fclose(cov_fp);
		fclose(missed_target_fp);
	}
}

void produceCaptureAllSitesReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char * chrom_id, User_Input *user_inputs, FILE *fh_all_sites, MYSQL *con, Regions_Skip_MySQL *all_site_reports) {
	uint32_t i=0;
	uint64_t cov_total=0;
	int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	//int32_t chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, all_site_reports);
	//if (chrom_idx > 25) return;

	//if (chrom_idx == -1 || chrom_idr == -1) return;
	if (chrom_idx == -1) return;

	for (i = begin; i < begin+length; i++) {
		cov_total += chrom_tracking->coverage[chrom_idx][i];
	}

	uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(length) + 0.5);
	fprintf(fh_all_sites, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"", chrom_tracking->chromosome_ids[chrom_idx], begin, begin+length-1, length, ave_coverage);

	/*int32_t prev_report_location = 0;
	if (all_site_reports) {
		int32_t report_location = binary_search(all_site_reports, begin, begin+length-1, chrom_idr, prev_report_location);
		if (report_location == -1) {
			fprintf(stderr, "Search return -1 for %"PRIu32" and %"PRIu32" and chromosome id is %s\n", begin, begin+length-1, chrom_id);
			return;
		}

		if (prev_report_location < report_location) prev_report_location = report_location;
		fprintf(fh_all_sites, "\t%s\n", all_site_reports->gene[chrom_idr][report_location]);
	} else {
		char *annotation = produceGeneAnnotations(begin, begin+length-1, chrom_tracking->chromosome_ids[chrom_idx], con);
		fprintf(fh_all_sites, "\t%s\n", annotation);
		if (annotation) free(annotation);
	}*/

	char *annotation = produceGeneAnnotations(begin, begin+length-1, chrom_tracking->chromosome_ids[chrom_idx], con);
    fprintf(fh_all_sites, "\t%s\n", annotation);
    if (annotation) free(annotation);
}

void writeAnnotations(char *chrom_id, Bed_Info *Ns_bed_info, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, MYSQL *con, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions, Regions_Skip_MySQL *all_site_reports) {
	// First, we need to find the index that is used to track current chromosome chrom_id
    int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);

	if(user_inputs->wgs_coverage) {
		// now need to report those regions with low or too high coverages
	    // NOTE: the bed format is different here, the end position is included!
		//
	    FILE *wgs_low_x_fp  = fopen(user_inputs->wgs_low_cov_file, "a");
    	FILE *wgs_high_x_fp = fopen(user_inputs->wgs_high_cov_file, "a");

	    writeLow_HighCoverageReport(0, chrom_tracking->chromosome_lengths[chrom_idx], chrom_tracking, chrom_id, user_inputs, wgs_low_x_fp, wgs_high_x_fp, con, inter_genic_regions, intronic_regions, exon_regions);

        fclose(wgs_low_x_fp);
        fclose(wgs_high_x_fp);
	}
}

uint32_t writeLow_HighCoverageReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char *chrom_id, User_Input *user_inputs,FILE *fh_low, FILE *fh_high, MYSQL *con, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions) {
	// for debugging
	//if (strcmp(chrom_tracking->chromosome_ids[chrom_idx], "1") == 0) return i;
	
	// First, we need to find the index that is used to track current chromosome chrom_id
    int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
    int32_t chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, intronic_regions);
    if (chrom_idx == -1 || chrom_idr == -1) return;
	
	int32_t index_inter_genic_location=-1, index_intronic_location=-1, exon_location=-1;
	int32_t prev_index_inter_genic_location=0, prev_index_intronic_location=0, prev_exon_location=0;
	char *intronic_info = calloc(1, sizeof(char));
	char *exon_info = calloc(1, sizeof(char));
	uint32_t i=0;

	for (i = begin; i < begin+length; i++) {
        uint32_t start=0, end=0;
        uint64_t cov_total=0;

        // for low coverage
        if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] < user_inputs->low_coverage_to_report) {
            start = i;

            while(i < begin+length && chrom_tracking->coverage[chrom_idx][i] < user_inputs->low_coverage_to_report) {
                cov_total += chrom_tracking->coverage[chrom_idx][i];
                i++;
            }
            end = i;
			if (start == end) continue;

            //float ave_coverage = (float)cov_total / (float)(end - start);
            uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);

			fprintf(fh_low, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage);

			// generate the annotation here
			if (user_inputs->annotation_on) {
				// For annotation, I only used normal chromosome IDs, anything outside here won't be used
				if ( strstr(chrom_tracking->chromosome_ids[chrom_idx], "GL000") ||
					 strstr(chrom_tracking->chromosome_ids[chrom_idx], "hs37d") || (chrom_idr >= 25)) {
					fprintf(fh_low, "\t\t.\t.\t.\t.\t.\t.\n");
				} else {
					// check if the position located in inter-genic regions
					if (index_inter_genic_location == -1) {
						index_inter_genic_location = checkInterGenicRegion(inter_genic_regions, start, end-1, chrom_idr, prev_index_inter_genic_location);
						if (index_inter_genic_location != -1) {
							prev_index_inter_genic_location = index_inter_genic_location;
							fprintf(fh_low, "\t\t.\t.\t.\t.\t.\t.\n");
						}
					} else {
						// for future round, we can skip the search, just check to see if current one is also part of the previous inter-genic region
						if (verifyIndex(inter_genic_regions, start, end-1, chrom_idr, index_inter_genic_location)) {
							fprintf(fh_low, "\t\t.\t.\t.\t.\t.\t.\n");
						} else {
							index_inter_genic_location = checkInterGenicRegion(inter_genic_regions, start, end-1, chrom_idr, prev_index_inter_genic_location);
							if (index_inter_genic_location != -1) {
								prev_index_inter_genic_location = index_inter_genic_location;
								fprintf(fh_low, "\t\t.\t.\t.\t.\t.\t.\n");
							} else {
								index_inter_genic_location = -1;
							}
						}
					}

					// check if the position located in intronic regions
					if (index_inter_genic_location == -1) {
						if (index_intronic_location == -1) {
						   	index_intronic_location = checkIntronicRegion(intronic_regions, start, end-1, chrom_idr, &intronic_info, prev_index_intronic_location);
							if (index_intronic_location != -1 && intronic_info) {
								fprintf(fh_low, "\t%s\t.\t.\t.\t.\n", intronic_info);
								prev_index_intronic_location = index_intronic_location;
								//printf("intronic_info %s\n", intronic_info);
								//free(intronic_info);
							}
						} else {
							if (verifyIndex(intronic_regions, start, end-1, chrom_idr, index_intronic_location)) {
								fprintf(fh_low, "\t%s\t.\t.\t.\t.\n", intronic_info);
								//free(intronic_info);
							} else {
								index_intronic_location = checkIntronicRegion(intronic_regions, start, end-1, chrom_idr, &intronic_info, prev_index_intronic_location);
								if (index_intronic_location != -1 && intronic_info) {
									prev_index_intronic_location = index_intronic_location;
									fprintf(fh_low, "\t%s\t.\t.\t.\t.\n", intronic_info);
								} else {
									index_intronic_location = -1;
									strcpy(intronic_info, "");
								}
							}
						}
					}

					// check if it locates at which exon location if both the above failed
					/*if (index_inter_genic_location == -1 && index_intronic_location == -1) {
						if (exon_location == -1) {
							exon_location = checkExonRegion(exon_regions, start, end-1, chrom_idr, &exon_info);
							if (exon_location != -1 && exon_info) {
								fprintf(fh_low, "\t%s\n", exon_info);
							}
						} else {
							if (verifyIndex(intronic_regions, start, end-1, chrom_idr, index_intronic_location)) {
								fprintf(fh_low, "\t%s\n", exon_info);
							} else {
								exon_location = checkExonRegion(exon_regions, start, end-1, chrom_idr, &exon_info);
								if (exon_location != -1 && exon_info) {
									fprintf(fh_low, "\t%s\n", exon_info);
								} else {
									exon_location = -1;
									strcpy(exon_info, "");
								}
							}
						}
					}*/

					//if (index_inter_genic_location == -1 && index_intronic_location == -1 && exon_location == -1) 
					if (index_inter_genic_location == -1 && index_intronic_location == -1) {
						char *annotation = produceGeneAnnotations(start, end-1, chrom_tracking->chromosome_ids[chrom_idx], con);
						fprintf(fh_low, "\t%s\n", annotation);
						if (annotation) free(annotation);
					}
				}
			} else {
		        fprintf(fh_low, "\n");
			}
        }
		//fflush(fh_low);

        // for high coverage
		start = 0;
		end = 0;
		cov_total = 0;
        if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] >= user_inputs->high_coverage_to_report) {
            start = i;

            while( i < begin+length && chrom_tracking->coverage[chrom_idx][i] >= user_inputs->high_coverage_to_report) {
				if ( (user_inputs->upper_bound_to_report == -1 ) || 
						(user_inputs->upper_bound_to_report > -1 && chrom_tracking->coverage[chrom_idx][i] < user_inputs->upper_bound_to_report)) {
					cov_total += chrom_tracking->coverage[chrom_idx][i];
					i++;
				} else {
					break;
				}
            }
            end = i;
			if (start == end) continue;

            //float ave_coverage = (float)cov_total / (float)(end - start);
            uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);

            //fprintf(fh_low, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%.2f\n", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage);
            fprintf(fh_high, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage);
        }
		//fflush(fh_high);
    }

	if (intronic_info) free(intronic_info); 
	if (exon_info) free(exon_info); 

	return i;
}

void addBaseStats(Stats_Info *stats_info, uint32_t cov_val, uint8_t target, uint8_t wgs) {
	// need to check and update max coverage information
	if (cov_val > stats_info->cov_stats->max_coverage) {
		stats_info->cov_stats->max_coverage = cov_val;
		stats_info->cov_stats->base_with_max_coverage = 1;
	} else if (cov_val == stats_info->cov_stats->max_coverage) {
		stats_info->cov_stats->base_with_max_coverage += 1;
	}

	if (cov_val < 0) {
        fprintf(stderr, "Coverage less than 0!!!!!!!\n");
    }

	// for histogram only
	uint32_t tmp_val = cov_val;
	if (tmp_val > 1000) tmp_val = 1000;
	if (target == 1) addValueToKhashBucket32(stats_info->target_cov_histogram, tmp_val, 1);
	if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_cov_histogram, tmp_val, 1);

    if (target == 1) stats_info->cov_stats->total_target_coverage += (uint64_t) cov_val;
    if (wgs == 1)    stats_info->cov_stats->total_genome_coverage += (uint64_t) cov_val;

    if (cov_val > 0) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 1, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 1, 1);
    } else if (cov_val == 0) {
		if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 0, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 0, 1);
	}

	if (cov_val >= 5) {
		if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 5, 1);
		if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 5, 1);
	}

	if (cov_val >= 10) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 10, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 10, 1);
    }

	if (cov_val >= 15) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 15, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 15, 1);
    }

	if (cov_val >= 20) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 20, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 20, 1);
    }

	if (cov_val >= 30) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 30, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 30, 1);
    }

	if (cov_val >= 40) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 40, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 40, 1);
    }

	if (cov_val >= 50) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 50, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 50, 1);
    }

	if (cov_val >= 60) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 60, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 60, 1);
    }

	if (cov_val >= 100) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 100, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 100, 1);
    }

	if (cov_val >= 500) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 500, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 500, 1);
    }

	if (cov_val >= 1000) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 1000, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 1000, 1);
    }

	if (target == 1) addValueToKhashBucket32(stats_info->target_coverage_for_median, cov_val, 1);
	if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_coverage_for_median, cov_val, 1);
}

void writeReport(Stats_Info *stats_info, User_Input *user_inputs) {
	if (TARGET_FILE_PROVIDED && stats_info->cov_stats->total_targeted_bases == 0) {
		fprintf(stderr, "Total targeted bases is zero.  This means that no read has aligned to a chromosome that contains a target.");
		fprintf(stderr, "No target matches a chromosome in the BAM, or something else went wrong.  Aborting.\n");
		exit(1);
	}

	if (stats_info->cov_stats->total_reads_aligned == 0) {
		fprintf(stderr, "No reads aligned. Aborting.\n");
		exit(1);
	}

	uint32_t non_duplicate_reads = stats_info->cov_stats->total_reads_aligned - stats_info->cov_stats->total_duplicate_reads;
	if (non_duplicate_reads == 0) {
        fprintf(stderr, "All reads are duplicates. Aborting.\n");
        exit(1);
    }

	if (TARGET_FILE_PROVIDED && stats_info->cov_stats->total_targets == 0) {
        //I don't think we should ever see this error, as its dealt with above.
        fprintf(stderr, "No target regions given.  Aborting.\n");
        exit(1);
    }

	//printf("Before Median calculation\n");
   	uint64_t sum=0;
	int32_t i=0;
	double average_coverage;
	int ret;
	uint16_t bins[13] = { 0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 100, 500, 1000 };
	khiter_t k_iter;
	char *file_name;

	if (user_inputs->wgs_coverage) {
		//Do not consider the Ns for Median calculation.
		uint64_t total_genome_non_Ns_bases = stats_info->cov_stats->total_genome_bases - stats_info->cov_stats->total_Ns_bases;
		for (k_iter=0; k_iter!=kh_end(stats_info->genome_coverage_for_median); k_iter++) {
			if (kh_exist(stats_info->genome_coverage_for_median, k_iter)) {

				if (sum >= (total_genome_non_Ns_bases/2)) {
					stats_info->cov_stats->median_genome_coverage = k_iter--;
					break;
				}else{
					sum += kh_value(stats_info->genome_coverage_for_median, k_iter);;
				}
			}
        }

		file_name = calloc(strlen(user_inputs->out_file) + 30, sizeof(char));
		strcpy(file_name, user_inputs->out_file);
		strcat(file_name, ".WGCoverageReport.csv");
		FILE *out_fp = fopen(file_name, "w");

		average_coverage = (double) stats_info->cov_stats->total_genome_coverage/ (double) total_genome_non_Ns_bases;
		outputGeneralInfo(out_fp, stats_info, average_coverage, 1);

		fprintf(out_fp, "Base Stats\n");
	    fprintf(out_fp, "Bases Targeted:,%"PRIu64"\n", total_genome_non_Ns_bases);

		for(i=0; i<13; i++) {
			int32_t val = getValueFromKhash32(stats_info->genome_base_with_N_coverage, bins[i]);
			if (i==0) { val -= stats_info->cov_stats->total_Ns_bases; }		// need to remove all Ns

			float percent = calculatePercentage(val, total_genome_non_Ns_bases);
			if (i==0) fprintf(out_fp, "Bases with %d coverage:,%"PRIu32",(%.2f%%)\n", bins[i], val, percent);
			if (i>0)  fprintf(out_fp, "Bases with %d+ coverage:,%"PRIu32",(%.2f%%)\n", bins[i], val, percent);
		}

		fprintf(out_fp, "Max Coverage:,%"PRIu32",(%"PRIu32")\n", stats_info->cov_stats->max_coverage, stats_info->cov_stats->base_with_max_coverage);
	    fprintf(out_fp, "\n");
		fprintf(out_fp, "Coverage Histogram for Whole Genome (may look weird if target regions overlap...)\n");

		for (i=1; i<=1000; i++) {
            k_iter = kh_put(m32, stats_info->genome_cov_histogram, i, &ret);
            if (ret == 1)
                kh_value(stats_info->genome_cov_histogram, k_iter) = 0;

			fprintf(out_fp, "%"PRIu32",", kh_key(stats_info->genome_cov_histogram, k_iter)-1);
        }
		fprintf(out_fp, "1000,\n");

		for (k_iter=kh_begin(stats_info->genome_cov_histogram); k_iter!=kh_end(stats_info->genome_cov_histogram); k_iter++) {
			if (kh_exist(stats_info->genome_cov_histogram, k_iter)) {
				fprintf(out_fp, "%"PRIu32",", kh_value(stats_info->genome_cov_histogram, k_iter));
			}
		}
	    fprintf(out_fp, "\n");

		fclose(out_fp);
		free(file_name);
	}

	// Now we need to process target information if target bed file is provided
	if (TARGET_FILE_PROVIDED) {
		// First we need to calculate the coverage for median, this is like N50 for sequencing
		sum = 0;
		for (k_iter=0; k_iter!=kh_end(stats_info->target_coverage_for_median); k_iter++) {
			if (kh_exist(stats_info->target_coverage_for_median, k_iter)) {

				// Divya: 29826 number of Ns in our VCrome + PKv2 design; will need to change if design changes
				//if (sum >= (stats_info->cov_stats->total_targeted_bases - 29826)/2) {
				if (sum >= stats_info->cov_stats->total_targeted_bases/2) {
					stats_info->cov_stats->median_target_coverage = k_iter--;
					break;
				} else {
					sum += kh_value(stats_info->target_coverage_for_median, k_iter);
				}
			}
        }

		average_coverage = (double)stats_info->cov_stats->total_target_coverage/(double)stats_info->cov_stats->total_targeted_bases;
		file_name = calloc(strlen(user_inputs->out_file) + 30, sizeof(char));
		strcpy(file_name, user_inputs->out_file);
		strcat(file_name, ".CoverageReport.csv");
		FILE *trt_fp = fopen(file_name, "w");	// trt_fp: target_fp

		//printf("Before output general information\n");
		outputGeneralInfo(trt_fp, stats_info, average_coverage, 2);
		//printf("After output general information\n");

		fprintf(trt_fp, "Target Stats\n");

		float percent = calculatePercentage(stats_info->cov_stats->hit_target_count, stats_info->cov_stats->total_targets);
        fprintf(trt_fp, "Targets Hit:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->hit_target_count, percent);

		percent = calculatePercentage(stats_info->cov_stats->hit_target_buffer_only_count, stats_info->cov_stats->total_targets);
        fprintf(trt_fp, "Target Buffers Hit:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->hit_target_buffer_only_count, percent);
        fprintf(trt_fp, "Total Targets:,%"PRIu32"\n", stats_info->cov_stats->total_targets);
        fprintf(trt_fp, "Non target regions with high coverage:,%"PRIu32"\n", stats_info->cov_stats->non_target_good_hits);
        fprintf(trt_fp, "Base Stats\n");
        fprintf(trt_fp, "Bases Targeted:,%"PRIu32"\n", stats_info->cov_stats->total_targeted_bases);
        fprintf(trt_fp, "Buffer Bases:,%"PRIu32"\n", stats_info->cov_stats->total_buffer_bases);

		for(i=0; i<13; i++) {
            uint32_t val = getValueFromKhash32(stats_info->targeted_base_with_N_coverage, bins[i]);

            float percent = calculatePercentage(val, stats_info->cov_stats->total_targeted_bases);
            if (i==0) fprintf(trt_fp, "Bases with %d coverage:,%"PRIu32",(%.2f%%)\n", bins[i], val, percent);
            if (i>0)  fprintf(trt_fp, "Bases with %d+ coverage:,%"PRIu32",(%.2f%%)\n", bins[i], val, percent);
        }

		//printf("After the base count for 1000\n");

		fprintf(trt_fp, "\n");
		fprintf(trt_fp, "Coverage Histogram (may look weird if target regions overlap...)\n");

		for (i=1; i<=1000; i++) {
			k_iter = kh_put(m32, stats_info->target_cov_histogram, i, &ret);
			if (ret == 1)
				kh_value(stats_info->target_cov_histogram, k_iter) = 0;

			fprintf(trt_fp, "%"PRIu32",", kh_key(stats_info->target_cov_histogram, k_iter)-1);
		}
    	fprintf(trt_fp, "1000,\n");

		for (k_iter=kh_begin(stats_info->target_cov_histogram); k_iter!=kh_end(stats_info->target_cov_histogram); k_iter++) {
			if (kh_exist(stats_info->target_cov_histogram, k_iter)) {
				fprintf(trt_fp, "%"PRIu32",", kh_value(stats_info->target_cov_histogram, k_iter));
			}
		}
		fprintf(trt_fp, "\n");
		printf("After histogram plotting\n");

		fprintf(trt_fp, "Target and region coverage plot\n");
		fprintf(trt_fp, "Position,5'count,3'count\n");
		for(i = 20; i <= PRIMER_SIZE; i+=20) {
        	fprintf(trt_fp, "%d,%"PRIu32",%"PRIu32"\n", i, stats_info->five_prime[PRIMER_SIZE-(i-1)-1], stats_info->three_prime[i-1]);
        }

		//fputs("%tar-Pos,count\n", trt_fp);
		fprintf(trt_fp, "target_Pos,count\n");
        for(i = 0; i < 101; i+=2) {
            fprintf(trt_fp, "%d,%"PRIu32"\n", i, stats_info->target_coverage[i]);
        }

		//printf("Before free the file name second time\n");
		fclose(trt_fp);
		free(file_name);
	}
}

void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage, uint8_t type) {
	fprintf(fp, "Version: %s\n", VERSION);
    if (type == 2) fprintf(fp, "BUFFER size:,%d\n", BUFFER);
    fprintf(fp, "Read Stats\n");
     
    fprintf(fp, "Total Reads Produced:,%"PRIu32"\n", stats_info->cov_stats->total_reads_produced);

    uint64_t yield = stats_info->cov_stats->read_length * (uint64_t) stats_info->cov_stats->total_reads_produced;
    fprintf(fp, "Total Yield Produced:,%d,%"PRIu64"\n", stats_info->cov_stats->read_length, yield);

    yield = stats_info->cov_stats->read_length * (uint64_t) (stats_info->cov_stats->total_reads_aligned - stats_info->cov_stats->total_duplicate_reads);
    fprintf(fp, "Total Unique Yield Produced:,%"PRIu64"\n", yield);

	float percent = calculatePercentage(stats_info->cov_stats->total_duplicate_reads,stats_info->cov_stats->total_reads_aligned);
    fprintf(fp, "Duplicate Reads:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->total_duplicate_reads, percent);

    percent = calculatePercentage(stats_info->cov_stats->total_reads_aligned,stats_info->cov_stats->total_reads_produced);
    fprintf(fp, "Total Reads Aligned:,%"PRIu32",(%.2f%%)", stats_info->cov_stats->total_reads_aligned, percent);

	fprintf(fp, ",reads paired:,%"PRIu32, stats_info->cov_stats->total_reads_paired);
    fprintf(fp, ",reads paired with mapped mates:,%"PRIu32"\n", stats_info->cov_stats->total_paired_reads_with_mapped_mates);

	if (type == 2) {
		percent = calculatePercentage(stats_info->cov_stats->in_buffer_read_hit_count, stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Aligned Reads On-Buffer:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->in_buffer_read_hit_count, percent);

		percent = calculatePercentage(stats_info->cov_stats->on_target_read_hit_count, stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Aligned Reads On-Target:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->on_target_read_hit_count, percent);
	}

    fprintf(fp, "Average Coverage:,-,(%.2f)\n", average_coverage);
    if (type == 1)
		fprintf(fp, "Median Coverage:,-,(%d)\n", stats_info->cov_stats->median_genome_coverage);

    if (type == 2)
		fprintf(fp, "Median Coverage:,-,(%d)\n", stats_info->cov_stats->median_target_coverage);

	if (type == 2) {
		uint32_t reads_hit_target_or_buffer = stats_info->cov_stats->on_target_read_hit_count + stats_info->cov_stats->in_buffer_read_hit_count;
		percent = calculatePercentage(reads_hit_target_or_buffer, stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Reads that hit target or buffer:,%"PRIu32",(%.2f%%)\n", reads_hit_target_or_buffer, percent);

		fprintf(fp, "Total Aligned Reads (expected):,%"PRIu32"\n", stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Total Aligned Reads (calculated):,%"PRIu32"\n", stats_info->cov_stats->on_target_read_hit_count + stats_info->cov_stats->in_buffer_read_hit_count + stats_info->cov_stats->off_target_read_hit_count);
	}
}

void produceOffTargetWigFile(Chromosome_Tracking *chrom_tracking, char *chrom_id, Bed_Info *target_bed_info, User_Input *user_inputs, Stats_Info *stats_info) {
	int32_t i, j;
	int32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	if (idx == -1) return;

	FILE *wp = fopen(user_inputs->wgs_wig_file, "a");
	fprintf(wp, "track type=wiggle_0 name=%s\n", user_inputs->bam_file);

	for(i=0; i<target_bed_info->size; i++) {

		// here we are going to erase everything that is on-target and in-buffer and set them to 0
		if (strcmp(chrom_id, target_bed_info->coords[i].chrom_id) == 0) {
			uint32_t start = target_bed_info->coords[i].start;
			uint32_t stop  = target_bed_info->coords[i].end;
			for(j=start-500; j<=stop+500; j++) {
				if (j < 0 || j > chrom_tracking->chromosome_lengths[idx])
					continue;

				chrom_tracking->coverage[idx][j] = 0;
			}
		}
	}

	// once captured area + BUFFER is initialized to 0, we want to check if there is any coverage > 20 in off-target regions
	for(i=0; i<chrom_tracking->chromosome_lengths[idx]; i++) {
		if (chrom_tracking->coverage[idx][i] > 20) {
			j = i;
			stats_info->cov_stats->non_target_good_hits += 1;

			// right side 
			while(i < chrom_tracking->chromosome_lengths[idx] && chrom_tracking->coverage[idx][i] > 0)
				i++;

			// left side
			while(j > 0 && chrom_tracking->coverage[idx][i] > 0)
				j--;

			// now write to the off target wig file
			fprintf(wp, "fixedStep chrom_id=%s start=%"PRId32" step=1\n", chrom_id, j);
			uint32_t h = j;
			for(h=j; h<i; h++)
				fprintf(wp,"%"PRIu32"\n", chrom_tracking->coverage[idx][h]);
		}
	}
	fclose(wp);
}
