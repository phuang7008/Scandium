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

void writeCoverage(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info,Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions) {

	// First, we need to find the index that is used to track current chromosome chrom_id
	int32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	if (idx == -1) return;

	// For whole genome (WGS) related outputs
	//
	uint32_t i=0;
	if(user_inputs->wgs_coverage) {

		if (user_inputs->Write_WGS) {
			FILE *wgs_cov_fasta_fp = fopen(user_inputs->wgs_cov_file, "a");
			printf("Whole Genome cov.fasta output for chrom id %s is on\n", chrom_id);

			// no need to add newline as the next line will take care of it!
			//
			fprintf(wgs_cov_fasta_fp, ">%s 1 %"PRIu32, chrom_id, chrom_tracking->chromosome_lengths[idx]-1);

			for (i = 0; i < chrom_tracking->chromosome_lengths[idx]; i++) {
				if(i%100==0) fputc('\n', wgs_cov_fasta_fp);
				fprintf(wgs_cov_fasta_fp, "%d ", chrom_tracking->coverage[idx][i]);

				addBaseStats(stats_info, (uint32_t) chrom_tracking->coverage[idx][i], 0, 1);
			}

			fputc('\n', wgs_cov_fasta_fp);
			fclose(wgs_cov_fasta_fp);
		} else {
			// no need to write WGS cov.fasta output
			//
			for (i = 0; i < chrom_tracking->chromosome_lengths[idx]; i++)
				addBaseStats(stats_info, (uint32_t) chrom_tracking->coverage[idx][i], 0, 1);
		}
    }

	// if the target bed file is available, we will need to handle it here and write the results to cov.fasta file
	//
	if (TARGET_FILE_PROVIDED) {
		FILE * cov_fp = fopen(user_inputs->capture_cov_file, "a");
		FILE * missed_target_fp = fopen(user_inputs->missed_targets_file, "a");

		// now need to report those Capture regions with low or too high coverages
		// NOTE: the bed format is different here, the end position is included!
		//
		FILE *capture_low_x_fp    = fopen(user_inputs->capture_low_cov_file, "a");
		FILE *capture_high_x_fp   = fopen(user_inputs->capture_high_cov_file, "a");
		FILE *capture_all_site_fp = fopen(user_inputs->capture_all_site_file, "a");
		FILE *capture_range_fp    = fopen(user_inputs->capture_range_file, "a");

		for(i = 0; i < target_info->size; i++) {

			if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
				continue;

			stats_info->cov_stats->total_targets++;
			int32_t start = target_info->coords[i].start;		// need to use signed value as it will get -1 later!!!
			int32_t end = target_info->coords[i].end;
			int length = end - start + 1;
			bool collect_target_cov = length > 99 ? true : false ;  // TODO: is it the min length of the exon? Check.

			//if (start == 89623861) {
			//	printf("Debug stopped!\n");
			//}

			int32_t j=0;
			if (collect_target_cov) {
				for(j = 0; j < PRIMER_SIZE; j++) {
					int32_t k = start - j;		// k could go negative, so it is a signed integer
					if ( k < 0 || (end + j) >= chrom_tracking->chromosome_lengths[idx])
						continue;

					stats_info->five_prime[j]  += chrom_tracking->coverage[idx][k];
					stats_info->three_prime[j] += chrom_tracking->coverage[idx][end+j];
				}
			}

			//In the original java code, the following is used for the pc and pc2 definition
			//short pc[101], pc2[101];
			//After discussion with Divya and Qiaoyan, here is what we understand.
			//pc and pc2 is used to setup the 'target_coverage' variable which is defined as an array of 101 in the origain java code
			//here I will use percentage_bin and percentage_count_per_bin
			//
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
				//
				fprintf(missed_target_fp, "%s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, start, end);
				bool hit = false;
				for (j = start - user_inputs->target_buffer_size; j < start && !hit; j++) {
					if (j >= chrom_tracking->chromosome_lengths[idx])
						continue;

					if ( chrom_tracking->coverage[idx][j] > 0 )
						hit = true;
				}

				if (hit)
					stats_info->cov_stats->hit_target_buffer_only_count += 1;
			}

			// For All Sites Report
			// the end position is not included based on the bed format
			//
			produceCaptureAllSitesReport(start, length-1, chrom_tracking, chrom_id, user_inputs, capture_all_site_fp, inter_genic_regions, intronic_regions, exon_regions);

			// For low coverage and high coverage Report
			//
			writeLow_HighCoverageReport(start, length, chrom_tracking, chrom_id, user_inputs, capture_low_x_fp, capture_high_x_fp, inter_genic_regions, intronic_regions, exon_regions, 2);

			// For range block coverage information for graphing
			//
			writeCoverageRanges(start, length, chrom_tracking, idx, user_inputs, capture_range_fp);
		}

		fclose(capture_low_x_fp);
		fclose(capture_high_x_fp);
		fclose(capture_all_site_fp);
		fclose(capture_range_fp);

		fclose(cov_fp);
		fclose(missed_target_fp);
	}
}

void produceCaptureAllSitesReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char * chrom_id, User_Input *user_inputs, FILE *fh_all_sites, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions) {
	uint32_t i=0;
	uint64_t cov_total=0;
	int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);

	//if (chrom_idx == -1 || chrom_idr == -1) return;
	if (chrom_idx == -1) return;

	for (i = begin; i < begin+length; i++) {
		cov_total += chrom_tracking->coverage[chrom_idx][i];
	}

	uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(length) + 0.5);
	//uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float) length);
	// Here I need to add 1 to the length as needed to compare with Java version of ExCID
	//
	fprintf(fh_all_sites, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\t", chrom_id, begin, begin+length, length+1, ave_coverage);

	if (user_inputs->annotation_on) {
		char *annotation = getRegionAnnotation(begin, begin+length-1, chrom_id, inter_genic_regions, intronic_regions, exon_regions, 2);
		if (annotation) {
			fprintf(fh_all_sites, "%s", annotation);
			free(annotation); 
			annotation = NULL; 
		} else {
			fprintf(stderr, "No annotation for %s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, begin, begin+length);
		}
	} else {
		fprintf(fh_all_sites, ".\t.\t.\t.\t.\t.\t.\n");
	}
}

// This function is used to write low coverage bases/regions and 
// high coverage bases/regions for WGS (whole genome) using its own thread.
//
void writeAnnotations(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions) {
	// First, we need to find the index that is used to track current chromosome chrom_id
	// 
    int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);

	if(user_inputs->wgs_coverage) {
	    // NOTE: the bed format is different here, the end position is included!
		//
	    FILE *wgs_low_x_fp  = fopen(user_inputs->wgs_low_cov_file, "a");
    	FILE *wgs_high_x_fp = fopen(user_inputs->wgs_high_cov_file, "a");

	    writeLow_HighCoverageReport(0, chrom_tracking->chromosome_lengths[chrom_idx], chrom_tracking, chrom_id, user_inputs, wgs_low_x_fp, wgs_high_x_fp, inter_genic_regions, intronic_regions, exon_regions, 1);

        fclose(wgs_low_x_fp);
        fclose(wgs_high_x_fp);
	}
}

uint32_t writeLow_HighCoverageReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char *chrom_id, User_Input *user_inputs,FILE *fh_low, FILE *fh_high, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions, uint8_t type) {
	// for debugging
	//if (strcmp(chrom_tracking->chromosome_ids[chrom_idx], "1") == 0) return i;
	
	// First, we need to find the index that is used to track current chromosome chrom_id
	int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	int32_t chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, intronic_regions);
	if (chrom_idx == -1 || chrom_idr == -1) return 0;
	
	uint32_t i=0;

	for (i = begin; i < begin+length; i++) {
		uint32_t start=0, end=0;
		uint64_t cov_total=0;

		// for low coverage checking and writing
		//
		if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] < user_inputs->low_coverage_to_report) {
			start = i;

			while(i < begin+length && chrom_tracking->coverage[chrom_idx][i] < user_inputs->low_coverage_to_report) {
				//if (start == 10568 || start == 11118) {
            	//	printf("beginning coverage is %"PRIu32"\n", chrom_tracking->coverage[chrom_idx][i]);
        		//}
				cov_total += chrom_tracking->coverage[chrom_idx][i];
				i++;
			}
			end = i;

			if (start < end) {
				//float ave_coverage = (float)cov_total / (float)(end - start);
				uint32_t ave_coverage = (uint32_t) (((float)cov_total / (float)(end - start)) + 0.5);
				fprintf(fh_low, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\t", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage);
			}

			// generate the annotation here
			//
			if (user_inputs->annotation_on) {
				//char *annotation = getRegionAnnotation(start, end-1, chrom_id, inter_genic_regions, intronic_regions, exon_regions, type);
				char *annotation = getRegionAnnotation(start, end-1, chrom_id, inter_genic_regions, intronic_regions, exon_regions, 2);
				if (annotation != NULL) {
					fprintf(fh_low, "%s", annotation);
					free(annotation);
					annotation = NULL;
				} else {
					fprintf(fh_low, "%s", "Annotation not available\n");
				}
			} else {
				fprintf(fh_low, ".\t.\t.\t.\t.\t.\t.\n");
			}
        }
		//fflush(fh_low);

        // for high coverage
		start = 0;
		end = 0;
		cov_total = 0;

		// now for anything above the High coverage
		//
        if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] >= user_inputs->high_coverage_to_report) {
            start = i;

            while( i < begin+length && chrom_tracking->coverage[chrom_idx][i] >= user_inputs->high_coverage_to_report) {
					cov_total += chrom_tracking->coverage[chrom_idx][i];
					i++;
            }
            end = i;

			if (start < end) {
				//float ave_coverage = (float)cov_total / (float)(end - start);
				uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);
				fprintf(fh_high, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\t", chrom_tracking->chromosome_ids[chrom_idx], start, end-1, end-start, ave_coverage);
			}

			// generate the annotation here
			//
            if (user_inputs->annotation_on) {
				//char *annotation = getRegionAnnotation(start, end-1, chrom_id, inter_genic_regions, intronic_regions, exon_regions, type);
				char *annotation = getRegionAnnotation(start, end-1, chrom_id, inter_genic_regions, intronic_regions, exon_regions, 2);
				if (annotation != NULL) {
					fprintf(fh_high, "%s", annotation);
					free(annotation);
					annotation = NULL;
				} else {
					fprintf(fh_high, "%s", "Annotation not available\n");
				}
            } else {
				fprintf(fh_high, ".\t.\t.\t.\t.\t.\t.\n");
			}
        }
    }

	return i;
}

// type 1 mean doing speed up seearch (as we need to remember the previous search index), while type 2 won't 
//
char * getRegionAnnotation(uint32_t start, uint32_t end, char *chrom_id, Regions_Skip_MySQL *inter_genic_regions, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions, uint8_t type) {
	char *annotation = calloc(50, sizeof(char));
	strcpy(annotation, ".");

	//if (start >= 1167104 && end <=1167198) {
	//	printf("Stop\n");
	//}

	int32_t index_intronic_location=-1, index_exon_location=-1;
	uint32_t tmp_loc_idx = 0;
	int32_t chrom_idr = -1;

	if (type == 1) {
		// to speed up the search, we need to check to see if the region is still part of the previous search lower bound
		// First check exon_region
		//
		chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, exon_regions);
		if (chrom_idr == exon_regions->prev_search_chrom_index && exon_regions->prev_search_loc_index > 0) {
			tmp_loc_idx = exon_regions->prev_search_loc_index;

			if (verifyIndex(exon_regions, start, end, chrom_idr, tmp_loc_idx)) {
				if (strlen(exon_regions->exon_info[chrom_idr][tmp_loc_idx]) > 1) {
					uint16_t str_len_needed = strlen(exon_regions->gene[chrom_idr][tmp_loc_idx]) + strlen(exon_regions->Synonymous[chrom_idr][tmp_loc_idx]) + strlen(exon_regions->prev_genes[chrom_idr][tmp_loc_idx]) + strlen(exon_regions->exon_info[chrom_idr][tmp_loc_idx]) + 50;
					char *tmp = realloc(annotation, str_len_needed * sizeof(char));
		        	if (!tmp) {
			        	fprintf(stderr, "Memory re-allocation for string failed in checkExonRegion\n");
				    	exit(1);
					}
					annotation = tmp;
			    	sprintf(annotation, "%s\t%s\t%s\t%s\n", exon_regions->gene[chrom_idr][tmp_loc_idx], exon_regions->prev_genes[chrom_idr][tmp_loc_idx], exon_regions->Synonymous[chrom_idr][tmp_loc_idx], exon_regions->exon_info[chrom_idr][tmp_loc_idx]);
					return annotation;
				}
			}
		}

		// now speed search for intronic regions
		chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, intronic_regions);
		if (chrom_idr == intronic_regions->prev_search_chrom_index && intronic_regions->prev_search_loc_index > 0) {
			tmp_loc_idx = intronic_regions->prev_search_loc_index;

			if (verifyIndex(intronic_regions, start, end, chrom_idr, tmp_loc_idx)) {
				uint16_t str_len_needed = strlen(intronic_regions->gene[chrom_idr][tmp_loc_idx]) + strlen(intronic_regions->Synonymous[chrom_idr][tmp_loc_idx]) + strlen(intronic_regions->prev_genes[chrom_idr][tmp_loc_idx]) + 50;
	            char *tmp = realloc(annotation, str_len_needed * sizeof(char));
		        if (!tmp) {
			        fprintf(stderr, "Memory re-allocation for string failed in checkExonRegion\n");
				    exit(1);
				}
	            annotation = tmp;
		        sprintf(annotation, "%s\t%s\t%s\t.\t.\t.\t.\n", intronic_regions->gene[chrom_idr][tmp_loc_idx], intronic_regions->prev_genes[chrom_idr][tmp_loc_idx], intronic_regions->Synonymous[chrom_idr][tmp_loc_idx]);
			    return annotation;
			}
		}
    }

	// Now, check if the region locates at exon area
    if (index_exon_location == -1) {
		chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, exon_regions);
		if (chrom_idr >= 0) {
			if (type == 1) {
				index_exon_location = checkExonRegion(exon_regions, start, end, chrom_idr, &annotation, exon_regions->prev_search_loc_index);
			} else {
				index_exon_location = checkExonRegion(exon_regions, start, end, chrom_idr, &annotation, 0);
			}

		    if (index_exon_location != -1) {
				if (type == 1) {
					exon_regions->prev_search_loc_index = index_exon_location;
					exon_regions->prev_search_chrom_index = chrom_idr;
				}

				strcat(annotation, "\n");
				return annotation;
			}
        }
    }

	//Next, check if the region locates at the intronic area
	if (index_exon_location == -1) {
		chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, intronic_regions);

        if (chrom_idr != -1 && index_intronic_location == -1) {
			if (type == 1) {
				index_intronic_location = checkIntronicRegion(intronic_regions, start, end, chrom_idr, &annotation, intronic_regions->prev_search_loc_index);
			} else {
				index_intronic_location = checkIntronicRegion(intronic_regions, start, end, chrom_idr, &annotation, 0);
			}
		
            if (index_intronic_location != -1) {
				if (type == 1) {
					intronic_regions->prev_search_loc_index   = index_intronic_location;
					intronic_regions->prev_search_chrom_index = chrom_idr;
				}

				strcat(annotation, "\t.\t.\t.\t.\n");
                return annotation;
            }
        }
    }

	// if the search comes here, it has to be inter-genic region. And there is not need to go further!
    strcpy(annotation, ".\t.\t.\t.\t.\t.\t.\n");
	return annotation;
}

// this is just a wrapper function to help run the inner function writeCoverageRanges()
//
void coverageRangeInfoForGraphing(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info) {
	// First, we need to find the index that is used to track current chromosome chrom_id
	// 
	int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	FILE *wgs_range_fp = fopen(user_inputs->wgs_range_file, "a");
	
	writeCoverageRanges(0, chrom_tracking->chromosome_lengths[chrom_idx], chrom_tracking, chrom_idx, user_inputs, wgs_range_fp);

	fclose(wgs_range_fp);
}

// here we stick to the bed format for the output
//
void writeCoverageRanges(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, uint16_t chrom_idx, User_Input *user_inputs, FILE *fh_range) {
	// NOTE: the bed format is different here, the end position is included!
	//
	uint32_t i=0;
	for (i=begin; i<begin+length; i++) {
		uint32_t start=0, end=0;
		uint64_t cov_total=0;
		
		// now handle within range block. Inclusive for both boundaries [lower_bound, upper_bound]
		//	lower_bound  ========================== upper_bound
		//	    the base to be checked  - 
		//
		if (i < begin+length && (chrom_tracking->coverage[chrom_idx][i] >= user_inputs->lower_bound)
				&& (chrom_tracking->coverage[chrom_idx][i] <= user_inputs->upper_bound) ) {
			start = i;

			while( i < begin+length && (chrom_tracking->coverage[chrom_idx][i] >= user_inputs->lower_bound)
					&& (chrom_tracking->coverage[chrom_idx][i] <= user_inputs->upper_bound) ) {
					cov_total += chrom_tracking->coverage[chrom_idx][i];
					i++;
			}
			end = i;

			if (start <= end) {
				uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);
				fprintf(fh_range, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);
			}
		}

		// now handle anything below the lower bound (here I will treat them as a single group, ExCID's way)
		//                   lower_bound ========================== upper_bound
		//  the base to be checked  -
		//
		start = 0;
		end = 0;
		cov_total = 0;

		if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] < user_inputs->lower_bound) {
			start = i;

			while (i < begin+length && chrom_tracking->coverage[chrom_idx][i] < user_inputs->lower_bound) {
				cov_total += chrom_tracking->coverage[chrom_idx][i];
				i++;
            }
            end = i;

			if (start < end) {
				uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);
				fprintf(fh_range, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);
			}
		}

		// now handle anything beyond the upper bound using gVCF approach
		// here we will use gVCF approach with the block max <= (100%+1) * min
		// the implementation is based on the group discussion 
		//	lower_bound ========================== upper_bound
		//							the base to be checked  -
		//
		start = 0;
		end = 0;
		cov_total = 0;

		if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] > user_inputs->upper_bound) {
			start = i;

			uint32_t min_cov=10000000, max_cov=0;

			while (i < begin+length && chrom_tracking->coverage[chrom_idx][i] > user_inputs->upper_bound) {
				// check using gVCF approach
				//
				if (min_cov > chrom_tracking->coverage[chrom_idx][i])
					min_cov = chrom_tracking->coverage[chrom_idx][i];

				if (max_cov < chrom_tracking->coverage[chrom_idx][i])
					max_cov = chrom_tracking->coverage[chrom_idx][i];

				// check if max_cov is more than user_inputs->target_buffer_size*100% of min
				//
				if ( min_cov <= max_cov && max_cov <= min_cov*user_inputs->target_buffer_size) {
					cov_total += chrom_tracking->coverage[chrom_idx][i];
					i++;
				} else {
					break;
				}
			}

			end=i;
			if (start < end) {
				uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);
				fprintf(fh_range, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);
			}
		}

		// Need to decrease i here as the loop will increase it again
		//
		i--;
	}
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

	if (cov_val >= 6) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 6, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 6, 1);
    }

	if (cov_val >= 10) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 10, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 10, 1);
    }

	if (cov_val >= 11) {
        if (target == 1) addValueToKhashBucket32(stats_info->targeted_base_with_N_coverage, 11, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->genome_base_with_N_coverage, 11, 1);
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
	uint16_t bins[15] = { 0, 1, 5, 6, 10, 11, 15, 20, 30, 40, 50, 60, 100, 500, 1000 };
	khiter_t k_iter;

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

		// open WGS coverage summary report file handle
		FILE *out_fp = fopen(user_inputs->wgs_cov_report, "a");

		average_coverage = (double) stats_info->cov_stats->total_genome_coverage/ (double) total_genome_non_Ns_bases;
		outputGeneralInfo(out_fp, stats_info, average_coverage, user_inputs, 1);

		fprintf(out_fp, "Base Stats\n");
	    fprintf(out_fp, "Bases Targeted:,%"PRIu64"\n", total_genome_non_Ns_bases);

		for(i=0; i<15; i++) {
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
	}

	// Now we need to process target information if target bed file is provided
	if (TARGET_FILE_PROVIDED) {
		// First we need to calculate the coverage for median, this is like N50 for sequencing
		sum = 0;
		for (k_iter=0; k_iter!=kh_end(stats_info->target_coverage_for_median); k_iter++) {
			if (kh_exist(stats_info->target_coverage_for_median, k_iter)) {

				// Divya: 29826 number of Ns in our VCrome + PKv2 design; will need to change if design changes
				//if (sum >= (stats_info->cov_stats->total_targeted_bases - 29826)/2) 
				if (sum >= stats_info->cov_stats->total_targeted_bases/2) {
					stats_info->cov_stats->median_target_coverage = k_iter--;
					break;
				} else {
					sum += kh_value(stats_info->target_coverage_for_median, k_iter);
				}
			}
        }

		average_coverage = (double)stats_info->cov_stats->total_target_coverage/(double)stats_info->cov_stats->total_targeted_bases;
		FILE *trt_fp = fopen(user_inputs->capture_cov_report, "a");	// trt_fp: target_fp

		//printf("Before output general information\n");
		outputGeneralInfo(trt_fp, stats_info, average_coverage, user_inputs, 2);
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

		for(i=0; i<15; i++) {
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
	}
}

void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage, User_Input *user_inputs, uint8_t type) {
    if (type == 2) fprintf(fp, "BUFFER size:,%d\n", user_inputs->target_buffer_size);
    fprintf(fp, "Read Stats\n");
     
    fprintf(fp, "Total Reads Produced:,%"PRIu32"\n", stats_info->cov_stats->total_reads_produced);

    uint64_t yield = stats_info->cov_stats->read_length * (uint64_t) stats_info->cov_stats->total_reads_produced;
    fprintf(fp, "Total Yield Produced:,%d,%"PRIu64"\n", stats_info->cov_stats->read_length, yield);

    yield = stats_info->cov_stats->read_length * (uint64_t) (stats_info->cov_stats->total_reads_aligned - stats_info->cov_stats->total_duplicate_reads);
    fprintf(fp, "Total Unique Yield Produced:,%"PRIu64"\n", yield);

	float percent = calculatePercentage(stats_info->cov_stats->total_duplicate_reads,stats_info->cov_stats->total_reads_aligned);
    fprintf(fp, "Duplicate Reads:,%"PRIu32",(%.2f%%)\n", stats_info->cov_stats->total_duplicate_reads, percent);

	//if (user_inputs->remove_supplementary_alignments) {
		percent = calculatePercentage(stats_info->cov_stats->total_supplementary_reads,stats_info->cov_stats->total_reads_aligned);
		fprintf(fp, "Supplementary Reads: %"PRIu32",(%.2f%%)\n", stats_info->cov_stats->total_supplementary_reads, percent);
	//}

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
				if (j < 0 || j >= chrom_tracking->chromosome_lengths[idx])
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

void calculateGenePercentageCoverage(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Databases *dbs, Low_Coverage_Genes *low_cov_genes) {
	// find out the index that is used to track current chromosome id
	int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	if (chrom_idx == -1) return;

	// print low_cov_genes for debugging
	//printLowCoverageGeneStructure(low_cov_genes);

	// if the target bed file is available, we will need to calculate percentage of gene bases that are covered
	uint32_t i, j;
	if (TARGET_FILE_PROVIDED) {

		for(i = 0; i < target_info->size; i++) {
			if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
				continue;

			uint32_t begin  = target_info->coords[i].start;
			uint32_t finish = target_info->coords[i].end;
			//int length = finish - begin + 1;

			for (j=begin; j<finish; j++) {

				// now we need to find out if this exon has low coverage bases and how many
				uint32_t start=0, end=0;

				// for low coverage
				if (j < finish && chrom_tracking->coverage[chrom_idx][j] < user_inputs->low_coverage_to_report) {
					start = j;

					while(j < finish && chrom_tracking->coverage[chrom_idx][j] < user_inputs->low_coverage_to_report) {
						j++;
					}
					end = j;
					if (start == end) continue;

					// now update the Low_Coverage_Genes variable low_cov_genes
					//
					produceGenePercentageCoverageInfo(start, end, chrom_id, low_cov_genes);
					//produceGenePercentageCoverageInfo(start, end-1, chrom_id, low_cov_genes);
				}
			}
		}
	}
}

void outputGenePercentageCoverage(char *chrom_id, Bed_Info *target_info, User_Input *user_inputs, Low_Coverage_Genes *low_cov_genes, Transcript_Coverage *transcript_cov, Databases *dbs) {
	
    // if the target bed file is available, we will need to calculate percentage of gene bases that are covered
    uint32_t i;

    if (TARGET_FILE_PROVIDED) {
		// open file handle for writing/appending
		FILE *gene_pct_fp = fopen(user_inputs->low_cov_gene_pct_file, "a");
		FILE *exon_pct_fp = fopen(user_inputs->low_cov_exon_pct_file, "a");

		// we need to sort the Gene_Coverage before we could use them as it will make it faster
		qsort(low_cov_genes->gene_coverage, low_cov_genes->total_size, sizeof(Gene_Coverage), &compare);
		//for(i = 0; i < low_cov_genes->total_size; i++) {
		//	printf("%s\t%"PRIu32"\n", low_cov_genes->gene_coverage[i].gene_name, low_cov_genes->gene_coverage[i].exon_start);
		//}

		// The following is to output the low coverage report for individual exon/cds (or target). Need to find out which one to use
		//
        for(i = 0; i < low_cov_genes->total_size; i++) {
			uint32_t cds_target_start = low_cov_genes->gene_coverage[i].cds_target_start;
			uint32_t cds_target_end   = low_cov_genes->gene_coverage[i].cds_target_end;

			// As we don't stick to the official bed format in the bed file, sometimes the cds_target_end == cds_target_start (such as, SNPs)
			// In this case, we get the divide by 0 error. To overcome this, we add 1 to the case where cds_target_end == cds_target_start
			//
			float pct = 0.000;
			if (cds_target_end == cds_target_start) {
				pct = (float) ((cds_target_end - cds_target_start + 1 - low_cov_genes->gene_coverage[i].num_of_low_cov_bases) * 100) / (float) (cds_target_end - cds_target_start + 1);
			} else {
				pct = (float) ((cds_target_end - cds_target_start - low_cov_genes->gene_coverage[i].num_of_low_cov_bases) * 100) / (float) (cds_target_end - cds_target_start);
			}

			// Note: there are several cases need to be handled.
			// 1). SNP: if it is for SNPs, the exon_id is set to -1 in the database, but the output will be '_'
			// 2). eMerge: three exons/CDSs are not designed but will be output into the final results, the exon_id is set to -exon_id (such as -13)
			// 3). regular genes output
			//
			if (low_cov_genes->gene_coverage[i].low_cov_regions) {
				if (low_cov_genes->gene_coverage[i].exon_id == -1) {
					fprintf(exon_pct_fp, "%s\t%s\t%s\t_\t%"PRIu32"\t%"PRIu32"\t%0.3f%%\t%s\n", chrom_id, low_cov_genes->gene_coverage[i].gene_symbol, low_cov_genes->gene_coverage[i].gene_name, cds_target_start, cds_target_end, pct, low_cov_genes->gene_coverage[i].low_cov_regions);
				} else if (low_cov_genes->gene_coverage[i].exon_id < -1 ) { 
					fprintf(exon_pct_fp, "%s\t%s\t%s\tcds_%d\t%"PRIu32"\t%"PRIu32"\t0.000%%\tNOT_Designed\n", chrom_id, low_cov_genes->gene_coverage[i].gene_symbol, low_cov_genes->gene_coverage[i].gene_name, low_cov_genes->gene_coverage[i].exon_id, cds_target_start, cds_target_end );	
				} else{	
					fprintf(exon_pct_fp, "%s\t%s\t%s\tcds_%"PRIu16"\t%"PRIu32"\t%"PRIu32"\t%0.3f%%\t%s\n", chrom_id, low_cov_genes->gene_coverage[i].gene_symbol, low_cov_genes->gene_coverage[i].gene_name, low_cov_genes->gene_coverage[i].exon_id, cds_target_start, cds_target_end, pct, low_cov_genes->gene_coverage[i].low_cov_regions);
				}

			} else {
				if (low_cov_genes->gene_coverage[i].exon_id == -1) {
					fprintf(exon_pct_fp, "%s\t%s\t%s\t_\t%"PRIu32"\t%"PRIu32"\t%0.3f%%\t%s\n", chrom_id, low_cov_genes->gene_coverage[i].gene_symbol, low_cov_genes->gene_coverage[i].gene_name, cds_target_start, cds_target_end, pct, ".");
				} else if (low_cov_genes->gene_coverage[i].exon_id < -1 ) {
					fprintf(exon_pct_fp, "%s\t%s\t%s\tcds_%d\t%"PRIu32"\t%"PRIu32"\t0.000%%\tNOT_Designed\n", chrom_id, low_cov_genes->gene_coverage[i].gene_symbol, low_cov_genes->gene_coverage[i].gene_name, low_cov_genes->gene_coverage[i].exon_id, cds_target_start, cds_target_end );	
				} else {
					fprintf(exon_pct_fp, "%s\t%s\t%s\tcds_%"PRIu16"\t%"PRIu32"\t%"PRIu32"\t%0.3f%%\t%s\n", chrom_id, low_cov_genes->gene_coverage[i].gene_symbol, low_cov_genes->gene_coverage[i].gene_name, low_cov_genes->gene_coverage[i].exon_id, cds_target_start, cds_target_end, pct, ".");
				}
			}
		}

		// The following is to output the low coverage report for each refseq/gene
		// Since multiple refseq could be start at the same position, we need to record previous gene/refseq name
		char *prev_gene_name = calloc(50, sizeof(char));
		char *prev_gene_symbol = NULL;
		strcpy(prev_gene_name, ",");

		// declare prev_cds_size and low base count here as they need to be accumulated
		uint32_t prev_cds_size=0, prev_gene_low_bases=0;
		
		// Some capture targarts come from the same exon, so we need to skip it when accumulate the refseq length
		uint16_t exon_count=0;
		float pct;

		uint32_t transcript_counter=0;

		for (i = 0; i < low_cov_genes->total_size; i++) {

			// Check to see if refseq transcript name has changed???
			//
			if (strcmp(prev_gene_name, low_cov_genes->gene_coverage[i].gene_name) != 0) {
				// begining of a new refseq
				// need to output current gene/refseq info if it is not the first item!
				//
				//if (strcmp(prev_gene_name, "NM_001077") == 0) {
				//	printf("stop!\n");
				//}

				if (i > 0 && prev_gene_symbol) {
					pct = (float) ((prev_cds_size - prev_gene_low_bases) * 100) / (float) (prev_cds_size);
					fprintf(gene_pct_fp, "%s\t%s\t%s\t%"PRIu32"\t%"PRIu16"\t%0.2f\n", chrom_id, prev_gene_symbol, prev_gene_name, prev_cds_size, exon_count, pct);

					// now we need to update the Transcript Coverage Percentage information
					transcript_cov->transcript_cov_pct[transcript_counter].gene_symbol = calloc(strlen(prev_gene_symbol)+1, sizeof(char));
					strcpy(transcript_cov->transcript_cov_pct[transcript_counter].gene_symbol, prev_gene_symbol);

					transcript_cov->transcript_cov_pct[transcript_counter].gene_name = calloc(strlen(prev_gene_name)+1, sizeof(char));
					strcpy(transcript_cov->transcript_cov_pct[transcript_counter].gene_name, prev_gene_name);

					transcript_cov->transcript_cov_pct[transcript_counter].gene_cov_percentage = pct;
					transcript_counter++;
				}

				// free 'previous' info for char*
				free(prev_gene_name);
				prev_gene_name=NULL;

				free(prev_gene_symbol);
				prev_gene_symbol=NULL;

				exon_count = 0;

				// now reset all the 'previous' information
				prev_gene_name = calloc(strlen(low_cov_genes->gene_coverage[i].gene_name)+1, sizeof(char));
				strcpy(prev_gene_name, low_cov_genes->gene_coverage[i].gene_name);

				prev_gene_symbol = calloc(strlen(low_cov_genes->gene_coverage[i].gene_symbol)+1, sizeof(char));
				strcpy(prev_gene_symbol, low_cov_genes->gene_coverage[i].gene_symbol);

				//prev_cds_size = low_cov_genes->gene_coverage[i].cds_target_end - low_cov_genes->gene_coverage[i].cds_target_start + 1;	
				prev_cds_size = low_cov_genes->gene_coverage[i].cds_length;	
				prev_gene_low_bases = low_cov_genes->gene_coverage[i].num_of_low_cov_bases;
				exon_count++;

			} else {
				// accumulate values
				//prev_cds_size += low_cov_genes->gene_coverage[i].cds_target_end - low_cov_genes->gene_coverage[i].cds_target_start + 1;
				prev_cds_size = low_cov_genes->gene_coverage[i].cds_length;
				prev_gene_low_bases += low_cov_genes->gene_coverage[i].num_of_low_cov_bases;
				exon_count++;
			}
		}

		// output the last one
		if (prev_gene_symbol) {
			pct = (float) ((prev_cds_size - prev_gene_low_bases) * 100) / (float) (prev_cds_size);
			fprintf(gene_pct_fp, "%s\t%s\t%s\t%"PRIu32"\t%"PRIu16"\t%0.2f\n", chrom_id, prev_gene_symbol, prev_gene_name, prev_cds_size, exon_count, pct);

			// now we need to update the Transcript Coverage Percentage information
            transcript_cov->transcript_cov_pct[transcript_counter].gene_symbol = calloc(strlen(prev_gene_symbol)+1, sizeof(char));
            strcpy(transcript_cov->transcript_cov_pct[transcript_counter].gene_symbol, prev_gene_symbol);

            transcript_cov->transcript_cov_pct[transcript_counter].gene_name = calloc(strlen(prev_gene_name)+1, sizeof(char));
            strcpy(transcript_cov->transcript_cov_pct[transcript_counter].gene_name, prev_gene_name);

            transcript_cov->transcript_cov_pct[transcript_counter].gene_cov_percentage = pct;
		}

		// now we need to write everything regarding transcript coverage percentage into a file
		writeTranscriptCoveragePercentage(transcript_cov, user_inputs);

		fclose(gene_pct_fp);
		fclose(exon_pct_fp);
	}
}

void writeTranscriptCoveragePercentage(Transcript_Coverage *transcript_cov,  User_Input *user_inputs) {
	if (transcript_cov->num_of_genes == 0)
		return;

	// we need to sort the transcripts on gene_symbol first
	qsort(transcript_cov->transcript_cov_pct, transcript_cov->num_of_genes, sizeof(Transcript_Coverage_Percentage), &compare2);
	
	FILE *pct_cov_fp  = fopen(user_inputs->low_cov_transcript_file, "a");

	char *prev_gene_symbol=calloc(50, sizeof(char));
	strcpy(prev_gene_symbol, "");

	float ave_cov_pct;
	uint32_t num_of_transcripts, i;

	for (i=0; i<transcript_cov->num_of_genes; i++) {

		if (strcmp(prev_gene_symbol, transcript_cov->transcript_cov_pct[i].gene_symbol) != 0) {
			if (strlen(prev_gene_symbol) != 0) {
				// first we need to finish the previous gene symbol/gene_name transcript
				ave_cov_pct = ave_cov_pct / (float) num_of_transcripts;
				fprintf(pct_cov_fp, "\t%0.2f\n", ave_cov_pct);
			}

			// now we need to report the current transcript
			fprintf(pct_cov_fp, "%s\t%s(%0.2f)", transcript_cov->transcript_cov_pct[i].gene_symbol, transcript_cov->transcript_cov_pct[i].gene_name, transcript_cov->transcript_cov_pct[i].gene_cov_percentage);

			strcpy(prev_gene_symbol, transcript_cov->transcript_cov_pct[i].gene_symbol);
			ave_cov_pct = transcript_cov->transcript_cov_pct[i].gene_cov_percentage;
			num_of_transcripts = 1;
		} else {
			fprintf(pct_cov_fp, "; %s(%0.2f)", transcript_cov->transcript_cov_pct[i].gene_name, transcript_cov->transcript_cov_pct[i].gene_cov_percentage);
			ave_cov_pct += transcript_cov->transcript_cov_pct[i].gene_cov_percentage;
			num_of_transcripts++;
		}
	}

	// finish the last one
	ave_cov_pct = ave_cov_pct / (float) num_of_transcripts;
	fprintf(pct_cov_fp, "\t%.2f\n", ave_cov_pct);

	fclose(pct_cov_fp);
}
