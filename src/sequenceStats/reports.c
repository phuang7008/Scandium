/*
 * =====================================================================================
 *
 *       Filename:  reports.c
 *
 *    Description:  the detailed implementation of report.h
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

void writeCoverage(char *chrom_id, Bed_Info **target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL **exon_regions) {

    // First, we need to find the index that is used to track current chromosome chrom_id
    int32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
    if (idx == -1) return;

    // For whole genome (WGS) related outputs
    //
    uint32_t i=0;
    if(user_inputs->wgs_coverage) {

        if (user_inputs->Write_WGS_cov_fasta) {
            FILE *wgs_cov_fasta_fp = fopen(user_inputs->wgs_cov_file, "a");
            printf("Whole Genome cov.fasta output for chromosome id %s is on\n", chrom_id);

            // no need to add newline as the next line will take care of it!
            //
            fprintf(wgs_cov_fasta_fp, ">%s 1 %"PRIu32, chrom_id, chrom_tracking->chromosome_lengths[idx]-1);

            for (i = 0; i < chrom_tracking->chromosome_lengths[idx]; i++) {
                if(i%100==0) fputc('\n', wgs_cov_fasta_fp);
                fprintf(wgs_cov_fasta_fp, "%d ", chrom_tracking->coverage[idx][i]);

                addBaseStats(stats_info, (uint32_t) chrom_tracking->coverage[idx][i], 0, 1, 0);
            }

            fputc('\n', wgs_cov_fasta_fp);
            fclose(wgs_cov_fasta_fp);
        } else {
            // no need to write WGS cov.fasta output
            //
            for (i = 0; i < chrom_tracking->chromosome_lengths[idx]; i++)
                addBaseStats(stats_info, (uint32_t) chrom_tracking->coverage[idx][i], 0, 1, 0);
        }
    }

    // if the target bed file is available, we will need to handle it here and write the results to cov.fasta file
    //
    if (TARGET_FILE_PROVIDED) {

        for (i=0; i<user_inputs->num_of_target_files; i++) {
            generateCaptureStats(chrom_id, target_info[i], chrom_tracking, user_inputs, stats_info, idx, i);

            // now need to report those Capture regions with low or too high coverages
            // NOTE: the bed format is different here, the end position is included!
            //
            produceReportsOnThresholds(chrom_id, target_info[i], chrom_tracking, user_inputs, intronic_regions, exon_regions, i);
        }
    }
}

void generateCaptureStats(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, int32_t chrom_idx, uint8_t target_file_index) {
    FILE * capture_cov_fp = NULL;
    if (user_inputs->Write_Capture_cov_fasta) capture_cov_fp = fopen(user_inputs->capture_cov_files[target_file_index], "a");

    uint32_t i;
    for(i = 0; i < target_info->size; i++) {

        if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
            continue;

        stats_info->capture_cov_stats[target_file_index]->total_targets++;
        int32_t start = target_info->coords[i].start;       // need to use signed value as it will get -1 later!!!
        int32_t end = target_info->coords[i].end;
        int length = end - start;
        bool collect_target_cov = length > 99 ? true : false ;  // TODO: is it the min length of the exon? Check.

        int32_t j=0;
        if (collect_target_cov) {
            for(j = 0; j < PRIMER_SIZE; j++) {
                int32_t k = start - j;      // k could go negative, so it is a signed integer
                if ( k < 0 || (uint32_t) (end + j) >= chrom_tracking->chromosome_lengths[chrom_idx])
                    continue;

                stats_info->capture_cov_stats[target_file_index]->five_prime[j]  += chrom_tracking->coverage[chrom_idx][k];
                stats_info->capture_cov_stats[target_file_index]->three_prime[j] += chrom_tracking->coverage[chrom_idx][end+j];
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
        if (user_inputs->Write_Capture_cov_fasta)
            fprintf(capture_cov_fp, ">%s %"PRIu32" %"PRIu32"\n", chrom_id, start, end);

        bool space_it = false;
        if(end - start > 10000) space_it = true;

        for(j = 0; j < length; j++) {
            if ((uint32_t) j+start >= chrom_tracking->chromosome_lengths[chrom_idx])
                continue;

            // enter a new line after every 100 bases
            if (user_inputs->Write_Capture_cov_fasta && space_it && j%100 == 0) fputc('\n', capture_cov_fp);

            uint32_t cov = chrom_tracking->coverage[chrom_idx][j+start];
            addBaseStats(stats_info, cov, 1, 0, target_file_index);

            if (!target_hit && (cov > 0))
                target_hit = true;

            // output to the cov.fasta file
            //
            if (user_inputs->Write_Capture_cov_fasta) fprintf(capture_cov_fp, "%d ", cov);

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
        if (user_inputs->Write_Capture_cov_fasta) fputc('\n', capture_cov_fp);

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
                stats_info->capture_cov_stats[target_file_index]->target_coverage[j] += percentage_bin[j];
        }

        if (target_hit) {
            stats_info->capture_cov_stats[target_file_index]->hit_target_count += 1;
        } else {
            // need to write to the missed target file
            //
            //fprintf(missed_target_fp, "%s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, start, end);
            bool hit = false;
            for (j = start - user_inputs->target_buffer_size; j < start && !hit; j++) {
                if ((uint32_t) j >= chrom_tracking->chromosome_lengths[chrom_idx])
                    continue;

                if ( chrom_tracking->coverage[chrom_idx][j] > 0 )
                    hit = true;
            }

            if (hit)
                stats_info->capture_cov_stats[target_file_index]->hit_target_buffer_only_count += 1;          
        }
    }
    
    if (user_inputs->Write_Capture_cov_fasta && capture_cov_fp != NULL) fclose(capture_cov_fp);
}

void produceReportsOnThresholds(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL **exon_regions, uint8_t target_file_index) {

    FILE *capture_low_x_fp  = NULL;
    FILE *capture_high_x_fp = NULL;
    FILE *capture_all_site_fp = NULL;

    capture_low_x_fp    = fopen(user_inputs->capture_low_cov_files[target_file_index], "a");
    capture_all_site_fp = fopen(user_inputs->capture_all_site_files[target_file_index], "a");

    if (user_inputs->above_10000_on) capture_high_x_fp = fopen(user_inputs->capture_low_cov_files[target_file_index], "a");

    uint32_t i;
    for(i = 0; i < target_info->size; i++) {
        if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
            continue;

        int32_t start = target_info->coords[i].start;       // need to use signed value as it will get -1 later!!!
        int32_t end = target_info->coords[i].end;
        int length = end - start;

        // For All Sites Report
        // the end position is not included based on the bed format
        //
        if (USER_DEFINED_DATABASE) {
            produceCaptureAllSitesReport(start, length, chrom_tracking, chrom_id, &capture_all_site_fp, NULL, exon_regions[target_file_index]);
        } else {
            produceCaptureAllSitesReport(start, length, chrom_tracking, chrom_id, &capture_all_site_fp, intronic_regions, exon_regions[target_file_index]);
        }

        // For low coverage and high coverage Reporf (user_inputs->Write_Capture_cov_fasta && capture_cov_fp != NULL) fclose(capture_cov_fp);t
        //
        if (USER_DEFINED_DATABASE) {
            writeLow_HighCoverageReport(start, length, chrom_tracking, chrom_id, user_inputs, capture_low_x_fp, capture_high_x_fp, NULL, exon_regions[target_file_index], 2);
        } else {
            writeLow_HighCoverageReport(start, length, chrom_tracking, chrom_id, user_inputs, capture_low_x_fp, capture_high_x_fp, intronic_regions, exon_regions[target_file_index], 2);
        }

        // For uniformity data coverage information for graphing
        //
        //writeCoverageRanges(start, length, chrom_tracking, idx, user_inputs, capture_uniformity_fp);
    }

    fclose(capture_low_x_fp);
    fclose(capture_all_site_fp);
    if (user_inputs->above_10000_on && capture_high_x_fp != NULL) fclose(capture_high_x_fp);
}


// this is for Capture only, so I don't have to check for the annotation_on flag as it is for WGS only
//
void produceCaptureAllSitesReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char * chrom_id, FILE **fh_all_sites, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions) {
    uint32_t i=0;
    uint64_t cov_total=0;
    int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);

    if (chrom_idx == -1) return;

    for (i = begin; i < begin+length; i++) {
        // check if it passes the end of the chromosome
        //
        if (i <= chrom_tracking->chromosome_lengths[chrom_idx]) {
            cov_total += chrom_tracking->coverage[chrom_idx][i];
        } else {
            //continue;
            break;
        }
    }

    uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(length) + 0.5);

    char *annotation = NULL;
    annotation = getRegionAnnotation(begin, begin+length, chrom_id, intronic_regions, exon_regions, 2);
    if (annotation == NULL) {
        annotation = calloc(30, sizeof(char));
        strcpy(annotation, ".\t.\t.\t.\t.\t.\t.\t.\n");
    }
    
    fprintf(*fh_all_sites, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\t%s", chrom_id, begin, begin+length, length, ave_coverage, annotation);
    //fflush(*fh_all_sites);

    if (annotation != NULL) {
        free(annotation);
        annotation = NULL;
    }
}

// This function is used to write low coverage bases/regions and 
// high coverage bases/regions for WGS (whole genome) using its own thread.
//
void writeAnnotations(char *chrom_id, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions) {
    // First, we need to find the index that is used to track current chromosome chrom_id
    // 
    int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);

    if(user_inputs->wgs_coverage) {
        // NOTE: the bed format is different here, the end position is included!
        //
        FILE *wgs_low_x_fp  = fopen(user_inputs->wgs_low_cov_file, "a");
        FILE *wgs_high_x_fp = NULL;
           if (user_inputs->above_10000_on) wgs_high_x_fp = fopen(user_inputs->wgs_high_cov_file, "a");

        writeLow_HighCoverageReport(0, chrom_tracking->chromosome_lengths[chrom_idx], chrom_tracking, chrom_id, user_inputs, wgs_low_x_fp, wgs_high_x_fp, intronic_regions, exon_regions, 1);

        fclose(wgs_low_x_fp);
        if (user_inputs->above_10000_on && wgs_high_x_fp != NULL) fclose(wgs_high_x_fp);
    }
}

// type is the mode of process: 1 for WGS, 2 for Capture
//
uint32_t writeLow_HighCoverageReport(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, char *chrom_id, User_Input *user_inputs, FILE *fh_low, FILE *fh_high, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions, uint8_t type) {
    // for debugging
    //if (strcmp(chrom_tracking->chromosome_ids[chrom_idx], "1") == 0) return i;
    
    // First, we need to find the index that is used to track current chromosome chrom_id
    int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
    int32_t chrom_idr = -1;
       if (exon_regions != NULL) {
        chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, exon_regions);
        //if (chrom_idr == -1) return 0;
    }

    if (chrom_idx == -1) return 0;
    
    uint32_t i=0;    // need to be signed as 'coverage' is a signed value

    for (i = begin; i < begin+length; i++) {
        int32_t start=0, end=0;
        uint64_t cov_total=0;

        // for low coverage checking and writing
        //
        if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] < (int32_t) user_inputs->low_coverage_to_report) {
            start = i;

            while(i < begin+length && chrom_tracking->coverage[chrom_idx][i] < (int32_t) user_inputs->low_coverage_to_report) {
                //if (start == 10568 || start == 11118) {
                //    printf("beginning coverage is %"PRIu32"\n", chrom_tracking->coverage[chrom_idx][i]);
                //}
                cov_total += chrom_tracking->coverage[chrom_idx][i];
                i++;
            }
            end = i;

            if (start < end) {
                //float ave_coverage = (float)cov_total / (float)(end - start);
                uint32_t ave_coverage = (uint32_t) (((float)cov_total / (float)(end - start)) + 0.5);
                fprintf(fh_low, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\t", chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);
            }

            // generate the annotation here
            //
            produceAnnotation(type, chrom_idr, chrom_id, start, end, user_inputs, fh_low, intronic_regions, exon_regions);
        }
        //fflush(fh_low);

        // For High coverage
        // Skip if fh_high is NULL (this is for USER_DEFINED_DATABASE)
        //
        if (fh_high == NULL)
            continue;

        start = 0;
        end = 0;
        cov_total = 0;

        // now for anything above the High coverage
        //
        if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] >= (int32_t) user_inputs->high_coverage_to_report) {
            start = i;

            while( i < begin+length && chrom_tracking->coverage[chrom_idx][i] >= (int32_t) user_inputs->high_coverage_to_report) {
                    cov_total += chrom_tracking->coverage[chrom_idx][i];
                    i++;
            }
            end = i;

            if (start < end) {
                //float ave_coverage = (float)cov_total / (float)(end - start);
                uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);
                fprintf(fh_high, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\t", chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);
            }

            // generate the annotation here
            //
            produceAnnotation(type, chrom_idr, chrom_id, start, end, user_inputs, fh_high, intronic_regions, exon_regions);
        }
    }

    return i;
}

// type: mode of processing. 1 for WGS, 2 for Capture
//
void produceAnnotation(uint8_t type, int32_t chrom_idr, char *chrom_id, uint32_t start, uint32_t end, User_Input *user_inputs, FILE *out_fh, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions) {

    bool good = false;

    if (type == 1) {
        // For WGS
        //
        if (user_inputs->wgs_annotation_on) {
            if (chrom_idr != -1 && exon_regions != NULL)
                good = fetchAnnotation(chrom_id, start, end, out_fh, intronic_regions, exon_regions);
        }
    } else {
        // Capture
        //
        if (chrom_idr != -1 && exon_regions != NULL) 
            good = fetchAnnotation(chrom_id, start, end, out_fh, intronic_regions, exon_regions);
    }

    if (!good)
        fprintf(out_fh, ".\t.\t.\t.\t.\t.\t.\t.\n");
}

bool fetchAnnotation(char *chrom_id, uint32_t start, uint32_t end, FILE *out_fh, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions) {
    char *annotation = getRegionAnnotation(start, end, chrom_id, intronic_regions, exon_regions, 2);
    if (annotation != NULL) { 
        fprintf(out_fh, "%s", annotation);
        free(annotation);
        annotation = NULL;
        return true;
    }

    return false;
}

// type 1 mean doing speed up seearch (as we need to remember the previous search index), while type 2 won't 
//
char * getRegionAnnotation(uint32_t start, uint32_t end, char *chrom_id, Regions_Skip_MySQL *intronic_regions, Regions_Skip_MySQL *exon_regions, uint8_t type) {
    char *annotation = NULL;
    //dynamicStringAllocation(".", &annotation);

    // just return ".\t.\t.\t.\t.\t.\t." if exon_regions is NULL
    //
    if (exon_regions == NULL) {
        dynamicStringExpansion(".\t.\t.\t.\t.\t.\t.\t.\n", &annotation);
        return annotation;
    }

    int32_t index_intronic_location=-1, index_exon_location=-1;
    uint32_t tmp_loc_idx = 0;
    int32_t chrom_idr = -1;

    if (type == 1) {
        // to speed up the search, we need to check to see if the region is still part of the previous search lower bound
        // First check exon_region
        //
        chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, exon_regions);
        if (chrom_idr >= 0 && (uint32_t) chrom_idr == exon_regions->prev_search_chrom_index && exon_regions->prev_search_loc_index > 0) {
            tmp_loc_idx = exon_regions->prev_search_loc_index;

            if (verifyIndex(exon_regions, start, end, chrom_idr, tmp_loc_idx)) {
                if (strlen(exon_regions->exon_info[chrom_idr][tmp_loc_idx]) > 1) {
                    uint16_t str_len_needed = strlen(exon_regions->gene[chrom_idr][tmp_loc_idx]) + strlen(exon_regions->Synonymous[chrom_idr][tmp_loc_idx]) + strlen(exon_regions->prev_genes[chrom_idr][tmp_loc_idx]) + strlen(exon_regions->exon_info[chrom_idr][tmp_loc_idx]) + 50;
                    char *tmp = realloc(annotation, str_len_needed * sizeof(char));
                    if (!tmp) {
                        fprintf(stderr, "ERROR: Memory re-allocation for string failed in checkExonRegion\n");
                        exit(EXIT_FAILURE);
                    }
                    annotation = tmp;
                    sprintf(annotation, "%s\t%s\t%s\t%s\n", exon_regions->gene[chrom_idr][tmp_loc_idx], exon_regions->prev_genes[chrom_idr][tmp_loc_idx], exon_regions->Synonymous[chrom_idr][tmp_loc_idx], exon_regions->exon_info[chrom_idr][tmp_loc_idx]);
                    return annotation;
                }
            }
        }

        chrom_idr = locateChromosomeIndexForRegionSkipMySQL(chrom_id, intronic_regions);
        if (chrom_idr >= 0 && (uint32_t) chrom_idr == intronic_regions->prev_search_chrom_index && intronic_regions->prev_search_loc_index > 0) {
            tmp_loc_idx = intronic_regions->prev_search_loc_index;

            if (verifyIndex(intronic_regions, start, end, chrom_idr, tmp_loc_idx)) {
                uint16_t str_len_needed = strlen(intronic_regions->gene[chrom_idr][tmp_loc_idx]) + strlen(intronic_regions->Synonymous[chrom_idr][tmp_loc_idx]) + strlen(intronic_regions->prev_genes[chrom_idr][tmp_loc_idx]) + 50;
                char *tmp = realloc(annotation, str_len_needed * sizeof(char));
                if (!tmp) {
                    fprintf(stderr, "ERROR: Memory re-allocation for string failed in checkExonRegion\n");
                    exit(EXIT_FAILURE);
                }
                annotation = tmp;
                sprintf(annotation, "%s\t%s\t%s\t.\t.\t.\t.\t.\n", intronic_regions->gene[chrom_idr][tmp_loc_idx], intronic_regions->prev_genes[chrom_idr][tmp_loc_idx], intronic_regions->Synonymous[chrom_idr][tmp_loc_idx]);
                return annotation;
            }
        }
    }

    // Now, check if the region locates at exon area (non-speed way)
    //
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

                dynamicStringExpansion("\n", &annotation);
                return annotation;
            }
        }
    }

    //Next, check if the region locates at the intronic area
    //
    if (index_exon_location == -1) {
        if (intronic_regions == NULL) {
            dynamicStringExpansion(".\t.\t.\t.\t.\t.\t.\t.\n", &annotation);
            return annotation;
        }

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

                dynamicStringExpansion("\t.\t.\t.\t.\t.\n", &annotation);
                return annotation;
            }
        }
    }

    // if the search comes here, it has to be inter-genic region. And there is not need to go further!
    dynamicStringExpansion(".\t.\t.\t.\t.\t.\t.\t.\n", &annotation);
    return annotation;
}

// this is just a wrapper function to help run the inner function writeCoverageRanges()
//
void coverageRangeInfoForGraphing(char *chrom_id, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs) {
    // First, we need to find the index that is used to track current chromosome chrom_id
    // 
    int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
    FILE *wgs_uniformity_fp = fopen(user_inputs->wgs_uniformity_file, "a");
    
    writeCoverageRanges(0, chrom_tracking->chromosome_lengths[chrom_idx], chrom_tracking, chrom_idx, user_inputs, wgs_uniformity_fp);

    fclose(wgs_uniformity_fp);
}

// here we stick to the bed format for the output
//
void writeCoverageRanges(uint32_t begin, uint32_t length, Chromosome_Tracking *chrom_tracking, uint16_t chrom_idx, User_Input *user_inputs, FILE *fh_uniformity) {
    // NOTE: the bed format is different here, the end position is included!
    //
    uint32_t i=0;
    for (i=begin; i<begin+length; i++) {
        uint32_t start=0, end=0;
        uint64_t cov_total=0;
        
        // now handle within uniformity range data. Inclusive for both boundaries [lower_bound, upper_bound]
        //    lower_bound  ========================== upper_bound
        //        the base to be checked  - 
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
                fprintf(fh_uniformity, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);
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
                fprintf(fh_uniformity, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);
            }
        }

        // now handle anything beyond the upper bound using gVCF approach
        // here we will use gVCF approach with the block max <= (100%+1) * min
        // the implementation is based on the group discussion 
        //    lower_bound ========================== upper_bound
        //                            the base to be checked  -
        //
        start = 0;
        end = 0;
        cov_total = 0;

        if (i < begin+length && chrom_tracking->coverage[chrom_idx][i] > user_inputs->upper_bound) {
            start = i;

            int32_t min_cov=10000000, max_cov=0;    // need to be signed as 'coverage' is a signed int

            while (i < begin+length && chrom_tracking->coverage[chrom_idx][i] > user_inputs->upper_bound) {
                // check using gVCF approach
                //
                if (min_cov > chrom_tracking->coverage[chrom_idx][i])
                    min_cov = chrom_tracking->coverage[chrom_idx][i];

                if (max_cov < chrom_tracking->coverage[chrom_idx][i])
                    max_cov = chrom_tracking->coverage[chrom_idx][i];

                // check if max_cov is more than user_inputs->gVCF_percentage*100% of min
                //
                if ( min_cov <= max_cov && max_cov <= min_cov*user_inputs->gVCF_percentage) {
                    cov_total += chrom_tracking->coverage[chrom_idx][i];
                    i++;
                } else {
                    break;
                }
            }

            end=i;
            if (start < end) {
                uint32_t ave_coverage = (uint32_t) ((float)cov_total / (float)(end - start) + 0.5);
                fprintf(fh_uniformity, "%s\t%"PRIu32"\t%"PRIu32"\t%d\t%"PRIu32"\n", chrom_tracking->chromosome_ids[chrom_idx], start, end, end-start, ave_coverage);
            }
        }

        // Need to decrease i here as the loop will increase it again
        //
        i--;
    }
}

void addBaseStats(Stats_Info *stats_info, uint32_t cov_val, uint8_t target, uint8_t wgs, uint8_t target_file_index) {
    // need to check and update max coverage information
    //
    if (wgs == 1) {
        if (cov_val > stats_info->wgs_cov_stats->wgs_max_coverage) {
            stats_info->wgs_cov_stats->wgs_max_coverage = cov_val;
            stats_info->wgs_cov_stats->base_with_wgs_max_coverage = 1;
        } else if (cov_val == stats_info->wgs_cov_stats->wgs_max_coverage) {
            stats_info->wgs_cov_stats->base_with_wgs_max_coverage += 1;
        }
    }

    if (target == 1) {
        if (cov_val > stats_info->capture_cov_stats[target_file_index]->target_max_coverage) {
            stats_info->capture_cov_stats[target_file_index]->target_max_coverage = cov_val;
            stats_info->capture_cov_stats[target_file_index]->base_with_target_max_coverage = 1;
        } else if (cov_val == stats_info->capture_cov_stats[target_file_index]->target_max_coverage) {
            stats_info->capture_cov_stats[target_file_index]->base_with_target_max_coverage += 1;
        }
    }

    // for histogram only
    uint32_t tmp_val = cov_val;
    if (tmp_val > 1000) tmp_val = 1000;
    if (target == 1) stats_info->capture_cov_stats[target_file_index]->target_cov_histogram[tmp_val]++;
    if (wgs == 1)    stats_info->wgs_cov_stats->genome_cov_histogram[tmp_val]++;

    if (target == 1) stats_info->capture_cov_stats[target_file_index]->total_target_coverage += (uint64_t) cov_val;
    if (wgs == 1)    stats_info->wgs_cov_stats->total_genome_coverage += (uint64_t) cov_val;

    if (cov_val > 0) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 1, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 1, 1);
    } else if (cov_val == 0) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 0, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 0, 1);
    }

    if (cov_val >= 5) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 5, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 5, 1);
    }

    if (cov_val >= 6) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 6, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 6, 1);
    }

    if (cov_val >= 10) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 10, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 10, 1);
    }

    if (cov_val >= 11) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 11, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 11, 1);
    }

    if (cov_val >= 15) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 15, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 15, 1);
    }

    if (cov_val >= 20) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 20, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 20, 1);
    }

    if (cov_val >= 30) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 30, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 30, 1);
    }

    if (cov_val >= 40) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 40, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 40, 1);
    }

    if (cov_val >= 50) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 50, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 50, 1);
    }

    if (cov_val >= 60) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 60, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 60, 1);
    }

    if (cov_val >= 70) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 70, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 70, 1);
    }

    if (cov_val >= 100) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 100, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 100, 1);
    }

    if (cov_val >= 500) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 500, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 500, 1);
    }

    if (cov_val >= 1000) {
        if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_base_with_N_coverage, 1000, 1);
        if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, 1000, 1);
    }

    if (target == 1) addValueToKhashBucket32(stats_info->capture_cov_stats[target_file_index]->target_coverage_for_median, cov_val, 1);
    if (wgs == 1)    addValueToKhashBucket32(stats_info->wgs_cov_stats->genome_coverage_for_median, cov_val, 1);
}

void writeWGSReports(Stats_Info *stats_info, User_Input *user_inputs) {
    if (stats_info->read_cov_stats->total_reads_aligned == 0) {
        fprintf(stderr, "ERROR: No reads aligned. Aborting.\n");
        exit(EXIT_FAILURE);
    }

    uint32_t non_duplicate_reads = stats_info->read_cov_stats->total_reads_aligned - stats_info->read_cov_stats->total_duplicate_reads;
    if (non_duplicate_reads == 0) {
        fprintf(stderr, "ERROR: All reads are duplicates. Aborting.\n");
        exit(EXIT_FAILURE);
    }

    uint64_t sum=0;
    int32_t i=0;
    double average_coverage=0.0;
    uint16_t bins[16] = { 0, 1, 5, 6, 10, 11, 15, 20, 30, 40, 50, 60, 70, 100, 500, 1000 };
    khiter_t k_iter;
    //uint32_t cov_bin_size = 20000;
    uint32_t coverage_bins[20000] = {0};    // initialize array contents to 0

    if (user_inputs->wgs_coverage) {

        //Do not consider the Ns for Median calculation.
        uint64_t total_genome_non_Ns_bases = stats_info->wgs_cov_stats->total_genome_bases - stats_info->wgs_cov_stats->total_Ns_bases;

        for (k_iter=0; k_iter!=kh_end(stats_info->wgs_cov_stats->genome_coverage_for_median); k_iter++) {
            if (kh_exist(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter)) {

                //fprintf(stderr, "%d\t%d\t%"PRIu32"\n", k_iter, kh_key(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter), kh_value(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter));
                if (kh_key(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter) < 20000)
                    coverage_bins[kh_key(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter)] =
                        kh_value(stats_info->wgs_cov_stats->genome_coverage_for_median, k_iter);
            }
        }

        for (i=0; i<20000; i++) {

            if (sum >= (total_genome_non_Ns_bases/2)) {
                stats_info->wgs_cov_stats->median_genome_coverage = i--;
                break;
            }else{
                sum += coverage_bins[i];
                if (i == 0)
                    sum -= stats_info->wgs_cov_stats->total_Ns_bases;
            }
        }

        // open WGS coverage summary report file handle
        FILE *out_fp = fopen(user_inputs->wgs_cov_report, "a");

        average_coverage = (double) stats_info->wgs_cov_stats->total_genome_coverage/ (double) total_genome_non_Ns_bases;
        outputGeneralInfo(out_fp, stats_info, average_coverage);
        
        fprintf(out_fp, "WGS Max_Coverage\t%"PRIu32"\n", stats_info->wgs_cov_stats->wgs_max_coverage);
        fprintf(out_fp, "Number_of_Bases_with_WGS_Max_Coverage\t%"PRIu32"\n", stats_info->wgs_cov_stats->base_with_wgs_max_coverage);

        fprintf(out_fp, "#Base_Stats\n");
        fprintf(out_fp, "Total_Genome_Base_Targeted_w/o_Ns\t%"PRIu64"\n", total_genome_non_Ns_bases);
        fprintf(out_fp, "Total_Ns_Bases_in_Ns_Regions\t%"PRIu32"\n", stats_info->wgs_cov_stats->total_Ns_bases);
        fprintf(out_fp, "Total_Mapped_Bases\t%"PRIu64"\n", stats_info->wgs_cov_stats->total_mapped_bases);
        fprintf(out_fp, "Total_Uniquely_Aligned_Bases\t%"PRIu64"\n", stats_info->wgs_cov_stats->total_uniquely_aligned_bases);

        float percent = calculatePercentage64(stats_info->wgs_cov_stats->total_overlapped_bases, stats_info->wgs_cov_stats->total_mapped_bases);
        fprintf(out_fp, "Total_Overlapped_Aligned_Bases\t%"PRIu32"\n", stats_info->wgs_cov_stats->total_overlapped_bases);
        fprintf(out_fp, "PCT_Overlapped_Aligned_Bases\t%0.2f%%\n", percent);

        percent = calculatePercentage64(stats_info->wgs_cov_stats->base_quality_20, stats_info->wgs_cov_stats->total_mapped_bases);
        fprintf(out_fp, "Aligned_Q20_Bases\t%"PRIu64"\n", stats_info->wgs_cov_stats->base_quality_20);
        fprintf(out_fp, "PCT_Aligned_Q20_Bases\t%0.2f%%\n", percent);

        percent = calculatePercentage64(stats_info->wgs_cov_stats->base_quality_30, stats_info->wgs_cov_stats->total_mapped_bases);
        fprintf(out_fp, "Aligned_Q30_Bases\t%"PRIu64"\n", stats_info->wgs_cov_stats->base_quality_30);
        fprintf(out_fp, "PCT_Aligned_Q30_Bases\t%0.2f%%\n", percent);

        fprintf(out_fp, "Median_Coverage\t%d\n", stats_info->wgs_cov_stats->median_genome_coverage);
        fprintf(out_fp, "Mode_Coverage_For_Uniformity\t%d\n", stats_info->wgs_cov_stats->mode);
        fprintf(out_fp, "Uniformity_Primary_Autosome_Only\t%.3f\n", stats_info->wgs_cov_stats->uniformity_metric_primary_autosome_only);
        //fprintf(out_fp, "Uniformity_Primary_Autosome_plus_X_and_Y\t%.3f\n", stats_info->wgs_cov_stats->uniformity_metric_all_primary);
        //fprintf(out_fp, "Uniformity_Primary_Autosome_plus_Alt_Decoys_HLA\t%.3f\n", stats_info->wgs_cov_stats->uniformity_metric_autosome_only);
        //fprintf(out_fp, "Uniformity_Primary_Autosome_plus_X_Y_Alt_Decoys_HLA_(Everything)\t%.3f\n", stats_info->wgs_cov_stats->uniformity_metric_all);

        for(i=0; i<16; i++) {
            uint32_t val = getValueFromKhash32(stats_info->wgs_cov_stats->genome_base_with_N_coverage, bins[i]);
            if (i==0) { val -= stats_info->wgs_cov_stats->total_Ns_bases; }        // need to remove all Ns

            percent = calculatePercentage32_64(val, total_genome_non_Ns_bases);
            fprintf(out_fp, "Bases_with_%dx_coverage\t%"PRIu32"\n", bins[i], val);
            fprintf(out_fp, "PCT_of_Bases_with_%dx_coverage\t%.2f%%\n", bins[i], percent);
        }

        fprintf(out_fp, "\n");
        fprintf(out_fp, "#Coverage_Frequency_Distribution_for_Whole_Genome\n");

        for (i=0; i<=1000; i++)
            fprintf(out_fp, "%d,", i);

        fprintf(out_fp, "\n");
        fprintf(out_fp, ">>");

        for (i=0; i<=1000; i++)
            fprintf(out_fp, "%"PRIu32",", stats_info->wgs_cov_stats->genome_cov_histogram[i]);

        fprintf(out_fp, "\n");
        fclose(out_fp);
    }
}

void writeCaptureReports(Stats_Info *stats_info, User_Input *user_inputs) {
    int32_t i, fidx;
    if (TARGET_FILE_PROVIDED) {
        for (fidx=0; fidx<user_inputs->num_of_target_files; fidx++) {
            if (stats_info->capture_cov_stats[fidx]->total_targeted_bases == 0) {
                fprintf(stderr, "ERROR: Total targeted bases is zero. Not Possible\n");
                fprintf(stderr, "No target matches a chromosome in the BAM, or something else went wrong.  Aborting.\n");
                exit(EXIT_FAILURE);
            }

            if (stats_info->capture_cov_stats[fidx]->total_targets == 0) {
                //I don't think we should ever see this error, as its dealt with above.
                fprintf(stderr, "ERROR: No target regions given.  Aborting.\n");
                exit(EXIT_FAILURE);
            }

            // First we need to calculate the coverage for median, this is like N50 for sequencing
            uint64_t sum = 0;
            khiter_t k_iter;
            for (k_iter=0; k_iter!=kh_end(stats_info->capture_cov_stats[fidx]->target_coverage_for_median); k_iter++) {
                if (kh_exist(stats_info->capture_cov_stats[fidx]->target_coverage_for_median, k_iter)) {

                    if (sum >= stats_info->capture_cov_stats[fidx]->total_targeted_bases/2) {
                        stats_info->capture_cov_stats[fidx]->median_target_coverage = k_iter--;
                        break;
                    } else {
                        sum += kh_value(stats_info->capture_cov_stats[fidx]->target_coverage_for_median, k_iter);
                    }
                }
            }

            double average_coverage = (double)stats_info->capture_cov_stats[fidx]->total_target_coverage/(double)stats_info->capture_cov_stats[fidx]->total_targeted_bases;
            FILE *trt_fp = fopen(user_inputs->capture_cov_reports[fidx], "a");    // trt_fp: target_fp

            fprintf(trt_fp, "BUFFER_size\t%d\n", user_inputs->target_buffer_size);
            //printf("Before output general information\n");
            outputGeneralInfo(trt_fp, stats_info, average_coverage);

            float percent = calculatePercentage32_64(stats_info->capture_cov_stats[fidx]->in_buffer_read_hit_count, stats_info->read_cov_stats->total_reads_aligned);
            fprintf(trt_fp, "Aligned_Reads_On-Buffer\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->in_buffer_read_hit_count);
            fprintf(trt_fp, "PCT_of_Aligned_Reads_On-Buffer_(agst_AR)\t%.2f%%\n", percent);

            percent = calculatePercentage32_64(stats_info->capture_cov_stats[fidx]->on_target_read_hit_count, stats_info->read_cov_stats->total_reads_aligned);
            if (user_inputs->remove_duplicate) {
                fprintf(trt_fp, "Aligned_Reads_On-Target_(Total_Usable_Reads)\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->on_target_read_hit_count);
                fprintf(trt_fp, "PCT_of_Aligned_Reads_On-Target_(agst_AR)\t%.2f%%\n", percent);
            } else {
                fprintf(trt_fp, "Aligned_Reads_On-Target_(Include_Duplicates)\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->on_target_read_hit_count);
                fprintf(trt_fp, "PCT_of_Aligned_Reads_On-Target_(agst_AR)\t%.2f%%\n", percent);
            }

            fprintf(trt_fp, "Median_Coverage\t%d\n", stats_info->capture_cov_stats[fidx]->median_target_coverage);

            fprintf(trt_fp, "Target_Max_Coverage\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->target_max_coverage);
            fprintf(trt_fp, "Number_of_Bases_with_Target_Max_Coverage\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->base_with_target_max_coverage);

            uint32_t reads_hit_target_or_buffer = stats_info->capture_cov_stats[fidx]->on_target_read_hit_count + stats_info->capture_cov_stats[fidx]->in_buffer_read_hit_count;
            percent = calculatePercentage32_64(reads_hit_target_or_buffer, stats_info->read_cov_stats->total_reads_aligned);
            fprintf(trt_fp, "Reads_that_hit_target_or_buffer\t%"PRIu32"\n", reads_hit_target_or_buffer);
            fprintf(trt_fp, "PCT_of_Reads_that_hit_target_or_buffer_(agst_AR)\t%.2f%%\n", percent);

            fprintf(trt_fp, "Total_Aligned_Reads_(expected)\t%"PRIu64"\n", stats_info->read_cov_stats->total_reads_aligned);
            fprintf(trt_fp, "Total_Aligned_Reads_(calculated)\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->on_target_read_hit_count + stats_info->capture_cov_stats[fidx]->in_buffer_read_hit_count + stats_info->capture_cov_stats[fidx]->off_target_read_hit_count);

            fprintf(trt_fp, "#Target_Stats\n");
            fprintf(trt_fp, "Total_Targets\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->total_targets);

            percent = calculatePercentage32(stats_info->capture_cov_stats[fidx]->hit_target_count, stats_info->capture_cov_stats[fidx]->total_targets);
            fprintf(trt_fp, "Total_Number_of_Targets_with_Coverage_Hit\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->hit_target_count);
            fprintf(trt_fp, "PCT_of_Total_Number_of_Targets_with_Coverage_Hit\t%.2f%%\n", percent);

            percent = calculatePercentage32(stats_info->capture_cov_stats[fidx]->hit_target_buffer_only_count, stats_info->capture_cov_stats[fidx]->total_targets);
            fprintf(trt_fp, "Number_of_Targets_without_Coverage_but_Buffers_have_Coverage_Hit\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->hit_target_buffer_only_count);
            fprintf(trt_fp, "PCT_of_Number_of_Targets_without_Coverage_but_Buffers_have_Coverage_Hit\t%.2f%%\n", percent);

            fprintf(trt_fp, "Non_target_regions_with_high_coverage\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->non_target_good_hits);

            fprintf(trt_fp, "#Base_Stats\n");
            fprintf(trt_fp, "Bases_Targeted\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->total_targeted_bases);
            fprintf(trt_fp, "Buffer_Bases\t%"PRIu32"\n", stats_info->capture_cov_stats[fidx]->total_buffer_bases);

            uint16_t bins[16] = { 0, 1, 5, 6, 10, 11, 15, 20, 30, 40, 50, 60, 70, 100, 500, 1000 };
            for(i=0; i<16; i++) {
                uint32_t val = getValueFromKhash32(stats_info->capture_cov_stats[fidx]->target_base_with_N_coverage, bins[i]);

                float percent = calculatePercentage32(val, stats_info->capture_cov_stats[fidx]->total_targeted_bases);
                fprintf(trt_fp, "Bases_with_%dx_coverage\t%"PRIu32"\n", bins[i], val);
                fprintf(trt_fp, "PCT_of_Bases_with_%dx_coverage\t%.2f%%\n", bins[i], percent);
            }

            //printf("After the base count for 1000\n");

            fprintf(trt_fp, "\n");
            fprintf(trt_fp, "#Coverage_Frequency_Distribution_For_Capture_Enrichment\n");

            for (i=0; i<=1000; i++)
                fprintf(trt_fp, "%d,", i);

            fprintf(trt_fp, "\n");
            fprintf(trt_fp, ">>");

            for (i=0; i<=1000; i++)
                fprintf(trt_fp, "%"PRIu32",", stats_info->capture_cov_stats[fidx]->target_cov_histogram[i]);

            fprintf(trt_fp, "\n");

            /*fprintf(trt_fp, "Target and region coverage plot\n");
            fprintf(trt_fp, "Position,5'count,3'count\n");
            for(i = 20; i <= PRIMER_SIZE; i+=20) {
                fprintf(trt_fp, "%d,%"PRIu32",%"PRIu32"\n", i, stats_info->five_prime[PRIMER_SIZE-(i-1)-1], stats_info->three_prime[i-1]);
            }

            //fputs("%tar-Pos,count\n", trt_fp);
            fprintf(trt_fp, "target_Pos,count\n");
            for(i = 0; i < 101; i+=2) {
                fprintf(trt_fp, "%d,%"PRIu32"\n", i, stats_info->target_coverage[i]);
            }*/

            fclose(trt_fp);
        }
    }
}

void outputGeneralInfo(FILE *fp, Stats_Info *stats_info, double average_coverage) {
    fprintf(fp, "#Read_Stats\n");
     
    fprintf(fp, "Total_Reads(TR)\t%"PRIu64"\n", stats_info->read_cov_stats->total_reads_produced);

    uint64_t yield = stats_info->read_cov_stats->read_length * (uint64_t) stats_info->read_cov_stats->total_reads_produced;
    fprintf(fp, "Sequenced_Read_Length\t%d\n", stats_info->read_cov_stats->read_length);
    fprintf(fp, "Total_Yield\t%"PRIu64"\n", yield);

    float percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_aligned,stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Aligned_Reads(AR)\t%"PRIu64"\n", stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_Reads_Aligned\t%.2f%%\n", percent);

    yield = stats_info->read_cov_stats->read_length * (uint64_t) (stats_info->read_cov_stats->total_reads_aligned - stats_info->read_cov_stats->total_duplicate_reads);
    uint64_t uniquely_aligned = stats_info->read_cov_stats->total_reads_aligned - stats_info->read_cov_stats->total_duplicate_reads;
    percent = calculatePercentage64(uniquely_aligned, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Uniquely_Aligned_Yield\t%"PRIu64"\n", yield);
    fprintf(fp, "Unique_Aligned_Reads\t%"PRIu64"\n", uniquely_aligned); 
    fprintf(fp, "PCT_of_Unique_Aligned_Reads_(agst_TR)\t%.2f%%\n", percent); 

    percent = calculatePercentage64(uniquely_aligned, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Unique_Aligned_Reads_(agst_AR)\t%.2f%%\n", percent); 

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_duplicate_reads, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Duplicate_Reads\t%"PRIu32"\n", stats_info->read_cov_stats->total_duplicate_reads);
    fprintf(fp, "PCT_of_Duplicate_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_duplicate_reads,stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Duplicate_Reads_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_supplementary_reads,stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Supplementary_Reads\t%"PRIu32"\n", stats_info->read_cov_stats->total_supplementary_reads);
    fprintf(fp, "PCT_of_Supplementary_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_supplementary_reads,stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Supplementary_Reads_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_paired, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Paired_READ\t%"PRIu64"\n", stats_info->read_cov_stats->total_reads_paired);
    fprintf(fp, "PCT_of_Paired_READ_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_paired, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Paired_READ_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_paired_reads_with_mapped_mates, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Paired_Reads_With_Mapped_Mates\t%"PRIu64"\n", stats_info->read_cov_stats->total_paired_reads_with_mapped_mates);
    fprintf(fp, "PCT_of_Paired_Reads_With_Mapped_Mates_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_paired_reads_with_mapped_mates, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Paired_Reads_With_Mapped_Mates_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_proper_paired, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Properly_Paired_Reads\t%"PRIu64"\n",stats_info->read_cov_stats->total_reads_proper_paired);
    fprintf(fp, "PCT_of_Properly_Paired_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage64(stats_info->read_cov_stats->total_reads_proper_paired, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Properly_Paired_Reads_(agst_AR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_chimeric_reads, stats_info->read_cov_stats->total_reads_produced);
    fprintf(fp, "Chimeric_Reads\t%"PRIu32"\n", stats_info->read_cov_stats->total_chimeric_reads);
    fprintf(fp, "PCT_of_Chimeric_Reads_(agst_TR)\t%.2f%%\n", percent);

    percent = calculatePercentage32_64(stats_info->read_cov_stats->total_chimeric_reads, stats_info->read_cov_stats->total_reads_aligned);
    fprintf(fp, "PCT_of_Chimeric_Reads_(agst_AR)\t%.2f%%\n", percent);

    fprintf(fp, "Average_Coverage\t%.2f\n", average_coverage);
}

void produceOffTargetWigFile(Chromosome_Tracking *chrom_tracking, char *chrom_id, Bed_Info *target_bed_info, User_Input *user_inputs, Stats_Info *stats_info, uint8_t target_file_index) {
    int32_t i, j;
    int32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
    if (idx == -1) return;

    FILE *wp = NULL;

    if (user_inputs->Write_WIG) {
        wp = fopen(user_inputs->capture_wig_files[target_file_index], "a");
        fprintf(wp, "track type=wiggle_0 name=%s\n", user_inputs->bam_file);
    }

    for(i=0; (uint32_t) i<target_bed_info->size; i++) {

        // here we are going to erase everything that is on-target and in-buffer and set them to 0
        if (strcmp(chrom_id, target_bed_info->coords[i].chrom_id) == 0) {
            uint32_t start = target_bed_info->coords[i].start;
            uint32_t stop  = target_bed_info->coords[i].end;
            for(j=start-500; (uint32_t) j<=stop+500; j++) {
                if (j < 0 || (uint32_t) j >= chrom_tracking->chromosome_lengths[idx])
                    continue;

                chrom_tracking->coverage[idx][j] = 0;
            }
        }
    }

    // once captured area + BUFFER is initialized to 0, we want to check if there is any coverage > 20 in off-target regions
    for(i=0; (uint32_t) i<chrom_tracking->chromosome_lengths[idx]; i++) {
        if (chrom_tracking->coverage[idx][i] > 20) {
            j = i;
            stats_info->capture_cov_stats[target_file_index]->non_target_good_hits += 1;

            // right side 
            while((uint32_t) i < chrom_tracking->chromosome_lengths[idx] && chrom_tracking->coverage[idx][i] > 0)
                i++;

            // left side
            while(j > 0 && chrom_tracking->coverage[idx][i] > 0)
                j--;

            // now write to the off target wig file
            if (user_inputs->Write_WIG) {
                fprintf(wp, "fixedStep chrom_id=%s start=%"PRId32" step=1\n", chrom_id, j);
                uint32_t h = j;
                for(h=j; h<(uint32_t)i; h++)
                    fprintf(wp,"%"PRIu32"\n", chrom_tracking->coverage[idx][h]);
            }
        }
    }
    if (wp != NULL) fclose(wp);
}

// Note: type 1 is for capture target, while type 2 is for user-defined-database
//
void calculateGenePercentageCoverage(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, khash_t(khStrLCG) *low_cov_gene_hash) {
    // find out the index that is used to track current chromosome id
    int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
    if (chrom_idx == -1) return;

    // print low_cov_genes for debugging
    //printLowCoverageGeneStructure(low_cov_genes);

    // if the target bed file is available, we will need to calculate percentage of gene bases that are covered
    uint32_t i, j;
    if (TARGET_FILE_PROVIDED || USER_DEFINED_DATABASE) {

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
                    produceGenePercentageCoverageInfo(start, end, low_cov_gene_hash);
                }
            }
        }
    }
}

// Note: type 1 is for capture target, while type 2 is for user-defined-database
// The function will loop through several data structure to find the information
// Here is the outline:
//
// gene_transcript[gene_symbol] -> a hash table of khash_t(khStrStrArray) with key as gene symbol and value as a stringArray variable
//                                -> [0] transcript_name0
//                                -> [1] transcript_name1
//                                -> ...
//
// transcript_hash[transcript_name] -> a hash table of khash_t(khStrLCG) with key as transcript_name and value as a Low_Coverage_Genes variable
//                              -> it is just like low_cov_gene_hash table, except everything is grouped by transcript names
//                              -> cds_1
//                              -> cds_2
//                              -> cds_2        we might have multiple same cds as it across multiple 1000 bases
//                              -> cds_3
//                              -> cds_4
//                              -> cds_4
//                              -> cds_4
//                              -> cds_5
//                              -> ...
//
// cds_hash[cds_target_start-cds_target_end] 
//            ==> a hash table of khash_t(khStrLCG) with key as cds tartget start and end with value as a Low_Coverage_Genes variable
//                              -> it is just like low_cov_gene_hash table, except everything is grouped by CDS target start positions
//                              -> low_cov_region[0]
//                              -> low_cov_region[1]
//                              -> low_cov_region[2]
//                              -> ...
//
//
void storeGenePercentageCoverage(char *chrom_id, User_Input *user_inputs, khash_t(khStrLCG) *transcript_hash, khash_t(khStrStrArray) *gene_transcripts, khash_t(khStrInt) *hgmd_transcripts, khash_t(khStrGTP) **gene_transcript_percentage_hash, uint8_t file_index) {
    // if the target bed file is available, we will need to calculate percentage of gene bases that are covered
    //
    if (TARGET_FILE_PROVIDED) {
        // open file handle for writing/appending
        //
        FILE *exon_pct_fp;
        FILE *tspt_pct_fp;

        exon_pct_fp = fopen(user_inputs->low_cov_exon_pct_files[file_index], "a");
        tspt_pct_fp = fopen(user_inputs->low_cov_transcript_files[file_index], "a");

        // loop through all genes in the gene_transcripts hash table
        //
        khiter_t iter_gt;
        for (iter_gt=kh_begin(gene_transcripts); iter_gt!=kh_end(gene_transcripts); ++iter_gt) {
            if (kh_exist(gene_transcripts, iter_gt)) {
                // one gene at a time
                //
                char gene_symbol[30];
                strcpy(gene_symbol, kh_key(gene_transcripts, iter_gt));

                // create a bucket for current gene in gene_transcript_percentage_hash
                //
                int absent;
                uint32_t i, j;
                khiter_t iter_gtp = kh_put(khStrGTP, gene_transcript_percentage_hash[file_index], gene_symbol, &absent);;
                if (absent) {
                    kh_key(gene_transcript_percentage_hash[file_index], iter_gtp) = strdup(gene_symbol);
                    kh_value(gene_transcript_percentage_hash[file_index], iter_gtp) = calloc(1, sizeof(Gene_Transcript_Percentage));
                    kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->capacity = 2;
                    kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->size = 0;
                    kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->has_HGMD = false;
                    kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->transcript_percentage = 
                        calloc(kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->capacity, sizeof(Transcript_Percentage));

                    // initialize Transcript_Percentage variables to NULL 
                    //
                    for (i=0; i<kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->capacity; i++) {
                        kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->transcript_percentage[i].transcript_name = NULL;
                        kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->transcript_percentage[i].percentage = 0.000;
                        kh_value(gene_transcript_percentage_hash[file_index], iter_gtp)->transcript_percentage[i].HGMD = false;
                    }
                }

                for (i=0; i<kh_value(gene_transcripts, iter_gt)->size; i++) {
                    // one transcript at a time
                    // check if current transcript name exists in transcript_hash table
                    //
                    char *transcript_name = calloc(strlen(kh_value(gene_transcripts, iter_gt)->theArray[i])+1, sizeof(char));
                    strcpy(transcript_name, kh_value(gene_transcripts, iter_gt)->theArray[i]);
                    khiter_t iter_tsh = kh_get(khStrLCG, transcript_hash, transcript_name);

                    if (iter_tsh == kh_end(transcript_hash)) {
                        // transcript name doesn't exist
                        //
                        continue;
                    } else {
                        // It is possible that we have duplicate entries, so we need to remove duplicate here first
                        // by re-group everything based on cds target start positions
                        //
                        khash_t(khStrLCG) *cds_hash = kh_init(khStrLCG);

                        for (j=0; j<kh_value(transcript_hash, iter_tsh)->total_size; j++) {
                            // check if current cds in the transcript_hash is targeted!
                            //
                            if (!kh_value(transcript_hash, iter_tsh)->gene_coverage[j].targeted)
                                continue;

                            // check if the current cds exists in cds_hash
                            //
                            char cur_cds_key_str[30];
                            sprintf(cur_cds_key_str, "%"PRIu32"-%"PRIu32, kh_value(transcript_hash, iter_tsh)->gene_coverage[j].cds_target_start, kh_value(transcript_hash, iter_tsh)->gene_coverage[j].cds_target_end);

                            // Note: we need to add all CDSs with the same start position here to form an array
                            // The array is needed because there might be multiple low coverage regions
                            //
                            //    cds_target_start ------------------------------------------------------ cds_target_end
                            //                     1000            2000            3000            4000
                            //                       ***        **            ****        **********
                            //                        bucket1            bucket2            bucket3
                            //
                            // As we can see that some of them have the same start and end position, 
                            // they will be recorded into different cds instances with different low coverage regions
                            // Therefore, we need to record all of them
                            //

                            // Initialize the bucket with the current key
                            //
                            lowCoverageGeneHashBucketKeyInit(cds_hash, cur_cds_key_str);

                            // Let's retrieve the bucket
                            //
                            khiter_t iter_cds = kh_get(khStrLCG, cds_hash, cur_cds_key_str);
                            if (iter_cds == kh_end(cds_hash)) {
                                // the bucket with the current key_str doesn't exists!
                                // This shouldn't happen. So let's exit()
                                //
                                fprintf(stderr, "ERROR: There is no CDS exist, which shouldn't happen, so exit()!\n");
                                exit(EXIT_FAILURE);
                            }

                            // add current cds in to the array
                            //
                            uint32_t idx = kh_value(cds_hash, iter_cds)->total_size;
                            copyGeneCoverageLowCovRegions(&kh_value(transcript_hash, iter_tsh)->gene_coverage[j],
                                    &kh_value(cds_hash, iter_cds)->gene_coverage[idx], false);
                            kh_value(cds_hash, iter_cds)->total_size++;
                        } // end current transcript name

                        // Everything is now grouped by cds target start position,
                        // Process one CDS at a time and record necessary information for current transcript
                        //
                        uint32_t cds_length=0, total_low_cov_bases=0, exon_count=0;

                        khiter_t it_cds;
                        for (it_cds=kh_begin(cds_hash); it_cds!=kh_end(cds_hash); ++it_cds) {
                            if (kh_exist(cds_hash, it_cds)) {
                                int32_t exon_id=0;        // signed int as it needs to compare with negative values
                                uint32_t cds_target_start=0, cds_target_end=0;

                                // process current CDS and record its low coverage regioins
                                //
                                khash_t(khStrInt) *low_cov_regions_hash = kh_init(khStrInt);
                                uint16_t num_of_low_cov_regions = 0;

                                uint16_t p, q;
                                bool first=true;
                                for (p=0; p<kh_value(cds_hash, it_cds)->total_size; p++) {
                                    // Create a tmp Gene_Coverage point as it will be easy to handle
                                    // record information for the current CDS
                                    //
                                    Gene_Coverage* gc = &kh_value(cds_hash, it_cds)->gene_coverage[p];
                                    exon_id = gc->exon_id;
                                    if (first) {
                                        cds_target_start = gc->cds_target_start;
                                        cds_target_end = gc->cds_target_end;
                                        //cds_length += cds_target_end - cds_target_start +1; //only for debug, will remove +1 later
                                        cds_length += cds_target_end - cds_target_start;
                                        exon_count++;
                                        first = false;
                                    }

                                    if (gc->low_cov_regions != NULL) {
                                        for (q = 0; q < gc->low_cov_regions->size; q++) {
                                            khiter_t it_lcr = kh_put(khStrInt, low_cov_regions_hash, gc->low_cov_regions->theArray[q], &absent);
                                            if (absent) {
                                                kh_key(low_cov_regions_hash, it_lcr) = strdup(gc->low_cov_regions->theArray[q]);
                                                num_of_low_cov_regions++;
                                            }
                                            kh_value(low_cov_regions_hash, it_lcr) = 1;
                                        }
                                    }
                                }

                                // Once we remove the duplicates through low_cov_regions_hash table
                                // we need to merge them if needed as shown by the following example:
                                // gc->low_cov_regions[0]->theArray[0] = "25316999-25317052"
                                // gc->low_cov_regions[0]->theArray[1] = "25316999-25317079"
                                //
                                char * low_cov_region_string=calloc(2, sizeof(char));
                                strcpy(low_cov_region_string, ".");
                                uint32_t cur_low_cov_count=0;

                                if (num_of_low_cov_regions > 0) {
                                    stringArray *mergedRegions = calloc(1, sizeof(stringArray));

                                    if (num_of_low_cov_regions == 1) {
                                        getOneLowCovRegions(low_cov_regions_hash, mergedRegions);
                                    } else {
                                        mergeLowCovRegions(low_cov_regions_hash, mergedRegions, num_of_low_cov_regions, cds_target_start, cds_target_end);
                                    }

                                    // we then walk through the merged low_cov_regions to obtain the necessary information
                                    //
                                    if (mergedRegions->size > 0) {
                                        //cur_low_cov_count = processLowCovRegionFromKhash(mergedRegions, &low_cov_region_string);
                                        cur_low_cov_count = processLowCovRegionFromStrArray(mergedRegions, &low_cov_region_string);
                                        total_low_cov_bases += cur_low_cov_count;
                                    }

                                    if (mergedRegions != NULL) {                                                                           
                                        stringArrayDestroy(mergedRegions);
                                        free(mergedRegions);
                                    }
                                }

                                if (cur_low_cov_count > (cds_target_end - cds_target_start)) {
                                    fprintf(stderr, "Something is wrong with this CDS, its start is %"PRIu32" and its end is %"PRIu32", while its low coverage base count is %"PRIu32"; gene:%s with transcript %s\n", cds_target_start, cds_target_end, cur_low_cov_count, gene_symbol, transcript_name);
                                    cur_low_cov_count = cds_target_end - cds_target_start;
                                }

                                // now we have everything and we should be able to do the calculation and output the current exon/cds result
                                //
                                float pct = (float) ((cds_target_end - cds_target_start - cur_low_cov_count) * 100) / (float) (cds_target_end - cds_target_start);

                                // Note: there are several cases need to be handled.
                                // 1). SNP: if it is for SNPs, the exon_id is set to -1 in the database, but the output will be '_'
                                // 2). eMerge: cds are not designed but will be output into the final results, the exon_id is set to -exon_id (such as -13)
                                // 3). regular genes output
                                //
                                if (exon_id == -1) {
                                    fprintf(exon_pct_fp, "%s\t%s\t%s\t_\t%"PRIu32"\t%"PRIu32"\t%0.3f%%\t%s\n", chrom_id, gene_symbol, transcript_name, cds_target_start, cds_target_end, pct, low_cov_region_string);
                                } else if (exon_id < -1 ) {
                                    fprintf(exon_pct_fp, "%s\t%s\t%s\tcds_%d\t%"PRIu32"\t%"PRIu32"\t0.000%%\tNOT_Designed\n", chrom_id, gene_symbol, transcript_name, exon_id, cds_target_start, cds_target_end);
                                } else{
                                    fprintf(exon_pct_fp, "%s\t%s\t%s\tcds_%"PRIu16"\t%"PRIu32"\t%"PRIu32"\t%0.3f%%\t%s\n", chrom_id, gene_symbol, transcript_name, exon_id, cds_target_start, cds_target_end, pct, low_cov_region_string);
                                    //fprintf(exon_pct_fp, "%s\t", gene_symbol);
                                    //fprintf(exon_pct_fp, "%s\t", transcript_name);
                                    //fprintf(exon_pct_fp, "%s\n", low_cov_region_string);
                                }

                                if (low_cov_region_string != NULL) {
                                    free(low_cov_region_string);
                                    low_cov_region_string = NULL;
                                }

                                // clean-up
                                //
                                cleanKhashStrInt(low_cov_regions_hash);

                            } // finishes current CDS group if it exists!

                        } // end for loop for kh_begin() and kh_end() for grouping everything based on CDS start position

                        // if there is no CDSs associated with this gene/transcript, don't output anything
                        //
                        if (cds_length == 0) {
                            // but we still need to clean cds_hash first, otherwise there will be memory leak
                            //
                            genePercentageCoverageDestroy(cds_hash);
                            continue;
                        }

                        // now need to output everything for current transcripts
                        // but first we need to double check with HGMD entries
                        //
                        float transcript_pct = (float) (cds_length - total_low_cov_bases)*100 / cds_length;

                        // output transcript information first
                        //
                        fprintf(tspt_pct_fp, "%s\t%s\t%s\t%"PRIu32"\t%"PRIu16"\t%0.3f%%\n", chrom_id, gene_symbol, transcript_name, cds_length, exon_count, transcript_pct);

                        // checking HGMD and if gene is in HGMD
                        // now need to check the transcript name
                        // Note: some genes have multiple transcript with the same name, I added -1, -2 suffix to the end of the transcript
                        // Therefore, the transcript_name key will be changed and not found.
                        // Need to do some modification here
                        //
                        if (HGMD_PROVIDED && hgmd_transcripts) {
                            char * tmp_t_name  = calloc(strlen(transcript_name)+1, sizeof(char));
                            strcpy(tmp_t_name, transcript_name);
                            char * orig_t_name = calloc(strlen(transcript_name)+1, sizeof(char));

                            char *ret = strstr(transcript_name, "-");

                            if (ret != NULL) {
                                // need to modify the transcript_name
                                //
                                char *savePtr, *tokPtr;
                                savePtr = tmp_t_name;
                                uint32_t i = 0;

                                while ((tokPtr = strtok_r(savePtr, "-", &savePtr))) {
                                    if (i == 0)
                                        strcpy(orig_t_name, tokPtr);
                                    i++;
                                }
                            } else {
                                strcpy(orig_t_name, transcript_name);
                            }

                            khiter_t trt_iter = kh_get(khStrInt, hgmd_transcripts, orig_t_name);
                            if (trt_iter != kh_end(hgmd_transcripts)) {
                                storeHGMD_TranscriptPercentageInfo(gene_transcript_percentage_hash[file_index], gene_symbol, transcript_name, transcript_pct, true);

                            } else {
                                storeHGMD_TranscriptPercentageInfo(gene_transcript_percentage_hash[file_index], gene_symbol, transcript_name, transcript_pct, false);
                            }

                            if (orig_t_name) free(orig_t_name);
                            if (tmp_t_name) free(tmp_t_name);
                        } else {
                            storeHGMD_TranscriptPercentageInfo(gene_transcript_percentage_hash[file_index], gene_symbol, transcript_name, transcript_pct, false);
                        }

                        // clean-up
                        //
                        genePercentageCoverageDestroy(cds_hash);

                    } // end of transcript_hash exist for current gene key
                    
                    if (transcript_name) free(transcript_name);

                } // end of all the transcripts for one gene (one transcript at a time)
            } // end if current gene kh_exist()

        } // end for loop kh_begin() and kh_end() for gene list (one gene at a time)

        fclose(exon_pct_fp);
        fclose(tspt_pct_fp);
    }
}

void storeHGMD_TranscriptPercentageInfo(khash_t(khStrGTP) *gene_transcript_percentage_hash, char *gene_symbol,  char* transcript_name, float percentage, bool HGMD_on) {
    // now store the information for Transcript_Percentage variable
    //
    khiter_t tmp_iter = kh_get(khStrGTP, gene_transcript_percentage_hash, gene_symbol);
    if (tmp_iter == kh_end(gene_transcript_percentage_hash)) {
        fprintf(stderr, "The key %s doesn't exist in gene_transcript_percentage_hash\n", gene_symbol);
    } else {
        // check to see if the bucket is big enough, otherwise, expand the memory
        //
        if (kh_value(gene_transcript_percentage_hash, tmp_iter)->capacity == kh_value(gene_transcript_percentage_hash, tmp_iter)->size) {
            kh_value(gene_transcript_percentage_hash, tmp_iter)->capacity = kh_value(gene_transcript_percentage_hash, tmp_iter)->capacity * 2;
            kh_value(gene_transcript_percentage_hash, tmp_iter)->transcript_percentage =
                realloc(kh_value(gene_transcript_percentage_hash, tmp_iter)->transcript_percentage, 
                        kh_value(gene_transcript_percentage_hash, tmp_iter)->capacity * sizeof(Transcript_Percentage));
            if (kh_value(gene_transcript_percentage_hash, tmp_iter)->transcript_percentage == NULL) {
                fprintf(stderr, "ERROR: Memory re-allocation failed at storeHGMD_TranscriptPercentageInfo\n");
                exit(EXIT_FAILURE);
            }
        }

        // do the insertion
        //
        kh_value(gene_transcript_percentage_hash, tmp_iter)->transcript_percentage[kh_value(gene_transcript_percentage_hash, tmp_iter)->size].transcript_name = calloc(strlen(transcript_name)+1, sizeof(char));
        strcpy(kh_value(gene_transcript_percentage_hash, tmp_iter)->transcript_percentage[kh_value(gene_transcript_percentage_hash, tmp_iter)->size].transcript_name, transcript_name);
        kh_value(gene_transcript_percentage_hash, tmp_iter)->transcript_percentage[kh_value(gene_transcript_percentage_hash, tmp_iter)->size].percentage = percentage;
        kh_value(gene_transcript_percentage_hash, tmp_iter)->transcript_percentage[kh_value(gene_transcript_percentage_hash, tmp_iter)->size].HGMD = HGMD_on;

        if (HGMD_on)
            kh_value(gene_transcript_percentage_hash, tmp_iter)->has_HGMD = HGMD_on;

        kh_value(gene_transcript_percentage_hash, tmp_iter)->size++;
    }
}

void outputGeneCoverage(khash_t(khStrGTP) *gene_transcript_percentage_hash, User_Input *user_inputs, uint8_t file_index) {
    FILE *gene_pct_fp;
    gene_pct_fp = fopen(user_inputs->low_cov_gene_pct_files[file_index], "a");

    uint16_t num_of_transcripts=0;
    float total_percentage=0.00;
    khiter_t iter;

    for (iter=kh_begin(gene_transcript_percentage_hash); iter!=kh_end(gene_transcript_percentage_hash); ++iter) {
        if (kh_exist(gene_transcript_percentage_hash, iter)) {
            if (kh_value(gene_transcript_percentage_hash, iter)->size == 0)
                continue;

            // output gene symbol => kh_key(gene_transcript_percentage_hash, iter) exist
            //
            fprintf(gene_pct_fp, "%s\t", kh_key(gene_transcript_percentage_hash, iter));

            int i;
            bool first=true;
            bool hgmd_flag=false;
            for (i=0; i<kh_value(gene_transcript_percentage_hash, iter)->size; i++) {
                if (HGMD_PROVIDED) {
                    // if has_HGMD is false, there is no HGMD transcripts, then output everything
                    // OR output HGMD only transcripts
                    //
                    if (!kh_value(gene_transcript_percentage_hash, iter)->has_HGMD ||
                            kh_value(gene_transcript_percentage_hash, iter)->transcript_percentage[i].HGMD ) {
                        if (first) {
                            first = false;
                        } else {
                            fprintf(gene_pct_fp, ";");
                        }

                        fprintf(gene_pct_fp, "%s(%.3f)", 
                            kh_value(gene_transcript_percentage_hash, iter)->transcript_percentage[i].transcript_name,
                            kh_value(gene_transcript_percentage_hash, iter)->transcript_percentage[i].percentage);

                        num_of_transcripts++;
                        total_percentage += kh_value(gene_transcript_percentage_hash, iter)->transcript_percentage[i].percentage;
                        
                    }

                    if (kh_value(gene_transcript_percentage_hash, iter)->has_HGMD)
                        hgmd_flag=true;

                } else {
                    // output everything no matter what
                    //
                    if (first) {
                        first = false;
                    } else {
                        fprintf(gene_pct_fp, ";");
                    }

                    fprintf(gene_pct_fp, "%s(%.3f)",
                        kh_value(gene_transcript_percentage_hash, iter)->transcript_percentage[i].transcript_name,
                        kh_value(gene_transcript_percentage_hash, iter)->transcript_percentage[i].percentage);

                    num_of_transcripts++;
                    total_percentage += kh_value(gene_transcript_percentage_hash, iter)->transcript_percentage[i].percentage;

                }
            } // end of transcripts for one gene

            if (num_of_transcripts == 0) continue;

            float ave_cov_pct_for_gene = total_percentage / (float) (num_of_transcripts);
            if (HGMD_PROVIDED) {
                if (hgmd_flag) {
                    fprintf(gene_pct_fp, "\t%.2f\tHGMD\n", ave_cov_pct_for_gene);
                } else {
                    fprintf(gene_pct_fp, "\t%.2f\n", ave_cov_pct_for_gene);
                }
            } else {
                fprintf(gene_pct_fp, "\t%.2f\n", ave_cov_pct_for_gene);
            }

            num_of_transcripts = 0;
            total_percentage = 0.0;
        } // end of khash keyed by gene

    }
    fclose(gene_pct_fp);
}
