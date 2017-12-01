/*
 * =====================================================================================
 *
 *       Filename:  annotation.c
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
#include <strings.h>
#include "annotation.h"
#include "utils.h"

// The followings are codes that related to gene annotation
//
void finish_with_error(MYSQL *con)
{
    fprintf(stderr, "%s\n", mysql_error(con));
    mysql_close(con);
    exit(1);
}

// make the MySQL connection and create the name of databases to be used by the program
//
void databaseSetup(Databases *dbs, User_Input *user_inputs) {
	dbs->con = mysql_init(NULL);
	if (dbs->con == NULL)
		finish_with_error(dbs->con);

	if (mysql_library_init(0, NULL, NULL)) {
		fprintf(stderr, "Could not initialize MySQL Library!");
		exit(1);
	}

	if (mysql_real_connect(dbs->con, "sug-esxa-db1", "phuang", "phuang", "GeneAnnotations", 0, NULL, 0) == NULL) {
		finish_with_error(dbs->con);
	}

	// once the MYSQL connection established, we need to find out which database to use here
	//
	dbs->db_introns = calloc(30, sizeof(char));

	if (user_inputs->annotation_type == 1) {
		// dynamic way,
		// use string case insensitive way
		//
		fprintf(stderr, "Dynamic generation of gene annotation is ON\n");

		dbs->db_dynamic = calloc(30, sizeof(char));
		dbs->db_annotation = calloc(30, sizeof(char));

		if (strcasecmp(user_inputs->database_version, "hg38") == 0) {
			strcpy(dbs->db_dynamic, "Gene_RefSeq_CDS38_orig");
			strcpy(dbs->db_introns, "Intron_Regions38");
			strcpy(dbs->db_annotation, "VCRomePKv2_Annotation38");

		} else {
			strcpy(dbs->db_dynamic, "Gene_RefSeq_CDS37_orig");
			strcpy(dbs->db_introns, "Intron_Regions37");
			strcpy(dbs->db_annotation, "VCRomePKv2_Annotation37");

		}

		fprintf(stderr, "The database for dynamic usage is: %s\n", dbs->db_dynamic);
		fprintf(stderr, "The dynamic database for gene annotation is %s\n", dbs->db_annotation);
	} else {
		// static way
		//
		fprintf(stderr, "Static use of prev-existing gene annotation database is ON\n");

		dbs->db_coords = calloc(30, sizeof(char));
		dbs->db_annotation = calloc(30, sizeof(char));

		if (strcasecmp(user_inputs->database_version, "hg38") == 0) {
			if (user_inputs->database_category == 1) {
				strcpy(dbs->db_coords, "Gene_RefSeq_CDS38_orig");
				strcpy(dbs->db_annotation, "VCRomePKv2_Annotation38");

			} else if (user_inputs->database_category == 2) {
				strcpy(dbs->db_coords, "eMerge_CDS38");
				strcpy(dbs->db_annotation, "eMerge_Annotation38");

			} else {
				strcpy(dbs->db_coords, "Right10K_CDS38");
				strcpy(dbs->db_annotation, "Right10K_Annotation38");

			}

			strcpy(dbs->db_introns, "Intron_Regions38");
		} else {
			if (user_inputs->database_category == 1) {
				strcpy(dbs->db_coords, "PKv2_CDS37");
				strcpy(dbs->db_annotation, "VCRomePKv2_Annotation37");

			} else if (user_inputs->database_category == 2) {
				strcpy(dbs->db_coords, "eMerge_CDS37");
				strcpy(dbs->db_annotation, "eMerge_Annotation37");
			} else {
				strcpy(dbs->db_coords, "Right10K_SNP37");
				strcpy(dbs->db_annotation, "Right10K_Annotation37");
			}
				
			strcpy(dbs->db_introns, "Intron_Regions37");
		}

		fprintf(stderr, "The static annotation database used is: %s\n", dbs->db_annotation);
		fprintf(stderr, "The static coordinate database (for CDS) used is: %s\n", dbs->db_coords);

	}

	fprintf(stderr, "The Intronic database is: %s\n\n", dbs->db_introns);
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions
//
uint32_t fetchTotalCount(uint8_t type, Databases *dbs, char *chrom_id, User_Input *user_inputs) {
	char *pre_sql = calloc(100, sizeof(char));

	if (type == 1) {
		sprintf(pre_sql, "SELECT count(*) from Intergenic_Regions WHERE chrom='%s'", chrom_id);
	} else if (type == 2) {
		sprintf(pre_sql, "SELECT count(*) from %s WHERE chrom='%s'", dbs->db_introns, chrom_id);
	} else if (type == 3) {
		sprintf(pre_sql, "SELECT count(*) from %s WHERE chrom='%s'", dbs->db_annotation, chrom_id);
	}
	//printf("%s\n", pre_sql);

	if (mysql_query(dbs->con,pre_sql))
        finish_with_error(dbs->con);

    dbs->mysql_results = mysql_store_result(dbs->con);
    if (dbs->mysql_results == NULL)
        finish_with_error(dbs->con);

    MYSQL_ROW row;
	uint32_t count=0;
    while ((row = mysql_fetch_row(dbs->mysql_results))) {
		if (row[0] > 0)
			count = (uint32_t) atol(row[0]);
	}

	if (dbs->mysql_results) { mysql_free_result(dbs->mysql_results); dbs->mysql_results = NULL; }

	free(pre_sql);
	pre_sql = NULL;

	return count;
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions
// this is used for general gene annotation, but not for the percentage coverage annotation
//
void regionsSkipMySQLInit(Databases *dbs, Regions_Skip_MySQL *regions_in, User_Input *user_inputs, uint8_t type) {
	// here we need to find out how many chromosomes we are dealing with
	char *sql = calloc(250, sizeof(char));
		
	// Here NOT LIKE will exclude non-primary chromosomes such as 19_GL949752v1_alt
	// For hg38, we need to include them in
	//
	if (type == 1) {
		//sprintf(sql, "SELECT DISTINCT chrom FROM Intergenic_Regions WHERE chrom NOT LIKE '%%\\_%%'");
		sprintf(sql, "SELECT DISTINCT chrom FROM Intergenic_Regions");
	} else if (type == 2) {
		//sprintf(sql, "SELECT DISTINCT chrom FROM %s WHERE chrom NOT LIKE '%%\\_%%'", dbs->db_introns);
		sprintf(sql, "SELECT DISTINCT chrom FROM %s", dbs->db_introns);
	} else {
		//sprintf(sql, "SELECT DISTINCT chrom FROM %s WHERE chrom NOT LIKE '%%\\_%%'", dbs->db_annotation);
		sprintf(sql, "SELECT DISTINCT chrom FROM %s", dbs->db_annotation);
	}

	if (mysql_query(dbs->con,sql))
        finish_with_error(dbs->con);

	dbs->mysql_results = mysql_store_result(dbs->con);
    if (dbs->mysql_results == NULL)
        finish_with_error(dbs->con);

	regions_in->chrom_list_size = mysql_num_rows(dbs->mysql_results);
	regions_in->chromosome_ids = calloc(regions_in->chrom_list_size, sizeof(char*));

	uint16_t chrom_idx=0;
    MYSQL_ROW row;
    while ((row = mysql_fetch_row(dbs->mysql_results))) {
		regions_in->chromosome_ids[chrom_idx] = calloc(strlen(row[0])+1, sizeof(char));
		strcpy(regions_in->chromosome_ids[chrom_idx], row[0]);
		chrom_idx++;
	}
	if (dbs->mysql_results) { mysql_free_result(dbs->mysql_results); dbs->mysql_results = NULL; }
	free(sql);
	sql = NULL;

	regions_in->starts = calloc(regions_in->chrom_list_size, sizeof(uint32_t*));
	regions_in->ends   = calloc(regions_in->chrom_list_size, sizeof(uint32_t*));
	regions_in->size_r = calloc(regions_in->chrom_list_size, sizeof(uint32_t));

	regions_in->prev_search_loc_index   = 0;
	regions_in->prev_search_chrom_index = 0;

	// for exon and intronic regions
	if (type > 1) {
        regions_in->gene = calloc(regions_in->chrom_list_size, sizeof(char**));
        if (!regions_in->gene) {
            fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Gene");
            exit(1);
        }

        regions_in->Synonymous = calloc(regions_in->chrom_list_size, sizeof(char**));
	    if (!regions_in->Synonymous) {
		    fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Synonymous");
			exit(1);
        }

	    regions_in->prev_genes = calloc(regions_in->chrom_list_size, sizeof(char**));
		if (!regions_in->prev_genes) {
			fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Prev_genes");
	        exit(1);
        }

		// for exon regions only
		if (type == 3) {
			regions_in->exon_info = calloc(regions_in->chrom_list_size, sizeof(char**));
	        if (!regions_in->exon_info) {
		        fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Gene");
			    exit(1);
			}
        }
	}

	uint32_t i, j;
	for (i=0; i<regions_in->chrom_list_size; i++) {
		regions_in->size_r[i] = fetchTotalCount(type, dbs, regions_in->chromosome_ids[i], user_inputs);
		regions_in->starts[i] = calloc(regions_in->size_r[i], sizeof(uint32_t));

		if (!regions_in->starts[i]) {
			fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Starts");
			exit(1);
		}

		regions_in->ends[i] = calloc(regions_in->size_r[i], sizeof(uint32_t));
		if (!regions_in->ends) {
			fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Ends");
			exit(1);
		}
	
		// for exon and intronic regions
		if (type > 1) {
			regions_in->gene[i] = calloc(regions_in->size_r[i], sizeof(char*));
            regions_in->Synonymous[i] = calloc(regions_in->size_r[i], sizeof(char*));
	        regions_in->prev_genes[i] = calloc(regions_in->size_r[i], sizeof(char*));

			for (j=0; j<regions_in->size_r[i]; j++) {
				regions_in->gene[i][j] = NULL;
				regions_in->Synonymous[i][j] = NULL;
				regions_in->prev_genes[i][j] = NULL;

				//regions_in->gene[i][j] = calloc(1, sizeof(char));
				//regions_in->Synonymous[i][j] = calloc(1, sizeof(char));
				//regions_in->prev_genes[i][j] = calloc(1, sizeof(char));
			}

			// for exon regions only
			if (type == 3) {
				regions_in->exon_info[i] = calloc(regions_in->size_r[i], sizeof(char*));

				for (j=0; j<regions_in->size_r[i]; j++) {
					regions_in->exon_info[i][j] = NULL;
					//regions_in->exon_info[i][j] = calloc(1, sizeof(char));
				}
			}
		}

		// now populate these regions without gene annotation
		//printf("current chromosome is %s for type %d and at %d\n", regions_in->chromosome_ids[i], type, i);
		populateStaticRegionsForOneChromOnly(regions_in, dbs, regions_in->chromosome_ids[i], i, user_inputs, type); 
	}
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions, type 4 for all_site_reports 
void regionsSkipMySQLDestroy(Regions_Skip_MySQL *regions_in, uint8_t type) {
	uint32_t i, j;

	for (i=0; i<regions_in->chrom_list_size; i++) {
		// for exon and intronic regions
		if (type > 1) {
			//for (j=0; j<regions_in->size_r[i]; j++) {
			//	if (regions_in->gene[i][j]) free(regions_in->gene[i][j]);
			//	if (regions_in->Synonymous[i][j]) free(regions_in->Synonymous[i][j]);
			//	if (regions_in->prev_genes[i][j]) free(regions_in->prev_genes[i][j]);
			//}

			if (regions_in->gene[i]) free(regions_in->gene[i]);
			if (regions_in->Synonymous[i]) free(regions_in->Synonymous[i]);
			if (regions_in->prev_genes[i]) free (regions_in->prev_genes[i]);

			// for exon regions only
			if (type == 3) {
				for (j=0; j<regions_in->size_r[i]; j++) {
					if (regions_in->exon_info[i][j]) free(regions_in->exon_info[i][j]);
				}
				if (regions_in->exon_info[i]) free(regions_in->exon_info[i]);
			}
		}

		if (regions_in->starts[i]) free(regions_in->starts[i]);
		if (regions_in->ends[i])   free(regions_in->ends[i]);
		if (regions_in->chromosome_ids[i])   free(regions_in->chromosome_ids[i]);
	}

	if (regions_in->starts) free(regions_in->starts);
	if (regions_in->ends)   free(regions_in->ends);
	if (regions_in->size_r) free(regions_in->size_r);
	if (regions_in->chromosome_ids)   free(regions_in->chromosome_ids);

	// for exon and intronic regions
	if (type > 1) {
		if (regions_in->gene) free(regions_in->gene);
		if (regions_in->Synonymous) free(regions_in->Synonymous);
		if (regions_in->prev_genes) free (regions_in->prev_genes);

		// for exon regions only
    	if (type == 3)
        	if (regions_in->exon_info) free(regions_in->exon_info);
	}

	if (regions_in) free(regions_in);
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions
// This is for the general gene annotation, not for the gene/transcript/CDS coverage percentage calculation
//
void populateStaticRegionsForOneChromOnly(Regions_Skip_MySQL *regions_in, Databases *dbs, char *chrom_id, uint32_t chrom_idx, User_Input *user_inputs, uint8_t type) {
	char *sql = calloc(250, sizeof(char));

	if (type == 1) {
		sprintf(sql, "SELECT start, end from Intergenic_Regions WHERE chrom='%s' ORDER BY start", chrom_id);
	} else if (type == 2) {
		sprintf(sql, "SELECT start, end, gene_symbol, Synonymous, prev_gene_symbol from %s WHERE chrom='%s' ORDER BY start", dbs->db_introns, chrom_id);
	} else {
		sprintf(sql, "SELECT start, end, gene_symbol, Synonymous, prev_gene_symbol, annotation from %s WHERE chrom='%s' ORDER BY start", dbs->db_annotation, chrom_id);
	}
	//printf("%s\n", sql);

	if (mysql_query(dbs->con,sql))
        finish_with_error(dbs->con);

    dbs->mysql_results = mysql_store_result(dbs->con);
    if (dbs->mysql_results == NULL)
        finish_with_error(dbs->con);

    MYSQL_ROW row;
    uint32_t count=0;
	while ((row = mysql_fetch_row(dbs->mysql_results))) {
        regions_in->starts[chrom_idx][count] = (uint32_t) atol(row[0]);
        regions_in->ends[chrom_idx][count] = (uint32_t) atol(row[1]);

		if (type > 1) {
			dynamicStringAllocation(row[2], &regions_in->gene[chrom_idx][count]) ;
			dynamicStringAllocation(row[3], &regions_in->Synonymous[chrom_idx][count]);
			dynamicStringAllocation(row[4], &regions_in->prev_genes[chrom_idx][count]);

			if (type == 3) {
				char *tmp = calloc(strlen(row[5]) + 50, sizeof(char));
				sprintf(tmp, "%s", row[5]);
				dynamicStringAllocation(tmp, &regions_in->exon_info[chrom_idx][count]);
				if (tmp) {
					free(tmp);
					tmp = NULL;
				}
			}
		}

		count++;
	}

	if (dbs->mysql_results) { mysql_free_result(dbs->mysql_results); dbs->mysql_results = NULL; }

    free(sql);
	sql = NULL;
}

// As the array is sorted, I am going to use binary search for find the location
int32_t checkInterGenicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t low_search_index) {
	int32_t found = binarySearch(regions_in, start, end, chrom_idx, low_search_index);

	return found;
}

// Note: for bed file format, the end position is not part of the region to be checked (so exclude end position)!
//
int32_t binarySearch(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t low_search_index) {
	int32_t low = low_search_index;
    int32_t high = regions_in->size_r[chrom_idx] - 1;
    int32_t middle = (low + high)/2;

    while (low <= high) {
		if (((regions_in->starts[chrom_idx][middle] <= start) && (start < regions_in->ends[chrom_idx][middle])) ||
				((regions_in->starts[chrom_idx][middle] <= end) && (end < regions_in->ends[chrom_idx][middle]))) {
			// take care of the cases: region.start |=========================| region.end
    		//                                         start\----------------------\end
    		//                          start\----------------------\end
			return middle;
		} else if ((start <= regions_in->starts[chrom_idx][middle]) && (regions_in->ends[chrom_idx][middle] <= end)) {
			// take care of the case: region.start |=========================| region.end
			//                       start\----------------------------------------------\end
            return middle;
        } else {
            if (regions_in->starts[chrom_idx][middle] < start) {
                low = middle + 1;
            } else {
                high = middle - 1;
            }
        }

        middle = (low + high)/2;
    }

    // outside the while loop means low > high
    //if (low > high)
    //  return -1;

    return -1;
}

// For this search I need to handle the following case:
// ================================== 500 ******************* 600 ==============================
//                  400 --------------------520
//                      420 ---- 450
// Because I sorted everything by the start position, I will miss the overlap region at 400-520
// because the code would exist at 420-450
// Therefore, I need to use the cds_start and cds_end position instead (not by exon)
//
int32_t binarySearchLowCoverage(Low_Coverage_Genes *low_cov_genes, uint32_t start, uint32_t end, uint32_t lower_bound) {
	int32_t low = lower_bound;
	int32_t high = low_cov_genes->total_size - 1;
	int32_t middle = (low + high)/2;

	while (low <= high) {
		if ( (low_cov_genes->gene_coverage[middle].cds_start <= start && start <= low_cov_genes->gene_coverage[middle].cds_end) ||
               (low_cov_genes->gene_coverage[middle].cds_start <= end && end <= low_cov_genes->gene_coverage[middle].cds_end) ) {
			return middle;
		} else if ((start <= low_cov_genes->gene_coverage[middle].cds_start) && (low_cov_genes->gene_coverage[middle].cds_end <= end)) {
			return middle;
		} else {
			if (low_cov_genes->gene_coverage[middle].cds_start <= start) {
				low = middle + 1;
			} else {
				high = middle - 1;
			}
		}

		middle = (low + high)/2;
	}

	// outside the while loop means low > high
	return -1;
}

// As the array is sorted, I am going to use binary search for find the location
//
int32_t checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, char **info_in_and_out, uint32_t low_search_index) {
	int32_t found = binarySearch(regions_in, start, end, chrom_idx, low_search_index);

	if (found != -1) {
		uint16_t orig_str_len = strlen(*info_in_and_out);
		uint16_t str_len_needed = strlen(regions_in->gene[chrom_idx][found]) + strlen(regions_in->Synonymous[chrom_idx][found]) + strlen(regions_in->prev_genes[chrom_idx][found]) + 50;
		
		if (str_len_needed > orig_str_len) {
			char *tmp = realloc(*info_in_and_out, str_len_needed);
			if (!tmp) {
				fprintf(stderr, "Memory re-allocation for string failed in checkIntronicRegion\n");
				exit(1);
			}

			*info_in_and_out = tmp;
		} 

		// for debugging only
		if (strcmp(regions_in->gene[chrom_idx][found], ".") == 0) 
			regions_in->gene[chrom_idx][found] = "";
		sprintf(*info_in_and_out, "%s\t%s\t%s", regions_in->gene[chrom_idx][found], regions_in->prev_genes[chrom_idx][found], regions_in->Synonymous[chrom_idx][found]);
	}

	return found;
}

int32_t checkExonRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, char **info_in_and_out, uint32_t low_search_index) {
	int32_t found = binarySearch(regions_in, start, end, chrom_idx, low_search_index);

    if (found != -1) {
        uint16_t orig_str_len = strlen(*info_in_and_out);
        uint16_t str_len_needed = strlen(regions_in->gene[chrom_idx][found]) + strlen(regions_in->Synonymous[chrom_idx][found]) + strlen(regions_in->prev_genes[chrom_idx][found]) + strlen(regions_in->exon_info[chrom_idx][found]) + 50;

        if (str_len_needed > orig_str_len) {
            char *tmp = realloc(*info_in_and_out, str_len_needed);
            if (!tmp) {
                fprintf(stderr, "Memory re-allocation for string failed in checkExonRegion\n");
				exit(1);
            }

            *info_in_and_out = tmp;
        }

        // for debugging only
        if (strcmp(regions_in->gene[chrom_idx][found], ".") == 0)
            regions_in->gene[chrom_idx][found] = "";

		if (strlen(regions_in->exon_info[chrom_idx][found]) > 1) {
			sprintf(*info_in_and_out, "%s\t%s\t%s\t%s", regions_in->gene[chrom_idx][found], regions_in->prev_genes[chrom_idx][found], regions_in->Synonymous[chrom_idx][found], regions_in->exon_info[chrom_idx][found]);
		}
    }

    return found;
}

bool verifyIndex(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t location_index) {
	if (location_index >= regions_in->size_r[chrom_idx]) return false;

	if (regions_in->starts[chrom_idx][location_index] <= start && regions_in->ends[chrom_idx][location_index] > start) {
		return true;
	} else if (regions_in->starts[chrom_idx][location_index] <= end && regions_in->ends[chrom_idx][location_index] > end) {
		return true;
	} else {
		return false;
	}
}

// This is used for gene percentage calculation for the target captured regions.
// Two types of Databases are used, one for static calculation 
// while another one is dynamic calculation
//
void genePercentageCoverageInit(Low_Coverage_Genes *refseq_cds_genes, Low_Coverage_Genes *low_cov_genes, char *chrom_id, Databases *dbs, User_Input *user_inputs) {
	// mysql to obtain total number of distinct gene_name for this specific chromosome
	//
	char *sql = calloc(350, sizeof(char));

	if (user_inputs->annotation_type == 1) {
		// dynamic calculation
		//
		sprintf(sql, "SELECT cds_target_start, cds_target_end, exon_id, exon_count, cds_start, cds_end, cds_length, gene_symbol, gene_name FROM %s WHERE chrom='%s' ORDER BY cds_start, cds_target_start, cds_target_end", dbs->db_dynamic, chrom_id);
	} else {
		// static way
		//
		sprintf(sql, "SELECT cds_target_start, cds_target_end, exon_id, exon_count, cds_start, cds_end, cds_length, gene_symbol, gene_name FROM %s WHERE chrom='%s' ORDER BY cds_start, cds_target_start, cds_target_end", dbs->db_coords, chrom_id);
	}

	if (mysql_query(dbs->con,sql))
		finish_with_error(dbs->con);

	dbs->mysql_results = mysql_store_result(dbs->con);
	if (dbs->mysql_results == NULL)
		finish_with_error(dbs->con);

	MYSQL_ROW row;
	uint32_t counter = 0;

	if (user_inputs->annotation_type == 1) {
		// Dynamic: Update Low_Coverage_Genes variable refseq_cds_genes
		// if I use the array notation, I have to use "[0]." way, if I use point way, I have to use '->' instead
		//
		refseq_cds_genes->total_size = mysql_num_rows(dbs->mysql_results);
		refseq_cds_genes->gene_coverage = calloc(refseq_cds_genes[0].total_size, sizeof(Gene_Coverage));

		while ((row = mysql_fetch_row(dbs->mysql_results))) {
			// update Gene_Coverage variable ==> refseq_cds_genes->gene_coverage
			//
			refseq_cds_genes->gene_coverage[counter].gene_symbol = calloc(strlen(row[7])+1, sizeof(char));
			strcpy(refseq_cds_genes->gene_coverage[counter].gene_symbol, row[7]);

		    refseq_cds_genes->gene_coverage[counter].gene_name = calloc(strlen(row[8])+1, sizeof(char));
			strcpy(refseq_cds_genes->gene_coverage[counter].gene_name, row[8]);

			refseq_cds_genes->gene_coverage[counter].cds_start = (uint32_t) strtol(row[4], NULL, 10);
			refseq_cds_genes->gene_coverage[counter].cds_end   = (uint32_t) strtol(row[5], NULL, 10);
			refseq_cds_genes->gene_coverage[counter].cds_length= (uint32_t) strtol(row[6], NULL, 10);

			refseq_cds_genes->gene_coverage[counter].cds_target_start = (uint32_t) strtol(row[0], NULL, 10);
			refseq_cds_genes->gene_coverage[counter].cds_target_end   = (uint32_t) strtol(row[1], NULL, 10);
			refseq_cds_genes->gene_coverage[counter].exon_id    = (int16_t)  strtol(row[2], NULL, 10);
			refseq_cds_genes->gene_coverage[counter].exon_count = (uint16_t) strtol(row[3], NULL, 10);

			refseq_cds_genes->gene_coverage[counter].low_cov_regions = NULL;
			counter++;
		}
	} else  {
		// Static: Update Low_Coverage_Genes variable low_cov_genes
		//
		low_cov_genes->total_size = mysql_num_rows(dbs->mysql_results);
		low_cov_genes->gene_coverage = calloc(low_cov_genes[0].total_size, sizeof(Gene_Coverage));
		
		 while ((row = mysql_fetch_row(dbs->mysql_results))) {
			 // update Gene_Coverage variable ==> low_cov_genes->gene_coverage
			 //
			 low_cov_genes->gene_coverage[counter].gene_symbol = calloc(strlen(row[7])+1, sizeof(char));
			 strcpy(low_cov_genes->gene_coverage[counter].gene_symbol, row[7]);

			 low_cov_genes->gene_coverage[counter].gene_name = calloc(strlen(row[8])+1, sizeof(char));
			 strcpy(low_cov_genes->gene_coverage[counter].gene_name, row[8]);

			 low_cov_genes->gene_coverage[counter].cds_start = (uint32_t) strtol(row[4], NULL, 10);
			 low_cov_genes->gene_coverage[counter].cds_end   = (uint32_t) strtol(row[5], NULL, 10);
			 low_cov_genes->gene_coverage[counter].cds_length= (uint32_t) strtol(row[6], NULL, 10);

			 low_cov_genes->gene_coverage[counter].cds_target_start = (uint32_t) strtol(row[0], NULL, 10);
			 low_cov_genes->gene_coverage[counter].cds_target_end   = (uint32_t) strtol(row[1], NULL, 10);
			 low_cov_genes->gene_coverage[counter].exon_id    = (uint16_t) strtol(row[2], NULL, 10);
			 low_cov_genes->gene_coverage[counter].exon_count = (uint16_t) strtol(row[3], NULL, 10);

			 if (low_cov_genes->gene_coverage[counter].exon_id <= -1 && strcasecmp(row[7], "SNP") != 0) {
				low_cov_genes->gene_coverage[counter].low_cov_regions = calloc(50, sizeof(char));
				sprintf(low_cov_genes->gene_coverage[counter].low_cov_regions, "%s-%s", row[0], row[1]);
				low_cov_genes->gene_coverage[counter].exon_id = low_cov_genes->gene_coverage[counter].exon_id  * -1;
				low_cov_genes->gene_coverage[counter].num_of_low_cov_bases = low_cov_genes->gene_coverage[counter].cds_target_end - low_cov_genes->gene_coverage[counter].cds_target_start + 1;
			 } else {
				low_cov_genes->gene_coverage[counter].low_cov_regions = NULL;
			 }
			 counter++;
		 }
	}

	// tidy up the memory
	//
	if (dbs->mysql_results) { mysql_free_result(dbs->mysql_results); dbs->mysql_results = NULL; }
	if (sql) { free(sql); sql = NULL; }
	printf("End refseq_cds_genes creation for chromosome %s\n", chrom_id);
}

void genePercentageCoverageDestroy(Low_Coverage_Genes *low_cov_genes) {
	uint32_t i;
	for (i=0; i<low_cov_genes->total_size; i++) {
		//printf("Cleaning Gene symbol %s \n", low_cov_genes->gene_coverage[i].gene_symbol);
		if (low_cov_genes->gene_coverage[i].gene_symbol) {
			free(low_cov_genes->gene_coverage[i].gene_symbol);
			low_cov_genes->gene_coverage[i].gene_symbol = NULL;
		}

		if (low_cov_genes->gene_coverage[i].gene_name) {
			free(low_cov_genes->gene_coverage[i].gene_name);
			low_cov_genes->gene_coverage[i].gene_name = NULL;
		}

		if (low_cov_genes->gene_coverage[i].low_cov_regions) {
			free(low_cov_genes->gene_coverage[i].low_cov_regions);
			low_cov_genes->gene_coverage[i].low_cov_regions=NULL;
		}
	}
	
	if (low_cov_genes->gene_coverage) {
		free(low_cov_genes->gene_coverage);
		low_cov_genes->gene_coverage = NULL;
	}

	if (low_cov_genes) {
		free(low_cov_genes);
		low_cov_genes = NULL;
	}
}

// This is used for the dynamic calculation
// do the real intersect between refseq cds bed regions with target bed regions
//
void intersectTargetsAndRefSeqCDS(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, Low_Coverage_Genes *refseq_cds_genes, Low_Coverage_Genes *low_cov_genes) {
	// find out the index that is used to track current chromosome id
	//
	int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	if (chrom_idx == -1) return;

	if (!refseq_cds_genes || !refseq_cds_genes->total_size || refseq_cds_genes->total_size == 0)
		return;

	if (TARGET_FILE_PROVIDED) {
		// need to allocate the memory for the low_cov_genes variable, 
		// the max size should be the same as refseq_cds_genes->total_size
		//
		low_cov_genes->total_size = refseq_cds_genes->total_size;
		low_cov_genes->gene_coverage = calloc(refseq_cds_genes->total_size, sizeof(Gene_Coverage));

		uint32_t counter=0;
		uint32_t i;

		for(i = 0; i < target_info->size; i++) {
			if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
				continue;

			uint32_t start_t  = target_info->coords[i].start;
			uint32_t end_t    = target_info->coords[i].end;

			if (start_t == 144615130) {
				printf("stop\n");
			}

			// now use binary search to find the index of the intersects
			//
			int32_t found = binarySearchLowCoverage(refseq_cds_genes, start_t, end_t, 0);
			if (found == -1) continue;

			// here we need to re-wound the found position by searching through backward
			// so that the start position on the new low_cov_genes will be in the sorted order
			// Note: we have to use signed integer here as the walk could generate negative values
			//
			int8_t decrease=0;
			int32_t j;

			for (j=found-1; j>=0; j--) {
				// here I need to use cds_coordinates for the checking. The reason is given at binarySearchLowCoverage()
				//                       start_t ================= end_t
				// <<<== start_c ----------- end_c
				//
				uint32_t end_c = refseq_cds_genes->gene_coverage[j].cds_end;
				if (end_c < start_t) break;

				decrease++;
			}

			found -= decrease;

			// if found, we have to walk through both end of the found position as there might be more to intersect
			// update the info based on the found region
			//

			// Walk through the right hand size of the refseq_cds_genes array
			//
			for (j=found; j<refseq_cds_genes->total_size; j++) {
				//  start_t ============= end_t
				//                      start_c ----------- end_c ==>>>
				//
				uint32_t start_c = refseq_cds_genes->gene_coverage[j].cds_start;
				if (end_t < start_c) break;

				// now check the cds_target_coordinates
				//
				uint32_t start_r = refseq_cds_genes->gene_coverage[j].cds_target_start;
				uint32_t end_r = refseq_cds_genes->gene_coverage[j].cds_target_end;

				//                       start_t ================= end_t
				// <<<== start_r ----------- end_r
				//
				if (end_r < start_t) continue;

				//  start_t ============= end_t
				//                      start_r ----------- end_r ==>>>
				//
				if (end_t < start_r) continue;

				//    start_t ========== end_t
				// start_r -------- end_r
				//         start_r ---------- end_r
				// if they intersect, update the low_cov_genes variable
				//
				uint32_t start_tmp = (start_t >= start_r) ? start_t : start_r;
				uint32_t end_tmp = (end_t <= end_r) ? end_t : end_r;

				// because the end position is not included, we need to skip this if start_tmp == end_tmp
				//
				if (end_tmp > start_tmp) {
					recordIntersectedRegions(refseq_cds_genes, low_cov_genes, j, counter, start_tmp, end_tmp);
					counter++;
				}
			}
		}

		// reset the size of the low_cov_genes
		//
		low_cov_genes->total_size = counter;
	}
}

// Record intersected regions
//
void recordIntersectedRegions(Low_Coverage_Genes *refseq_cds_genes, Low_Coverage_Genes *low_cov_genes, int32_t refseq_found_index, uint32_t target_counter, uint32_t start, uint32_t end) {

	low_cov_genes->gene_coverage[target_counter].cds_target_start = start;
	low_cov_genes->gene_coverage[target_counter].cds_target_end   = end;

	low_cov_genes->gene_coverage[target_counter].gene_symbol = calloc(strlen(refseq_cds_genes->gene_coverage[refseq_found_index].gene_symbol)+1, sizeof(char));
	strcpy(low_cov_genes->gene_coverage[target_counter].gene_symbol, refseq_cds_genes->gene_coverage[refseq_found_index].gene_symbol);
	low_cov_genes->gene_coverage[target_counter].gene_name = calloc(strlen(refseq_cds_genes->gene_coverage[refseq_found_index].gene_name)+1, sizeof(char));
	strcpy(low_cov_genes->gene_coverage[target_counter].gene_name, refseq_cds_genes->gene_coverage[refseq_found_index].gene_name);

	low_cov_genes->gene_coverage[target_counter].cds_start = refseq_cds_genes->gene_coverage[refseq_found_index].cds_start;
	low_cov_genes->gene_coverage[target_counter].cds_end   = refseq_cds_genes->gene_coverage[refseq_found_index].cds_end;
	low_cov_genes->gene_coverage[target_counter].cds_length= refseq_cds_genes->gene_coverage[refseq_found_index].cds_length;
	low_cov_genes->gene_coverage[target_counter].exon_id    = refseq_cds_genes->gene_coverage[refseq_found_index].exon_id;
	low_cov_genes->gene_coverage[target_counter].exon_count = refseq_cds_genes->gene_coverage[refseq_found_index].exon_count;
	low_cov_genes->gene_coverage[target_counter].low_cov_regions = NULL;
}

// Here we will initialize the data structure that is used to store information related to transcript coverage percentage
// Two types of Databases are used, one for static calculation 
// while another one is dynamic calculation
//
void transcriptPercentageCoverageInit(char* chrom_id, Transcript_Coverage *transcript_cov, Low_Coverage_Genes *low_cov_genes, User_Input *user_inputs, Databases *dbs) {
	if (!low_cov_genes || !low_cov_genes->total_size || low_cov_genes->total_size == 0)
		return;

	// here is the static way
	//
	if (user_inputs->annotation_type == 2) {
		// mysql to obtain total number of distinct gene_name for this specific chromosome
		//
		char *sql = calloc(200, sizeof(char));
		sprintf(sql, "SELECT COUNT(distinct gene_name) FROM %s WHERE chrom='%s'", dbs->db_coords, chrom_id);

		if (mysql_query(dbs->con,sql))
			finish_with_error(dbs->con);

		dbs->mysql_results = mysql_store_result(dbs->con);
		if (dbs->mysql_results == NULL)
			finish_with_error(dbs->con);

		MYSQL_ROW row;
		while ((row = mysql_fetch_row(dbs->mysql_results))) {
			transcript_cov->num_of_genes = (uint32_t) strtol(row[0], NULL, 10);
			break;
		}

		// Update Transcript_Coverage_Percentage variable xcript_cov_pct
		// if I use the array notation, I have to use "[idx]." way, if I use point way, I have to use '->' instead
		//
		transcript_cov->transcript_cov_pct = calloc(transcript_cov->num_of_genes, sizeof(Transcript_Coverage_Percentage));
		uint32_t i;
		for (i=0; i<transcript_cov->num_of_genes; i++) {
			transcript_cov->transcript_cov_pct[i].gene_symbol = NULL;
			transcript_cov->transcript_cov_pct[i].gene_name = NULL;
			transcript_cov->transcript_cov_pct[i].gene_cov_percentage = 0.0;
		}

		// tidy up the memory
		//
		if (dbs->mysql_results) { mysql_free_result(dbs->mysql_results); dbs->mysql_results = NULL; }
		free(sql);
		sql=NULL;
	} else {
		// First, we need to find out the number of distinct gene_names
		//
		char** refseq_array_tmp = calloc(low_cov_genes->total_size, sizeof(char*));

		// if possible, we don't want to store the same gene_name more than once
		//
		char *prev_gene_name = calloc(50, sizeof(char));
		strcpy(prev_gene_name, ".");

		uint32_t i;
		uint32_t counter=0;

		for (i=0; i<low_cov_genes->total_size; i++) {
			if (strcmp(low_cov_genes->gene_coverage[i].gene_name, prev_gene_name) != 0) {
				refseq_array_tmp[counter] = calloc(strlen(low_cov_genes->gene_coverage[i].gene_name), sizeof(char));
				strcpy(refseq_array_tmp[counter], low_cov_genes->gene_coverage[i].gene_name);
				strcpy(prev_gene_name, low_cov_genes->gene_coverage[i].gene_name);
				counter++;
			}
		}

		// Now we know the size, declare a new variable with the exact same size
		//
		char** refseq_array = calloc(counter, sizeof(char*));
		for (i=0; i<counter; i++) {
			refseq_array[i] = calloc(strlen(refseq_array_tmp[i]), sizeof(char));
			strcpy(refseq_array[i], refseq_array_tmp[i]);
		}

		// print gene_name string array first before sorting:
		//
		//print_string_array(refseq_array, low_cov_genes->total_size);
     
		// Secondly, quick sort the gene_name string array
		//
		qsort(refseq_array, counter, sizeof(char *), compare3);
	
		// print sorted gene_name string array after sorting
		//
		//print_string_array(refseq_array, low_cov_genes->total_size);

		// Next, we need to eliminate the duplicates to obtain unique gene_names
		//
		strcpy(prev_gene_name, refseq_array[0]);
		uint32_t num_of_uniq_gene_names=1;

		for (i=1; i<counter; i++) {
			if (strcmp(prev_gene_name, refseq_array[i])!=0) {
				// different refseq names, updating...
				//
				num_of_uniq_gene_names++;
				strcpy(prev_gene_name, refseq_array[i]);
			}
		}

		printf("Number of unique refseq names is %d\n", num_of_uniq_gene_names);
		transcript_cov->num_of_genes = num_of_uniq_gene_names;

		// Update Transcript_Coverage_Percentage variable xcript_cov_pct
		// if I use the array notation, I have to use "[idx]." way, if I use point way, I have to use '->' instead
		//
	   	transcript_cov->transcript_cov_pct = calloc(transcript_cov->num_of_genes, sizeof(Transcript_Coverage_Percentage));
		for (i=0; i<transcript_cov->num_of_genes; i++) {
			transcript_cov->transcript_cov_pct[i].gene_symbol = NULL;
			transcript_cov->transcript_cov_pct[i].gene_name = NULL;
			transcript_cov->transcript_cov_pct[i].gene_cov_percentage = 0.0;
		}

		// tidy up the memory
		//
		free(prev_gene_name);
		prev_gene_name=NULL;

		for (i=0; i<low_cov_genes->total_size; i++) {
			free(refseq_array_tmp[i]);
			refseq_array_tmp[i]=NULL;
		}

		for (i=0; i<counter; i++) {
			free(refseq_array[i]);
			refseq_array[i]=NULL;
		}
	}
}

void transcriptPercentageCoverageDestroy(Transcript_Coverage *transcript_cov) {
	uint32_t i;
	for (i=0; i<transcript_cov->num_of_genes; i++) {
		if (transcript_cov->transcript_cov_pct[i].gene_symbol) {
			free(transcript_cov->transcript_cov_pct[i].gene_symbol);
			transcript_cov->transcript_cov_pct[i].gene_symbol=NULL;
		}

		if (transcript_cov->transcript_cov_pct[i].gene_name) {
			free(transcript_cov->transcript_cov_pct[i].gene_name);
			transcript_cov->transcript_cov_pct[i].gene_name=NULL;
		}
	}

	if (transcript_cov->transcript_cov_pct) {
		free(transcript_cov->transcript_cov_pct);
		transcript_cov->transcript_cov_pct=NULL;
	}

	if (transcript_cov) {
		free(transcript_cov);
		transcript_cov=NULL;
	}
}

// Here we are not going to go through the database, instead, we are going to use the Low_Coverage_Genes stored on the heap!!!
void produceGenePercentageCoverageInfo(uint32_t start_in, uint32_t stop_in, char *chrom_id, Low_Coverage_Genes *low_cov_genes) {
	// print low_cov_genes for debugging
	//printLowCoverageGeneStructure(low_cov_genes);
	//
	//if (start_in == 89623861) {
	//	printf("stop\n");
	//}
	
	int32_t found = binarySearchLowCoverage(low_cov_genes, start_in, stop_in, 0);

	if (found == -1) return;

	// if found, we have to roll back to find the start position, so that everything will be in the sorted order by start position
	//
	uint8_t decrease=0;
	int32_t i;

	for (i=found-1; i>=0; i--) {
		// here I need to use cds_start and cds_end
		//
		if ((low_cov_genes->gene_coverage[i].cds_start <= start_in && start_in < low_cov_genes->gene_coverage[i].cds_end) ||
               (low_cov_genes->gene_coverage[i].cds_start <= stop_in && stop_in < low_cov_genes->gene_coverage[i].cds_end) ||
               (start_in <= low_cov_genes->gene_coverage[i].cds_start && low_cov_genes->gene_coverage[i].cds_end <= stop_in)) {
			decrease++;
        } else {
            break;
        }
    }

	found -= decrease;

	// walk through the low_cov_genes->gene_coverage on the right hand size
	//
	for (i=found; i<low_cov_genes->total_size; i++) {
		// we will exit only when the cds_start is larger than stop_in
		//
		if (low_cov_genes->gene_coverage[i].cds_start > stop_in)
			break;

		//printf("index %"PRIu32"\tsymbol %s and gene to match %s\n", i, low_cov_genes->gene_coverage[i].gene_symbol, row[4]);
		if ((low_cov_genes->gene_coverage[i].cds_target_start <= start_in && start_in < low_cov_genes->gene_coverage[i].cds_target_end) ||
			   (low_cov_genes->gene_coverage[i].cds_target_start <= stop_in && stop_in < low_cov_genes->gene_coverage[i].cds_target_end) || 
			   (start_in <= low_cov_genes->gene_coverage[i].cds_target_start && low_cov_genes->gene_coverage[i].cds_target_end <= stop_in))	{
			processExonArrays(low_cov_genes, i, start_in, stop_in);
		} else {
			continue;
		}
	}
}

// This will calculate the number of overlap count to find out the true low coverage region
// here start is the low coverage region start
// end is the low coverage region end
//
void processExonArrays(Low_Coverage_Genes *low_cov_genes, uint32_t refseq_exon_index, uint32_t start, uint32_t end) {
	uint32_t cds_target_start = low_cov_genes->gene_coverage[refseq_exon_index].cds_target_start;
	uint32_t cds_target_end   = low_cov_genes->gene_coverage[refseq_exon_index].cds_target_end;
	uint32_t report_string_size = 50;

	// now we need to allocate memory space for low_cov_regions string for exon percentage report
	//
	if (!low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions) {
		low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions = calloc(report_string_size, sizeof(char));
		strcpy(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions, "");
	} else {
		report_string_size += strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions);
		char *tmp = realloc(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions, report_string_size);
		if (!tmp) {
			fprintf(stderr, "Memory re-allocation failed at processExonArrays\n");
			exit(1);
		}
		low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions = tmp;
	}

	// To begin with, we presume that there is no low coverage bases at all. ie, everything is high coverage
	//
	if (start <= cds_target_start && cds_target_end <= end) {
		//               cds_target_start ============== cds_target_end
		//  start ---------------------------------------------------------- end
		//
        //low_cov_genes->gene_coverage[refseq_exon_index].num_of_low_cov_bases += cds_target_end - cds_target_start + 1;
        low_cov_genes->gene_coverage[refseq_exon_index].num_of_low_cov_bases += cds_target_end - cds_target_start;

		if (strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions) == 0) {
			sprintf(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions, "%"PRIu32"-%"PRIu32, cds_target_start, cds_target_end);
		} else {
			sprintf(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions + strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions), ";%"PRIu32"-%"PRIu32, cds_target_start, cds_target_end);
		}
    } else if (cds_target_start <= start && end <= cds_target_end) {
		//  cds_target_start =================== cds_target_end
		//               start ------------- end
		//
        //low_cov_genes->gene_coverage[refseq_exon_index].num_of_low_cov_bases += end - start + 1; 
        low_cov_genes->gene_coverage[refseq_exon_index].num_of_low_cov_bases += end - start; 

		if (strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions) == 0) {
            sprintf(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions, "%"PRIu32"-%"PRIu32, start, end);
        } else {
            sprintf(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions + strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions), ";%"PRIu32"-%"PRIu32, start, end);
        }

    } else if (cds_target_start <= end && end <= cds_target_end) {
		//             cds_target_start ============== cds_target_end
		//  start ----------------------------- end
		//
        //low_cov_genes->gene_coverage[refseq_exon_index].num_of_low_cov_bases += end - cds_target_start + 1;
        low_cov_genes->gene_coverage[refseq_exon_index].num_of_low_cov_bases += end - cds_target_start;

		if (strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions) == 0) {
            sprintf(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions, "%"PRIu32"-%"PRIu32, cds_target_start, end);
        } else {
            sprintf(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions + strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions), ";%"PRIu32"-%"PRIu32, cds_target_start, end);
        }

    } else if (cds_target_start <= start && start <= cds_target_end) {
		// cds_target_start ================= cds_target_end
		//                     start ----------------------- end
		//
        //low_cov_genes->gene_coverage[refseq_exon_index].num_of_low_cov_bases += cds_target_end - start + 1;
        low_cov_genes->gene_coverage[refseq_exon_index].num_of_low_cov_bases += cds_target_end - start;

		if (strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions) == 0) {
            sprintf(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions, "%"PRIu32"-%"PRIu32, start, cds_target_end);
        } else {
            sprintf(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions + strlen(low_cov_genes->gene_coverage[refseq_exon_index].low_cov_regions), ";%"PRIu32"-%"PRIu32, start, cds_target_end);
        }
    }
}
