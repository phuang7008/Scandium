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
#include "annotation.h"
#include "utils.h"

// The followings are codes that related to gene annotation
void finish_with_error(MYSQL *con)
{
    fprintf(stderr, "%s\n", mysql_error(con));
    mysql_close(con);
    exit(1);
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions
uint32_t fetchTotalCount(uint8_t type, MYSQL *con, char *chrom_id) {
	char *pre_sql = calloc(100, sizeof(char));
	if (type == 1) {
		sprintf(pre_sql, "SELECT count(*) from Intergenic_Regions WHERE chrom='%s'", chrom_id);
	} else if (type == 2) {
		sprintf(pre_sql, "SELECT count(*) from Intron_Regions WHERE chrom='%s'", chrom_id);
	} else if (type == 3) {
		sprintf(pre_sql, "SELECT count(*) from Exon_Regions WHERE chrom='%s'", chrom_id);
	}
	//printf("%s\n", pre_sql);

	if (mysql_query(con,pre_sql))
        finish_with_error(con);

    MYSQL_RES *result = mysql_store_result(con);
    if (result == NULL)
        finish_with_error(con);

    MYSQL_ROW row;
	uint32_t count=0;
    while ((row = mysql_fetch_row(result))) {
		if (row[0] > 0)
			count = (uint32_t) atol(row[0]);
	}

	if (result) { mysql_free_result(result); result = NULL; }

	free(pre_sql);
	pre_sql = NULL;

	return count;
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions, type 4 for all_site_reports 
void regionsSkipMySQLInit(MYSQL *con, Regions_Skip_MySQL *regions_in, uint8_t type) {
	// here we need to find out how many chromosomes we are dealing with
	char *sql = calloc(250, sizeof(char));
    if (type == 1) {
        sprintf(sql, "SELECT DISTINCT chrom FROM Intergenic_Regions WHERE chrom NOT LIKE '%%\\_%%'");
    } else if (type == 2) {
        sprintf(sql, "SELECT DISTINCT chrom FROM Intron_Regions WHERE chrom NOT LIKE '%%\\_%%'");
    } else {
        sprintf(sql, "SELECT DISTINCT chrom FROM Exon_Regions WHERE chrom NOT LIKE '%%\\_%%'");
    }

	if (mysql_query(con,sql))
        finish_with_error(con);

	MYSQL_RES *result = mysql_store_result(con);
    if (result == NULL)
        finish_with_error(con);

	regions_in->chrom_list_size = mysql_num_rows(result);
	regions_in->chromosome_ids = calloc(regions_in->chrom_list_size, sizeof(char*));

	uint16_t chrom_idx=0;
    MYSQL_ROW row;
    while ((row = mysql_fetch_row(result))) {
		regions_in->chromosome_ids[chrom_idx] = calloc(strlen(row[0])+1, sizeof(char));
		strcpy(regions_in->chromosome_ids[chrom_idx], row[0]);
		chrom_idx++;
	}
	if (result) { mysql_free_result(result); result = NULL; }
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
		regions_in->size_r[i] = fetchTotalCount(type, con, regions_in->chromosome_ids[i]);
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
		populateStaticRegionsForOneChromOnly(regions_in, con, regions_in->chromosome_ids[i], i, type); 
	}
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions, type 4 for all_site_reports 
void regionsSkipMySQLDestroy(Regions_Skip_MySQL *regions_in, uint8_t type) {
	uint32_t i, j;

	for (i=0; i<regions_in->chrom_list_size; i++) {
		// for exon and intronic regions
		if (type > 1) {
			for (j=0; j<regions_in->size_r[i]; j++) {
				if (regions_in->gene[i][j]) free(regions_in->gene[i][j]);
				if (regions_in->Synonymous[i][j]) free(regions_in->Synonymous[i][j]);
				if (regions_in->prev_genes[i][j]) free(regions_in->prev_genes[i][j]);
			}

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

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions, type 4 for all_site_reports 
void populateStaticRegionsForOneChromOnly(Regions_Skip_MySQL *regions_in, MYSQL *con, char *chrom_id, uint32_t chrom_idx, uint8_t type) {
	char *sql = calloc(250, sizeof(char));
	if (type == 1) {
		sprintf(sql, "SELECT start, end from Intergenic_Regions WHERE chrom='%s' ORDER BY start", chrom_id);
	} else if (type == 2) {
		sprintf(sql, "SELECT start, end, gene_symbol, Synonymous, prev_gene_symbol from Intron_Regions WHERE chrom='%s' ORDER BY start", chrom_id);
	} else {
		sprintf(sql, "SELECT start, end, gene_symbol, Synonymous, prev_gene_symbol, annotation from Exon_Regions WHERE chrom='%s' ORDER BY start", chrom_id);
	}
	//printf("%s\n", sql);

	if (mysql_query(con,sql))
        finish_with_error(con);

    MYSQL_RES *result = mysql_store_result(con);
    if (result == NULL)
        finish_with_error(con);

    MYSQL_ROW row;
    uint32_t count=0;
	while ((row = mysql_fetch_row(result))) {
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

	if (result) { mysql_free_result(result); result = NULL; }

    free(sql);
	sql = NULL;
}

// As the array is sorted, I am going to use binary search for find the location
int32_t checkInterGenicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t low_search_index) {
	int32_t found = binary_search(regions_in, start, chrom_idx, low_search_index);

	if (found == -1 && end > start) 
		found = binary_search(regions_in, end, chrom_idx, low_search_index);

	return found;
}

// Note: for bed file format, the end position is not part of the region to be checked!
int32_t binary_search(Regions_Skip_MySQL *regions_in, uint32_t pos, uint32_t chrom_idx, uint32_t low_search_index) {
	int32_t low = low_search_index;
	int32_t high = regions_in->size_r[chrom_idx] - 1;
	int32_t middle = (low + high)/2;

	while (low <= high) {
		if (regions_in->starts[chrom_idx][middle] <= pos && pos < regions_in->ends[chrom_idx][middle]) {
			return middle;
		} else {
			if (regions_in->starts[chrom_idx][middle] < pos) {
				low = middle + 1;
			} else {
				high = middle - 1;
			}
		}

		middle = (low + high)/2;
	}

	// outside the while loop means low > high
	//if (low > high)
	//	return -1;

	return -1;
}

// As the array is sorted, I am going to use binary search for find the location
int32_t checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, char **info_in_and_out, uint32_t low_search_index) {
	int32_t found = binary_search(regions_in, start, chrom_idx, low_search_index);
	if (found == -1 && end > start) {
		found = binary_search(regions_in, end, chrom_idx, low_search_index);
	}

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
	//printf("chromsome index is %"PRIu32"\n", chrom_idx);
	int32_t found = binary_search(regions_in, start, chrom_idx, low_search_index);

	if (found == -1 && end > start) {
		found = binary_search(regions_in, end, chrom_idx, low_search_index);
	}

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

void genePercentageCoverageInit(Low_Coverage_Genes *low_cov_genes, char *chrom_id, MYSQL *con) {
	// first find out how many distinct gene_symbol this chromosome has 
	// and for each gene_symbol, how many distinct refseq_name it has
	char *sql = calloc(250, sizeof(char));
	sprintf(sql, "SELECT gene_symbol, COUNT(DISTINCT refseq_name) FROM Gene_Exons WHERE chrom='%s' GROUP BY gene_symbol ORDER BY gene_symbol, refseq_name", chrom_id);

	if (mysql_query(con,sql))
		finish_with_error(con);

	MYSQL_RES *result = mysql_store_result(con);
	if (result == NULL)
	    finish_with_error(con);

	// Update Low_Coverage_Genes variable low_cov_genes
	// allocate memory space and set member values for chrom_id, gene_pct_coverage and num_of_gene_symbol
	// if I use the array notation, I have to use "[0]." way, if I use point way, I have to use '->' instead
	low_cov_genes->num_of_gene_symbol = mysql_num_rows(result);
	low_cov_genes->chrom_id = calloc(strlen(chrom_id) + 1, sizeof(char));
	strcpy(low_cov_genes->chrom_id, chrom_id);
	low_cov_genes->gene_pct_coverage = calloc(low_cov_genes[0].num_of_gene_symbol, sizeof(Gene_Percentage_Coverage));

	MYSQL_ROW row;
    uint32_t counter = 0;
	while ((row = mysql_fetch_row(result))) {
		// update Gene_Percentage_Coverage variable ==> low_cov_genes->gene_pct_coverage
        low_cov_genes->gene_pct_coverage[counter].gene_symbol = calloc(strlen(row[0]) + 1, sizeof(char));
		strcpy(low_cov_genes->gene_pct_coverage[counter].gene_symbol, row[0]);

        low_cov_genes->gene_pct_coverage[counter].num_of_refseq_names = (uint16_t) strtol(row[1], NULL, 10);
		low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper = calloc((uint16_t) strtol(row[1], NULL, 10), sizeof(Exon_Wrapper));

		counter++;
    }

	// here we are going to query one gene symbol at a time
	for (counter=0; counter<low_cov_genes->num_of_gene_symbol; counter++) {
		if (counter == 771) {
			printf("stop\n");
		}
		// Need to clean MYSQL result here otherwise, valgrind will give possible memory leak error
		if (result) { 
			mysql_free_result(result); 
			result = NULL; 
		}
		strcpy(sql, "");

		// now for new query to retrieve exon info based on gene symbol 
		sprintf(sql, "SELECT start, end, exon_count, exon_id, refseq_name FROM Gene_Exons WHERE chrom='%s' and gene_symbol='%s' ORDER BY refseq_name, exon_id", chrom_id, low_cov_genes->gene_pct_coverage[counter].gene_symbol);

		if (mysql_query(con,sql))
			finish_with_error(con);

		MYSQL_RES *result = mysql_store_result(con);
		if (result == NULL)
			finish_with_error(con);

		char *prev_refseq_name = calloc(50, sizeof(char));
		strcpy(prev_refseq_name, "");

		int16_t refseq_index=-1;
		while ((row = mysql_fetch_row(result))) {
			uint16_t exon_count = (uint16_t) strtol(row[2], NULL, 10);
			uint16_t exon_id = (uint16_t) strtol(row[3], NULL, 10);			// or using atol(row[3])

			// check if current refseq name equal to the previous refseq name
			if (strcmp(row[4], prev_refseq_name) == 0) {
				low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons[exon_id].exon_id = exon_id;
                low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons[exon_id].exon_start = (uint32_t) strtol(row[0], NULL, 10);
                low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons[exon_id].exon_end   = (uint32_t) strtol(row[1], NULL, 10);
				//printf("Same refseq %s with id %d out of %d\n", row[4], exon_id, exon_count);
			} else {
				refseq_index++;
				//printf("%s counter is %d and refseq index is %d with name  %s to hold %d\n",low_cov_genes->gene_pct_coverage[counter].gene_symbol, counter, refseq_index, row[4], low_cov_genes->gene_pct_coverage[counter].num_of_refseq_names);

				// update Exon_Wrapper variable ==> low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper
				low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].num_of_exons = exon_count;
				low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].refseq_name = calloc(strlen(row[4]) + 1, sizeof(char));
				strcpy(low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].refseq_name, row[4]);
				strcpy(prev_refseq_name, row[4]);

				// update the Exon_Details variable ==> low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons
				low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons = calloc(exon_count, sizeof(Exon_Details));
				low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons[exon_id].exon_id = exon_id;
				low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons[exon_id].exon_start = (uint32_t) strtol(row[0], NULL, 10);
				low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons[exon_id].exon_end   = (uint32_t) strtol(row[1], NULL, 10);
				low_cov_genes->gene_pct_coverage[counter].refseq_exon_wrapper[refseq_index].exons[exon_id].num_of_low_cov_bases = 0;
			}
		}
		if (prev_refseq_name) { free(prev_refseq_name); prev_refseq_name = NULL; }
	}
	printf("End for chromosome %s\n", chrom_id);
	if (sql) { free(sql); sql = NULL; }
}

void genePercentageCoverageDestroy(Low_Coverage_Genes *low_cov_genes, char *chrom_id) {
	if (low_cov_genes->chrom_id) {
		free(low_cov_genes->chrom_id);
		low_cov_genes->chrom_id = NULL;
	}

	uint32_t i, j;
	for (i=0; i<low_cov_genes->num_of_gene_symbol; i++) {
		//printf("Cleaning Gene symbol %s \n", low_cov_genes->gene_pct_coverage[i].gene_symbol);
		if (low_cov_genes->gene_pct_coverage[i].gene_symbol) {
			free(low_cov_genes->gene_pct_coverage[i].gene_symbol);
			low_cov_genes->gene_pct_coverage[i].gene_symbol = NULL;
		}

		for (j=0; j<low_cov_genes->gene_pct_coverage[i].num_of_refseq_names; j++) {
			//printf("Cleaning refseq %s\n", low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].refseq_name);
			if (low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].refseq_name) {
				free(low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].refseq_name);
				low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].refseq_name = NULL;
			}

			if (low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].exons) {
				free(low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].exons);
				low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].exons = NULL;
			}
		}
	}
	
	if (low_cov_genes->gene_pct_coverage) {
		free(low_cov_genes->gene_pct_coverage);
		low_cov_genes->gene_pct_coverage = NULL;
	}

	if (low_cov_genes) {
		free(low_cov_genes);
		low_cov_genes = NULL;
	}
}

void produceGenePercentageCoverageInfo(uint32_t start_in, uint32_t stop_in, char *chrom_id, MYSQL *con, Low_Coverage_Genes *low_cov_genes) {
	// print low_cov_genes for debugging
	//printLowCoverageGeneStructure(low_cov_genes);

	char *sql   = calloc(350, sizeof(char)); 
	sprintf(sql,  "SELECT DISTINCT start, end, exon_count, exon_id, gene_symbol, refseq_name FROM Gene_Exons WHERE  chrom='%s' AND ((start <= %"PRIu32" AND %"PRIu32" <= end) OR (start <= %"PRIu32" AND %"PRIu32" <= end))", chrom_id, start_in, start_in, stop_in, stop_in);

	if (mysql_query(con,sql))
	finish_with_error(con);

    MYSQL_RES *result = mysql_store_result(con);
	if (result == NULL)
		finish_with_error(con);

	uint32_t i, gene_symbol_index;

	MYSQL_ROW row;
	uint16_t exon_count=0;

	while ((row = mysql_fetch_row(result))) {
		// find the index for current refseq index based on the gene symbol retrieved
		uint32_t refseq_exon_index;
		for (i=0; i<low_cov_genes->num_of_gene_symbol; i++) {
			//printf("index %"PRIu32"\tsymbol %s and gene to match %s\n", i, low_cov_genes->gene_pct_coverage[i].gene_symbol, row[4]);
			if ((strcmp(low_cov_genes->gene_pct_coverage[i].gene_symbol, row[4]) == 0)) {
				gene_symbol_index = i;
		
				uint16_t j;
				for (j=0; j<low_cov_genes->gene_pct_coverage[i].num_of_refseq_names; j++) {
					if (strcmp(low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].refseq_name, row[5]) == 0) {
						refseq_exon_index = j;
						exon_count = low_cov_genes->gene_pct_coverage[i].refseq_exon_wrapper[j].num_of_exons;
						processExonArrays(exon_count, low_cov_genes, gene_symbol_index, refseq_exon_index, start_in, stop_in);
					}
				}
				break;
			}
		}
	}

	if (result) { mysql_free_result(result); result=NULL; }
	//if (prev_sql) { free(prev_sql); prev_sql=NULL; }
	//if (mid_sql)  { free(mid_sql); mid_sql=NULL; }
	if (sql) { free(sql); sql=NULL; }
}

// This will return the number of overlap count
void processExonArrays(uint16_t exon_count, Low_Coverage_Genes *low_cov_genes, uint32_t gene_symbol_index, uint16_t refseq_exon_index, uint32_t start, uint32_t end) {
	uint32_t i;

	for (i=0; i<exon_count; i++) {
		uint32_t db_exon_start = low_cov_genes->gene_pct_coverage[gene_symbol_index].refseq_exon_wrapper[refseq_exon_index].exons[i].exon_start;
		uint32_t db_exon_end   = low_cov_genes->gene_pct_coverage[gene_symbol_index].refseq_exon_wrapper[refseq_exon_index].exons[i].exon_end;

		if (start <= db_exon_start && db_exon_end <= end) {
			low_cov_genes->gene_pct_coverage[gene_symbol_index].refseq_exon_wrapper[refseq_exon_index].exons[i].num_of_low_cov_bases += 
				db_exon_end - db_exon_start + 1;

		} else if (db_exon_start <= start && end <= db_exon_end) {
			low_cov_genes->gene_pct_coverage[gene_symbol_index].refseq_exon_wrapper[refseq_exon_index].exons[i].num_of_low_cov_bases += end - start + 1;

		} else if (db_exon_start <= end && end <= db_exon_end) {
			low_cov_genes->gene_pct_coverage[gene_symbol_index].refseq_exon_wrapper[refseq_exon_index].exons[i].num_of_low_cov_bases += 
				end - db_exon_start + 1;

		} else if (db_exon_start <= start && start <= db_exon_end) {
			low_cov_genes->gene_pct_coverage[gene_symbol_index].refseq_exon_wrapper[refseq_exon_index].exons[i].num_of_low_cov_bases += 
			   db_exon_end - start + 1;

		}
	}
}
