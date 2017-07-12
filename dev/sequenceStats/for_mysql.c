/*
 * =====================================================================================
 *
 *       Filename:  for_mysql.c
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
	free(pre_sql);
	pre_sql = NULL;

	return count;
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions, type 4 for all_site_reports 
void regionsSkipMySQLInit(MYSQL *con, Regions_Skip_MySQL *regions_in, bam_hdr_t *header, uint8_t type) {
	uint32_t i, j;

	regions_in->starts = calloc(25, sizeof(uint32_t*));
	regions_in->ends   = calloc(25, sizeof(uint32_t*));
	regions_in->size_r = calloc(25, sizeof(uint32_t));
	regions_in->chromosome_ids = calloc(25, sizeof(char*));

	// for exon and intronic regions
	if (type > 1) {
        regions_in->gene = calloc(25, sizeof(char**));
        if (!regions_in->gene) {
            fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Gene");
            exit(1);
        }

        regions_in->Synonymous = calloc(25, sizeof(char**));
	    if (!regions_in->Synonymous) {
		    fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Synonymous");
			exit(1);
        }

	    regions_in->prev_genes = calloc(25, sizeof(char**));
		if (!regions_in->prev_genes) {
			fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Prev_genes");
	        exit(1);
        }

		// for exon regions only
		if (type == 3) {
			regions_in->exon_info = calloc(25, sizeof(char**));
	        if (!regions_in->exon_info) {
		        fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Gene");
			    exit(1);
			}
        }
	}

	for(i=0; i<header->n_targets; i++) {
		if (i >= 25) break;

		regions_in->size_r[i] = fetchTotalCount(type, con, header->target_name[i]);
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

		regions_in->chromosome_ids[i] = calloc(strlen(header->target_name[i])+1, sizeof(char));
		strcpy(regions_in->chromosome_ids[i], header->target_name[i]);
	
		// for exon and intronic regions
		if (type > 1) {
			regions_in->gene[i] = calloc(regions_in->size_r[i], sizeof(char*));
            regions_in->Synonymous[i] = calloc(regions_in->size_r[i], sizeof(char*));
	        regions_in->prev_genes[i] = calloc(regions_in->size_r[i], sizeof(char*));

			for (j=0; j<regions_in->size_r[i]; j++) {
				//regions_in->gene[i][j] = NULL;
				regions_in->gene[i][j] = calloc(1, sizeof(char));

				regions_in->Synonymous[i][j] = calloc(1, sizeof(char));
				//regions_in->Synonymous[i][j] = NULL;
				regions_in->prev_genes[i][j] = calloc(1, sizeof(char));
				//regions_in->prev_genes[i][j] = NULL;
			}

			// for exon regions only
			if (type == 3) {
				regions_in->exon_info[i] = calloc(regions_in->size_r[i], sizeof(char*));

				for (j=0; j<regions_in->size_r[i]; j++) {
					//regions_in->exon_info[i][j] = NULL;
					regions_in->exon_info[i][j] = calloc(1, sizeof(char));
				}
			}
		}

		// now populate these regions without gene annotation
		populateStaticRegionsForOneChromOnly(regions_in, con, header->target_name[i], i, type); 
	}
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions, type 4 for all_site_reports 
void regionsSkipMySQLDestroy(Regions_Skip_MySQL *regions_in, uint8_t type) {
	uint32_t i, j;

	for (i=0; i<25; i++) {
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
void populateStaticRegionsForOneChromOnly(Regions_Skip_MySQL *regions_in, MYSQL *con, char *chrom_id, uint32_t index, uint8_t type) {
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
        regions_in->starts[index][count] = (uint32_t) atol(row[0]);
        regions_in->ends[index][count] = (uint32_t) atol(row[1]);

		if (type > 1) {
			regions_in->gene[index][count] = dynamicStringAllocation(row[2], regions_in->gene[index][count]) ;
			regions_in->Synonymous[index][count] = dynamicStringAllocation(row[3], regions_in->Synonymous[index][count]);
			regions_in->prev_genes[index][count] = dynamicStringAllocation(row[4], regions_in->prev_genes[index][count]);

			if (type == 3) {
				char *tmp = calloc(strlen(row[5]) + 50, sizeof(char));
				sprintf(tmp, "%s", row[5]);
				regions_in->exon_info[index][count] = dynamicStringAllocation(tmp, regions_in->exon_info[index][count]);
				if (tmp) {
					free(tmp);
					tmp = NULL;
				}
			}
		}

		count++;
	}
    free(sql);
	sql = NULL;
}

// As the array is sorted, I am going to use binary search for find the location
// here index is the chromosome index
int32_t checkInterGenicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, uint32_t low_index) {
	int32_t found = binary_search(regions_in, start, index, low_index);

	if (found == -1 && end > start) 
		found = binary_search(regions_in, end, index, low_index);

	return found;
}

// here index is the chromosome index
// Note: for bed file format, the end position is not part of the region to be checked!
int32_t binary_search(Regions_Skip_MySQL *regions_in, uint32_t pos, uint32_t index, uint32_t low_index) {
	int32_t low = low_index;
	int32_t high = regions_in->size_r[index] - 1;
	int32_t middle = (low + high)/2;

	while (low <= high) {
		if (regions_in->starts[index][middle] <= pos && pos < regions_in->ends[index][middle]) {
			return middle;
		} else {
			if (regions_in->starts[index][middle] < pos) {
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
// here index is the chromosome index
int32_t checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, char **info_in_and_out, uint32_t low_index) {
	int32_t found = binary_search(regions_in, start, index, low_index);
	if (found == -1 && end > start) {
		found = binary_search(regions_in, end, index, low_index);
	}

	if (found != -1) {
		uint16_t orig_str_len = strlen(*info_in_and_out);
		uint16_t str_len_needed = strlen(regions_in->gene[index][found]) + strlen(regions_in->Synonymous[index][found]) + strlen(regions_in->prev_genes[index][found]) + 50;
		
		if (str_len_needed > orig_str_len) {
			char *tmp = realloc(*info_in_and_out, str_len_needed);
			if (!tmp) {
				fprintf(stderr, "Memory re-allocation for string failed in checkIntronicRegion\n");
				exit(1);
			}

			*info_in_and_out = tmp;
		} 

		// for debugging only
		if (strcmp(regions_in->gene[index][found], ".") == 0) 
			regions_in->gene[index][found] = "";
		sprintf(*info_in_and_out, "%s\t%s\t%s", regions_in->gene[index][found], regions_in->prev_genes[index][found], regions_in->Synonymous[index][found]);
	}

	return found;
}

int32_t checkExonRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, char **info_in_and_out, uint32_t low_index) {
	int32_t found = binary_search(regions_in, start, index, low_index);

	if (found == -1 && end > start) {
		found = binary_search(regions_in, end, index, low_index);
	}

    if (found != -1) {
        uint16_t orig_str_len = strlen(*info_in_and_out);
        uint16_t str_len_needed = strlen(regions_in->gene[index][found]) + strlen(regions_in->Synonymous[index][found]) + strlen(regions_in->prev_genes[index][found]) + strlen(regions_in->exon_info[index][found]) + 50;

        if (str_len_needed > orig_str_len) {
            char *tmp = realloc(*info_in_and_out, str_len_needed);
            if (!tmp) {
                fprintf(stderr, "Memory re-allocation for string failed in checkExonRegion\n");
				exit(1);
            }

            *info_in_and_out = tmp;
        }

        // for debugging only
        if (strcmp(regions_in->gene[index][found], ".") == 0)
            regions_in->gene[index][found] = "";

		if (strlen(regions_in->exon_info[index][found]) > 1) {
			sprintf(*info_in_and_out, "%s\t%s\t%s\t%s\t.\t.\t.", regions_in->gene[index][found], regions_in->prev_genes[index][found], regions_in->Synonymous[index][found], regions_in->exon_info[index][found]);
		}
    }

    return found;
}

bool verifyIndex(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t location_index) {
	if (location_index >= regions_in->size_r[chrom_idx]) return false;

	if (regions_in->starts[chrom_idx][location_index] <= start && regions_in->ends[chrom_idx][location_index] >= end) {
		return true;
	} else {
		return false;
	}
}
