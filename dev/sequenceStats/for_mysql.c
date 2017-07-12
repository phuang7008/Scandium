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
#include "for_mysql.h"
#include "utils.h"

// The followings are codes that related to gene annotation
void finish_with_error(MYSQL *con)
{
    fprintf(stderr, "%s\n", mysql_error(con));
    mysql_close(con);
    exit(1);
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions, type 4 for all_site_reports 
uint32_t fetchTotalCount(uint8_t type, MYSQL *con, char *chrom_id) {
	char *pre_sql = calloc(100, sizeof(char));
	if (type == 1) {
		sprintf(pre_sql, "SELECT count(*) from Regions_With_No_Annotation WHERE chrom='%s'", chrom_id);
	} else if (type == 2) {
		sprintf(pre_sql, "SELECT count(*) from Intron_Regions WHERE chrom='%s'", chrom_id);
	} else if (type == 3) {
		sprintf(pre_sql, "SELECT count(*) from Unique_Gene_Exons WHERE chrom='%s'", chrom_id);
	} else {
		sprintf(pre_sql, "SELECT count(*) from All_Site_Reports WHERE chrom='%s'", chrom_id);
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

		if (type < 4) {
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

			if (type < 4) {
	            regions_in->Synonymous[i] = calloc(regions_in->size_r[i], sizeof(char*));
		        regions_in->prev_genes[i] = calloc(regions_in->size_r[i], sizeof(char*));
			}

			for (j=0; j<regions_in->size_r[i]; j++) {
				regions_in->gene[i][j] = NULL;
				//regions_in->gene[i][j] = calloc(500, sizeof(char));

				if (type < 4) {
					//regions_in->Synonymous[i][j] = calloc(250, sizeof(char));
					regions_in->Synonymous[i][j] = NULL;
				    //regions_in->prev_genes[i][j] = calloc(250, sizeof(char));
				    regions_in->prev_genes[i][j] = NULL;
				}
			}

			// for exon regions only
			if (type == 3) {
				regions_in->exon_info[i] = calloc(regions_in->size_r[i], sizeof(char*));

				for (j=0; j<regions_in->size_r[i]; j++) {
					regions_in->exon_info[i][j] = NULL;
					//regions_in->exon_info[i][j] = calloc(500, sizeof(char));
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
				free(regions_in->gene[i][j]);
				if (type < 4) {
					free(regions_in->Synonymous[i][j]);
					free(regions_in->prev_genes[i][j]);
				}
			}

			if (regions_in->gene[i]) free(regions_in->gene[i]);

			if (type < 4) {
				if (regions_in->Synonymous[i]) free(regions_in->Synonymous[i]);
				if (regions_in->prev_genes[i]) free (regions_in->prev_genes[i]);
			}

			// for exon regions only
			if (type == 3) {
				for (j=0; j<regions_in->size_r[i]; j++) {
					free(regions_in->exon_info[i][j]);
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
	}

	if (type == 2 || type == 3) {
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
		sprintf(sql, "SELECT start, end from Regions_With_No_Annotation WHERE chrom='%s' ORDER BY start", chrom_id);
	} else if (type == 2) {
		sprintf(sql, "SELECT start, end, gene_symbol, Synonymous, prev_gene_symbol from Intron_Regions WHERE chrom='%s' ORDER BY start", chrom_id);
	} else if (type == 3) {
		sprintf(sql, "SELECT start, end, gene_symbol, Synonymous, prev_gene_symbol, source_name, exon_id from Unique_Gene_Exons WHERE chrom='%s' ORDER BY start", chrom_id);
	} else {
		sprintf(sql, "SELECT start, end, annotation from All_Site_Reports WHERE chrom='%s' ORDER BY start", chrom_id);
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
		}

		if (type == 2 || type == 3) {
			regions_in->Synonymous[index][count] = dynamicStringAllocation(row[3], regions_in->Synonymous[index][count]);
			regions_in->prev_genes[index][count] = dynamicStringAllocation(row[4], regions_in->prev_genes[index][count]);

			if (type == 3) {
				uint16_t exon_id = (uint16_t) atoi(row[6]);
				char *tmp = calloc(strlen(row[5]) + CHAR_BIT *sizeof(uint16_t) + 1, sizeof(char));
				sprintf(tmp, "%s_exon_%"PRIu16"", row[5], exon_id);
				regions_in->exon_info[index][count] = dynamicStringAllocation(tmp, regions_in->exon_info[index][count]);
				if (tmp) free(tmp);
			}
		}

		count++;
	}
    free(sql);
}

// As the array is sorted, I am going to use binary search for find the location
// here index is the chromosome index
int32_t checkInterGenicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, uint32_t low_index) {
	return binary_search(regions_in, start, end, index, low_index);
}

// here index is the chromosome index
int32_t binary_search(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, uint32_t low_index) {
	int32_t low = low_index;
	int32_t high = regions_in->size_r[index] - 1;
	int32_t middle = (low + high)/2;

	while (low <= high) {
		if (regions_in->starts[index][middle] <= start && end <= regions_in->ends[index][middle]) {
			return middle;
		} else {
			if (regions_in->starts[index][middle] < start) {
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
	int32_t found = binary_search(regions_in, start, end, index, low_index);

	if (found != -1) {
		uint16_t orig_str_len = strlen(*info_in_and_out);
		uint16_t str_len_needed = strlen(regions_in->gene[index][found]) + strlen(regions_in->Synonymous[index][found]) + strlen(regions_in->prev_genes[index][found]) + 10;
		
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
	int32_t found = binary_search(regions_in, start, end, index, low_index);

    if (found != -1) {
        uint16_t orig_str_len = strlen(*info_in_and_out);
        uint16_t str_len_needed = strlen(regions_in->gene[index][found]) + strlen(regions_in->Synonymous[index][found]) + strlen(regions_in->prev_genes[index][found]) + strlen(regions_in->exon_info[index][found]) + 25;

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
			if (strstr(regions_in->exon_info[index][found], "NR_") || strstr(regions_in->exon_info[index][found], "NM_")) {
				sprintf(*info_in_and_out, "%s\t%s\t%s\t%s\t.\t.\t.", regions_in->gene[index][found], regions_in->prev_genes[index][found], regions_in->Synonymous[index][found], regions_in->exon_info[index][found]);
			} else if (strstr(regions_in->exon_info[index][found], "CCDS")) {
				sprintf(*info_in_and_out, "%s\t%s\t%s\t.\t%s\t.\t.", regions_in->gene[index][found], regions_in->prev_genes[index][found], regions_in->Synonymous[index][found], regions_in->exon_info[index][found]);
			} else if (strstr(regions_in->exon_info[index][found], "OTTHUM")) {
				sprintf(*info_in_and_out, "%s\t%s\t%s\t.\t.\t%s\t.", regions_in->gene[index][found], regions_in->prev_genes[index][found], regions_in->Synonymous[index][found], regions_in->exon_info[index][found]);
			} else {
				sprintf(*info_in_and_out, "%s\t%s\t%s\t.\t.\t.\t%s", regions_in->gene[index][found], regions_in->prev_genes[index][found], regions_in->Synonymous[index][found], regions_in->exon_info[index][found]);
			}
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

void fromStringToIntArray(char *str_in, uint32_t *array_in) {
    uint16_t counter = 0;
    char *pch;
    pch = strtok(str_in, ",|\"");
    while (pch != NULL) {
        if (strcmp(pch, "\"") != 0) {
            array_in[counter] = (uint32_t) atol(pch);
            //printf("value is %s and counter is %d\n", pch, counter);
            counter++;
        }
        pch = strtok(NULL, ",|");
    }
}

void processExonArrays(uint16_t exon_count, uint32_t *exon_starts, uint32_t *exon_ends, char *gene_name, uint32_t pos, khash_t(str) *ret_hash) {
    char string_to_add[100]="";
    uint16_t i = 0;
	int ret;
	khiter_t k_iter;

	// Create an instance of Temp_Coverage_Array and initialize it, 
	// Here we have to use Temp_Coverage_Array instead of uint32_t because it seems that I could have only one khash_t(str) type declared
	//
    //Temp_Coverage_Array *temp_cov_array = calloc(1, sizeof(Temp_Coverage_Array));
    //temp_cov_array->size = 0;

    for (i=0; i<exon_count; i++) {
        if (exon_starts[i] <= pos && pos <= exon_ends[i]) {
            //printf("Found with exon start at %d and end at %d with iteration of %d\n", exonStarts[i], exonEnds[i], i);
            sprintf(string_to_add, "%s_exon_%d", gene_name, i);

			// create the key
			k_iter = kh_put(str, ret_hash, string_to_add, &ret);
			if (ret)
				kh_key(ret_hash, k_iter) = strdup(string_to_add);
            //kh_value(coverage_hash, k_iter) = temp_cov_array;
        }
		strcpy(string_to_add, "");
    }
}

char* combinedEachAnnotation(khash_t(str) *hash_in) {
    khiter_t k_iter;
    int flag = 0, str_total_size=200;
	char *ret_string;
	ret_string = calloc(str_total_size, sizeof(char));
    strcpy(ret_string, ".");

    for (k_iter = kh_begin(hash_in); k_iter != kh_end(hash_in); ++k_iter) {
        if (kh_exist(hash_in, k_iter) && strlen(kh_key(hash_in, k_iter)) > 0) {
            if (flag == 0) {
                strcpy(ret_string, kh_key(hash_in, k_iter));
                flag++;
            } else {
				// check to see if the ret_string can hold the s_length
            	if (strlen(ret_string) >= (str_total_size - 30)) {
                	// need to dynamically allocate the memory
	                char *tmp=NULL;
					str_total_size *= 2;
					//printf("Dynamic allocation memory with size %d and string is %s\n", str_total_size, ret_string);
    	            tmp = realloc(ret_string, str_total_size);
        	        if (!tmp) {
						free(ret_string);
						ret_string = NULL;
                    	fprintf(stderr, "String realloc() failed! Exiting...");
	                    //exit(1);
    	            }
            	    ret_string = tmp;
				}
                strcat(ret_string, ",");
                strcat(ret_string, kh_key(hash_in, k_iter));
            }
        }
    }

	return ret_string;
}

// here hash_in refers to the refseq_hash, or ccds_hash or vega_hash or miRNA_hash 
void processingMySQL(MYSQL *con, char *sql, uint32_t pos_start, uint32_t pos_end, char *gene, khash_t(str) *prev_gene, khash_t(str) *Synonymous, khash_t(str) *hash_in) {
    if (mysql_query(con,sql))
        finish_with_error(con);

    MYSQL_RES *result = mysql_store_result(con);
    if (result == NULL)
        finish_with_error(con);

    //int num_fields = mysql_num_fields(result);
    //printf("Total number of fields is %d\n", num_fields);

    MYSQL_ROW row;
    int absent;
    uint16_t exon_count=0;
    khiter_t Synonymous_iter, prev_gene_iter;

    while ((row = mysql_fetch_row(result))) {
        // here I need to locate which exon it is part of and process them accordingly
        exon_count = (uint16_t) atoi(row[3]);
        uint32_t exon_starts[exon_count], exon_ends[exon_count];
        fromStringToIntArray(row[4], exon_starts);
        fromStringToIntArray(row[5], exon_ends);

        processExonArrays(exon_count, exon_starts, exon_ends, row[0], pos_start, hash_in);
		if (pos_start != pos_end)
	        processExonArrays(exon_count, exon_starts, exon_ends, row[0], pos_end, hash_in);

		//id_iter_hash = kh_put(m32, id_list, atoi(row[1]), &absent);

        // for gene row[6], Synonymous row[7], and prev_gene row[8]
        if (row[6] && strlen(row[6]) > 0 && strcmp(row[6], "NULL") != 0) {
            if (strlen(gene) == 0 || strcmp(gene, ".") == 0) {
                strcpy(gene, row[6]);
            } else {
                if (strcmp(gene, row[6]) != 0 && strcmp(gene, ".") != 0) {
                    Synonymous_iter = kh_put(str, Synonymous, row[6], &absent);
					if (absent)
						kh_key(Synonymous, Synonymous_iter) = strdup(row[6]);
                }
            }
        }

		/*
		if (row[7] && strlen(row[7]) > 0 && strcmp(row[7], "NULL") != 0) {
            Synonymous_iter = kh_put(str, Synonymous, row[7], &absent);
            if (absent)
				kh_key(Synonymous, Synonymous_iter) = strdup(row[7]);
		}*/

        if (row[8] && strlen(row[8]) > 0 && strcmp(row[8], "NULL") != 0) {
            prev_gene_iter = kh_put(str, prev_gene, row[8], &absent);
			if (absent)
				kh_key(prev_gene, prev_gene_iter) = strdup(row[8]);
        }
	}

	// Need to clean it here otherwise, valgrind will give possible memory leak error
	if (result) mysql_free_result(result);
}

char * produceGeneAnnotations(uint32_t start_in, uint32_t stop_in, char *chrom_id, MYSQL *con) {
    char *prev_sql  = " start, end, exon_count, exon_starts, exon_ends, gene_symbol, alias_gene_symbol, prev_gene_symbol ";
	char *mid_sql   = calloc(150, sizeof(char)); 
    sprintf(mid_sql,  " chrom='%s' AND ((start <= %"PRIu32" AND %"PRIu32" <= end) OR (start <= %"PRIu32" AND %"PRIu32" <= end))", chrom_id, start_in, start_in, stop_in, stop_in);

	// for debugging
	/*
	//if (strcmp(chrom_id, "1") == 0)
	//	return;

	//if (start_in < 70276350) {
	if (start_in < 53676957) {
		//printf("%s\n", sql);
		return;
	} else {
		printf("%s\n", mid_sql);
	}
	if (start_in == 1246709) {
		printf("%s\n", mid_sql);
	}*/

    char gene[100]="";
	khash_t(str) *refseq     = kh_init(str);     // hash_table using string as key
	khash_t(str) *ccds       = kh_init(str);     // hash_table using string as key
	khash_t(str) *vega       = kh_init(str);     // hash_table using string as key
	khash_t(str) *Synonymous = kh_init(str);     // hash_table using string as key
	khash_t(str) *prev_gene  = kh_init(str);     // hash_table using string as key
	khash_t(str) *miRNA      = kh_init(str);     // hash_table using string as key

	// for refseq
	char *sql = calloc(strlen(prev_sql) + strlen(mid_sql) + 75, sizeof(char));
    sprintf(sql, "SELECT DISTINCT refseq_name, %s FROM RefSeq_annotation WHERE %s ", prev_sql, mid_sql);
    //printf("%s\n", sql);
    processingMySQL(con, sql, start_in, stop_in, gene, prev_gene, Synonymous, refseq);
	memset(sql,0,strlen(sql));	 

    // for ccds
    sprintf(sql, "SELECT DISTINCT ccds_name, %s FROM CCDS_annotation WHERE %s ", prev_sql, mid_sql);
    //printf("%s\n", sql);
    processingMySQL(con, sql, start_in, stop_in, gene, prev_gene, Synonymous, ccds);
	memset(sql,0,strlen(sql));	

    // for vega
    sprintf(sql, "SELECT DISTINCT vega_name, %s FROM VEGA_annotation WHERE %s ", prev_sql, mid_sql);
    //printf("%s\n", sql);
    processingMySQL(con, sql, start_in, stop_in, gene, prev_gene, Synonymous, vega);
	memset(sql,0,strlen(sql));	 

    // now for miRNA
    sprintf(sql, "SELECT DISTINCT miRNA_name, %s FROM miRNA_annotation WHERE %s ", prev_sql, mid_sql);
    //printf("%s\n", sql);
    processingMySQL(con, sql, start_in, stop_in, gene, prev_gene, Synonymous, miRNA);
	memset(sql,0,strlen(sql));	 

    // now we need to combine everything together
	char *refseq_str = combinedEachAnnotation(refseq);
	char *ccds_str   = combinedEachAnnotation(ccds);
	char *vega_str   = combinedEachAnnotation(vega);
	char *miRNA_str  = combinedEachAnnotation(miRNA);
	char *prev_gene_str  = combinedEachAnnotation(prev_gene);
	char *Synonymous_str = combinedEachAnnotation(Synonymous);

	if (strlen(gene) == 0)
		strcpy(gene, ".");
	uint16_t annotation_size = strlen(gene) + strlen(prev_gene_str) + strlen(Synonymous_str) + strlen(refseq_str) + strlen(ccds_str) 
								+ strlen(vega_str) + strlen(miRNA_str);
	char *annotation = calloc(annotation_size + 50, sizeof(char));

	sprintf(annotation, "%s\t%s\t%s\t%s\t%s\t%s\t%s", gene, prev_gene_str, Synonymous_str, refseq_str, ccds_str, vega_str, miRNA_str);

	cleanKhashStr(refseq, 2);
	cleanKhashStr(ccds, 2);
	cleanKhashStr(vega, 2);
	cleanKhashStr(miRNA, 2);
	cleanKhashStr(Synonymous, 2);
	cleanKhashStr(prev_gene, 2);

	//printf("%s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, start_in, stop_in);
	//fflush(stdout);

	// properly empty/reset the strings I declared
	if (sql) free(sql);
	//free(prev_sql);
	if (mid_sql) free(mid_sql);

	memset(gene,0,sizeof(gene));
	if (refseq_str) free(refseq_str);
	if (prev_gene_str) free(prev_gene_str);
	if (Synonymous_str) free(Synonymous_str);
	if (ccds_str) free(ccds_str);
	if (vega_str) free(vega_str);
	if (miRNA_str) free(miRNA_str);
	//memset(,0,sizeof());
	//printf("return produceGeneAnnations %d\n", start_in);
	
	return annotation;
}
