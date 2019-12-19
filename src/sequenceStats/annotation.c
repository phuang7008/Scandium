/*
 * =====================================================================================
 *
 *       Filename:  annotation.c
 *
 *    Description:  the detailed implementation of the annotation header file
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
    exit(EXIT_FAILURE);
}

// make the MySQL connection and create the name of databases to be used by the program
//
void databaseSetup(Databases *dbs, User_Input *user_inputs) {
	dbs->con = mysql_init(NULL);
	if (dbs->con == NULL)
		finish_with_error(dbs->con);

	if (mysql_library_init(0, NULL, NULL)) {
		fprintf(stderr, "Could not initialize MySQL Library!\n");
		exit(EXIT_FAILURE);
	}

	if (user_inputs->user_name == NULL) {
		fprintf(stderr, "Please enter your MySQL DB login user name and password!\n");
		exit(EXIT_FAILURE);
	}

	if (user_inputs->passwd == NULL) {
		fprintf(stderr, "Please enter your MySQL DB login user name and password!\n");
		exit(EXIT_FAILURE);
	}

	if (mysql_real_connect(dbs->con, "sug-esxa-db1", user_inputs->user_name, user_inputs->passwd, "GeneAnnotations", 0, NULL, 0) == NULL) {
		finish_with_error(dbs->con);
	}

	// once the MYSQL connection established, we need to find out which database to use here
	//
	dbs->db_introns = calloc(50, sizeof(char));
	dbs->db_coords  = calloc(50, sizeof(char));
	dbs->db_annotation = calloc(50, sizeof(char));
	dbs->db_hgmd = calloc(50, sizeof(char));

	if (strcasecmp(user_inputs->database_version, "hg38") == 0) {
		strcpy(dbs->db_coords, "Gene_RefSeq_CDS_Coords38");
		strcpy(dbs->db_introns, "Intron_Regions38");
		strcpy(dbs->db_annotation, "Gene_Annotations38");
		strcpy(dbs->db_hgmd, "HGMD38");

	} else {
		strcpy(dbs->db_coords, "Gene_RefSeq_CDS_Coords37");
		strcpy(dbs->db_introns, "Intron_Regions37");
		strcpy(dbs->db_annotation, "Gene_Annotations37");
		strcpy(dbs->db_hgmd, "HGMD37");

	}

	fprintf(stderr, "The database for gene/transcript/exon percentage calculation is: %s\n", dbs->db_coords);
	fprintf(stderr, "The database for gene annotation (partitioned) is %s\n", dbs->db_annotation);
	fprintf(stderr, "The Intronic database is: %s\n", dbs->db_introns);
	if (HGMD_PROVIDED)
		fprintf(stderr, "The HGMD database is: %s\n", dbs->db_hgmd);
	fprintf(stderr, "\n");
}

void databaseCleanUp(Databases *dbs) {
	if (dbs->db_coords != NULL)
		free(dbs->db_coords);

	if (dbs->db_annotation != NULL)
		free(dbs->db_annotation);

	if (dbs->db_introns != NULL)
		free(dbs->db_introns);

	if (dbs->db_hgmd != NULL)
		free(dbs->db_hgmd);
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions
//
uint32_t fetchTotalCount(uint8_t type, Databases *dbs, char *chrom_id) {
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
		if (row[0] != 0)    // row[0] is char*, a pointer, need to use pointer comparison
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
void regionsSkipMySQLInit(Databases *dbs, Regions_Skip_MySQL *regions_in, uint8_t type) {
	// here we need to find out how many chromosomes we are dealing with
	char *sql = calloc(250, sizeof(char));
		
	// Here 'NOT LIKE' will exclude non-primary chromosomes such as 19_GL949752v1_alt
	// For hg38, however, we need to include them in. Thus, I remove them.
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

	// clean up the tmp variables
	//
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
            exit(EXIT_FAILURE);
        }

        regions_in->Synonymous = calloc(regions_in->chrom_list_size, sizeof(char**));
	    if (!regions_in->Synonymous) {
		    fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Synonymous");
			exit(EXIT_FAILURE);
        }

	    regions_in->prev_genes = calloc(regions_in->chrom_list_size, sizeof(char**));
		if (!regions_in->prev_genes) {
			fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Prev_genes");
	        exit(EXIT_FAILURE);
        }

		// for exon regions only
		if (type == 3) {
			regions_in->exon_info = calloc(regions_in->chrom_list_size, sizeof(char**));
	        if (!regions_in->exon_info) {
		        fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Gene");
			    exit(EXIT_FAILURE);
			}
        }
	}

	uint32_t i, j;
	for (i=0; i<regions_in->chrom_list_size; i++) {
		regions_in->size_r[i] = fetchTotalCount(type, dbs, regions_in->chromosome_ids[i]);
		regions_in->starts[i] = calloc(regions_in->size_r[i], sizeof(uint32_t));

		if (!regions_in->starts[i]) {
			fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Starts");
			exit(EXIT_FAILURE);
		}

		regions_in->ends[i] = calloc(regions_in->size_r[i], sizeof(uint32_t));
		if (!regions_in->ends) {
			fprintf(stderr, "Memory allocation for %s failed\n", "Regions With Less Annotation Ends");
			exit(EXIT_FAILURE);
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

		// now populate these regions with gene annotation
		//
		populateStaticRegionsForOneChromOnly(regions_in, dbs, regions_in->chromosome_ids[i], i, type); 
	}
}

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions, type 4 for all_site_reports 
void regionsSkipMySQLDestroy(Regions_Skip_MySQL *regions_in, uint8_t type) {
	uint32_t i, j;

	//printf("current chromsome id is %"PRIu32"\n", regions_in->chrom_list_size);
	for (i=0; i<regions_in->chrom_list_size; i++) {
		// for exon and intronic regions
		//
		if (type > 1) {
			for (j=0; j<regions_in->size_r[i]; j++) {
				if (regions_in->gene[i][j] != NULL) free(regions_in->gene[i][j]);
				if (regions_in->Synonymous[i][j] != NULL) free(regions_in->Synonymous[i][j]);
				if (regions_in->prev_genes[i][j] != NULL) free(regions_in->prev_genes[i][j]);
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

// type 1 for inter-genic regions, type 2 for intronic-regions, type 3 for exon regions
// This is for the general gene annotation, not for the gene/transcript/CDS coverage percentage calculation
//
void populateStaticRegionsForOneChromOnly(Regions_Skip_MySQL *regions_in, Databases *dbs, char *chrom_id, uint32_t chrom_idx, uint8_t type) {
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
        regions_in->starts[chrom_idx][count] = (uint32_t) strtol(row[0], NULL, 10);
        regions_in->ends[chrom_idx][count]   = (uint32_t) strtol(row[1], NULL, 10);

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

// As the array is sorted, I am going to use binary search for find the location
//
int32_t checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, char **info_in_and_out, uint32_t low_search_index) {
	int32_t found = binarySearch(regions_in, start, end, chrom_idx, low_search_index);

	if (found != -1) {
		uint32_t orig_str_len = strlen(*info_in_and_out);
		uint32_t str_len_needed = strlen(regions_in->gene[chrom_idx][found]) + strlen(regions_in->Synonymous[chrom_idx][found]) + strlen(regions_in->prev_genes[chrom_idx][found]) + 50;
		
		if (str_len_needed > orig_str_len) {
			char *tmp = realloc(*info_in_and_out, str_len_needed * sizeof(char));
			if (!tmp) {
				fprintf(stderr, "Memory re-allocation for string failed in checkIntronicRegion\n");
				exit(EXIT_FAILURE);
			}

			*info_in_and_out = tmp;
		} 

		sprintf(*info_in_and_out, "%s\t%s\t%s", regions_in->gene[chrom_idx][found], regions_in->prev_genes[chrom_idx][found], regions_in->Synonymous[chrom_idx][found]);
	}

	return found;
}

// this function will try to find all exons that overlap with a low coverage region and annotate them accordingly
//
int32_t checkExonRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, char **info_in_and_out, uint32_t low_search_index) {
	int32_t found = binarySearch(regions_in, start, end, chrom_idx, low_search_index);

	//if (start == 81003) {
	//	printf("stop\n");
	//}

    if (found != -1) {
		// define a variable that hold the annotation array with initial size of 5
		// the size of annotation array might differ, so we will have to increase the size dynamically later
		// Note: the size information are stored inside the annotation wrapper
		//
		Annotation_Wrapper *annotation_wrapper = calloc(1, sizeof(Annotation_Wrapper));
		annotation_wrapper->annotations = calloc(5, sizeof(Annotation));
		annotation_wrapper->allocated_size = 5;
		annotation_wrapper->real_size = 0;

		// Now copy the current exon info into the Annotation Array (annotations)
		//
		copyAnnotationDetails(annotation_wrapper, regions_in, chrom_idx, found);

		// Next, need to loop through left hand side to see if any left exons are part of low coverage regions
		// need to use SIGNED integer as it might go to negative and cause undefined behavior if using UNSIGNED!
		//
		int32_t i;
		int32_t j=0;
		for (i=found-1; i>=0; i--) {
			// the following is added to prevent annotation getting too long and out of control!
			//
			if (j>200) break;

			// check if the left exon overlaps with the low coverage region (lcr)
			//         lcr start   =======================================   lcr end
			//                                     left exon start  ------------   left exon end
			//					left exon start  ------ left exon end
			// left exon start  ------ left exon end
			//
			if ( (start <= regions_in->starts[chrom_idx][i] && regions_in->starts[chrom_idx][i] <= end) || 
				   (start <= regions_in->ends[chrom_idx][i] && regions_in->ends[chrom_idx][i] <= end) ) {
				copyAnnotationDetails(annotation_wrapper, regions_in, chrom_idx, i);
			} else {
				// because the records are partitioned, so we should break it here!
				//
				break;
			}
			j++;
		}

		// Here we also need to loop through the right hand side as well!
		// The checking criteria is the same as the schema used the above
		//
		j=0;
		for (i=found+1; (uint32_t) i<regions_in->size_r[chrom_idx]; i++) {
			// the following is added to prevent annotation getting too long and out of control!
			//
			if (j>200) break;

			if ( (start <= regions_in->starts[chrom_idx][i] && regions_in->starts[chrom_idx][i] <= end) ||
				   (start <= regions_in->ends[chrom_idx][i] && regions_in->ends[chrom_idx][i] <= end) ) {
				copyAnnotationDetails(annotation_wrapper, regions_in, chrom_idx, i);
			} else {
				break;
			}
			j++;
		}

		// Now need to walk through the Annoration array and put everything together
		//
		combineAllExonAnnotations(annotation_wrapper, info_in_and_out);

		// clean up
		//
		annotationWrapperDestroy(annotation_wrapper);
    }

    return found;
}

// The real copy of annotations for one exon only
//
void copyAnnotationDetails(Annotation_Wrapper *annotation_wrapper, Regions_Skip_MySQL *regions_in, uint32_t chrom_idx, uint32_t found_loc) {
	// need to dynamically expand the annotation array
	//
	if (annotation_wrapper->allocated_size <= annotation_wrapper->real_size) {
		annotation_wrapper->allocated_size = annotation_wrapper->allocated_size * 2;
		Annotation *tmp = realloc(annotation_wrapper->annotations, annotation_wrapper->allocated_size*sizeof(Annotation));
		if (!tmp) {
			fprintf(stderr, "Memory re-allocation for the struct Annotation failed in checkExonRegion\n");
			exit(EXIT_FAILURE);
		}
		annotation_wrapper->annotations = tmp;
	}

	// first store the gene info from the exon info found
	//
	annotation_wrapper->annotations[annotation_wrapper->real_size].gene = calloc(strlen(regions_in->gene[chrom_idx][found_loc])+1, sizeof(char));
	strcpy(annotation_wrapper->annotations[annotation_wrapper->real_size].gene, regions_in->gene[chrom_idx][found_loc]);

	// store Synonymous info
	//
	annotation_wrapper->annotations[annotation_wrapper->real_size].Synonymous = calloc(strlen(regions_in->Synonymous[chrom_idx][found_loc])+1, sizeof(char));
	strcpy(annotation_wrapper->annotations[annotation_wrapper->real_size].Synonymous, regions_in->Synonymous[chrom_idx][found_loc]);

	// store prev_gene info (the previous gene information will not be available)
	//
	//annotation_wrapper[annotation_wrapper->real_size].prev_genes = calloc(strlen(regions_in->prev_genes[chrom_idx][found_loc]), sizeof(char*));
	//strcpy(annotation_wrapper[annotation_wrapper->real_size].prev_genes, regions_in->prev_genes[chrom_idx][found_loc]);

	// store exon_info
	//
	annotation_wrapper->annotations[annotation_wrapper->real_size].exon_info = calloc(strlen(regions_in->exon_info[chrom_idx][found_loc])+1, sizeof(char));
	strcpy(annotation_wrapper->annotations[annotation_wrapper->real_size].exon_info, regions_in->exon_info[chrom_idx][found_loc]);

	annotation_wrapper->real_size++;
}

// To combine all the annotations together
//
void combineAllExonAnnotations(Annotation_Wrapper *annotation_wrapper, char **info_in_and_out) {
	uint32_t orig_str_len = strlen(*info_in_and_out);
	uint32_t str_len_needed = 0;

	if (annotation_wrapper->real_size == 1) {
		// just output the results
		//
		//str_len_needed = strlen(annotations[0].gene) + strlen(annotations[0].prev_genes) + strlen(annotations[0].Synonymous) + strlen(annotations[0].exon_info) + 10;	// 10 is added for extra spacing
		str_len_needed = strlen(annotation_wrapper->annotations[0].gene) + strlen(annotation_wrapper->annotations[0].Synonymous) + strlen(annotation_wrapper->annotations[0].exon_info) + 20;	// 20 is added for extra spacing

		if (str_len_needed > orig_str_len) {
			*info_in_and_out = realloc(*info_in_and_out, str_len_needed * sizeof(char));
			if (*info_in_and_out == NULL) {
				fprintf(stderr, "Memory re-allocation for string failed in checkExonRegion\n");
				exit(EXIT_FAILURE);
			}
		}

		sprintf(*info_in_and_out, "%s\t%s\t%s\t%s", annotation_wrapper->annotations[0].gene, annotation_wrapper->annotations[0].Synonymous, ".", annotation_wrapper->annotations[0].exon_info);
	} else {
		uint16_t i, j;
		khiter_t iter;		// khash iterator
		int absent = 0;

		// First, Handle gene list
		//
		khash_t(khStrInt) *gene_symbols = kh_init(khStrInt);

		// Note, here I used khash table to eliminate duplicates
		//
		for (i=0; i<annotation_wrapper->real_size; i++) {
			iter = kh_put(khStrInt, gene_symbols, annotation_wrapper->annotations[i].gene, &absent);
			if (absent) {
				kh_key(gene_symbols, iter)   = strdup(annotation_wrapper->annotations[i].gene);
				kh_value(gene_symbols, iter) = 0;
			}
			kh_value(gene_symbols, iter)++;
		}

		size_t total_gene_symbol_size = 0;
		uint16_t num_of_genes=0;
		for (iter = kh_begin(gene_symbols); iter != kh_end(gene_symbols); ++iter) {
			if (kh_exist(gene_symbols, iter)) {
				num_of_genes++;
				total_gene_symbol_size += strlen(kh_key(gene_symbols, iter)) + 5;
			}
		}
		//printf("The number of genes are %"PRIu16"\n", num_of_genes);

		// save everything into an array for later comparison with Synonymous
		//
		char **gene_arrays = NULL;
		if (num_of_genes > 0)
			gene_arrays = calloc(num_of_genes, sizeof(char*));

		i=0;
		for (iter = kh_begin(gene_symbols); iter != kh_end(gene_symbols); ++iter) {
			if (kh_exist(gene_symbols, iter)) {
				gene_arrays[i] = calloc(strlen(kh_key(gene_symbols, iter))+1, sizeof(char));
				strcpy(gene_arrays[i], kh_key(gene_symbols, iter));
				i++;
			}
		}

		// get the gene list into a string, need to leave some rooms for the separator ;
		//
		char *gene_list = calloc(total_gene_symbol_size, sizeof(char));

		// need to get rid of '.' if the size is >= 2
		// Note, here I need to do strcpy before I do the strcat
		//
		if (num_of_genes == 0) {
			strcpy(gene_list, ".");
		} else if (num_of_genes == 1) {
			strcpy(gene_list, gene_arrays[0]);
		} else {
			uint8_t flag = 0;
			for(i=0; i<num_of_genes; i++) {
				if (strcmp(gene_arrays[i], ".") != 0) {
					if (flag>0) {
						strcat(gene_list, ";");
						strcat(gene_list, gene_arrays[i]);
					} else {
						strcpy(gene_list, gene_arrays[i]);
					}
					flag++;
				}
			}
		}

		// Next, the Synonymous
		// First, need to find out how many Synonymous
		// Since the strtok_r() will destroy the original string, we need to make a copy of it
		// or store them somewhere for later usage
		//
		khash_t(khStrInt) *Synonymous = kh_init(khStrInt);
		char *savePtr;

		for (i=0; i<annotation_wrapper->real_size; i++) {
			char *copy_info = calloc(strlen(annotation_wrapper->annotations[i].Synonymous)+2, sizeof(char));
			strcpy(copy_info, annotation_wrapper->annotations[i].Synonymous);
			savePtr = copy_info;

			char *tokPtr;

			while ((tokPtr = strtok_r(savePtr, ";", &savePtr))) {
				iter = kh_put(khStrInt, Synonymous, tokPtr, &absent);
				if (absent) {
					kh_key(Synonymous, iter)   = strdup(tokPtr);
					kh_value(Synonymous, iter) = 0;
				}
				kh_value(Synonymous, iter)++;
			}

			if (copy_info != NULL) free(copy_info);
		}

		size_t total_Synonymous_size = 0;
		uint16_t num_of_Synonymous = 0;
		for (iter = kh_begin(Synonymous); iter != kh_end(Synonymous); ++iter) {
			if (kh_exist(Synonymous, iter)) {
				total_Synonymous_size += strlen(kh_key(Synonymous, iter)) + 5;
				num_of_Synonymous++;
			}
		}

		// save all Synonymous into a string
		// but we also need to get rid of "." if the size >= 2
		//
		char * Synonymous_list = calloc(total_Synonymous_size, sizeof(char));

		uint8_t flag = 0;
		if (num_of_Synonymous == 0) {
			strcpy(Synonymous_list, ".");
		} else {
			for (iter = kh_begin(Synonymous); iter != kh_end(Synonymous); ++iter) {
				if (kh_exist(Synonymous, iter)) {
					// check if it is also part of gene_symbol list
					//
					bool skip = false;
					for (j=0; j<num_of_genes; j++) {
						if ( (strcmp(gene_arrays[j], kh_key(Synonymous, iter)) == 0) ||
							 (strcmp(kh_key(Synonymous, iter), ".") == 0) )	{
							skip = true;
							break;
						}
					}

					if (!skip) {
						if (flag>0) {
							strcat(Synonymous_list, ";");
							strcat(Synonymous_list, kh_key(Synonymous, iter));
						} else {
							// this if else is needed because you need to do strcpy before strcat
							//
							strcpy(Synonymous_list, kh_key(Synonymous, iter));
						}
						flag++;
					}
				}
			}
		}

		if (flag == 0)
			strcpy(Synonymous_list, ".");

		if (gene_arrays != NULL) {
			for (i=0; i<num_of_genes; i++) {
				if (gene_arrays[i] != NULL) {
					free(gene_arrays[i]);
				}
			}

			free(gene_arrays);
		}

		// Now handle exon_annotations
		// Need to handle the SNP or Pseudo-gene annotations as well!
		// there are 6 sources of annotations in five categories: 
		// 0:RefSeq, 1:CCDS, 2:(VEGA for hg37, while Gencode for hg38), 3:miRNA and 4:everything else (including SNP and Pseudo-genes)
		// But will only have 4 place holders as VEGA is not available in hg38 and Gencode is not availabe for hg37
		// we need to handle them separately
		//
		khash_t(khStrInt) **exon_annotations = calloc(5, sizeof(khash_t(khStrInt) *));
		size_t size_ann=5;
		for (i=0; i<size_ann; i++) {
			exon_annotations[i] = kh_init(khStrInt);
		}

		for (i=0; i<annotation_wrapper->real_size; i++) {
			// make a copy and pass the copied one to the function
			//
			char *copy_info = calloc(strlen(annotation_wrapper->annotations[i].exon_info)+2, sizeof(char));
			strcpy(copy_info, annotation_wrapper->annotations[i].exon_info);
			savePtr = copy_info;

			char *tmp_token;
			uint32_t j=0;

			while ((tmp_token = strtok_r(savePtr, "\t", &savePtr))) {
				splitStringToKhash(tmp_token, exon_annotations, j);
				j++;
			}

			if (copy_info != NULL) free(copy_info);
		}

		size_t total_string_length=0;		// the string length could get quite long ...
		for (i=0; i<size_ann; i++) {
			for (iter = kh_begin(exon_annotations[i]); iter != kh_end(exon_annotations[i]); ++iter) {
				if (kh_exist(exon_annotations[i], iter)) {
					total_string_length += strlen(kh_key(exon_annotations[i], iter)) + 5;
				}
			}
		}

		// save all the exon information into a string
		//
		char *exon_list = calloc(total_string_length+50, sizeof(char));
		exon_list[0] = '\0';	// first character is now the null terminator, we can use strcat directly
		for (i=0; i<size_ann; i++) {
			// flag used to remove un-neccessary ";" in the output
			//
			uint8_t format_flag=0;

			for (iter = kh_begin(exon_annotations[i]); iter != kh_end(exon_annotations[i]); ++iter) {
				if (kh_exist(exon_annotations[i], iter)) {

					if (strcmp(kh_key(exon_annotations[i], iter), ".") == 0)
						continue;

					if (format_flag > 0)
						strcat(exon_list, ";");

					// we can use strcat directly as we set exon_list[0] = '\0' already
					//
					strcat(exon_list, kh_key(exon_annotations[i], iter));
					if (format_flag == 0) format_flag++;
				}
			}
			if (format_flag == 0)
				strcat(exon_list, ".");

			if (i<(size_ann-1))
				strcat(exon_list, "\t");
		}

		// final string to be output
		//
		str_len_needed = strlen(gene_list) + strlen(Synonymous_list) + strlen(exon_list) + 20;   // 20 is added for extra spacing
		if (str_len_needed > orig_str_len) {
			*info_in_and_out = realloc(*info_in_and_out, str_len_needed * sizeof(char));
			if (*info_in_and_out == NULL) {
				fprintf(stderr, "Memory re-allocation for string failed in checkExonRegion\n");
				exit(EXIT_FAILURE);
			}
		}

		sprintf(*info_in_and_out, "%s\t%s\t%s\t%s", gene_list, Synonymous_list, ".", exon_list);

		//fprintf(stderr, "Multiple genes output: %s\n", gene_list);

		// clean-up
		//
		cleanKhashStrInt(gene_symbols);
		cleanKhashStrInt(Synonymous);
		for (i=0; i<size_ann; i++)
			cleanKhashStrInt(exon_annotations[i]);

		if (exon_annotations != NULL) {
			free(exon_annotations);
			exon_annotations=NULL;
		}

		//if (genes != NULL) free(genes);
		//if (Synonymous != NULL) free(Synonymous);
		//if (exon_annotations != NULL) free(exon_annotations);

		if (exon_list != NULL) free(exon_list);
		if (gene_list != NULL) free(gene_list);
		if (Synonymous_list != NULL) free(Synonymous_list);
	}
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
// Two types of Databases are used, one uses MySQL, while another one uses user-defined database
// Here is the schema:
//
// low_cov_gene_hash: key every 1000 base such as 123,000 124,000 125,000 etc. 
//					 -> total_size
//					 -> capacity
//					 -> gene_coverage [an array of the Gene_Coverage structure type]
//									 -> [0] gene_symbol
//									 -> [0] transcript_name
//									 -> [0] cds_target_start
//									 -> [0] cds_target_end
//									 -> [0] exon_id
//									 -> [0] exon_count
//									 -> [0] ...
//									 -> [0] low_cov_regions [an array of the stringArray structure type]
//														   -> size
//														   -> capacity
//														   -> theArray [an array of string: char*]
//																	  ->[0] 235856732-235856734
//																	  ->[1] 235856749-235856786
//																	  ......	
//									 -> [1] gene_symbol
//									 -> [1] transcript_name
//									 ......
//
void genePercentageCoverageInit(khash_t(khStrLCG) *low_cov_gene_hash, char *chrom_id, Databases *dbs, khash_t(khStrStrArray) *gene_transcripts) {
	// mysql to obtain total number of distinct transcript_name for this specific chromosome
	//
	char *sql = calloc(350, sizeof(char));

	sprintf(sql, "SELECT cds_target_start, cds_target_end, exon_id, exon_count, cds_start, cds_end, cds_length, gene_symbol, transcript_name FROM %s WHERE chrom='%s' ORDER BY cds_start, cds_end", dbs->db_coords, chrom_id);

	if (mysql_query(dbs->con,sql))
		finish_with_error(dbs->con);

	dbs->mysql_results = mysql_store_result(dbs->con);
	if (dbs->mysql_results == NULL)
		finish_with_error(dbs->con);

	MYSQL_ROW row;

	// the following hash is used to check if we have already seen current transcript
	// if so, we don't need to add it to the gene_transcripts hash table
	//
	khash_t(khStrInt) *seen_transcript = kh_init(khStrInt);

	while ((row = mysql_fetch_row(dbs->mysql_results))) {
		// need to setup the gene-transcripts hash table here as we only want to handle this one for each row
		// row[7] is gene_symbol while row[8] is transcript name (ie, transcript_name)
		//
		addToGeneTranscriptKhashTable(row[7], row[8], gene_transcripts, seen_transcript);

		// now build the hash table every 1000 bases interval on current chromosome if any cds is present
		// need to find which key to use
		//
		uint32_t t_start = (uint32_t) strtol(row[0], NULL, 10);
		uint32_t t_end   = (uint32_t) strtol(row[1], NULL, 10);
		char *key_str = calloc(30, sizeof(char));
		uint32_t current_key = getHashKey(t_start);

		while (current_key < t_end) {
			sprintf(key_str, "%"PRIu32, current_key);

			// initialize the bucket with the current key_str
			//
			lowCoverageGeneHashBucketKeyInit(low_cov_gene_hash, key_str);

			// update Gene_Coverage variable ==> low_cov_gene_hash->gene_coverage
			// use a tmp Gene_Coverage for easy handling
			//
			khiter_t iter = kh_get(khStrLCG, low_cov_gene_hash, key_str);
			if (iter == kh_end(low_cov_gene_hash)) {  
				// iter will be equal to kh_end if key not present
				// However, this shouldn't happen as the previous lowCoverageGeneHashBucketKeyInit()
				// should generate the key bucket. So if it happens, exit!
				//
				fprintf(stderr, "Memory re-allocation failed at genePercentageCoverageInit()!\n");
				exit(EXIT_FAILURE);
			}

			uint32_t idx = kh_value(low_cov_gene_hash, iter)->total_size;
			Gene_Coverage *gc = &kh_value(low_cov_gene_hash, iter)->gene_coverage[idx];

			gc->gene_symbol = calloc(strlen(row[7])+1, sizeof(char));
			strcpy(gc->gene_symbol, row[7]);

			gc->transcript_name = calloc(strlen(row[8])+1, sizeof(char));
			strcpy(gc->transcript_name, row[8]);

			gc->cds_start = (uint32_t) strtol(row[4], NULL, 10);
			gc->cds_end   = (uint32_t) strtol(row[5], NULL, 10);
			gc->cds_length= (uint32_t) strtol(row[6], NULL, 10);

			gc->cds_target_start = (uint32_t) strtol(row[0], NULL, 10);
			gc->cds_target_end   = (uint32_t) strtol(row[1], NULL, 10);
			gc->exon_id    = (int16_t)  strtol(row[2], NULL, 10);
			gc->exon_count = (uint16_t) strtol(row[3], NULL, 10);

			gc->low_cov_regions = NULL;
			gc->targeted = false;

			kh_value(low_cov_gene_hash, iter)->total_size++;

			current_key += 1000;
		}
		free(key_str);
	}

	// tidy up the memory
	//
	if (dbs->mysql_results) { mysql_free_result(dbs->mysql_results); dbs->mysql_results = NULL; }
	if (sql) { free(sql); sql = NULL; }

	cleanKhashStrInt(seen_transcript);
}

void genePercentageCoverageDestroy(khash_t(khStrLCG) *low_cov_gene_hash) {
	uint32_t i;
	khiter_t iter;
	for (iter = kh_begin(low_cov_gene_hash); iter != kh_end(low_cov_gene_hash); ++iter) {
		if (kh_exist(low_cov_gene_hash, iter)) {
			for (i=0; i<kh_value(low_cov_gene_hash, iter)->total_size; i++) {
				if (kh_value(low_cov_gene_hash, iter)->gene_coverage[i].gene_symbol) {
					free(kh_value(low_cov_gene_hash, iter)->gene_coverage[i].gene_symbol);
					kh_value(low_cov_gene_hash, iter)->gene_coverage[i].gene_symbol = NULL;
				}

				if (kh_value(low_cov_gene_hash, iter)->gene_coverage[i].transcript_name) {
					free(kh_value(low_cov_gene_hash, iter)->gene_coverage[i].transcript_name);
					kh_value(low_cov_gene_hash, iter)->gene_coverage[i].transcript_name = NULL;
				}

				if (kh_value(low_cov_gene_hash, iter)->gene_coverage[i].low_cov_regions != NULL) {
					stringArrayDestroy(kh_value(low_cov_gene_hash, iter)->gene_coverage[i].low_cov_regions);
					free(kh_value(low_cov_gene_hash, iter)->gene_coverage[i].low_cov_regions);
					kh_value(low_cov_gene_hash, iter)->gene_coverage[i].low_cov_regions = NULL;
				}
			}
	
			if (kh_value(low_cov_gene_hash, iter)->gene_coverage) {
				free(kh_value(low_cov_gene_hash, iter)->gene_coverage);
				kh_value(low_cov_gene_hash, iter)->gene_coverage = NULL;
			}

			// now delete the key is the key exists
			//
			if (kh_key(low_cov_gene_hash, iter) != NULL)
				free((char *) kh_key(low_cov_gene_hash, iter));

			// now delete the value
			//
			if (kh_value(low_cov_gene_hash, iter) != NULL)
				free(kh_value(low_cov_gene_hash, iter));
		}
	}

	if (low_cov_gene_hash) {
		kh_destroy(khStrLCG, low_cov_gene_hash);
		low_cov_gene_hash = NULL;
	}
}

// do the real intersect between refseq cds bed regions with target bed regions
// it will re-adjust the start and end position based on the intersection.
// Here is one case:
//	bucket X:  60,377,000 ---------------------------------------------------------60,378,000
//	low coverage regions:    ******										********
//	CDS regions:			 ======										========
//							  cds3										  cds2
//
// Therefore, one bucket could contains more than one cds regions
//
void intersectTargetsAndRefSeqCDS(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, khash_t(khStrLCG) *low_cov_gene_hash) {
	// find out the index that is used to track current chromosome id
	//
	int32_t chrom_idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	if (chrom_idx == -1) return;

	if (low_cov_gene_hash == NULL) return;

	if (TARGET_FILE_PROVIDED) {

		uint32_t i,j;
		for(i = 0; i < target_info->size; i++) {
			if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
				continue;

			uint32_t start_t  = target_info->coords[i].start;
			uint32_t end_t    = target_info->coords[i].end;

			// now do the intersection
			//
			uint32_t current_key = getHashKey(start_t);
			while (current_key < end_t) {
				char key_str[20];
				sprintf(key_str, "%"PRIu32, current_key);

				// Need to find out if low_cov_gene_hash has such key entry
				//
				khiter_t iter = kh_get(khStrLCG, low_cov_gene_hash, key_str);

				if (iter == kh_end(low_cov_gene_hash)) {  // iter will be equal to kh_end if key not present
					current_key += 1000;
					continue;
				} else {
					// walk through the current key array
					//
					for (j=0; j<kh_val(low_cov_gene_hash,iter)->total_size; j++) {
						// now check the cds_target_coordinates
						//
						Gene_Coverage *gc = &kh_val(low_cov_gene_hash,iter)->gene_coverage[j];
						uint32_t start_ct = gc->cds_target_start;
						uint32_t end_ct   = gc->cds_target_end;

						// Here is the break out condition:
						//		start_t ============= end_t
						//			                start_ct ----------- end_ct ==>>>
						//
						if (end_t < start_ct) continue;

						//							start_t ================= end_t
						// <<<== start_ct ----------- end_ct
						//
						if (end_ct < start_t) continue;

						//    start_t ========== end_t
						// start_ct -------- end_ct
						//         start_ct ---------- end_ct
						// if they intersect, update the low_cov_gene_hash variable
						//
						uint32_t start_tmp = (start_t >= start_ct) ? start_t : start_ct;
						uint32_t end_tmp = (end_t <= end_ct) ? end_t : end_ct;

						// because the end position is not included, we need to skip this if start_tmp == end_tmp
						// Now record intersection regions
						//
						if (end_tmp > start_tmp) {
							gc->cds_target_start = start_tmp;
							gc->cds_target_end   = end_tmp;
							gc->targeted = true;
						}
					}
				}

				current_key += 1000;
			}
		}
	}
}

// Two types of Databases are used, one uses MySQL database, while another one uses user-defined database
// This is used to re-organize everything based on transcript name
//
// The original hash table (low_cov_gene_hash):
// key: every 1000
//			6,037,000 ---------------------------------------------- 6,038,000 ---------------------------------------- 6,039,000
//	low_cov_regions		****	**			*************************			*		****		*	
//			CDSs		===========			   ==========		=========================				============
//	Gene_Coverage		   cds1					  cds2					cds3								cds4
//
// For the new hash table (transcript_hash):
// key: transcript_name
//			NM_000069 => cds1
//			NM_000069 => cds2
//			NM_000069 => cds3
//			...
//			NM_000016 => cds0
//			NM_000016 => cds1
//			NM_000016 => cds2
//			...
//
void transcriptPercentageCoverageInit(khash_t(khStrLCG) *transcript_hash, khash_t(khStrLCG) *low_cov_gene_hash) {
	uint32_t i;
	khiter_t iter_ts, iter_lcg;

	// because at this time, all the low coverage regions will be counted already
	// iterate through the low coverage gene hash and re-group them based on transcript name
	//
	for (iter_lcg = kh_begin(low_cov_gene_hash); iter_lcg != kh_end(low_cov_gene_hash); ++iter_lcg) {
		if (kh_exist(low_cov_gene_hash, iter_lcg)) {
			// loop through the array based on current key such as (57,000 or 326,000 etc.)
			//
			for (i=0; i<kh_value(low_cov_gene_hash, iter_lcg)->total_size; i++) {
				// Initialize the bucket with the current key
				//
				lowCoverageGeneHashBucketKeyInit(transcript_hash, kh_value(low_cov_gene_hash, iter_lcg)->gene_coverage[i].transcript_name);

				// need to find out if the key for the current transcript name exists
				//
				iter_ts = kh_get(khStrLCG, transcript_hash, kh_value(low_cov_gene_hash, iter_lcg)->gene_coverage[i].transcript_name);
				if (iter_ts == kh_end(transcript_hash)) {
					// the key doesn't exists!
					// This shouldn't happen. So we exit right away!
					//
					fprintf(stderr, "key doesn't exist at the transcriptPercentageCoverageInit()!\n");
					exit(EXIT_FAILURE);
				}

				// now add the cds into the array
				//
				uint32_t idx = kh_value(transcript_hash, iter_ts)->total_size;
				copyGeneCoverageLowCovRegions(&kh_value(low_cov_gene_hash, iter_lcg)->gene_coverage[i], 
						&kh_value(transcript_hash, iter_ts)->gene_coverage[idx], false);
				kh_value(transcript_hash, iter_ts)->total_size++;
			}
		}
	}
}

// Here we are not going to go through the database, instead, we are going to use the Low_Coverage_Genes stored on the heap!!!
//
void produceGenePercentageCoverageInfo(uint32_t start_in, uint32_t stop_in, khash_t(khStrLCG) *low_cov_gene_hash) {
	khiter_t iter;
	uint32_t current_key = getHashKey(start_in);

	while (current_key < stop_in) {
		char key_str[20];
		sprintf(key_str, "%"PRIu32, current_key);

		// locate the iterator from the low_cov_gene_hash table through key
		//
		iter = kh_get(khStrLCG, low_cov_gene_hash, key_str);
		if (iter == kh_end(low_cov_gene_hash)) {  
			// iter will be equal to kh_end if key not present
			//
			current_key += 1000;
			continue;
		} else {
			// FOUND!
			// Even we found the bucket, it doesn't mean there will be an intersection with low coverage region!
			// Need to check each region one by one and find if they truly intersect with the low coverage region!
			// Now loop through the array pointed by this key
			//
			uint32_t i;
			for (i=0; i<kh_value(low_cov_gene_hash, iter)->total_size; i++) {
				Gene_Coverage *gc = &kh_value(low_cov_gene_hash, iter)->gene_coverage[i];

				uint32_t cds_target_start = gc->cds_target_start;
				uint32_t cds_target_end   = gc->cds_target_end;

				// check if they intersect
				//
				if (start_in <= cds_target_start && cds_target_end <= stop_in) {
					//               cds_target_start ============== cds_target_end
					//  start_in -------------------------------------------------------- stop_in
					//
					// kh_value(low_cov_gene_hash, iter)->gene_coverage[i].num_of_low_cov_bases += cds_target_end - cds_target_start;
					processExonArrays(gc, cds_target_start, cds_target_end);
				} else if (cds_target_start <= start_in && stop_in <= cds_target_end) {
					//  cds_target_start =================== cds_target_end
					//          start_in ------------- orig_end
					//
					processExonArrays(gc, start_in, stop_in);
				} else if (cds_target_start <= stop_in && stop_in <= cds_target_end) {
					//             cds_target_start ============== cds_target_end
					//  start_in ----------------------------- stop_in
					//
					processExonArrays(gc, cds_target_start, stop_in);
				} else if (cds_target_start <= start_in && start_in <= cds_target_end) {
					// cds_target_start ================= cds_target_end
					//                     start_in ----------------------- stop_in
					//
					processExonArrays(gc, start_in, cds_target_end);
				} else {
					// not interested!
					continue;
				}
			}
		}

		current_key += 1000;
	}
}

// This will calculate the number of overlap count to find out the true low coverage region
// here start is the low coverage region start
// end is the low coverage region end
//
void processExonArrays(Gene_Coverage *gc, uint32_t start, uint32_t end) {
	if (start == end)
		return;

	// let's first find the string size to use
	//
	char tmp_string[50];
	sprintf(tmp_string, "%"PRIu32, start);
	uint32_t report_string_size = strlen(tmp_string) + 1;

	sprintf(tmp_string, "%"PRIu32, start);
	report_string_size += strlen(tmp_string) + 1;

	int i;

	// now we need to allocate memory space for low_cov_regions string array for exon percentage report
	//
	if (gc->low_cov_regions == NULL) {
		gc->low_cov_regions = calloc(1, sizeof(stringArray));
		gc->low_cov_regions->capacity = 2;
		gc->low_cov_regions->theArray = calloc(gc->low_cov_regions->capacity, sizeof(char*));
		gc->low_cov_regions->size = 0;

		for (i=0; i<gc->low_cov_regions->capacity; i++)
			gc->low_cov_regions->theArray[i]=NULL;
	} 
	
	if (gc->low_cov_regions->capacity == gc->low_cov_regions->size) {
		gc->low_cov_regions->capacity *= 2;
		gc->low_cov_regions->theArray = 
			realloc(gc->low_cov_regions->theArray, gc->low_cov_regions->capacity * sizeof(char*));
		if (gc->low_cov_regions->theArray == NULL) {
			fprintf(stderr, "Memory re-allocation failed at processExonArrays\n");
			exit(EXIT_FAILURE);
		}

		for (i=gc->low_cov_regions->size; i<gc->low_cov_regions->capacity; i++)
			gc->low_cov_regions->theArray[i]=NULL;
	}

	// Now append the low coverage region to the array
	//
	uint16_t idx = gc->low_cov_regions->size;
	gc->low_cov_regions->theArray[idx] = calloc(report_string_size, sizeof(char));
	sprintf(gc->low_cov_regions->theArray[idx], "%"PRIu32"-%"PRIu32, start, end);
	gc->low_cov_regions->size++;
}
