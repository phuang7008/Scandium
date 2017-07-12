/*
 * =====================================================================================
 *
 *       Filename:  annotation.h
 *
 *    Description:  The header file for sequence stats analysis
 *
 *        Version:  1.0
 *        Created:  02/24/2017 03:47:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine
 *
 * =====================================================================================
 */

#ifndef ANNOTATION_H
#define ANNOTATION_H

//#include <my_global.h>
//#include <mysql.h>
//#include <stdbool.h>
#include "htslib/sam.h"
#include "terms.h"

/*
 * produce the error message for possible MySQL queries/executions
 * @param con: the MySQL connection object/handler
 */
void finish_with_error(MYSQL *con);

/*
 * The values stored in the MySQL database are strings. We need to extract them and store them into an INT array
 * @param str_in: the string that contains all the starts OR ends
 * @param array_in: the integer array to store all the exon starts OR ends
 */
void fromStringToIntArray(char *str_in, uint32_t *array_in);

/*
 * this function is used to process the INT exon array and find the intercepted regions
 * @param exon_count: the number of exons need to be known before calling this function
 * @param exon_starts: the start positions for every exon region
 * @param exon_ends:   the end positions for every exon region
 * @param gene_name: the name of the gene from MySQL query
 * @param pos: the position to be intercepted by exons
 * @param ret_val: the returned hash table with string as key
 */
void processExonArrays(uint16_t exon_count, Low_Coverage_Genes *low_cov_genes, uint32_t LCG_array_index, uint16_t refseq_index, uint32_t start, uint32_t end);

/**
 * the function will combine all the strings stored as key and formatted them for output
 * @param hash_in: khash_t(str) object
 * @return ret_string: the string to store the formatted results for output
 */
char* combinedEachAnnotation(khash_t(str) *hash_in);

uint32_t fetchTotalCount(uint8_t type, MYSQL *con, char *chrom_id);

void regionsSkipMySQLInit(MYSQL *con, Regions_Skip_MySQL *regions_in, uint8_t type);

void regionsSkipMySQLDestroy(Regions_Skip_MySQL *regions_in, uint8_t type);

void populateStaticRegionsForOneChromOnly(Regions_Skip_MySQL *regions_in, MYSQL *con, char *chrom_id, uint32_t index, uint8_t type);

int32_t checkInterGenicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, uint32_t low_index);

//char * checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index);
int32_t checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, char **info_in_and_out, uint32_t low_index);

int32_t checkExonRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, char **info_in_and_out, uint32_t low_index);

int32_t binary_search(Regions_Skip_MySQL *regions_in, uint32_t pos, uint32_t index, uint32_t low_index);

bool verifyIndex(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t location_index);


void fromStringToIntArray(char *str_in, uint32_t *array_in);

void processingMySQL(MYSQL *con, char *sql, uint32_t pos_start, uint32_t pos_end, Low_Coverage_Genes *low_cov_genes, uint32_t LCG_array_index, uint8_t type);

/**
 * this is used to initialize the gene percentage coverage structure
 * @param low_cov_genes, the variable to hold the low coverage genes info
 * @param size_in, the initial size of the Low_Coverage_Genes, the size would be expanded dynamically
 */
void genePercentageCoverageInit(Low_Coverage_Genes *low_cov_genes, char *chrom_id, MYSQL *con);

void genePercentageCoverageDestroy(Low_Coverage_Genes *low_cov_genes, char *chrom_id);

void produceGenePercentageCoverageInfo(uint32_t start_in, uint32_t stop_in, char *chrom_id, MYSQL *con, Low_Coverage_Genes *low_cov_genes);

/*
 * produce the gene annotation for the capture region
 * @param start_in: the start position for the capture region
 * @param stop_in: the stop position for the capture region
 * @param con: the MySQL connection object/handler
 */
char* produceGeneAnnotations(uint32_t start_in, uint32_t stop_in, char *chrom_id, MYSQL *con, omp_lock_t *query_lock);

void getGeneExonDetails(MYSQL *con, uint32_t start, uint32_t end, char *chrom_id, Gene_Exon_Array *gene_exon_array);

#endif // ANNOTATION_H
