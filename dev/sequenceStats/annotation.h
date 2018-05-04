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
 * It will make the database connection and create db names to be used by the program
 * @param dbs, the structure that contains the information related to the MySQL 
 * @user_inputs, the information about all user's inputs
 */
void databaseSetup(Databases *dbs, User_Input *user_inputs);

void databaseCleanUp(Databases *dbs);

/*
 * The values stored in the MySQL database are strings. We need to extract them and store them into an INT array
 * @param str_in: the string that contains all the starts OR ends
 * @param array_in: the integer array to store all the exon starts OR ends
 */
void fromStringToIntArray(char *str_in, uint32_t *array_in);

/*
 * this function is used to process the INT exon array and find the intercepted regions
 * @param exon_count: the number of exons need to be known before calling this function
 * @param start: the start positions for low coverage region
 * @param end:   the end positions for low coverage region
 * @low_cov_gene_hash: the structure that holds all the refseq information of one chromosome
 * @refseq_exon_index: the found index that points to the refseq regions intercept with low coverage region
 */
//void processExonArrays(Low_Coverage_Genes *low_cov_gene_hash, uint32_t refseq_exon_index, uint32_t start, uint32_t end);
void processExonArrays(Gene_Coverage *gene_coverage, uint32_t start, uint32_t end);

/**
 * the function will combine all the strings stored as key and formatted them for output
 * @param hash_in: khash_t(str) object
 * @return ret_string: the string to store the formatted results for output
 */
char* combinedEachAnnotation(khash_t(str) *hash_in);

uint32_t fetchTotalCount(uint8_t type, Databases *dbs, char *chrom_id, User_Input *user_inputs);

void regionsSkipMySQLInit(Databases *dbs, Regions_Skip_MySQL *regions_in, User_Input *user_inputs, uint8_t type);

void regionsSkipMySQLDestroy(Regions_Skip_MySQL *regions_in, uint8_t type);

void populateStaticRegionsForOneChromOnly(Regions_Skip_MySQL *regions_in, Databases *dbs, char *chrom_id, uint32_t index, User_Input *user_inputs, uint8_t type);

//char * checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index);
int32_t checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, char **info_in_and_out, uint32_t low_index);

int32_t checkExonRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, char **info_in_and_out, uint32_t low_index);

/**
 * The function is used to record one exon info that overlaps with the low coverage region (regions_in)
 * @param annotations, an array of Annotation that stores all of the exon info that overlap with low coverage region
 * @param regions_in, the low coverage region being checked and annotated
 * @param chrom_idx, the current chromosome index
 * @param found_loc, the found exon index location
 */
void copyAnnotationDetails(Annotation_Wrapper *annotation_wrapper, Regions_Skip_MySQL *regions_in, uint32_t chrom_idx, uint32_t found_loc);

void combineAllExonAnnotations(Annotation_Wrapper *annotation_wrapper, char **info_in_and_out);

int32_t binarySearch(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t low_search_index);

bool verifyIndex(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t location_index);

/**
 * this is used to initialize the gene percentage coverage structure
 * @param low_cov_gene_hash, the variable to hold the official refseq genes' cds info
 * @param chrom_id, the chromosome id we are going to handle
 * @param dbs, the structure contain MySQL and Database information
 * @param user_inputs, contains all the original user's inputs
 */
void genePercentageCoverageInit(khash_t(khStrLCG) *low_cov_gene_hash, char *chrom_id, Databases *dbs, User_Input *user_inputs, khash_t(khStrStrArray) *gene_transcripts);

/**
 * This is used to clean up the Low_Coverage_Genes variables
 * @param Low_Coverage_Genes variable
 */
//void genePercentageCoverageDestroy(Low_Coverage_Genes *low_cov_gene_hash);
void genePercentageCoverageDestroy(khash_t(khStrLCG) *low_cov_gene_hash);

/**
 * This is used to intersect the refseq cds bed regions with target bed regions
 * @param chrom_id, the chromosome id we are going to handle
 * @param target_info, the target bed regions that user provided
 * @param low_cov_gene_hash, official refseq gene's cds information
 * @param low_cov_gene_hash, the target regions that intersect with low_cov_gene_hash
 */
void intersectTargetsAndRefSeqCDS(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, khash_t(khStrLCG) *low_cov_gene_hash);

void produceGenePercentageCoverageInfo(uint32_t start_in, uint32_t stop_in, char *chrom_id, khash_t(khStrLCG) *low_cov_gene_hash);

/*
 * produce the gene annotation for the capture region
 * @param start_in: the start position for the capture region
 * @param stop_in: the stop position for the capture region
 * @param dbs, the structure contain MySQL and Database information
 */
//char* produceGeneAnnotations(uint32_t start_in, uint32_t stop_in, char *chrom_id, Databases *dbs, omp_lock_t *query_lock);

void transcriptPercentageCoverageInit(char* chrom_id, khash_t(khStrLCG) *transcript_hash, khash_t(khStrLCG) *low_cov_gene_hash);

/**
 * record the intersected regions between the target bed regions and the official refseq CDS regions
 * @param low_cov_gene_hash, the official refseq CDS regions
 * @param low_cov_gene_hash, the target bed regions (at this moment, we presume all regions are of high quality)
 * @param refseq_found_index, the index on the refseq CDS region that intersects with the target bed region
 * @param target_counter, the counter index on the target region that intersects with the low_cov_gene_hash
 * @param start, the start position for the intersection
 * @param end, the end position for the intersection
 */
void recordIntersectedRegions(Low_Coverage_Genes *low_cov_gene_hash, int32_t refseq_found_index, uint32_t target_counter, uint32_t start, uint32_t end);

#endif // ANNOTATION_H
