/*
 * =====================================================================================
 *
 *       Filename:  annotation.h
 *
 *    Description:  The header file for the gene annotation related functionalities
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

/*
 * a helper funciton that is used to clean-up the Databases variable
 * @param dbs: a Databases variable to be cleaned
 */
void databaseCleanUp(Databases *dbs);

/*
 * The values stored in the MySQL database are strings. We need to extract them and store them into an INT array
 * @param str_in: the string that contains all the starts OR ends
 * @param array_in: the integer array to store all the exon starts OR ends
 */
void fromStringToIntArray(char *str_in, uint32_t *array_in);

/*
 * this function is used to process the INT exon array and find the intercepted regions
 * @param gene_coverage: a Gene_Coverage variable that contains all the gene coverage information
 * @param start: the start positions for low coverage region that will be handled
 * @param end:   the end positions for low coverage region that will be handled
 */
void processExonArrays(Gene_Coverage *gene_coverage, uint32_t start, uint32_t end);

/**
 * the function will combine all the strings (stored as key) and formatted them for output
 * @param hash_in: khash_t(str) object
 * @return ret_string: the string to store the formatted results for output
 */
char* combinedEachAnnotation(khash_t(str) *hash_in);

/*
 * get total number of rows feteched from the MySQL database
 * @param dbs: a Databases variable that contains all the database connection information
 * @param chrom_id: only fetch those on current chromosome
 * @param user_inputs: it contains the database version used
 */
uint32_t fetchTotalCount(uint8_t type, Databases *dbs, char *chrom_id);

/*
 * if we query MySQL db for every region, it will be too slow.
 * Here we fetch everything in for a single query to speed things up
 * @param dbs: a Databases variable that contains all the database connection information
 * @param regions_in: used to store all the fetched regions
 * @param user_inputs: it contains the database version used
 * @param type: type of regions to fetch. 1 => intergenic, 2 => intron, 3 => exon
 */
void regionsSkipMySQLInit(Databases *dbs, Regions_Skip_MySQL *regions_in, uint8_t type);

/*
 * a helper function that is used to clean-up all fetched regions from MySQL db
 * @param regions_in: regions to be cleaned
 * @param type: type of regions to fetch. 1 => intergenic, 2 => intron, 3 => exon
 */
void regionsSkipMySQLDestroy(Regions_Skip_MySQL *regions_in, uint8_t type);

/*
 * this is the real function that is used to read each rows from MySQL query results
 * @param regions_in: used to store all the fetched regions
 * @param dbs: a Databases variable that contains all the database connection information
 * @param chrom_id: current chromosome id to be handled
 * @param index: the array index on the regions_in
 * @param type: type of regions to fetch. 1 => intergenic, 2 => intron, 3 => exon
 */
void populateStaticRegionsForOneChromOnly(Regions_Skip_MySQL *regions_in, Databases *dbs, char *chrom_id, uint32_t index, uint8_t type);

/*
 * this is used to check if the region (start -- end) is overlapped with intronic regions
 * @param regions_in: used to store all the fetched intronic regions
 * @param start: start position to be checked
 * @param end: end position to be checked
 * @param index: the array index for the regions_in
 * @param info_in_and_out: the output string to be stored
 * @param low_index: the lower bound to start the search
 * @return return 1 if it is found or -1 if it is not
 */
int32_t checkIntronicRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, char **info_in_and_out, uint32_t low_index);

/*
 * this is used to check if the region (start -- end) is overlapped with exonic regions
 * @param regions_in: used to store all the fetched exonic regions
 * @param start: start position to be checked 
 * @param end: end position to be checked
 * @param index: the array index for the regions_in
 * @param info_in_and_out: the output string to be stored 
 * @param low_index: the lower bound to start the search 
 * @return return 1 if it is found or -1 if it is not  
 */
int32_t checkExonRegion(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t index, char **info_in_and_out, uint32_t low_index);

/**
 * The function is used to record one exon info that overlaps with the low coverage region (regions_in)
 * @param annotations, an array of Annotation that stores all of the exon info that overlap with low coverage region
 * @param regions_in, the official annotation region being checked against
 * @param chrom_idx, the current chromosome index
 * @param found_loc, the found exon index location
 */
void copyAnnotationDetails(Annotation_Wrapper *annotation_wrapper, Regions_Skip_MySQL *regions_in, uint32_t chrom_idx, uint32_t found_loc);

/*
 * combines all the exon annotation for output
 * @param annotation_wrapper: a variable that is used to store all gene annotation information
 * @param info_in_and_out: the combined exon annotation for output
 */
void combineAllExonAnnotations(Annotation_Wrapper *annotation_wrapper, char **info_in_and_out);

/*
 * It uses binary search to see if the region (start -- end) overlapped with anything regions of interest
 * @param regions_in, the official annotation region being checked against
 * @param start: start position to be checked
 * @param end: end position to be checked
 * @param chrom_idx: the array index for the regions_in
 * @param low_search_index: the lower bound to start the search
 */
int32_t binarySearch(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t low_search_index);

/*
 * To check if a region (start -- end) is within the location_index before doing anything else
 * @param regions_in, the official annotation region being checked against
 * @param start: start position to be checked
 * @param end: end position to be checked
 * @param chrom_idx: the array index for the regions_in
 * @param low_search_index: the lower bound to start the search
 */
bool verifyIndex(Regions_Skip_MySQL *regions_in, uint32_t start, uint32_t end, uint32_t chrom_idx, uint32_t location_index);

/**
 * this is used to initialize the gene percentage coverage structure
 * @param low_cov_gene_hash, the variable to hold the official refseq genes' cds info also the low coverage regions
 * @param chrom_id, the chromosome id we are going to handle
 * @param dbs, the structure contain MySQL and Database information
 * @param user_inputs, contains all the original user's inputs
 * @param gene_transcripts: a hash variable that is to be initialized
 */
void genePercentageCoverageInit(khash_t(khStrLCG) *low_cov_gene_hash, char *chrom_id, Databases *dbs, khash_t(khStrStrArray) *gene_transcripts);

/**
 * This is used to clean up the Low_Coverage_Genes variables
 * @param low_cov_gene_hash: a khash_t variable using string as key and Low Coverage Gene as value
 */
void genePercentageCoverageDestroy(khash_t(khStrLCG) *low_cov_gene_hash);

/**
 * This is used to intersect the refseq cds bed regions with target bed regions
 * @param chrom_id: the chromosome id we are going to handle
 * @param target_info, the target bed regions that user provided
 * @param chrom_tracking: a variable to track chromosome processed
 * @param low_cov_gene_hash, contains all the low coverage gene information including low coverage regions
 */
void intersectTargetsAndRefSeqCDS(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, khash_t(khStrLCG) *low_cov_gene_hash);

/**
 * Calculate the gene/transcript coverage percentage
 * @param start_in, the start position for the intersection
 * @param stop_in, the end position for the intersection
 * @param chrom_id: the chromosome id we are going to handle
 * @param low_cov_gene_hash: a khash_t(khStrLCG) variable that is used to store all the low coverage information
 */
void produceGenePercentageCoverageInfo(uint32_t start_in, uint32_t stop_in, khash_t(khStrLCG) *low_cov_gene_hash);

/*
 * it is used to initialized the transcript_hash variable
 * @param transcript_hash: a khash_t(khStrLCG) variable to be initialized
 * @param low_cov_gene_hash: a low_cov_gene_hash variable that contains all the low coverage information
 */
void transcriptPercentageCoverageInit(khash_t(khStrLCG) *transcript_hash, khash_t(khStrLCG) *low_cov_gene_hash);

/**
 * record the intersected regions between the target bed regions and the official refseq CDS regions
 * @param low_cov_gene_hash, a Low_Coverage_Genes variable that is used to store all the low coverage regions 
 * @param refseq_found_index, the index on the refseq CDS region that intersects with the target bed region
 * @param target_counter, the counter index on the target region that intersects with the low_cov_gene_hash
 * @param start, the start position for the intersection
 * @param end, the end position for the intersection
 */
void recordIntersectedRegions(Low_Coverage_Genes *low_cov_gene_hash, int32_t refseq_found_index, uint32_t target_counter, uint32_t start, uint32_t end);

#endif // ANNOTATION_H
