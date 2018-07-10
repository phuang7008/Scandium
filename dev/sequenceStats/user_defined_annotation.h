/*
 * =====================================================================================
 *
 *       Filename:  user_defined_annotation.h
 *
 *    Description:  This is used to handle the annotation database provided by users
 *
 *        Version:  1.0
 *        Created:  02/12/2018 03:38:01 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */

#ifndef USER_DEFINED_ANNOTATION_H
#define USER_DEFINED_ANNOTATION_H

#include "terms.h"
#include "utils.h"
#include "reports.h"

/*
 * this function is used to pre-process the user-defined database before store anything to the exon_regions
 * @param user_inputs, that provides the necessary file name for the user defined database
 * @param udd_wrapper, a wrapper that contains general information such as the number of lines, number of chromosomes etc.
 */
void getUserDefinedDatabaseInfo(User_Input *user_inputs, User_Defined_Database_Wrapper *udd_wrapper, khash_t(khStrInt) *cds_lengths, khash_t(khStrInt) *cds_counts, khash_t(khStrInt) *user_defined_targets);

/*
 * if user specifies his/her own database, we will use this database only
 * @param user_inputs, that provides the necessary file name for the user defined database
 * @param exon_regions, will store the information for all user-defined exons
 * @param udd_wrapper, a User_Defined_Database_Wrapper that contains the cds count information for each chromsome
 * @param raw_user_defined_database, store the raw information from user-defined database file
 */
void processUserDefinedDatabase(User_Input *user_inputs, Regions_Skip_MySQL *exon_regions, User_Defined_Database_Wrapper *udd_wrapper, Raw_User_Defined_Database * raw_user_defined_database, khash_t(khStrInt) *cds_lengths, khash_t(khStrInt) *cds_counts);

/*
 * Initialize the coverage regions that user-defined (need to be analyzed)
 *
 */
void userDefinedGeneCoverageInit(khash_t(khStrLCG) *user_cds_gene_hash, char *chrom_id, Raw_User_Defined_Database *raw_user_defined_database, khash_t(khStrStrArray) *gene_transcripts);

void recordUserDefinedTargets(khash_t(khStrInt) *user_defined_targets, Bed_Info *user_defined_bed_info);

void writeCoverageForUserDefinedDB(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Regions_Skip_MySQL *exon_regions, Regions_Skip_MySQL *intron_regions);

void cleanUserDefinedDatabase(User_Defined_Database_Wrapper *udd_wrapper);

void cleanRawUserDefinedDatabase(Raw_User_Defined_Database * raw_user_defined_database);

#endif //USER_DEFINED_ANNOTATION_H
