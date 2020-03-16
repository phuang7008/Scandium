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
 *         Author:  Peiming (Peter) Huang, (phuang@bcm.edu)
 *        Company:  Baylor College of Medicine, Houston, TX 77030
 *
 * =====================================================================================
 */

#ifndef USER_DEFINED_ANNOTATION_H
#define USER_DEFINED_ANNOTATION_H

#include "terms.h"
#include "user_inputs.h"
#include "utils.h"
#include "reports.h"

/*
 * This function is used to validate the user input annotation file format and provide useful information if needed
 * @param user_inputs: the inputs that user provides
 */
void checkAnnotationFormat(User_Input *user_inputs);

/*
 * this function is used to pre-process the user-defined database before store anything to the exon_regions
 * @param user_inputs: it provides the necessary file name for the user defined database
 * @param udd_wrapper: a wrapper that contains general information such as the number of lines, number of chromosomes etc.
 * @param cds_lengths: a khash_t(khStrInt) variable that is used to store total cds length info for a transcript
 * @param cds_counts:  a khash_t(khStrInt) variable that is used to store number of used cds info for a transcript
 * @param user_defined_targets: a khash_t(khStrInt) variable that is used to store all targets from a target bed file
 *
 */
void getUserDefinedDatabaseInfo(User_Input *user_inputs, User_Defined_Database_Wrapper *udd_wrapper, khash_t(khStrInt) *cds_lengths, khash_t(khStrInt) *cds_counts, khash_t(khStrInt) *user_defined_targets);

/*
 * if user specifies his/her own database, we will use this database only
 * @param user_inputs: it provides the necessary file name for the user defined database
 * @param exon_regions: it stores the information for all user-defined exons
 * @param udd_wrapper: a User_Defined_Database_Wrapper variable that contains the cds information for each chromosome
 * @param raw_user_defined_database: stores the raw information from user-defined database file
 * @param cds_lengths: a khash_t(khStrInt) variable that is used to store total cds length info for a transcript
 * @param cds_counts:  a khash_t(khStrInt) variable that is used to store number of used cds info for a transcript
 */
void processUserDefinedDatabase(User_Input *user_inputs, Regions_Skip_MySQL *exon_regions, User_Defined_Database_Wrapper *udd_wrapper, Raw_User_Defined_Database * raw_user_defined_database, khash_t(khStrInt) *cds_lengths, khash_t(khStrInt) *cds_counts);

/*
 * Initialize the coverage regions to be analyzed
 * @param user_cds_gene_hash: a khash_t(khStrLCG) variable that is used to store all the genes' cds regions
 * @param chrom_id: current chromosome id to be handled
 * @param raw_user_defined_database: stores the raw information from user-defined database file
 * @param gene_transcripts: a khash_t(khStrStrArray) variable that stores all the transcripts for a gene
 */
void userDefinedGeneCoverageInit(khash_t(khStrLCG) *user_cds_gene_hash, char *chrom_id, Raw_User_Defined_Database *raw_user_defined_database, khash_t(khStrStrArray) *gene_transcripts);

/*
 * it will store all the target defined by end users and store them in user_defined_targets
 * @param user_defined_targets: a khash_t(khStrInt) variable that is used to store user defined target
 * @param user_defined_bed_info: the source of user defined target in bed format
 */
void recordUserDefinedTargets(khash_t(khStrInt) *user_defined_targets, Bed_Info *user_defined_bed_info);

void writeCoverageForUserDefinedDB(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Regions_Skip_MySQL *exon_regions, Regions_Skip_MySQL *intron_regions);

/*
 * A helper function used to clean-up user defined database
 * @param udd_wrapper: a User_Defined_Database_Wrapper variable that contains the cds information for each chromosome
 */
void cleanUserDefinedDatabase(User_Defined_Database_Wrapper *udd_wrapper);

/*
 * A helper function used to clean-up Raw user defined database
 * @param raw_user_defined_database: stores all raw user defined database as is
 */
void cleanRawUserDefinedDatabase(Raw_User_Defined_Database * raw_user_defined_database);

#endif //USER_DEFINED_ANNOTATION_H
