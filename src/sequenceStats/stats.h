/*
 * =====================================================================================
 *
 *       Filename:  stats.h
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

#ifndef STATS_H
#define STATS_H

//#include <my_global.h>
//#include <mysql.h>
//#include <stdbool.h>
#include "htslib/sam.h"
#include "terms.h"
#include "annotation.h"

/**
 * This function is used to process one chunk of aligned reads from bam file and process it using current thread
 * @param user_inputs: variable that contains all the user input info
 * @param rec: the individual alignment record to be processed
 * @param tmp_stats_info: variable used to store all the statistical information regarding bases and reads
 * @param header: bam/cram file header information that has chromosome id and length info 
 * @param chrom_tracking: the variable to store all the base coverage count information
 * @param chrom_index: current chromosome id to be handled
 * @param target_buffer_status: it contains the target and buffer info for capture sequencing
 *
 */
void processCurrentRecord(User_Input *user_inputs, bam1_t *rec, bam_hdr_t *header, Stats_Info *tmp_stats_info, Chromosome_Tracking *chrom_tracking, uint32_t chrom_index, Target_Buffer_Status *target_buffer_status);

/**
 * This function is used to process individual aligned read and put the results into a hash table where the key is the chromosome location
 * while the value is the count
 * @param user_inputs: variable that contains all the user input info
 * @param tmp_stats_info: variable used to store all the statistical information regarding bases and reads
 * @param rec: the individual alignment record to be processed
 * @param chrom_tracking: the variable to store all the base coverage count information
 * @param chrom_index: current chromosome index to be used
 * @param target_buffer_status: it contains the target and buffer info for capture sequencing
 */
void processRecord(User_Input *user_inputs, Stats_Info *tmp_stats_info, bam1_t *rec, Chromosome_Tracking *chrom_tracking, uint32_t chrom_index, Target_Buffer_Status *target_buffer_status);

//void setTargetBufferFlags(Target_Buffer_Status *target_buffer_status, uint8_t *on_target, uint8_t *on_buffer, uint32_t chrom_idx, uint32_t pos_idx);
void setTargetBufferFlags(Target_Buffer_Status *target_buffer_status, bool *on_target, bool *on_buffer, uint32_t chrom_idx, uint32_t pos_idx);

bool getOverlapInfo(User_Input *user_inputs, Stats_Info *stats_info, bam1_t *rec, uint32_t *m_pos_r_end);

#endif // STATS_H
