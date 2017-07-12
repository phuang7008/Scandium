/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  using openmp for parallel processing to gather sequencing statistics
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include "htslib/sam.h"

#include "terms.h"
#include "stats.h"
#include "utils.h"
#include "targets.h"

int main(int argc, char *argv[]) {
	//fprintf(stderr, "Starting ... %ld\n", time(NULL));

	// get user input options and then processing it accordingly
	User_Input *user_inputs = userInputInit();
	processUserOptions(user_inputs, argc, argv);
	//printf("User option for wgs is: %d\n", user_inputs->wgs_coverage);

	// now need to setup a Stats_Info variable to track various statistical information
    Stats_Info *stats_info = statsInfoInit();

	// now for the bam/cram file
    samFile *sfd = sam_open(user_inputs->bam_file, "r");
    if (sfd == 0) {
        fprintf(stderr, "Can not open file %s\n", argv[1]);
        return -1;
    }

    // use sam_hdr_read to process both bam and cram header
    bam_hdr_t *header=NULL;
    if ((header = sam_hdr_read(sfd)) == 0) return -1;

    fetchTotalGenomeBases(header, stats_info);

	// setup a tracking variable to track chromosome working status 
    Chromosome_Tracking *chrom_tracking = chromosomeTrackingInit(header);

	uint32_t i;
	Bed_Info *target_bed_info=NULL, *Ns_bed_info=NULL;				//store only bed info chrom_id, start, stop
	Target_Buffer_Status *target_buffer_status = NULL;

    if (TARGET_FILE_PROVIDED || N_FILE_PROVIDED) {
		target_buffer_status = calloc(header->n_targets, sizeof(Target_Buffer_Status));
	    for(i=0; i<header->n_targets; i++) {
		    strcpy(target_buffer_status[i].chrom_id, header->target_name[i]);
			target_buffer_status[i].size = header->target_len[i] + 1;
			target_buffer_status[i].index = -1;
	        target_buffer_status[i].status_array = calloc(header->target_len[i] + 1, sizeof(uint8_t));
		}
	}

	if (N_FILE_PROVIDED) {      // the file that contains regions of Ns in the reference genome
        Ns_bed_info = calloc(1, sizeof(Bed_Info));
        processBedFiles(user_inputs->n_file, Ns_bed_info, stats_info, target_buffer_status, header, 2);
    }
    //outputForDebugging(Ns_bed_info);
    //fflush(stdout);
    //return 0;

	if (TARGET_FILE_PROVIDED) {
        target_bed_info = calloc(1, sizeof(Bed_Info));
        processBedFiles(user_inputs->target_file, target_bed_info, stats_info, target_buffer_status, header, 1);
    }
    //outputForDebugging(target_bed_info);
    //fflush(stdout);
    //return 0;

	// can't set to be static as openmp won't be able to handle it
	uint32_t total_chunk_of_reads = 1000000;
	if (user_inputs->wgs_coverage) {
		if (user_inputs->num_of_threads == 1) {
			total_chunk_of_reads = 45000000;
		} else if (user_inputs->num_of_threads == 2) {
            total_chunk_of_reads = 25000000;
		} else if (user_inputs->num_of_threads == 4) {
            total_chunk_of_reads = 12000000;
		} else if (user_inputs->num_of_threads == 6) {
			total_chunk_of_reads = 8500000;
		} else if (user_inputs->num_of_threads == 8) {
            total_chunk_of_reads = 5000000;
		} else { 
            total_chunk_of_reads = 4000000;
		}
	}

	// try to allocate the bam1_t array here for each thread, so they don't have to create and delete the array at each loop
	Read_Buffer *read_buff = calloc(user_inputs->num_of_threads, sizeof(Read_Buffer));
	for (i=0; i<user_inputs->num_of_threads; i++) {
		read_buff[i].chunk_of_reads = calloc(total_chunk_of_reads, sizeof(bam1_t));
		read_buff[i].size = total_chunk_of_reads;
	}

	fflush(stdout);

	// now let's do the parallelism
	while(chrom_tracking->more_to_read) {
#pragma omp parallel num_threads(user_inputs->num_of_threads)
		{
			int num_of_threads = omp_get_num_threads();	
			int thread_id = omp_get_thread_num();
			readBufferInit(&read_buff[thread_id]);
			uint32_t num_records = 0;
			//printf("Before File Reading: number of threads is %d and current thread id is %d\n", num_of_threads, thread_id);

#pragma omp critical 
			{
				// this part of the code need to run atomically, that is only one thread should allow access to read
				num_records = readBam(sfd, header, chrom_tracking, &read_buff[thread_id]);
				read_buff[thread_id].size = num_records;
			}
			//printf("First Critical position for thread %d\n", thread_id);

			if (num_records == 0) {
				printf("No more to read for thread %d for a total of %d threads!!!!!!!!!!!!\n", thread_id, num_of_threads);
				readBufferDestroy(&read_buff[thread_id]);
				if (chrom_tracking->more_to_read) chrom_tracking->more_to_read = false;
			}

			Coverage_Stats *cov_stats = coverageStatsInit();
			khash_t(str) *coverage_hash = kh_init(str);		// hash_table using string as key

			// can not use the else {} with the previous if {} block, otherwise the barrier will wait forever!
			if (num_records > 0) {
				printf("After file reading: Thread %d performed %d iterations of the loop.\n", thread_id, num_records);

				processBamChunk(user_inputs, cov_stats, coverage_hash, header, &read_buff[thread_id], target_buffer_status, thread_id);

				// release the allocated chunk of buffer for alignment reads after they have been processed!
			    printf("cleaning the read buffer hash for thread %d...\n\n", thread_id);
				readBufferDestroy(&read_buff[thread_id]);
			}
				//printf("Before Second Critical position for thread %d\n", thread_id);
#pragma omp critical
			if (num_records > 0) {
				//printf("Before writing to the coverage array for Thread %d\n", thread_id);
				combineThreadResults(chrom_tracking, coverage_hash, header);
				combineCoverageStats(stats_info, cov_stats);
				//printf("after writing to the coverage array for Thread %d\n", thread_id);

				cleanKhashStr(coverage_hash);
				free(cov_stats);

				// since all reads have been process, we need to set the status for all chromosomes (especially the last one to 2)
	        	if (!chrom_tracking->more_to_read) {
	            	for(i=0; i<header->n_targets; i++) {
        	        	if (chrom_tracking->chromosome_status[i] == 1)
                        	chrom_tracking->chromosome_status[i] = 2;
					}
            	}
        	}
			//printf("After Second Critical position for thread %d\n", thread_id);

// setup a barrier here and wait for every one of them to reach this point!
#pragma omp barrier 

			// check to see if any of the chromosomes has finished. If so, write the results out
#pragma omp critical
			{
				if (num_records > 0) {
					for (i=0; i<header->n_targets; i++) {
						if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
							printf("Chromosome id %s in thread id %d has finished processing, now dumping\n", chrom_tracking->chromosome_ids[i], thread_id);
							writeCoverage(chrom_tracking->chromosome_ids[i], Ns_bed_info, target_bed_info, chrom_tracking, user_inputs, stats_info);
						
							// now write the off targets into a wig file
							if (TARGET_FILE_PROVIDED && user_inputs->Write_WIG)
								produceOffTargetWigFile(chrom_tracking, chrom_tracking->chromosome_ids[i], target_bed_info, user_inputs, stats_info);

							// clean up the array allocated
							//printf("free the big array for current chromosome %s\n", chrom_tracking->chromosome_ids[i]);
							if (chrom_tracking->coverage[i]) {
								free(chrom_tracking->coverage[i]);
								chrom_tracking->coverage[i] = NULL;
							}
							chrom_tracking->chromosome_status[i] = 3;
						}
					}
				}
			}
			//if (chrom_tracking->chromosome_status[0] == 3) chrom_tracking->more_to_read = false;
			printf("End of while loop before flush for thread %d\n", thread_id);
			fflush(stdout);
		}
	}

	sam_close(sfd);
	bam_hdr_destroy(header);

	//printf("Outside the threads\n");

	// Now need to write the report
	writeReport(stats_info, user_inputs);
	//printf("After the write report\n");

	chromosomeTrackingDestroy(chrom_tracking);
	//printf("After the chromosome tracking destroy\n");

	if (target_bed_info)
		cleanBedInfo(target_bed_info);
	//printf("After the target destroy\n");

	if (Ns_bed_info)
		cleanBedInfo(Ns_bed_info);

	//printf("Before buffer hash\n");
	if (target_buffer_status) {
		for (i=0; i<25; i++)
			free(target_buffer_status[i].status_array);
		free(target_buffer_status);
	}

	//printf("Before stats_info destroy\n");

	if (stats_info)
		statsInfoDestroy(stats_info);
	//printf("after stats_info destroy\n");

	if (chrom_tracking)
		free(chrom_tracking);

	if (read_buff) {
		for (i=0; i<user_inputs->num_of_threads; i++) {
			if (read_buff[i].chunk_of_reads) 
				free(read_buff[i].chunk_of_reads);
    	}
		free(read_buff);
	}

	userInputDestroy(user_inputs);

	return 0;
}
