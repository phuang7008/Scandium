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
#include "annotation.h"
#include "reports.h"

int main(int argc, char *argv[]) {
	//fprintf(stderr, "Starting ... %ld\n", time(NULL));

	// get user input options and then processing it accordingly
	// and output all the options to the end user
	//
	User_Input *user_inputs = userInputInit();
	processUserOptions(user_inputs, argc, argv);
	outputUserInputOptions(user_inputs);

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

	// Target_Buffer_Status need to be set even though there is no target or Ns region specified
	// It is because one of the method processRecord() need chromosome lengths information to be set
	//
	Target_Buffer_Status *target_buffer_status = calloc(header->n_targets, sizeof(Target_Buffer_Status));
	for(i=0; i<header->n_targets; i++) {
		strcpy(target_buffer_status[i].chrom_id, removeChr(header->target_name[i]));
		target_buffer_status[i].size = header->target_len[i] + 1;
		target_buffer_status[i].index = -1;
		target_buffer_status[i].status_array = calloc(header->target_len[i] + 1, sizeof(uint8_t));
		target_buffer_status[i].num_of_chromosomes = header->n_targets;
	}

	if (N_FILE_PROVIDED) {      // the file that contains regions of Ns in the reference genome
        Ns_bed_info = calloc(1, sizeof(Bed_Info));
        processBedFiles(user_inputs, Ns_bed_info, stats_info, target_buffer_status, header, 2);
    }
    //outputForDebugging(Ns_bed_info);
    //fflush(stdout);
    //return 0;

	if (TARGET_FILE_PROVIDED) {
        target_bed_info = calloc(1, sizeof(Bed_Info));
        processBedFiles(user_inputs, target_bed_info, stats_info, target_buffer_status, header, 1);
    }
    //outputForDebugging(target_bed_info);
    //fflush(stdout);
    //return 0;

	// for MySQL connection
	//
	Databases *dbs = calloc(1, sizeof(Databases));

	// For quick processing, we can't query MySQL for every single gene/exon/cds regions
	// Instead, we will fetch them only once and store them in the regions defined below
	//
    Regions_Skip_MySQL *inter_genic_regions = calloc(1, sizeof(Regions_Skip_MySQL));
    Regions_Skip_MySQL *intronic_regions    = calloc(1, sizeof(Regions_Skip_MySQL));
    Regions_Skip_MySQL *exon_regions        = calloc(1, sizeof(Regions_Skip_MySQL));

	if (user_inputs->annotation_on) {
		// Initialize all annotation related variables
		//
		databaseSetup(dbs, user_inputs);

		regionsSkipMySQLInit(dbs, inter_genic_regions, user_inputs, 1);
		regionsSkipMySQLInit(dbs, intronic_regions, user_inputs, 2);
		regionsSkipMySQLInit(dbs, exon_regions, user_inputs, 3);
	}

	// can't set to be static as openmp won't be able to handle it
	// check the bam/cram file size first
	//
	uint64_t input_bam_file_size = check_file_size(user_inputs->bam_file);

	uint32_t total_chunk_of_reads = 500000;		// for small bam/cram file
    //total_chunk_of_reads = 3000000;			// Good for 3 threads with 16gb of memory
    //total_chunk_of_reads = 2200000;			// Good for 3 threads with 9gb  of memory
	if (input_bam_file_size > 5000000000)		// anything > 5Gb
		total_chunk_of_reads = 1400000;			// Good for 3 threads with 9gb  of memory

	// try to allocate the bam1_t array here for each thread, so they don't have to create and delete the array at each loop
	Read_Buffer *read_buff = calloc(user_inputs->num_of_threads, sizeof(Read_Buffer));
	for (i=0; i<user_inputs->num_of_threads; i++) {
		read_buff[i].chunk_of_reads = calloc(total_chunk_of_reads, sizeof(bam1_t));
		read_buff[i].size = total_chunk_of_reads;
	}

	fflush(stdout);

	// set random seed (should only be called ONCE)
	//
	if(user_inputs->percentage < 1.0) {
		srand((uint32_t)time(NULL));    // set random seed
	}

	// now let's do the parallelism
    while(chrom_tracking->more_to_read) {
#pragma omp parallel shared(read_buff, chrom_tracking) num_threads(user_inputs->num_of_threads)
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
          //readBufferDestroy(&read_buff[thread_id]);
          if (chrom_tracking->more_to_read) chrom_tracking->more_to_read = false;
        }

        Coverage_Stats *cov_stats = coverageStatsInit();
        khash_t(str) *coverage_hash = kh_init(str);		// hash_table using string as key

        // can not use the else {} with the previous if {} block, otherwise the barrier will wait forever!
        if (num_records > 0) {
          printf("After file reading: Thread %d performed %d iterations of the loop.\n", thread_id, num_records);

          processBamChunk(user_inputs, cov_stats, coverage_hash, header, &read_buff[thread_id], target_buffer_status, thread_id);
        }

        // release the allocated chunk of buffer for alignment reads after they have been processed!
        printf("cleaning the read buffer hash for thread %d...\n\n", thread_id);
        readBufferDestroy(&read_buff[thread_id]);

#pragma omp critical
        {
          if (num_records > 0) {
            combineThreadResults(chrom_tracking, coverage_hash, header);
            combineCoverageStats(stats_info, cov_stats);

            // since all reads have been process, we need to set the status for all chromosomes (especially the last one to 2)
            if (!chrom_tracking->more_to_read) {
              for(i=0; i<header->n_targets; i++) {
                if (chrom_tracking->chromosome_status[i] == 1)
                  chrom_tracking->chromosome_status[i] = 2;
              }
            }
          }
        }

        cleanKhashStr(coverage_hash, 1);
        free(cov_stats);
		cov_stats = NULL;	// to eleminate dangling pointer

// setup a barrier here and wait for every one of them to reach this point!
#pragma omp barrier 

#pragma omp single
        {
          if (num_records > 0) {
            for (i=0; i<header->n_targets; i++) {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                // check to see if any of the chromosomes has finished. If so, write the results out
                // for the whole genome, we need to use the file that contains regions of all Ns in the reference
                // As they will be not used, so we are going to set the count info in these regions to 0
                if (N_FILE_PROVIDED)
                  zeroAllNsRegions(chrom_tracking->chromosome_ids[i], Ns_bed_info, chrom_tracking, target_buffer_status);
              }
            }
          }
		}

#pragma omp barrier

		i = 0;
		while (i<header->n_targets) {
          if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
              printf("Chromosome id %s in thread id %d has finished processing, now dumping\n", chrom_tracking->chromosome_ids[i], thread_id);
            }

#pragma omp sections
          {
		    //printf("Inside sections with chromosome id %s for thread id %d\n", chrom_tracking->chromosome_ids[i], thread_id);
#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                printf("Thread %d is now producing coverage information for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);
                writeCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info, inter_genic_regions, intronic_regions, exon_regions);
						
                // now write the off target regions with high coverage into a wig file if the Write_WIG flag is set
				//
                if (TARGET_FILE_PROVIDED && user_inputs->Write_WIG)
                  produceOffTargetWigFile(chrom_tracking, chrom_tracking->chromosome_ids[i], target_bed_info, user_inputs, stats_info);
              }
            }

#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                //if (user_inputs->annotation_on) {
                  printf("Thread %d is now writing annotation for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);
                  writeAnnotations(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info, inter_genic_regions, intronic_regions, exon_regions);
                //}

				// if user specifies the range information (usually for graphing purpose), need to handle it here
				//
				printf("Thread %d is now writing coverage range info for graphing for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);
				coverageRangeInfoForGraphing(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info);
              }
            }

#pragma omp section
			{
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {

                if (user_inputs->annotation_on) {
                  printf("Thread %d is now working on exon percentage calculation for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);

				  // For calculating the percentage of gene bases with low coverge for capture only
				  // we need to allocate memories for both refseq_cds_genes and low_cov_genes
				  // and use the intersect regions between refseq_cds_genes for official annotation and low_cov_genes for targets
				  //
				  Low_Coverage_Genes  *refseq_cds_genes  = calloc(1, sizeof(Low_Coverage_Genes));
				  Low_Coverage_Genes  *low_cov_genes     = calloc(1, sizeof(Low_Coverage_Genes));
				  Transcript_Coverage *transcript_cov    = calloc(1, sizeof(Transcript_Coverage));

				  genePercentageCoverageInit(refseq_cds_genes, low_cov_genes, chrom_tracking->chromosome_ids[i], dbs, user_inputs);
				  if (user_inputs->annotation_type == 1)	// dynamic only
					intersectTargetsAndRefSeqCDS(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, refseq_cds_genes, low_cov_genes);
				  transcriptPercentageCoverageInit(chrom_tracking->chromosome_ids[i], transcript_cov, low_cov_genes, user_inputs, dbs);

				  calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info, dbs, low_cov_genes);
				  outputGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, user_inputs, low_cov_genes, transcript_cov, dbs);

				  // clean-up the memory space
				  genePercentageCoverageDestroy(low_cov_genes);
				  genePercentageCoverageDestroy(refseq_cds_genes);
				  transcriptPercentageCoverageDestroy(transcript_cov);
                }
              }
			}
          }

		  //printf("Waiting for other thread to completely here for thread %d\n", thread_id);
#pragma omp barrier 

#pragma omp single
          {
            if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
		      printf("\n");
              // clean up the array allocated
              if (chrom_tracking->coverage[i]) {
                free(chrom_tracking->coverage[i]);
                chrom_tracking->coverage[i] = NULL;
              }
              chrom_tracking->chromosome_status[i] = 3;
            }
		    i++;
          }
		}

#pragma omp barrier 
        printf("End of while loop before flush for thread %d\n", thread_id);
      }
	  printf("\n");
      fflush(stdout);
    }

	sam_close(sfd);

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

	if (target_buffer_status) {
		for (i=0; i<header->n_targets; i++)
			free(target_buffer_status[i].status_array);
		free(target_buffer_status);
	}

	if (stats_info)
		statsInfoDestroy(stats_info);
	//printf("after stats_info destroy\n");

	if (chrom_tracking)
		free(chrom_tracking);

	if (read_buff) {
		for (i=0; i<user_inputs->num_of_threads; i++) {
			if (read_buff[i].chunk_of_reads) {
				free(read_buff[i].chunk_of_reads);
			}
    	}
		free(read_buff);
	}

	// do some cleanup
	regionsSkipMySQLDestroy(inter_genic_regions, 1);
	regionsSkipMySQLDestroy(intronic_regions, 2);
	regionsSkipMySQLDestroy(exon_regions, 3);

	userInputDestroy(user_inputs);
	bam_hdr_destroy(header);

	// MYSQL clean-up
	mysql_close(dbs->con);
	mysql_library_end();

	return 0;
}
