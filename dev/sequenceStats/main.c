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
	printf("Processing bam file\n");
	fprintf(stderr, "Starting ... %ld\n", time(NULL));

	// get user input options and then processing it accordingly
	User_Input *user_inputs = userInputInit();
	processUserOptions(user_inputs, argc, argv);
	printf("User option for target file is: %s\n", user_inputs->target_file);

	// now need to setup a Stats_Info variable to track various statistical information
    Stats_Info *stats_info = calloc(1, sizeof(Stats_Info));
    statsInfoInit(stats_info);

	int i;
	Bed_Info *target_bed_info, *Ns_bed_info;				//store only bed info chrom_id, start, stop
	khash_t(str) *target_buffer_hash, *Ns_buffer_hash;		//store all the bed info in hash for quick lookup 'T':target, 'B':buffer, 'N':Ns

	if (TARGET_FILE_PROVIDED) {
		target_bed_info = calloc(1, sizeof(Bed_Info));
		target_buffer_hash = kh_init(str);
		processBedFiles(user_inputs->target_file, target_bed_info, target_buffer_hash, stats_info, 1);
	}
	//outputForDebugging(target_bed_info, target_buffer_hash);
	//fflush(stdout);
	//return 0;

	if (N_FILE_PROVIDED) {		// the file that contains regions of Ns in the reference genome
		Ns_bed_info = calloc(1, sizeof(Bed_Info));
		Ns_buffer_hash = kh_init(str);
		processBedFiles(user_inputs->n_file, Ns_bed_info, Ns_buffer_hash, stats_info, 2);
	}

	// now for the bam/cram file
	samFile *sfd = sam_open(user_inputs->bam_file, "r");
	if (sfd == 0) {
		fprintf(stderr, "Can not open file %s\n", argv[1]);
		return -1;
	}

	// use sam_hdr_read to process both bam and cram header
	bam_hdr_t *header=NULL;
	if ((header = sam_hdr_read(sfd)) == 0) return -1;
	
	int thread_id=0;
	bool more_to_read=true;
	uint32_t total_chunk_of_reads = 500000;	// can't set to be static as openmp won't be able to handle it

	// setup a tracking variable to track chromosome working status 
	Chromosome_Tracking *chrom_tracking = calloc(1, sizeof(Chromosome_Tracking));
	chrom_tracking->number_tracked = 0;

	printf ("user input number of threads is %d\n", user_inputs->num_of_threads);
	printf("before threading\n\n");
	fflush(stdout);

	// now let's do the parallelism
	while(more_to_read) {
#pragma omp parallel firstprivate(thread_id, total_chunk_of_reads) shared(more_to_read) num_threads(user_inputs->num_of_threads)
		{
			int num_of_threads = omp_get_num_threads();	
			thread_id = omp_get_thread_num();
			printf("Before File Reading: number of threads is %d and current thread id is %d\n",num_of_threads, thread_id);

			Read_Buffer *read_buff;
			read_buff = malloc(sizeof(Read_Buffer));
			read_buff->chunk_of_reads = calloc(total_chunk_of_reads, sizeof(bam1_t));
			read_buff->size = total_chunk_of_reads;
			readBufferInit(read_buff);
			uint32_t num_records = 0;

#pragma omp critical 
			{
				// this part of the code need to run atomically, that is only one thread should allow access to read
				num_records = readBam(sfd, header, more_to_read, read_buff);
			}

#pragma omp flush(more_to_read)
			// flush the shared variable more_to_read to ensure other threads will observe the updated value
			if (num_records == 0) {
				printf("No more to read!!!!!!!!!!!!\n");
				more_to_read = false;
			}
			
			thread_id = omp_get_thread_num();
			printf("After file reading: Thread %d performed %d iterations of the loop.\n", thread_id, num_records);

			khash_t(str) *coverage_hash = kh_init(str);		// hash_table using string as key
			processBamChunk(user_inputs, stats_info, coverage_hash, header, read_buff);

			// release the allocated chunk of buffer for alignment reads after they have been processed!
            printf("cleaning the read buffer hash...\n\n");
            readBufferDestroy(read_buff);
            free(read_buff);

#pragma omp critical
			{
				printf("Before writing to the coverage array for Thread %d\n", thread_id);
				combineThreadResults(chrom_tracking, coverage_hash, header);
				printf("after writing to the coverage array for Thread %d\n", thread_id);

				cleanKhashStr(coverage_hash);
			}

			// setup a barrier here and wait for every one of them to reach this point!
#pragma omp barrier

			// check to see if any of the chromosomes has finished. If so, write the results out
			//int j;
			for (i=0; i<chrom_tracking->number_tracked; i++) {
                if ( chrom_tracking->chromosome_status[i] == 2) {
#pragma omp critical
					{
						printf("Chromosome id %s in thread id %d has finished processing, now dumping\n", chrom_tracking->chromosome_ids[i], thread_id);
						writeCoverage(chrom_tracking->chromosome_ids[i], NULL, target_bed_info, chrom_tracking, user_inputs, stats_info);
						
						// clean up the array allocated
						printf("free the big array for current chromosome %s\n", chrom_tracking->chromosome_ids[i]);
						if (chrom_tracking->coverage[i]) {
							//free(chrom_tracking->coverage[i]);
							chrom_tracking->coverage[i] = realloc(chrom_tracking->coverage[i], sizeof(uint16_t));
						}
						chrom_tracking->chromosome_status[i] = 3;
					}
                }
            }

			// release the allocated chunk of buffer for alignment reads after they have been processed!
			//printf("cleaning the read buffer hash...\n\n");
			//readBufferDestroy(read_buff);
			//free(read_buff);
			fflush(stdout);
		}
	}

	// Now need to write the report
	writeReport(stats_info, "the_output.txt");

	sam_close(sfd);

	bam_hdr_destroy(header);
	userInputDestroy(user_inputs);
	chromosomeTrackingDestroy(chrom_tracking);

	if (target_bed_info)
		cleanBedInfo(target_bed_info);

	if (Ns_bed_info)
		cleanBedInfo(Ns_bed_info);

	if (target_buffer_hash) 
		cleanKhashStr(target_buffer_hash);

	if (Ns_buffer_hash)
		cleanKhashStr(Ns_buffer_hash);

	statsInfoDestroy(stats_info);

	if (user_inputs)
		free(user_inputs);

	return 0;
}
