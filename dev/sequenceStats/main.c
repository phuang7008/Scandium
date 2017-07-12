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

	// now process the target bed file if it is available
	// First, let's get the total number of lines(items or count) within the target file
    int target_line_count = getTargetCount(user_inputs->target_file);

    // Now initialize the storage for arrays that are used to store target coordinates
    // Do need to remember to free these allocated memories at the end of the program!!!
    Target_Coords *target_coords = (Target_Coords*)  calloc(target_line_count, sizeof(Target_Coords));

    // Finally, load target file again and store the target information (starts, stops and chromosome ids)
    loadTargets(user_inputs->target_file, target_coords);

    // Output stored coordinates for verification
    /*for (i=0; i<target_line_count; i++) {
        printf("%s\t%ld\t%ld\n", targets[i].chr, targets[i].start, targets[i].end);
    }*/

	// Now we are going to generate target-buffer lookup table for all the loaded targets
    // we will store targets and buffers information based on chromosome ID
    // we will have 22 + X + Y = 24 chromosomes, Here X=23 and Y will be 24
	int i;
    khash_t(m32) *target_buffer_hash[24];
    for (i = 0; i < 26; i++) {
        target_buffer_hash[i] = kh_init(m32);
    }
    //generateTargetBufferLookupTable(target_coords, target_line_count, &target_buffer_hash);


	samFile *sfd = sam_open("NA19431-ny-r1-NA19431-o-2.recal.cram", "r");
	//samFile *sfd = sam_open("/stornext/snfs5/next-gen/scratch/phuang/dev/excid/test.bam", "r");
	if (sfd == 0) {
		fprintf(stderr, "Can not open file %s\n", argv[1]);
		return -1;
	}

	// use sam_hdr_read to process both bam and cram header
	bam_hdr_t *header=NULL;
	if ((header = sam_hdr_read(sfd)) == 0) return -1;
	//if ((header = bam_hdr_read(sfd->fp.bgzf)) == 0) return -1;

	int thread_id=0;
	bool more_to_read=true;
	uint32_t total_chunk_of_reads = 1000000;	// can't set to be static as openmp won't be able to handle it

	// setup a tracking variable to track chromosome working status 
	Chromosome_Tracking *chrom_tracking = calloc(1, sizeof(Chromosome_Tracking));
	chrom_tracking->number_tracked = 0;

	// initialize current working chromosome id
    strcpy(CURRENT_CHROMOSOME_ID, "nothing");

	NUM_OF_THREADS = 8;
	// now let's do the parallelism
	while(more_to_read) {
#pragma omp parallel firstprivate(thread_id, total_chunk_of_reads) shared(more_to_read) num_threads(NUM_OF_THREADS)
		{
			int num_of_threads = omp_get_num_threads();	
			thread_id = omp_get_thread_num();
			printf("Before File Reading: number of threads is %d and current thread id is %d\n",num_of_threads, thread_id);

			Read_Buffer *read_buff;
			read_buff = malloc(sizeof(Read_Buffer));
			read_buff->chunk_of_reads = calloc(total_chunk_of_reads, sizeof(bam1_t));
			read_buff->size = total_chunk_of_reads;
			//bam1_t *read_buff[total_chunk_of_reads];
			read_buff_init(read_buff);
			uint32_t num_records = 0;

#pragma omp critical 
			{
				// this part of the code need to run atomically, that is only one thread should allow access to read
				num_records = read_bam(sfd, header, more_to_read, read_buff);
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
			process_chunk_of_bam(thread_id, coverage_hash, header, read_buff);

#pragma omp critical
			{
				printf("Before writing to the coverage array for Thread %d\n", thread_id);
				combine_thread_results(chrom_tracking, coverage_hash, header);
				printf("after writing to the coverage array for Thread %d\n", thread_id);

			}

			// setup a barrier here and wait for every one of them to reach this point!
#pragma omp barrier

			// check to see if any of the chromosomes has finished. If so, write the results out
			int j;
			for (i=0; i<chrom_tracking->number_tracked; i++) {
                if ( chrom_tracking->chromosome_status[i] == 2) {
#pragma omp critical
					{
						printf(">%s\n", chrom_tracking->chromosome_ids[i]);
						for (j=1; j<chrom_tracking->chromosome_lengths[i]; j++) {
							printf("%d\t", chrom_tracking->coverage[i][j]);
						}
					}
					// clean up the array allocated
					free(chrom_tracking->coverage[i]);
                }
            }

			// release the allocated chunk of buffer for alignment reads after they have been processed!
			read_buff_destroy(read_buff);
		}
	}

	bam_hdr_destroy(header);
	userInputDestroy(user_inputs);
	chromosome_tracking_destroy(chrom_tracking);

	return 0;
}
