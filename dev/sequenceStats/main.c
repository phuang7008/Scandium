/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  using openmp for parallel processing sequencing data
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
#include "user_defined_annotation.h"

// for khash: key -> string (char*);    value -> string (char*)
//
int khStrStr = 30;

// for khash: key -> string (char*);	value -> string array (stringArray*)
//
int khStrStrArray = 31;

// For khash: key -> string (char*);	value -> int
//
int khStrInt = 34;

// For khash: key ->string (char*);		value -> Low_Coverage_Genes*
//
int khStrLCG = 33;

// For khash: key ->string (char*);		value -> Transcript_Percentage*
//
int khStrGTP = 35;

int main(int argc, char *argv[]) {
	//fprintf(stderr, "Starting ... %ld\n", time(NULL));

	// get user input options and then processing it accordingly
	// and output all the options to the end user
	//
	User_Input *user_inputs = userInputInit();
	processUserOptions(user_inputs, argc, argv);
	outputUserInputOptions(user_inputs);

	// now need to setup a Stats_Info variable to track various statistical information
	//
    Stats_Info *stats_info = calloc(1, sizeof(Stats_Info));
	statsInfoInit(stats_info);

	// now for the bam/cram file open it for read
	//
    samFile *sfd = sam_open(user_inputs->bam_file, "r");
    if (sfd == 0) {
        fprintf(stderr, "Can not open file %s\n", argv[1]);
        return -1;
    }

    // use sam_hdr_read to process both bam and cram header
	//
    bam_hdr_t *header=NULL;
    if ((header = sam_hdr_read(sfd)) == 0) return -1;

	// setup a tracking variable to track chromosome working status 
	//
    Chromosome_Tracking *chrom_tracking = calloc(1, sizeof(Chromosome_Tracking));
	
	// Target_Buffer_Status need to be set even though there is no target or Ns region specified
	// It is because one of the method processRecord() need chromosome lengths information to be set
	//
	Target_Buffer_Status *target_buffer_status = NULL;
	
	// setup a variable to store chromosomes that specified by the user
	//
	khash_t(khStrInt) *wanted_chromosome_hash = kh_init(khStrInt);
	uint32_t i=0;
	uint32_t num_of_chroms = 0;
	
	// because hash keys are not in order, therefore, I need to store the chromosome ids in an array
	// to make them the same order as those in bam/cram file
	//
	if (user_inputs->chromosome_bed_file != NULL) {
		loadWantedChromosomes(wanted_chromosome_hash, user_inputs, stats_info);
		num_of_chroms = chromosomeTrackingInit2(wanted_chromosome_hash, chrom_tracking, header);

		// here we need to verify if the chromosome naming convention matches 
		// between the bam/cram file and the chromosome bed file specified by the end user
		//
		checkNamingConvention(header, wanted_chromosome_hash);

		target_buffer_status = calloc(num_of_chroms, sizeof(Target_Buffer_Status));
		TargetBufferStatusInit2(target_buffer_status, wanted_chromosome_hash, num_of_chroms);
	} else {
		loadGenomeInfoFromBamHeader(wanted_chromosome_hash, header, stats_info, user_inputs);
		num_of_chroms = header->n_targets;
		chromosomeTrackingInit1(num_of_chroms, chrom_tracking, wanted_chromosome_hash, header);

		target_buffer_status = calloc(num_of_chroms, sizeof(Target_Buffer_Status));
		TargetBufferStatusInit(target_buffer_status, header);
	}

	fprintf(stderr, "The total genome bases is %"PRIu32"\n", stats_info->cov_stats->total_genome_bases);

	// For target bed file and Ns region bed file
	//
	Bed_Info *target_bed_info=NULL, *Ns_bed_info=NULL;				

	if (N_FILE_PROVIDED) {      // the file that contains regions of Ns in the reference genome
        Ns_bed_info = calloc(1, sizeof(Bed_Info));
        processBedFiles(user_inputs, Ns_bed_info, stats_info, target_buffer_status, header, wanted_chromosome_hash, 2);
		fprintf(stderr, "The Ns base is %"PRIu32"\n", stats_info->cov_stats->total_Ns_bases);
    }

	if (TARGET_FILE_PROVIDED) {
        target_bed_info = calloc(1, sizeof(Bed_Info));
        processBedFiles(user_inputs, target_bed_info, stats_info, target_buffer_status, header, wanted_chromosome_hash, 1);
    }

	// there are two types of annotations available, MySQL or User_Defined_Database
	// for MySQL connection
	//
	Databases *dbs = NULL;
	khash_t(khStrInt) *hgmd_genes = NULL;
	khash_t(khStrInt) *hgmd_transcripts = NULL;
	
	// For quick processing, it would be very slow if we query MySQL for every single gene/exon/cds regions
	// Instead, we will fetch them only once and store them in the regions defined below
	// Note: We only setup MySQL database connection if TARGET_FILE_PROVIDED is on
	// For capture targets, we always use MySQL for the annotation and coverage percentage calculation
	//
	Regions_Skip_MySQL *intronic_regions    = NULL;
	Regions_Skip_MySQL *exon_regions        = NULL;

	// for user-defined database
	//
	User_Defined_Database_Wrapper *udd_wrapper = NULL;

	// We use Raw_User_Defined_Database to store everything from user_defined_database
	// so that we don't have to process the file over and over
	//
	Raw_User_Defined_Database * raw_user_defined_database = NULL;

	// If user provides user-defined-database, we need to record them here
	// Note, as there are many duplicates, we need to move all duplicates
	//
	khash_t(khStrInt) *user_defined_targets=NULL;
	Bed_Info *user_defined_bed_info=NULL;

	// Store cds length and cds count informtion from user-defined-database
	//
	khash_t(khStrInt) *cds_lengths=NULL;
	khash_t(khStrInt) *cds_counts =NULL;

	if (user_inputs->annotation_on && TARGET_FILE_PROVIDED) {
		exon_regions = calloc(1, sizeof(Regions_Skip_MySQL));

		if (USER_DEFINED_DATABASE) {
			udd_wrapper = calloc(1, sizeof(User_Defined_Database_Wrapper));
			raw_user_defined_database = calloc(1, sizeof(Raw_User_Defined_Database));

			user_defined_targets = kh_init(khStrInt);
			cds_lengths = kh_init(khStrInt);
			cds_counts  = kh_init(khStrInt);

			getUserDefinedDatabaseInfo(user_inputs, udd_wrapper, cds_lengths, cds_counts, user_defined_targets);
			processUserDefinedDatabase(user_inputs, exon_regions, udd_wrapper, raw_user_defined_database, cds_lengths, cds_counts);

		} else {
			// Initialize all annotation related variables
			//
			dbs = calloc(1, sizeof(Databases));
			databaseSetup(dbs, user_inputs);

			intronic_regions    = calloc(1, sizeof(Regions_Skip_MySQL));

			regionsSkipMySQLInit(dbs, intronic_regions, user_inputs, 2);
			regionsSkipMySQLInit(dbs, exon_regions, user_inputs, 3);
		}

		if (HGMD_PROVIDED && hgmd_genes == NULL) {
			if (dbs == NULL) {
				dbs = calloc(1, sizeof(Databases));
				databaseSetup(dbs, user_inputs);
			}

			hgmd_genes = kh_init(khStrInt);
			hgmd_transcripts = kh_init(khStrInt);
			recordHGMD(dbs, user_inputs, hgmd_genes, hgmd_transcripts);
		}
	}

	// initialize the Gene_Transcript_Percentage variable
	//
	khash_t(khStrGTP) *gene_transcript_percentage_hash = NULL;
	if (TARGET_FILE_PROVIDED)
		gene_transcript_percentage_hash = kh_init(khStrGTP);

	// can't set to be static as openmp won't be able to handle it
	// check the bam/cram file size first
	//
	uint64_t input_bam_file_size = check_file_size(user_inputs->bam_file);

	uint32_t total_chunk_of_reads = 300000;		// for small bam/cram file
	if (input_bam_file_size > 5000000000)		// anything > 5Gb
		total_chunk_of_reads = 1400000;			// Good for 3 threads with 9gb  of memory

	// try to allocate the bam1_t array here for each thread, so they don't have to create and delete the array at each loop
	// Here there are two layers of read_buff are created through calloc(), which need to be freed at the end for two levels
	//
	uint32_t j;
	Read_Buffer *read_buff = calloc(user_inputs->num_of_threads, sizeof(Read_Buffer));
	for (i=0; i<user_inputs->num_of_threads; i++) {
		read_buff[i].chunk_of_reads = calloc(total_chunk_of_reads, sizeof(bam1_t*));
		for (j=0; j<total_chunk_of_reads; j++)
			read_buff[i].chunk_of_reads[j]=NULL;

		read_buff[i].capacity = total_chunk_of_reads;
		read_buff[i].size = 0;
	}

	fflush(stdout);

	// set random seed (should only be called ONCE)
	//
	if(user_inputs->percentage < 1.0)
		srand((uint32_t)time(NULL));    // set random seed

	// now let's do the parallelism
	//
	uint64_t total_reads=0;
    while(chrom_tracking->more_to_read) {
#pragma omp parallel shared(read_buff, chrom_tracking) num_threads(user_inputs->num_of_threads)
      {
        //int num_of_threads = omp_get_num_threads();	
        int thread_id = omp_get_thread_num();
        readBufferInit(&read_buff[thread_id]);		// third level of memory at created through bam_init1()
        uint32_t num_records = 0;
        //printf("Before File Reading: number of threads is %d and current thread id is %d\n", num_of_threads, thread_id);

#pragma omp critical 
        {
          // this part of the code need to run atomically, that is only one thread should allow access to read
          num_records = readBam(sfd, header, chrom_tracking, &read_buff[thread_id]);
          read_buff[thread_id].size = num_records;
		  total_reads += num_records;
        }

        Coverage_Stats *cov_stats = calloc(1, sizeof(Coverage_Stats)); 
		coverageStatsInit(cov_stats);
        khash_t(str) *coverage_hash = kh_init(str);		// hash_table using string as key

        // can not use the else {} with the previous if {} block, otherwise the barrier will wait forever!
		//
        if (num_records > 0) {
          printf("Reading: %"PRIu32" records\t\tTotal: %"PRIu64"\t\tThread id: %d.\n", num_records, total_reads, thread_id);

          processBamChunk(user_inputs, cov_stats, coverage_hash, header, &read_buff[thread_id], target_buffer_status, thread_id, wanted_chromosome_hash, chrom_tracking);

		}

		// release the allocated chunk of buffer for aligned reads after they have been processed!
		//
        //printf("cleaning the read buffer hash for thread %d...\n\n", thread_id);
        readBufferDestroy(&read_buff[thread_id]);		// third level of memory allocation is destroyed and freed here

        if (num_records == 0) {
          printf("No more to read for thread %d !!!!!!!!!!!!\n", thread_id);
          if (chrom_tracking->more_to_read) chrom_tracking->more_to_read = false;
		}

#pragma omp critical
        {
          if (num_records > 0) {
            combineThreadResults(chrom_tracking, coverage_hash);
            combineCoverageStats(stats_info, cov_stats);

            // since all reads have been processed for current chromosome, we need to set the status to 2
			// 
            if (!chrom_tracking->more_to_read) {
              for(i=0; i<num_of_chroms; i++) {
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
            for (i=0; i<num_of_chroms; i++) {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                // check to see if any of the chromosomes has finished. If so, write the results out
                // for the whole genome, we need to use the file that contains regions of all Ns in the reference
                // As they will be not used, so we are going to set the count info in these regions to 0
				//
                if (N_FILE_PROVIDED)
                  zeroAllNsRegions(chrom_tracking->chromosome_ids[i], Ns_bed_info, chrom_tracking, target_buffer_status);
              }
            }
          }
		}

#pragma omp barrier

		i = 0;
		while (i<num_of_chroms) {
          //if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
          //  printf("Chromosome id %s in thread id %d has finished processing, now dumping\n", chrom_tracking->chromosome_ids[i], thread_id);
          //}

#pragma omp sections
          {
#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                printf("Thread %d is now producing coverage information for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);
						
				// The following will be produced no matter whether annotation_on is set or not
				//
                writeCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info, intronic_regions, exon_regions);

                // now write the off target regions with high coverage into a wig file if the Write_WIG flag is set
				//
                if (TARGET_FILE_PROVIDED && user_inputs->Write_WIG)
                  produceOffTargetWigFile(chrom_tracking, chrom_tracking->chromosome_ids[i], target_bed_info, user_inputs, stats_info);
              }
            }

#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
				// For Whole Genome Annotations (use MySQL for the annotation)
				// if the annotation is not on, it will just output . . . . . . . )
				//
                if (user_inputs->wgs_coverage) {
                  printf("Thread %d is now writing WGS annotation for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);
                  writeAnnotations(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info, intronic_regions, exon_regions);

				  // if user specifies the range information (usually for graphing purpose), need to handle it here
				  //
				  printf("Thread %d is now writing coverage uniformity data for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);
				  coverageRangeInfoForGraphing(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info);
                }
              }
            }

#pragma omp section
			{
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {

                if (user_inputs->annotation_on && TARGET_FILE_PROVIDED) {
                  printf("Thread %d is now calculating gene/transcript/cds coverage percentage for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[i]);

                  khash_t(khStrLCG) *transcript_hash = kh_init(khStrLCG);
                  khash_t(khStrStrArray) *gene_transcripts = kh_init(khStrStrArray);

				  if (USER_DEFINED_DATABASE) {
					khash_t(khStrLCG) *user_defined_cds_gene_hash = kh_init(khStrLCG);

					userDefinedGeneCoverageInit(user_defined_cds_gene_hash, chrom_tracking->chromosome_ids[i], raw_user_defined_database, gene_transcripts);
					calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info, user_defined_cds_gene_hash);
					transcriptPercentageCoverageInit(chrom_tracking->chromosome_ids[i], transcript_hash, user_defined_cds_gene_hash);
					storeGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, user_inputs, transcript_hash, gene_transcripts, hgmd_genes, hgmd_transcripts, gene_transcript_percentage_hash);

					// clean-up
					//
					genePercentageCoverageDestroy(user_defined_cds_gene_hash);
					user_defined_cds_gene_hash=NULL;
				  } else {
				    // Allocate memories for the low_cov_gene_hash, which get all its content from the MySQL database (official RefSeq DB)
				    // For calculating the percentage of gene bases with low coverge for capture only
				    // and use the intersect regions between refseq_cds_genes for official annotation and low_cov_genes for targets
					//
					khash_t(khStrLCG) *low_cov_gene_hash = kh_init(khStrLCG);

				    genePercentageCoverageInit(low_cov_gene_hash, chrom_tracking->chromosome_ids[i], dbs, user_inputs, gene_transcripts);
					intersectTargetsAndRefSeqCDS(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, low_cov_gene_hash);

					calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, chrom_tracking, user_inputs, stats_info, low_cov_gene_hash);
				    transcriptPercentageCoverageInit(chrom_tracking->chromosome_ids[i], transcript_hash, low_cov_gene_hash);

					storeGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, user_inputs, transcript_hash, gene_transcripts, hgmd_genes, hgmd_transcripts, gene_transcript_percentage_hash);

				    // clean-up the memory space
				    //
				    genePercentageCoverageDestroy(low_cov_gene_hash);
					low_cov_gene_hash=NULL;
                  }

                  // More clean-ups
				  //
                  genePercentageCoverageDestroy(transcript_hash);
				  cleanKhashStrStrArray(gene_transcripts);
				  transcript_hash=NULL;
				  gene_transcripts=NULL;
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
        //printf("End of while loop before flush for thread %d\n", thread_id);
      }
	  printf("\n");
      fflush(stdout);
    }

	sam_close(sfd);

	/* calculate the uniformity metric*/
	khash_t(m32) *cov_freq_dist = kh_init(m32);
	if (user_inputs->wgs_coverage) {
		// Autosomes only including alt decoys, but without X and Y chromosomes
		//
		//calculateUniformityMetrics(stats_info, user_inputs, wanted_chromosome_hash, 1, 0);		

		// Primary autosomes only without X, Y, alt and decoys!
		//
		calculateUniformityMetrics(stats_info, user_inputs, wanted_chromosome_hash, cov_freq_dist, 1, 1);

		// All (including X and Y chromosomes), also include alt, decoys
		//
		//calculateUniformityMetrics(stats_info, user_inputs, wanted_chromosome_hash, 0, 0);

		// All Primaries (including X and Y chromosomes), BUT without alt, decoy
		//
		//calculateUniformityMetrics(stats_info, user_inputs, wanted_chromosome_hash, cov_freq_dist, 0, 1);
	}

	// Now need to write the report
	//
	writeReport(stats_info, user_inputs);

	if (user_inputs->wgs_coverage)
		outputFreqDistribution(user_inputs, cov_freq_dist);

	// output gene converage
	//
	if (TARGET_FILE_PROVIDED)
		outputGeneCoverage(gene_transcript_percentage_hash, user_inputs);

	/* clean-up everything*/

	cleanKhashInt(cov_freq_dist);

	chromosomeTrackingDestroy(chrom_tracking);

	if (target_bed_info != NULL)
		cleanBedInfo(target_bed_info);

	if (Ns_bed_info != NULL)
		cleanBedInfo(Ns_bed_info);

	if (user_defined_bed_info != NULL)
		cleanBedInfo(user_defined_bed_info);

	if (target_buffer_status) {
		for (i=0; i<num_of_chroms; i++) {
			free(target_buffer_status[i].status_array);
		}
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
				read_buff[i].chunk_of_reads=NULL;
			}
    	}
		free(read_buff);
	}

	if (intronic_regions != NULL) regionsSkipMySQLDestroy(intronic_regions, 2);
	if (exon_regions != NULL) regionsSkipMySQLDestroy(exon_regions, 3);

	userInputDestroy(user_inputs);
	bam_hdr_destroy(header);

	// MYSQL clean-up
	//
	if (dbs != NULL) {
		mysql_close(dbs->con);
		mysql_library_end();

		databaseCleanUp(dbs);
		free(dbs);
	}

	if (wanted_chromosome_hash != NULL)
		cleanKhashStrInt(wanted_chromosome_hash);

	if (hgmd_genes != NULL)
		cleanKhashStrInt(hgmd_genes);

	if (hgmd_transcripts != NULL)
		cleanKhashStrInt(hgmd_transcripts);

	if (gene_transcript_percentage_hash != NULL)
		cleanGeneTranscriptPercentage(gene_transcript_percentage_hash);

	if (udd_wrapper) 
		cleanUserDefinedDatabase(udd_wrapper);

	if (raw_user_defined_database) 
		cleanRawUserDefinedDatabase(raw_user_defined_database);

	// now need to clean-up the khash_t variables                                                     
	//
	if (cds_lengths) 
		cleanKhashStrInt(cds_lengths);                                                                    

	if (cds_counts) 
		cleanKhashStrInt(cds_counts);

	if (user_defined_targets) 
		cleanKhashStrInt(user_defined_targets);

	printf("Program Finished Successfully!\n\n");

	return 0;
}
