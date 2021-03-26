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

#include "annotation.h"
#include "reports.h"
#include "stats.h"
#include "targets.h"
#include "terms.h"
#include "user_defined_annotation.h"
#include "user_inputs.h"
#include "utils.h"

// for khash: key -> string (char*);    value -> string (char*)
//
int khStrStr = 30;

// for khash: key -> string (char*);    value -> string array (stringArray*)
//
int khStrStrArray = 31;

// For khash: key -> string (char*);    value -> int
//
int khStrInt = 34;

// For khash: key ->string (char*);        value -> Low_Coverage_Genes*
//
int khStrLCG = 33;

// For khash: key ->string (char*);        value -> Transcript_Percentage*
//
int khStrGTP = 35;

// for khash: key -> string (char*); value -> Regions_Skip_MySQL
//
//int khStrRSM = 36;

int main(int argc, char *argv[]) {
    //fprintf(stderr, "Starting ... %ld\n", time(NULL));

    // get user input options and then processing it accordingly
    // and output all the options to the end user
    //
    User_Input *user_inputs = userInputInit();
    processUserOptions(user_inputs, argc, argv);
    //outputUserInputOptions(user_inputs);

    // now need to setup a Stats_Info variable to track various statistical information
    //
    Stats_Info *stats_info = calloc(1, sizeof(Stats_Info));
    statsInfoInit(stats_info, user_inputs);

    // now for the bam/cram file open it for read
    //
    samFile *sfd = sam_open(user_inputs->bam_file, "r");
    if (sfd == 0) {
        fprintf(stderr, "Cannot open file \n%s\n", user_inputs->bam_file);
        return -1;
    }

    // Set the reference if it is the cram file
    //
    char * fn_ref = 0;
    if (user_inputs->reference_file) {
        fn_ref = getReferenceFaiPath(user_inputs->reference_file);

        if (hts_set_fai_filename(sfd, fn_ref) != 0) {
            fprintf(stderr, "Failed to use reference at hts_set_fai_filename() \"%s\".\n", fn_ref);
            if (fn_ref) free(fn_ref);
            return -1;
        }
    } else {
        if ( sfd->is_cram || sfd->format.format == cram ) {
            fprintf(stderr, "Please provide the reference sequences for the input CRAM file \n%s\n", user_inputs->bam_file);
            return -1;
        }
    }

    // check if bam file named as .cram and cram file named as .bam
    //
    checkFileExtension(user_inputs, sfd);

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
    
    // because hash keys are not in order, therefore, I need to store the chromosome ids in an array
    // to make them the same order as those in bam/cram file
    //
    if (user_inputs->chromosome_bed_file != NULL) {
        loadWantedChromosomes(wanted_chromosome_hash, user_inputs, stats_info);
        chromosomeTrackingInit2(wanted_chromosome_hash, chrom_tracking, header);

        // here we need to verify if the chromosome naming convention matches 
        // between the bam/cram file and the chromosome bed file specified by the end user
        //
        checkNamingConvention(header, wanted_chromosome_hash);

        target_buffer_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(Target_Buffer_Status));
        TargetBufferStatusInit2(target_buffer_status, wanted_chromosome_hash);
    } else {
        loadGenomeInfoFromBamHeader(wanted_chromosome_hash, header, stats_info, user_inputs);
        chrom_tracking->number_of_chromosomes = header->n_targets;
        chromosomeTrackingInit1(chrom_tracking, wanted_chromosome_hash, header);

        target_buffer_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(Target_Buffer_Status));
        TargetBufferStatusInit(target_buffer_status, header);
    }

    fprintf(stderr, "The total genome bases is %"PRIu64"\n", stats_info->wgs_cov_stats->total_genome_bases);

    // For target bed file and Ns region bed file
    //
    Bed_Info **target_bed_info=NULL;
    Bed_Info *Ns_bed_info=NULL;                

    if (N_FILE_PROVIDED) {      // the file that contains regions of Ns in the reference genome
        Ns_bed_info = calloc(1, sizeof(Bed_Info));
        processBedFiles(user_inputs, Ns_bed_info, stats_info, target_buffer_status, 
                wanted_chromosome_hash, user_inputs->n_file,chrom_tracking->number_of_chromosomes,0, 2);
        fprintf(stderr, "The Ns base is %"PRIu32"\n", stats_info->wgs_cov_stats->total_Ns_bases);
    }

    if (TARGET_FILE_PROVIDED) {
        target_bed_info = calloc(user_inputs->num_of_target_files, sizeof(Bed_Info*));
        for (i=0; i<user_inputs->num_of_target_files; i++) {
            target_bed_info[i] = calloc(1, sizeof(Bed_Info));
            processBedFiles(user_inputs, target_bed_info[i], stats_info, target_buffer_status, 
                    wanted_chromosome_hash, user_inputs->target_files[i], chrom_tracking->number_of_chromosomes, i, 1);
        }
    }

    // there are two types of annotations available, MySQL or User_Defined_Database
    // For quick processing, it would be very slow if we query MySQL for every single gene/exon/cds regions
    // Instead, we will fetch them only once and store them in the regions defined below
    // Note: We only setup MySQL database connection if WGS annotation_on is true
    // For capture targets (-t is on), we always use user_defined_annotation database files if they are provided.
    //      However, if user doesn't provide user_defined_annotatioin files, we will useMySQL for the annotation
    //
    Regions_Skip_MySQL **exon_regions    = NULL;

    // We use Raw_User_Defined_Database to store everything from user_defined_database
    // so that we don't have to process the file over and over
    //
    Raw_User_Defined_Database **raw_user_defined_databases = NULL;

    // If user provides user-defined-database, we need to record them here
    // Note, as there are many duplicates, we need to move all duplicates
    //
    khash_t(khStrInt) *user_defined_targets=NULL;
    Bed_Info *user_defined_bed_info=NULL;

    if (TARGET_FILE_PROVIDED) {

        if (USER_DEFINED_DATABASE) {
            exon_regions = calloc(user_inputs->num_of_target_files, sizeof(Regions_Skip_MySQL*));
            raw_user_defined_databases = calloc(user_inputs->num_of_target_files, sizeof(Raw_User_Defined_Database*));

            for (i=0; i<user_inputs->num_of_annotation_files; i++) {
                // for user-defined database
                //
                User_Defined_Database_Wrapper *udd_wrapper = NULL;

                // Store cds length and cds count informtion from user-defined-database
                //
                khash_t(khStrInt) *cds_lengths=NULL;
                khash_t(khStrInt) *cds_counts =NULL;

                // need to check if the annotation format is correct!
                // Stop is the format is wrong!
                //
                checkAnnotationFormat(user_inputs);

                udd_wrapper = calloc(1, sizeof(User_Defined_Database_Wrapper));

                user_defined_targets = kh_init(khStrInt);
                cds_lengths = kh_init(khStrInt);
                cds_counts  = kh_init(khStrInt);

                getUserDefinedDatabaseInfo(user_inputs, udd_wrapper, cds_lengths, cds_counts, user_defined_targets, i);

                exon_regions[i] = calloc(1, sizeof(Regions_Skip_MySQL));
                raw_user_defined_databases[i] = calloc(1, sizeof(Raw_User_Defined_Database));
                processUserDefinedDatabase(user_inputs, exon_regions[i], udd_wrapper, raw_user_defined_databases[i], cds_lengths, cds_counts, i);

                if (udd_wrapper)
                    cleanUserDefinedDatabase(udd_wrapper);

                if (cds_lengths)
                    cleanKhashStrInt(cds_lengths);

                if (cds_counts)
                    cleanKhashStrInt(cds_counts);

                if (user_defined_targets)
                    cleanKhashStrInt(user_defined_targets);
            }
        }
    }

    // for MySQL connection
    //
    Databases *dbs = NULL;
    khash_t(khStrInt) *hgmd_genes = NULL;
    khash_t(khStrInt) *hgmd_transcripts  = NULL;
    Regions_Skip_MySQL *intronic_regions = NULL;
    Regions_Skip_MySQL* db_exon_regions  = NULL;

    if (user_inputs->wgs_annotation_on) {
        // Initialize all annotation related variables
        //
        if (dbs == NULL) setupMySQLDB(&dbs, user_inputs);

        intronic_regions = calloc(1, sizeof(Regions_Skip_MySQL));
        db_exon_regions  = calloc(1, sizeof(Regions_Skip_MySQL));

        regionsSkipMySQLInit(dbs, intronic_regions, 2);
        regionsSkipMySQLInit(dbs, db_exon_regions, 3);
    }

    if (HGMD_PROVIDED && hgmd_genes == NULL) {
        if (dbs == NULL) setupMySQLDB(&dbs, user_inputs);

        hgmd_genes = kh_init(khStrInt);
        hgmd_transcripts = kh_init(khStrInt);
        recordHGMD(dbs, hgmd_genes, hgmd_transcripts);
    }

    // initialize the Gene_Transcript_Percentage variable
    //
    khash_t(khStrGTP) **gene_transcript_percentage_hash = NULL;
    if (TARGET_FILE_PROVIDED) {
        gene_transcript_percentage_hash = calloc(user_inputs->num_of_target_files, sizeof(khash_t(khStrGTP)*));
        for (i=0; i<user_inputs->num_of_target_files; i++)
            gene_transcript_percentage_hash[i] = kh_init(khStrGTP);
        
        // here we need to check if we have more capture files than annotation files
        //
        uint8_t num_of_additional_annotations = user_inputs->num_of_target_files > user_inputs->num_of_annotation_files;
        if (num_of_additional_annotations > 0) {
            if (exon_regions == NULL) {
                exon_regions = calloc(user_inputs->num_of_target_files, sizeof(Regions_Skip_MySQL*));
                raw_user_defined_databases = calloc(user_inputs->num_of_target_files, sizeof(Raw_User_Defined_Database*));
            }

            if (db_exon_regions == NULL) {
                intronic_regions = calloc(1, sizeof(Regions_Skip_MySQL));
                db_exon_regions  = calloc(1, sizeof(Regions_Skip_MySQL));

                if (dbs == NULL) setupMySQLDB(&dbs, user_inputs);
                regionsSkipMySQLInit(dbs, intronic_regions, 2);
                regionsSkipMySQLInit(dbs, db_exon_regions, 3);
            }

            for (i=user_inputs->num_of_annotation_files; i<user_inputs->num_of_target_files; i++) {
                exon_regions[i] = db_exon_regions;
            }
        }
    }

    setupOutputReportFiles(user_inputs);

    // can't set to be static as openmp won't be able to handle it
    // check the bam/cram file size first
    //
    uint64_t input_bam_file_size = check_file_size(user_inputs->bam_file);

    uint32_t total_chunk_of_reads = 200000;        // for small bam/cram file
    if (input_bam_file_size > 5000000000)        // anything > 5Gb
        total_chunk_of_reads = 1000000;            // Good for 3 threads with 9gb  of memory

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
        readBufferInit(&read_buff[thread_id]);      // third level of memory at created through bam_init1()
        uint32_t num_records = 0;
        //printf("Before File Reading: number of threads is %d and current thread id is %d\n", num_of_threads, thread_id);

#pragma omp critical 
        {
          // this part of the code need to run atomically, that is only one thread should allow access to read
          num_records = readBam(sfd, header, chrom_tracking, &read_buff[thread_id]);
          read_buff[thread_id].size = num_records;
          total_reads += num_records;
        }

        Stats_Info *tmp_stats_info = calloc(1, sizeof(Stats_Info));
        statsInfoInit(tmp_stats_info, user_inputs);
        //Coverage_Stats *cov_stats = calloc(1, sizeof(Coverage_Stats)); 
        //coverageStatsInit(cov_stats);
        khash_t(str) *coverage_hash = kh_init(str);     // hash_table using string as key

        // can not use the else {} with the previous if {} block, otherwise the barrier will wait forever!
        //
        if (num_records > 0) {
          printf("Reading: %"PRIu32" records\t\tTotal: %"PRIu64"\t\tThread id: %d.\n", num_records, total_reads, thread_id);

          processBamChunk(user_inputs, tmp_stats_info, coverage_hash, header, &read_buff[thread_id], 
                  target_buffer_status, thread_id, wanted_chromosome_hash, chrom_tracking->number_of_chromosomes);

        }

        // release the allocated chunk of buffer for aligned reads after they have been processed!
        //
        //printf("cleaning the read buffer hash for thread %d...\n\n", thread_id);
        readBufferDestroy(&read_buff[thread_id]);       // third level of memory allocation is destroyed and freed here

        if (num_records == 0) {
          printf("No more to read for thread %d !!!!!!!!!!!!\n", thread_id);
          if (chrom_tracking->more_to_read) chrom_tracking->more_to_read = false;
        }

#pragma omp critical
        {
          if (num_records > 0) {
            combineThreadResults(chrom_tracking, coverage_hash);
            combineCoverageStats(stats_info, tmp_stats_info, user_inputs);

            // if all reads have been processed for the entire file, we need to set the status to 2 for all
            // 
            if (!chrom_tracking->more_to_read) {
              for(i=0; i<(uint32_t)chrom_tracking->number_of_chromosomes; i++) {
                if (chrom_tracking->chromosome_status[i] == 1)
                  chrom_tracking->chromosome_status[i] = 2;
              }
            }
          }
        }

        cleanKhashStr(coverage_hash, 1);
        statsInfoDestroy(tmp_stats_info, user_inputs);
        tmp_stats_info = NULL;   // to eleminate dangling pointer

// setup a barrier here and wait for every one of them to reach this point!
#pragma omp barrier 

#pragma omp single
        {
          if (num_records > 0) {
            for (i=0; i<(uint32_t)chrom_tracking->number_of_chromosomes; i++) {
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
        while (i<(uint32_t)chrom_tracking->number_of_chromosomes) {

#pragma omp sections
          {
            if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2)
              printf("\n");

#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                printf("Thread %d is now producing coverage information for chromosome %s\n", 
                        thread_id, chrom_tracking->chromosome_ids[i]);
                        
                // The following will be produced no matter whether annotation_on is set or not
                //
                writeCoverage(chrom_tracking->chromosome_ids[i], target_bed_info, 
                        chrom_tracking, user_inputs, stats_info, intronic_regions, exon_regions);

                // now write the off target regions with high coverage into a wig file if the Write_WIG flag is set
                //
                //if (TARGET_FILE_PROVIDED && user_inputs->Write_WIG)
                if (TARGET_FILE_PROVIDED) {
                  int j;
                  for (j=0; j<user_inputs->num_of_target_files; j++) {
                    produceOffTargetWigFile(chrom_tracking, chrom_tracking->chromosome_ids[i], 
                            target_bed_info[j], user_inputs, stats_info, j);
                  }
                }
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
                  writeAnnotations(chrom_tracking->chromosome_ids[i], chrom_tracking, user_inputs, intronic_regions, db_exon_regions);

                  // if user specifies the range information (usually for graphing purpose), need to handle it here
                  //
                  printf("Thread %d is now writing coverage uniformity data for chromosome %s\n", 
                          thread_id, chrom_tracking->chromosome_ids[i]);
                  coverageRangeInfoForGraphing(chrom_tracking->chromosome_ids[i], chrom_tracking, user_inputs);
                }
              }
            }

#pragma omp section
            {
              if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {

                if (TARGET_FILE_PROVIDED) {
                  printf("Thread %d is now calculating gene/transcript/cds coverage percentage for chromosome %s\n", 
                          thread_id, chrom_tracking->chromosome_ids[i]);

                  // process one annotation file at a time
                  //
                  int x;
                  for (x=0; x<user_inputs->num_of_target_files; x++) {
                    khash_t(khStrLCG) *transcript_hash = kh_init(khStrLCG);
                    khash_t(khStrStrArray) *gene_transcripts = kh_init(khStrStrArray);

                    if (USER_DEFINED_DATABASE) {
                        khash_t(khStrLCG) *user_defined_cds_gene_hash = kh_init(khStrLCG);

                        if (raw_user_defined_databases[x]) {
                          userDefinedGeneCoverageInit(user_defined_cds_gene_hash, 
                                  chrom_tracking->chromosome_ids[i], raw_user_defined_databases[x], gene_transcripts);
                        } else {
                          genePercentageCoverageInit(user_defined_cds_gene_hash, 
                                  chrom_tracking->chromosome_ids[i], dbs, gene_transcripts);
                          intersectTargetsAndRefSeqCDS(chrom_tracking->chromosome_ids[i], 
                                  target_bed_info[x], chrom_tracking, user_defined_cds_gene_hash);
                        }

                        calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info[x], 
                                chrom_tracking, user_inputs, user_defined_cds_gene_hash);
                        transcriptPercentageCoverageInit(transcript_hash, user_defined_cds_gene_hash);
                        storeGenePercentageCoverage(chrom_tracking->chromosome_ids[i], user_inputs, 
                                transcript_hash, gene_transcripts, hgmd_transcripts, gene_transcript_percentage_hash, x);

                        // clean-up
                        genePercentageCoverageDestroy(user_defined_cds_gene_hash);
                        user_defined_cds_gene_hash=NULL;
                    } else {
                        // Allocate memories for the low_cov_gene_hash, which get all its content from the MySQL database (official RefSeq DB)
                        // For calculating the percentage of gene bases with low coverge for capture only
                        // and use the intersect regions between refseq_cds_genes for official annotation and low_cov_genes for targets
                        //
                        khash_t(khStrLCG) *low_cov_gene_hash = kh_init(khStrLCG);

                        genePercentageCoverageInit(low_cov_gene_hash, chrom_tracking->chromosome_ids[i], dbs, gene_transcripts);
                        intersectTargetsAndRefSeqCDS(chrom_tracking->chromosome_ids[i], target_bed_info[x], 
                                chrom_tracking, low_cov_gene_hash);

                        calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[i], target_bed_info[x], 
                                chrom_tracking, user_inputs, low_cov_gene_hash);
                        transcriptPercentageCoverageInit(transcript_hash, low_cov_gene_hash);

                        storeGenePercentageCoverage(chrom_tracking->chromosome_ids[i], user_inputs, 
                                transcript_hash, gene_transcripts, hgmd_transcripts, gene_transcript_percentage_hash, x);

                        // clean-up the memory space
                        genePercentageCoverageDestroy(low_cov_gene_hash);
                        low_cov_gene_hash=NULL;
                    }

                    // More clean-ups
                    genePercentageCoverageDestroy(transcript_hash);
                    cleanKhashStrStrArray(gene_transcripts);
                    transcript_hash=NULL;
                    gene_transcripts=NULL;
                  }
                }
              }
            }
          }

#pragma omp barrier 

#pragma omp single
          {
            if ( chrom_tracking->chromosome_ids[i] && chrom_tracking->chromosome_status[i] == 2) {
                printf("\n");

              // clean up the array allocated
              //
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
      }     // End of while loop before flush for thread

      //printf("\n");
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
        // need to mock the cov_freq_dist so that they won't affect cov_freq_dist calculation
        //
        khash_t(m32) *cov_freq_dist_XY = kh_init(m32);
        calculateUniformityMetrics(stats_info, user_inputs, wanted_chromosome_hash, cov_freq_dist_XY, 0, 1);
        cleanKhashInt(cov_freq_dist_XY);
    }

    // Now need to write the report
    //
    writeWGSReports(stats_info, user_inputs);
    writeCaptureReports(stats_info, user_inputs);

    if (user_inputs->wgs_coverage)
        outputFreqDistribution(user_inputs, cov_freq_dist);

    // output gene converage
    //
    if (TARGET_FILE_PROVIDED) {
        for (i=0; i<user_inputs->num_of_target_files; i++)
            outputGeneCoverage(gene_transcript_percentage_hash[i], user_inputs, i);
    }

    /* clean-up everything*/

    cleanKhashInt(cov_freq_dist);

    TargetBufferStatusDestroy(target_buffer_status, chrom_tracking->number_of_chromosomes);

    chromosomeTrackingDestroy(chrom_tracking);

    for (i=0; i<user_inputs->num_of_target_files; i++) {
        if (target_bed_info[i] != NULL)
            cleanBedInfo(target_bed_info[i]);
    }
    free(target_bed_info);

    if (Ns_bed_info != NULL)
        cleanBedInfo(Ns_bed_info);

    if (user_defined_bed_info != NULL)
        cleanBedInfo(user_defined_bed_info);

    if (stats_info)
        statsInfoDestroy(stats_info, user_inputs);

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

    if (gene_transcript_percentage_hash != NULL) {
        for (i=0; i<user_inputs->num_of_target_files; i++)
            cleanGeneTranscriptPercentage(gene_transcript_percentage_hash[i]);
        free(gene_transcript_percentage_hash);
    }

    if (raw_user_defined_databases) {
        for (i=0; i<user_inputs->num_of_target_files; i++)
            if (raw_user_defined_databases[i])
                cleanRawUserDefinedDatabase(raw_user_defined_databases[i]);
        free(raw_user_defined_databases);
    }

    // now need to clean-up the khash_t variables                                                     
    //
    if (exon_regions != NULL) {
        // here I use num_of_annotation_files to delete any exon_regions
        // associated with user defined annotation database
        // However, for the exon_regions links to the db_exon_regions from
        // the MySQL database, we only need to delete once 
        // It will be handled after this loop stated 8 lines below
        //
        for (i=0; i<user_inputs->num_of_annotation_files; i++)
            if (exon_regions[i] != NULL)
                regionsSkipMySQLDestroy(exon_regions[i], 3);
        free(exon_regions);
    }

    if (intronic_regions != NULL) regionsSkipMySQLDestroy(intronic_regions, 2);
    if (db_exon_regions != NULL) regionsSkipMySQLDestroy(db_exon_regions, 3);

    if (fn_ref) free(fn_ref);

    userInputDestroy(user_inputs);

    printf("Program Finished Successfully!\n\n");

    return 0;
}
