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

#include "data_structure.h"
#include "coverage_tracking.h"
#include "utility.h"

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

    // now for the bam/cram file open it for read for multi-threading
    //
    int t;
    samFile **sfh = calloc(user_inputs->num_of_threads, sizeof(samFile*));
    for (t=0; t<user_inputs->num_of_threads; t++) {
        sfh[t] = sam_open(user_inputs->bam_file, "r");

        if (sfh[t] == 0) {
            fprintf(stderr, "ERROR: Cannot open file \n%s\n", user_inputs->bam_file);
            return -1;
        }
    }

    // since we are going to handle one chromosome per thread, we need to get the index file
    //
    hts_idx_t **sfh_idx = calloc(user_inputs->num_of_threads, sizeof(hts_idx_t*));
    for (t=0; t<user_inputs->num_of_threads; t++) {
        sfh_idx[t] = sam_index_load(sfh[t], user_inputs->bam_file);
        if (sfh_idx[t] == NULL) {
            fprintf(stderr, "ERROR: Can't locate the index file\n");
            return -1;
        }
    }

    // Set the reference if it is the cram file
    //
    char * fn_ref = 0;
    if (user_inputs->reference_file) {
        fn_ref = getReferenceFaiPath(user_inputs->reference_file);

        if (hts_set_fai_filename(sfh[0], fn_ref) != 0) {
            fprintf(stderr, "Failed to use reference at hts_set_fai_filename() \"%s\".\n", fn_ref);
            if (fn_ref) free(fn_ref);
            return -1;
        }
    } else {
        if ( sfh[0]->is_cram || sfh[0]->format.format == cram ) {
            fprintf(stderr, "Please provide the reference sequences for the input CRAM file \n%s\n", user_inputs->bam_file);
            return -1;
        }
    }

    // check if bam file named as .cram and cram file named as .bam
    //
    checkFileExtension(user_inputs->bam_file, sfh[0]);

    // use sam_hdr_read to process both bam and cram header
    //
    bam_hdr_t **headers = calloc(user_inputs->num_of_threads, sizeof(bam_hdr_t*));
    for (t=0; t<user_inputs->num_of_threads; t++) {
        if ((headers[t] = sam_hdr_read(sfh[t])) == 0) return -1;
    }

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
        stats_info->wgs_cov_stats->total_genome_bases = 
            loadWantedChromosomes(wanted_chromosome_hash, user_inputs->database_version, user_inputs->chromosome_bed_file);
        chromosomeTrackingInit2(wanted_chromosome_hash, chrom_tracking, headers[0]);

        // here we need to verify if the chromosome naming convention matches 
        // between the bam/cram file and the chromosome bed file specified by the end user
        //
        checkNamingConvention(headers[0], wanted_chromosome_hash);

        target_buffer_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(Target_Buffer_Status));
        TargetBufferStatusInit2(target_buffer_status, wanted_chromosome_hash);
    } else {
        stats_info->wgs_cov_stats->total_genome_bases =
            loadGenomeInfoFromBamHeader(wanted_chromosome_hash, headers[0], user_inputs->database_version);
        chrom_tracking->number_of_chromosomes = headers[0]->n_targets;
        chromosomeTrackingInit1(chrom_tracking, wanted_chromosome_hash, headers[0]);

        target_buffer_status = calloc(chrom_tracking->number_of_chromosomes, sizeof(Target_Buffer_Status));
        TargetBufferStatusInit(target_buffer_status, headers[0]);
    }

    fprintf(stderr, "The total genome bases is %"PRIu64"\n", stats_info->wgs_cov_stats->total_genome_bases);

    // setup the stats_info for each individual chromosome for each threads and initialize them
    //
    Stats_Info **stats_info_per_chr = calloc(chrom_tracking->number_of_chromosomes, sizeof(Stats_Info*));
    uint32_t chrom_index = 0;
    for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index) {
        stats_info_per_chr[chrom_index] = calloc(1, sizeof(Stats_Info));
        statsInfoInit(stats_info_per_chr[chrom_index], user_inputs);
    }

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

    fflush(stdout);

    // set random seed (should only be called ONCE)
    //
    if(user_inputs->percentage < 1.0)
        srand((uint32_t)time(NULL));    // set random seed

    // first we need to handle unmapped reads and all other mapped reads that we are not intested in 
    // as these reads will be used for final stats calculation
    //
    hts_itr_t *iter = sam_itr_querys(sfh_idx[0], headers[0], "*");
    bam1_t *b = bam_init1();
    while (sam_itr_next(sfh[0], iter, b) >= 0)
        processCurrentRecord(user_inputs, b, headers[0], stats_info, chrom_tracking, 0, target_buffer_status);


    // now let's do the parallelism
    //
#pragma omp parallel shared(chrom_tracking) num_threads(user_inputs->num_of_threads)
    {
#pragma omp single
      {
        for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index)
        {
            // need to check if we need to process current chromosome
            // 
            if (wanted_chromosome_hash != NULL) {
                khiter_t iter = kh_get(khStrInt, wanted_chromosome_hash, chrom_tracking->chromosome_ids[chrom_index]);
                if (iter == kh_end(wanted_chromosome_hash)) {
                    // skip it if not needed
                    //
                    continue;
                }
            }
#pragma omp task
            {
                //int num_of_threads = omp_get_num_threads();
                int thread_id = omp_get_thread_num();
                printf("Current thread id: %d\n", thread_id);

                // get the iterator for the current chromosome
                //
                hts_itr_t *iter = sam_itr_querys(sfh_idx[thread_id], headers[thread_id], chrom_tracking->chromosome_ids[chrom_index]);
                if (iter == NULL) {
                    fprintf(stderr, "ERROR: iterator creation failed: chr %s\n", chrom_tracking->chromosome_ids[chrom_index]);
                    exit(EXIT_FAILURE);
                }

                chromosomeTrackingUpdate(chrom_tracking, chrom_tracking->chromosome_lengths[chrom_index], chrom_index);

                bam1_t *b = bam_init1();
                while (sam_itr_next(sfh[thread_id], iter, b) >= 0)
                    processCurrentRecord(user_inputs, b, headers[thread_id], stats_info_per_chr[chrom_index], chrom_tracking, chrom_index, target_buffer_status);

                bam_destroy1(b);
                hts_itr_destroy(iter);

                if (N_FILE_PROVIDED)
                    zeroAllNsRegions(chrom_tracking->chromosome_ids[chrom_index], Ns_bed_info, chrom_tracking, target_buffer_status, 0);

                printf("Thread %d is now producing coverage information for chromosome %s\n", 
                        thread_id, chrom_tracking->chromosome_ids[chrom_index]);
                        
#pragma omp critical
              {
                // The following will be produced no matter whether annotation_on is set or not
                //
                writeCoverage(chrom_tracking->chromosome_ids[chrom_index], target_bed_info, 
                        chrom_tracking, user_inputs, stats_info, intronic_regions, exon_regions);

                // now write the off target regions with high coverage into a wig file if the Write_WIG flag is set
                //
                if (TARGET_FILE_PROVIDED)
                    produceOffTargetWigFile(chrom_tracking, chrom_tracking->chromosome_ids[chrom_index], 
                            target_bed_info[chrom_index], user_inputs, stats_info, chrom_index);

                // For Whole Genome Annotations (use MySQL for the annotation)
                // if the annotation is not on, it will just output . . . . . . . )
                //
                if (user_inputs->wgs_coverage) {
                    printf("Thread %d is now writing WGS annotation for chromosome %s\n", thread_id, chrom_tracking->chromosome_ids[chrom_index]);
                    writeAnnotations(chrom_tracking->chromosome_ids[chrom_index], chrom_tracking, user_inputs, intronic_regions, db_exon_regions);

                    // if user specifies the range information (usually for graphing purpose), need to handle it here
                    //
                    printf("Thread %d is now writing coverage uniformity data for chromosome %s\n", 
                              thread_id, chrom_tracking->chromosome_ids[chrom_index]);
                    coverageRangeInfoForGraphing(chrom_tracking->chromosome_ids[chrom_index], chrom_tracking, user_inputs);
                }

                if (TARGET_FILE_PROVIDED) {
                    printf("Thread %d is now calculating gene/transcript/cds coverage percentage for chromosome %s\n", 
                              thread_id, chrom_tracking->chromosome_ids[chrom_index]);

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
                                      chrom_tracking->chromosome_ids[chrom_index], raw_user_defined_databases[x], gene_transcripts);
                            } else {
                                genePercentageCoverageInit(user_defined_cds_gene_hash, 
                                      chrom_tracking->chromosome_ids[chrom_index], dbs, gene_transcripts);
                                intersectTargetsAndRefSeqCDS(chrom_tracking->chromosome_ids[chrom_index], 
                                      target_bed_info[x], chrom_tracking, user_defined_cds_gene_hash);
                            }

                            calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[chrom_index], target_bed_info[x], 
                                    chrom_tracking, user_inputs, user_defined_cds_gene_hash);
                            transcriptPercentageCoverageInit(transcript_hash, user_defined_cds_gene_hash);
                            storeGenePercentageCoverage(chrom_tracking->chromosome_ids[chrom_index], user_inputs, 
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

                            genePercentageCoverageInit(low_cov_gene_hash, chrom_tracking->chromosome_ids[chrom_index], dbs, gene_transcripts);
                            intersectTargetsAndRefSeqCDS(chrom_tracking->chromosome_ids[chrom_index], target_bed_info[x], 
                                    chrom_tracking, low_cov_gene_hash);

                            calculateGenePercentageCoverage(chrom_tracking->chromosome_ids[chrom_index], target_bed_info[x], 
                                    chrom_tracking, user_inputs, low_cov_gene_hash);
                            transcriptPercentageCoverageInit(transcript_hash, low_cov_gene_hash);

                            storeGenePercentageCoverage(chrom_tracking->chromosome_ids[chrom_index], user_inputs, 
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

            // clean-up chrom_tracking->coverage[chrom_index], which store the base coverage in array
            //
            if (chrom_tracking->coverage[chrom_index]) {
                free(chrom_tracking->coverage[chrom_index]);
                chrom_tracking->coverage[chrom_index] = NULL;
            }
        }

#pragma omp taskwait
        // now need to combine all the stats_info for the final results
        //
#pragma omp critical
        {
            for (chrom_index=0; chrom_index<chrom_tracking->number_of_chromosomes; ++chrom_index) {
                fprintf(stderr, "\nFor chrom id %d and chrom %s\t", chrom_index, chrom_tracking->chromosome_ids[chrom_index]);
                combineCoverageStats(stats_info, stats_info_per_chr[chrom_index], user_inputs);

                // clean up the array allocated
                //
                if (stats_info_per_chr[chrom_index])
                    statsInfoDestroy(stats_info_per_chr[chrom_index], user_inputs);
            }
        }

      }     // End of single

      fflush(stdout);
    }   // end parallel loop

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

    for (t=0; t<user_inputs->num_of_threads; t++) {
        sam_close(sfh[t]);
        bam_hdr_destroy(headers[t]);
        hts_idx_destroy(sfh_idx[t]);
    }

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
