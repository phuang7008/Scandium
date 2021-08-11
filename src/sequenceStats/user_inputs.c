/*
 * =====================================================================================
 *
 *        Filename:        user_inputs.c
 *
 *        Description:    The implementation file for the user inputs
 *
 *      Version:        1.0
 *      Created:        03/16/2020 04:45:04 PM
 *      Revision:        none
 *      Compiler:        gcc
 *
 *      Author:            Peiming (Peter) Huang, phuang@bcm.edu
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>        // for file access() and getopt()
#include <dirent.h>        // for checking output directory
#include <libgen.h>        // for function basename()
#include "terms.h"
#include "utils.h"
#include "user_inputs.h"
#include "annotation.h"

// define (initialize) global variables declared at the terms.h file
//
bool N_FILE_PROVIDED = false;
bool HGMD_PROVIDED   = false;
bool TARGET_FILE_PROVIDED  = false;
bool USER_DEFINED_DATABASE = false;


void recordHGMD(Databases *dbs, khash_t(khStrInt) *hgmd_genes, khash_t(khStrInt) *hgmd_transcripts) {
    // query from the MySQL database
    //
    char *sql = calloc(100, sizeof(char));
    sprintf(sql, "SELECT gene_symbol, transcript FROM %s", dbs->db_hgmd);

    if (mysql_query(dbs->con,sql))
        finish_with_error(dbs->con);

    dbs->mysql_results = mysql_store_result(dbs->con);
    if (dbs->mysql_results == NULL)
        finish_with_error(dbs->con);

    MYSQL_ROW row;
    while ((row = mysql_fetch_row(dbs->mysql_results))) {
        int absent;

        // check to see if the current gene exists
        //
        khiter_t iter = kh_put(khStrInt, hgmd_genes, row[0], &absent);
        if (absent) {
            kh_key(hgmd_genes, iter) = strdup(row[0]);
        }
        kh_value(hgmd_genes, iter) = 1;

        // check to see if current transcripts exists
        //
        int i=0;
        char *savePtr = row[1];
        char *tokPtr;
        while ((tokPtr = strtok_r(savePtr, ".", &savePtr))) {
            if (i==0) {
                iter = kh_put(khStrInt, hgmd_transcripts, tokPtr, &absent);
                if (absent) {
                    kh_key(hgmd_transcripts, iter) = strdup(tokPtr);
                }
                kh_value(hgmd_transcripts, iter) = 1;
            }
            i++;
        }
    }

    if (dbs->mysql_results) { 
        mysql_free_result(dbs->mysql_results);
        dbs->mysql_results = NULL;
    }
    free(sql);
}

void annotationWrapperDestroy(Annotation_Wrapper *annotation_wrapper) {
    uint16_t i;

    // Free the inner most part first
    //
    for (i=0; i<annotation_wrapper->real_size; i++) {
        if (annotation_wrapper->annotations[i].gene != NULL)
            free(annotation_wrapper->annotations[i].gene);

        if (annotation_wrapper->annotations[i].Synonymous != NULL)
            free(annotation_wrapper->annotations[i].Synonymous);

        //free(annotation_wrapper->annotations[i].prev_genes);

        if (annotation_wrapper->annotations[i].exon_info != NULL)
            free(annotation_wrapper->annotations[i].exon_info);
    }

    // free the next level
    //
    if (annotation_wrapper->annotations != NULL)
        free(annotation_wrapper->annotations);

    if (annotation_wrapper != NULL)
        free(annotation_wrapper);
}

// print out the help information to show the general usage of the package
//
void usage() {
    printf("Version %s\n\n", VERSION_ );
    printf("Usage:  scandium -input_bam bam/cram -output_dir output_directory [options ...]\n");
    printf("Note:   this is a multi-threading program. Each thread needs 4Gb of memory. So please allocate them accordingly!\n");
    printf("\tfor example: 3 threads would use 12Gb of memory, while 4 threads would need 16Gb of memory, etc.\n\n");
    printf("Mandatory:\n");
    printf("--input_bam     -i  BAM/CRAM alignment file (multiple files are not allowed!).\n");
    printf("                    It Is Mandatory\n");
    printf("--output_dir    -o  output directory. It Is Mandatory\n");
    printf("--reference     -R  the file path of the reference sequence. \n");
    printf("                    It is Mandatory for CRAM files\n\n");

    printf("The Followings Are Optional:\n");
    printf("--min_base_qual -b  minimal base quality\n");
    printf("                    to filter out any bases with base quality less than b. Default 0\n");
    printf("--min_map_qual  -m  minimal mapping quality\n");
    printf("                    to filter out any reads with mapping quality less than m. Default 0\n");
    printf("--f1 --f2 ... --f8  user defined annotation file (one tag per annotation file).\n");
    printf("                    Note, the maximum number of annotatioin files allowed is 8\n");
    printf("--t1 --t2 ... --t8  capture target file. (one tag per capture file)\n");
    printf("                    Note: the maximum number of capture files allowed is 8\n");
    printf("--gvcf_block    -g  the percentage used for gVCF blocking: Default 10 for 1000%%>\n");
    printf("                    The value of gvcf_block should be larger than or equal to 1\n");
    printf("--peak_size     -k  number of points around peak (eg, Mode) area for the area \n");
    printf("                    under histogram calculation (for WGS Uniformity only)\n");
    printf("                    Dynamically Selected Based on Average Coverage of the Sample\n");
    printf("--Ns_regions    -n  file name that contains regions of Ns in the reference genome in bed format\n");
    printf("--percentage    -p  the percentage (fraction) of reads used for this analysis. Default 1.0 (ie, 100%%)\n");
    printf("--chr_list      -r  file name that contains chromosomes and their regions \n");
    printf("                    need to be processed in bed format. Default: Not Provided\n");
    printf("--buffer        -B  the Buffer size immediate adjacent to a target region. Default: 100\n");
    printf("--DB_version    -D  the version of human genome database (either hg19 [or hg37], or hg38).\n");
    printf("                    Default: hg19/hg37>\n");
    printf("--threshold_high -H the high coverage threshold/cutoff value.\n");
    printf("                    Any coverages larger than or equal to it will be outputted. Default=10000\n");
    printf("--threshold_low -L  the low coverage threshold/cutoff value.\n");
    printf("                    Any coverages smaller than it will be outputted. Default: 20\n");
    printf("--password      -P  the MySQL DB login user's Password\n");
    printf("--threads       -T  the number of threads \n");
    printf("                    (Note: when used with HPC's msub, make sure that the number of\n"); 
    printf("                    processors:ppn matches to number of threads). Default 3\n");
    printf("--username      -U  the MySQL DB login User name\n");

    printf("\nThe Followings generate block regions used for data smoothing in Coverage Uniformity Analysis\n");
    printf("--lower_bound   -l  the lower bound for the block region output. Default: 1\n");
    printf("--upper_bound   -u  the upper bound for the block region output. Default: 150\n\n");

    printf("The Followings Are Flags\n");
    printf("--annotation    -a  Specify this flag only when you want to Turn ON the annotation for WGS only.\n");
    printf("                    Default: Annotation is OFF>\n");
    printf("--capture_depth -C  To produce Capture coverage depth (cov.fasta) file that contains \n");
    printf("                    coverage count info for each targeted base. Default: OFF\n");
    printf("--duplicate     -d  Specify this flag only when you want to keep Duplicates reads.\n");
    printf("                    Default: Remove Duplicate is ON\n");
    printf("--supplemental  -s  Remove Supplementary alignments and DO NOT use them for statistics. Default: off\n");
    printf("--wgs           -w  conducting whole genome coverage analysis. Default: off\n");
    printf("--wig_output    -G  Write/Dump the WIG formatted file. Default: off\n");
    printf("--hgmd          -M  Use HGMD annotation. Default: off\n");
    printf("--overlap       -O  Remove Overlapping Bases to avoid double counting. Default: off\n");
    printf("--high_cov_out  -V  Output regions with high coverage (used with -H: default 10000). Default: off\n");
    printf("--wgs_depth     -W  Write/Dump the WGS base coverage depth into Coverage.fasta file \n");
    printf("                    (both -w and -W needed). Default: off\n");
    printf("--help          -h  Print this help/usage message\n");
}

// Get command line arguments in and check the sanity of user inputs 
//
void processUserOptions(User_Input *user_inputs, int argc, char *argv[]) {
    int arg, i;
    bool input_error_flag=false;
    bool flag_float=true;
    khash_t(khStrStr) *capture_files = kh_init(khStrStr);
    khash_t(khStrStr) *annotation_files = kh_init(khStrStr);

    // Flag set by '--verbose'
    //static int verbose_flag;

    while (1) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose",  no_argument,  &verbose_flag,  9},
            //{"brief",    no_argument,  &verbose_flag,  0},
            /* These options don't set a flag. We distinguish them by their indices. */
            {"f1",      required_argument,  0,  'e'},
            {"f2",      required_argument,  0,  'j'},
            {"f3",      required_argument,  0,  'x'},
            {"f4",      required_argument,  0,  'z'},
            {"f5",      required_argument,  0,  'E'},
            {"f6",      required_argument,  0,  'J'},
            {"f7",      required_argument,  0,  'X'},
            {"f8",      required_argument,  0,  'Z'},
            {"t1",      required_argument,  0,  '1'},
            {"t2",      required_argument,  0,  '2'},
            {"t3",      required_argument,  0,  '3'},
            {"t4",      required_argument,  0,  '4'},
            {"t5",      required_argument,  0,  '5'},
            {"t6",      required_argument,  0,  '6'},
            {"t7",      required_argument,  0,  '7'},
            {"t8",      required_argument,  0,  '8'},
            {"buffer",          required_argument,  0,  'B'},
            {"chr_list",        required_argument,  0,  'r'},
            {"DB_version",      required_argument,  0,  'D'},
            {"input_bam",       required_argument,  0,  'i'},
            {"gvcf_block",      required_argument,  0,  'g'},
            {"threshold_high",  required_argument,  0,  'H'},
            {"threshold_low",   required_argument,  0,  'L'},
            {"lower_bound",     required_argument,  0,  'l'},
            {"upper_bound",     required_argument,  0,  'u'},
            {"min_base_qual",   required_argument,  0,  'b'},
            {"min_map_qual",    required_argument,  0,  'm'},
            {"output_dir",      required_argument,  0,  'o'},
            {"reference",       required_argument,  0,  'R'},
            {"percentage",      required_argument,  0,  'p'},
            {"peak_size",       required_argument,  0,  'k'},
            {"password",        required_argument,  0,  'P'},
            {"username",        required_argument,  0,  'U'},
            {"Ns_regions",      required_argument,  0,  'n'},
            {"threads",         required_argument,  0,  'T'},
            {"annotation",          no_argument,  0,  'a'},
            {"capture_depth",       no_argument,  0,  'C'},
            {"duplicate",           no_argument,  0,  'd'},
            {"help",                no_argument,  0,  'h'},
            {"hgmd",                no_argument,  0,  'M'},
            {"overlap",             no_argument,  0,  'O'},
            {"supplemental",        no_argument,  0,  's'},
            {"high_cov_out",        no_argument,  0,  'V'},
            {"wig_output",          no_argument,  0,  'G'},
            {"wgs",                 no_argument,  0,  'w'},
            {"wgs_depth",           no_argument,  0,  'W'},
            {0,  0,  0,  0},
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        arg = getopt_long_only (argc, argv, 
                    "ab:B:CdD:f:g:GH:i:k:L:l:m:Mn:No:Op:P:r:R:st:T:u:U:VwWy:h01:2:3:4:5:6:7:8:9e:E:j:J:x:X:z:Z:", 
                    long_options, &option_index);

        /* Detect the end of the options. */
        if (arg == -1) break;

        //printf("User options for %c is %s\n", arg, optarg);
        switch(arg) {
            case '1':
                addValueToKhashBucketStrStr(capture_files, "1", optarg);
                user_inputs->num_of_target_files++;
                TARGET_FILE_PROVIDED = true;
                break;
            case '2':
                addValueToKhashBucketStrStr(capture_files, "2", optarg);
                user_inputs->num_of_target_files++;
                TARGET_FILE_PROVIDED = true;
                break;
            case '3':
                addValueToKhashBucketStrStr(capture_files, "3", optarg);
                user_inputs->num_of_target_files++;
                TARGET_FILE_PROVIDED = true;
                break;
            case '4':
                addValueToKhashBucketStrStr(capture_files, "4", optarg);
                user_inputs->num_of_target_files++;
                TARGET_FILE_PROVIDED = true;
                break;
            case '5':
                addValueToKhashBucketStrStr(capture_files, "5", optarg);
                user_inputs->num_of_target_files++;
                TARGET_FILE_PROVIDED = true;
                break;
            case '6':
                addValueToKhashBucketStrStr(capture_files, "6", optarg);
                user_inputs->num_of_target_files++;
                TARGET_FILE_PROVIDED = true;
                break;
            case '7':
                addValueToKhashBucketStrStr(capture_files, "7", optarg);
                user_inputs->num_of_target_files++;
                TARGET_FILE_PROVIDED = true;
                break;
            case '8':
                addValueToKhashBucketStrStr(capture_files, "8", optarg);
                user_inputs->num_of_target_files++;
                TARGET_FILE_PROVIDED = true;
                break;
            case 'e':
                addValueToKhashBucketStrStr(annotation_files, "1", optarg);
                user_inputs->num_of_annotation_files++;
                USER_DEFINED_DATABASE = true;
                break;
            case 'j':
                addValueToKhashBucketStrStr(annotation_files, "2", optarg);
                user_inputs->num_of_annotation_files++;
                USER_DEFINED_DATABASE = true;
                break;
            case 'x':
                addValueToKhashBucketStrStr(annotation_files, "3", optarg);
                user_inputs->num_of_annotation_files++;
                USER_DEFINED_DATABASE = true;
                break;
            case 'z':
                addValueToKhashBucketStrStr(annotation_files, "4", optarg);
                user_inputs->num_of_annotation_files++;
                USER_DEFINED_DATABASE = true;
                break;
            case 'E':
                addValueToKhashBucketStrStr(annotation_files, "5", optarg);
                user_inputs->num_of_annotation_files++;
                USER_DEFINED_DATABASE = true;
                break;
            case 'J':
                addValueToKhashBucketStrStr(annotation_files, "6", optarg);
                user_inputs->num_of_annotation_files++;
                USER_DEFINED_DATABASE = true;
                break;
            case 'X':
                addValueToKhashBucketStrStr(annotation_files, "7", optarg);
                user_inputs->num_of_annotation_files++;
                USER_DEFINED_DATABASE = true;
                break;
            case 'Z':
                addValueToKhashBucketStrStr(annotation_files, "8", optarg);
                user_inputs->num_of_annotation_files++;
                USER_DEFINED_DATABASE = true;
                break;
            case 'a':
                user_inputs->wgs_annotation_on = true; break;
            case 'b':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered base quality filter score %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->min_base_quality = atoi(optarg);
                break;
            case 'B': user_inputs->target_buffer_size = (uint16_t) strtol(optarg, NULL, 10); break;
            case 'C': user_inputs->Write_Capture_cov_fasta = true; break;
            case 'd': user_inputs->remove_duplicate = false; break;
            case 'D': 
                strcpy(user_inputs->database_version, optarg); 

                // change all to lower case
                for(i = 0; user_inputs->database_version[i]; i++){
                    user_inputs->database_version[i] = tolower(user_inputs->database_version[i]);
                }

                if (strcmp(user_inputs->database_version, "hg19") == 0)
                    strcpy(user_inputs->database_version, "hg37");

                break;
            case 'g': 
                user_inputs->gVCF_percentage = (uint16_t) strtol(optarg, NULL, 10);
                if (user_inputs->gVCF_percentage < 1) {
                    fprintf(stderr, "ERROR: The gVCF block size %s should be larger than or equal to 1\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'G': user_inputs->Write_WIG = true; break;
            case 'h': usage(); exit(EXIT_FAILURE);
            case 'H':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered High coverage cutoff value %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->high_coverage_to_report = (uint32_t) strtol(optarg, NULL, 10);
                break;
            case 'i':
                user_inputs->bam_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));

                strcpy(user_inputs->bam_file, optarg);
                break;
            case 'k':
                user_inputs->size_of_peak_area = (uint8_t) strtol(optarg, NULL, 10);
                user_inputs->user_set_peak_size_on = true; break;
            case 'L':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered Lower coverage cutoff value %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->low_coverage_to_report = (uint16_t) strtol(optarg, NULL, 10);
                break;
            case 'l':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered lower_bound value %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->lower_bound = (uint16_t) strtol(optarg, NULL, 10);
                break;
            case 'm':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered map quality filter score %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->min_map_quality = atoi(optarg);
                break;
            case 'M':
                HGMD_PROVIDED = true;
                break;
            case 'N':
                user_inputs->non_MC_tag_ON = true;
                break;
            case 'n':
                N_FILE_PROVIDED = true;
                user_inputs->n_file = malloc(strlen(optarg)+1 * sizeof(char));
                strcpy(user_inputs->n_file, optarg);
                break;
            case 'o':
                user_inputs->output_dir = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->output_dir, optarg);
                break;
            case 'O':
                user_inputs->excluding_overlapping_bases = true;
                break;
            case 'p': 
                   flag_float = isFloat(optarg, &(user_inputs->percentage)); 
                if (!flag_float) {
                    fprintf(stderr, "ERROR: Entered percentage value %s is not a float decimal number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'P':
                user_inputs->passwd = malloc(strlen(optarg)+1 * sizeof(char));
                strcpy(user_inputs->passwd, optarg);
                break;
            case 'r': 
                user_inputs->chromosome_bed_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->chromosome_bed_file, optarg);
                break;
            case 'R':
                user_inputs->reference_file = (char *) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(user_inputs->reference_file, optarg);
                break;
            case 's': user_inputs->remove_supplementary_alignments = true; break;
            case 'T':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered number of threads %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->num_of_threads = atoi(optarg);
                break;
            case 'u':
                if (!isNumber(optarg)) {
                    fprintf (stderr, "ERROR: Entered upper_bound value %s is not a number\n", optarg);
                    exit(EXIT_FAILURE);
                }
                user_inputs->upper_bound = (uint16_t) strtol(optarg, NULL, 10);
                break;
            case 'U':
                user_inputs->user_name = malloc(strlen(optarg)+1 * sizeof(char));
                strcpy(user_inputs->user_name, optarg);
                break;
            case 'V': user_inputs->above_10000_on = true; break;
            case 'w': user_inputs->wgs_coverage = true; break;
            case 'W': user_inputs->Write_WGS_cov_fasta = true; break;
            case '?':   // "ab:B:CdD:f:g:GH:i:k:L:l:m:Mn:No:Op:P:r:R:st:T:u:U:VwWy:h01:2:3:4:5:6:7:8:9e:E:j:J:x:X:z:Z:"
                if (   optopt == 'b' || optopt == 'B' || optopt == 'D' || optopt == 'g' || optopt == 'H'
                    || optopt == 'k' || optopt == 'i' || optopt == 'L' || optopt == 'l' || optopt == 'm'
                    || optopt == 'n' || optopt == 'o' || optopt == 'p' || optopt == 'P' 
                    || optopt == 'r' || optopt == 't' || optopt == 'T' || optopt == 'u' || optopt == 'U'
                    || optopt == '1' || optopt == '2' || optopt == '3' || optopt == '4' || optopt == '5' 
                    || optopt == '6' || optopt == '7' || optopt == '8' || optopt == '1' || optopt == 'e'
                    || optopt == 'E' || optopt == 'j' || optopt == 'J' || optopt == 'x' || optopt == 'X' 
                    || optopt == 'z' || optopt == 'Z' )
                    fprintf(stderr, "ERROR: Option -%c requires an argument.\n", optopt);
                //else if (optopt == 0)
                //    fprintf (stderr, "===>Unknown option `-%d'.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "ERROR: You have entered an Unknown option. See the above error message\n");
                else
                    fprintf (stderr, "ERROR: You have entered an unknown option. See the above error message\n");
                exit(EXIT_FAILURE);
                break;
            default: 
                fprintf(stderr, "ERROR: Non-option argument %c\n", optopt); input_error_flag=true; break;
                exit(EXIT_FAILURE);
        }
    }

    outputUserInputOptions(user_inputs, capture_files, annotation_files);
    checkInputCaptureAndAnnotationFiles(user_inputs);
    formTargetAnnotationFileArray(capture_files, annotation_files, user_inputs);

    // don't proceed if the user doesn't specify either -t or -w or both
    //
    if (!user_inputs->wgs_coverage && !TARGET_FILE_PROVIDED) {
        fprintf(stderr, "ERROR: You specify neither --t (for Capture Project) nor --wgs (for WGS analysis)\n");
        fprintf(stderr, "\tPlease specify either --t or --wgs or both before proceeding. Thanks!\n");
        input_error_flag=true;
    }

    // check the mandatory arguments (will turn this on for the final test/run)
    if (user_inputs->bam_file == NULL) {
        fprintf(stderr, "ERROR: --input_bam (or -i)\toption is mandatory!\n");
        input_error_flag=true;
    }

    // check database version
    if ( (strcmp(user_inputs->database_version, "hg19") != 0) && (strcmp(user_inputs->database_version, "hg37") != 0)
            && strcmp(user_inputs->database_version, "hg38") != 0) {
        fprintf(stderr, "ERROR: -D\toption is not correct! It should be either hg19 or hg37 or hg38! All in lower case, please! Thanks!\n");
        input_error_flag=true;
    }
    
    if (user_inputs->output_dir == NULL) {
        fprintf(stderr, "ERROR: --output_dir (or -o)\toption is mandatory!\n");
        input_error_flag=true;
    } else {
        // check to see if the directory exist!
        DIR* dir = opendir(user_inputs->output_dir);
        if (dir) {
            /* Directory exists */
            closedir(dir);
        } else if (ENOENT == errno) {
            fprintf(stderr, "ERROR: The output directory \n%s\n doesn't exist! \n", user_inputs->output_dir);
            fprintf(stderr, "Please double check the output directory and try again. Thanks!!\n");
            input_error_flag=true;
        } else {
            /* opendir() failed for some other reason, such as permission */
            fprintf(stderr, "ERROR: Can't open the output directory \n%s\n", user_inputs->output_dir);
            fprintf(stderr, "Please check to see if the permission is set correctly. Thanks!\n");
            input_error_flag=true;
        }
    }

    if (user_inputs->upper_bound < user_inputs->lower_bound) {
        fprintf(stderr, "ERROR: The value for --upper_bound (or -u) (%d) should be larger than the value for --lower_bound (or -l) (%d) option \n", user_inputs->upper_bound, user_inputs->lower_bound);
        input_error_flag=true;
    }

    // Need to check out that all files user provided exist before proceeding
    if (user_inputs->bam_file && !checkFile(user_inputs->bam_file)) input_error_flag=true;
    if (N_FILE_PROVIDED && !checkFile(user_inputs->n_file)) input_error_flag=true;
    if (TARGET_FILE_PROVIDED) {
        for (i=0; i<user_inputs->num_of_target_files; i++) {
            if (!checkFile(user_inputs->target_files[i])) {
                fprintf(stderr, "ERROR: Capture file \n%s\n doesn't exist\n", user_inputs->target_files[i]);
                input_error_flag=true;
            }
        }
    }

    for (i=0; i<user_inputs->num_of_annotation_files; i++) {
        if (!checkFile(user_inputs->user_defined_annotation_files[i])) {
            fprintf(stderr, "ERROR: Annotation file \n%s\n doesn't exist\n", user_inputs->user_defined_annotation_files[i]);
            input_error_flag=true;
        }
    }

    if (input_error_flag) {
        //usage();
        fprintf(stderr, "Please use --help (or -h) for all Scandium options\n");
        exit(EXIT_FAILURE);
    }

    cleanKhashStrStr(annotation_files);
    cleanKhashStrStr(capture_files);
}

void setupOutputReportFiles(User_Input *user_inputs) {
    // need to get the basename from BAM/CRAM filename
    char *tmp_basename = basename(user_inputs->bam_file);
    if (!tmp_basename || strlen(tmp_basename) == 0) {
        fprintf(stderr, "ERROR: Something went wrong for extracting the basename from the input BAM/CRAM file\n");
        fprintf(stderr, "Please use --help (or -h) for all Scandium options\n");
        exit(EXIT_FAILURE);
    }

    char string_to_add[350];

    // For all Capture related output files
    //
    int p = 0;
    if (TARGET_FILE_PROVIDED) {
        checkRepeatedCaptureFiles(user_inputs);

        // produce the base filename information
        //
        getBaseFilenameWithoutExtension(user_inputs, 1);

        user_inputs->capture_all_site_files = calloc(user_inputs->num_of_target_files, sizeof(char*));
        for (p=0; p<user_inputs->num_of_target_files; p++) {
            sprintf(string_to_add, ".%s.Capture_AllSites_REPORT.txt", user_inputs->target_file_basenames[p]);
            createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_all_site_files[p], string_to_add, VERSION_);
            writeHeaderLine(user_inputs->capture_all_site_files[p], user_inputs, p, 1);
        }

        // For capture coverage summary report
        //
        user_inputs->capture_cov_reports = calloc(user_inputs->num_of_target_files, sizeof(char*));
        for (p=0; p<user_inputs->num_of_target_files; p++) {
            sprintf(string_to_add, ".%s.Capture_Coverage_Summary_Report.txt", user_inputs->target_file_basenames[p]);
            createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_cov_reports[p], string_to_add, VERSION_);
        }

        // for cov.fasta file name
        //
        if (user_inputs->Write_Capture_cov_fasta) {
            user_inputs->capture_cov_bedfiles = calloc(user_inputs->num_of_target_files, sizeof(char*));
            for (p=0; p<user_inputs->num_of_target_files; p++) {
                sprintf(string_to_add, ".%s.Capture_cov.bed", user_inputs->target_file_basenames[p]);
                createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_cov_bedfiles[p], string_to_add, VERSION_);
            }

            user_inputs->capture_cov_files = calloc(user_inputs->num_of_target_files, sizeof(char*));
            for (p=0; p<user_inputs->num_of_target_files; p++) {
                sprintf(string_to_add, ".%s.Capture_cov.fasta", user_inputs->target_file_basenames[p]);
                createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_cov_files[p], string_to_add, VERSION_);
            }
        }

        // for off target good hit wig.fasta file name
        //
        if (user_inputs->Write_WIG) {
            user_inputs->capture_wig_files = calloc(user_inputs->num_of_target_files, sizeof(char*));
            for (p=0; p<user_inputs->num_of_target_files; p++) {
                sprintf(string_to_add, ".%s.Capture_off_target_good_hits.wig.fasta", user_inputs->target_file_basenames[p]);
                createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_wig_files[p], string_to_add, VERSION_);
            }
        }

        // output low coverage regions for target (capture)
        //
        user_inputs->capture_low_cov_files = calloc(user_inputs->num_of_target_files, sizeof(char*));
        for (p=0; p<user_inputs->num_of_target_files; p++) {
            sprintf(string_to_add, ".%s.Capture_below%dx_REPORT.txt", user_inputs->target_file_basenames[p], user_inputs->low_coverage_to_report);
            createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_low_cov_files[p], string_to_add, VERSION_);
            writeHeaderLine(user_inputs->capture_low_cov_files[p], user_inputs, p, 1);
        }

        // output too high coverage regions for target (capture)
        //
        if (user_inputs->above_10000_on) {
            user_inputs->capture_high_cov_files = calloc(user_inputs->num_of_target_files, sizeof(char*));
            for (p=0; p<user_inputs->num_of_target_files; p++) {
                sprintf(string_to_add, ".%s.Capture_above%dx_REPORT.txt", user_inputs->target_file_basenames[p], user_inputs->high_coverage_to_report);
                createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->capture_high_cov_files[p], string_to_add, VERSION_);
                writeHeaderLine(user_inputs->capture_high_cov_files[p], user_inputs, p, 1);
            }
        }

        // for low coverage gene/exon/transcript reports
        // allocate the memory
        //
        user_inputs->low_cov_gene_pct_files = calloc(user_inputs->num_of_target_files, sizeof(char*));
        user_inputs->low_cov_exon_pct_files = calloc(user_inputs->num_of_target_files, sizeof(char*));
        user_inputs->low_cov_transcript_files = calloc(user_inputs->num_of_target_files, sizeof(char*));

        for (p=0; p<user_inputs->num_of_target_files; p++) {
            sprintf(string_to_add, ".%s.Capture_Gene_pct.txt", user_inputs->target_file_basenames[p]);
            createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_gene_pct_files[p], string_to_add, VERSION_);
            writeHeaderLine(user_inputs->low_cov_gene_pct_files[p], user_inputs, p, 3);

            sprintf(string_to_add, ".%s.Capture_CDS_pct.txt", user_inputs->target_file_basenames[p]);
            createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_exon_pct_files[p], string_to_add, VERSION_);
            writeHeaderLine(user_inputs->low_cov_exon_pct_files[p], user_inputs, p, 4);

            sprintf(string_to_add, ".%s.Capture_Transcript_pct.txt", user_inputs->target_file_basenames[p]);
            createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->low_cov_transcript_files[p], string_to_add, VERSION_);
            writeHeaderLine(user_inputs->low_cov_transcript_files[p], user_inputs, p, 5);
        }
    }

    // For whole Genome report
    //
    if (user_inputs->wgs_coverage) {
        // output WGS coverage summary report
        sprintf(string_to_add, ".WGS_Coverage_Summary_Report.txt");
        createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_cov_report, string_to_add, VERSION_);

        // output low coverage regions for WGS
        sprintf(string_to_add, ".WGS_below%dx_REPORT.txt", user_inputs->low_coverage_to_report);
        createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_low_cov_file, string_to_add, VERSION_);
        writeHeaderLine(user_inputs->wgs_low_cov_file, user_inputs, 20, 1);

        // output too high coverage regions for the whole genome
        //
        if (user_inputs->above_10000_on) {
            sprintf(string_to_add, ".WGS_above%dx_REPORT.txt", user_inputs->high_coverage_to_report);
            createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_high_cov_file, string_to_add, VERSION_);
            writeHeaderLine(user_inputs->wgs_high_cov_file, user_inputs, 20, 1);
        }

        // output the uniformity data file for Uniformity Analysis
        //
        sprintf(string_to_add, ".WGS_uniformity_REPORT.txt");
        createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_uniformity_file, string_to_add, VERSION_);

        // for whole genome (wgs) file name
        if (user_inputs->Write_WGS_cov_fasta) {
            createFileName(user_inputs->output_dir, tmp_basename, &user_inputs->wgs_cov_file, ".WGS_cov.fasta", VERSION_);
            //printf("Create wgs file name %s\n", user_inputs->wgs_file);
        }
    }

    // KEEP the following Please!
    // free the memory
    // However, here is the quote from the basename() official site
    // Both dirname() and basename() return pointers to null-terminated strings. (Do not pass these pointers to free(3))
    //if (tmp_basename) {
    //    free(tmp_basename);
    //    tmp_basename=NULL;
    //}

    // string_to_add is declared at the stack, so no need to free it!
}

void checkInputCaptureAndAnnotationFiles(User_Input *user_inputs) {
    if (user_inputs->num_of_annotation_files > user_inputs->num_of_target_files) {
        fprintf(stderr, "ERROR: You have entered more annotation files than capture files\n");
        fprintf(stderr, "Please ensure the annotation file number < the capture file number!\n");
        exit(EXIT_FAILURE);
    }
}

// for handling multiple user defined database annotation files and multiple capture target bed files
// type 1: target files     2: annotation files
//
void formTargetAnnotationFileArray(khash_t(khStrStr) *capture_files, khash_t(khStrStr) *annotation_files, User_Input *user_inputs) {
    user_inputs->target_files = calloc(user_inputs->num_of_target_files, sizeof(char*));
    user_inputs->user_defined_annotation_files = calloc(user_inputs->num_of_annotation_files, sizeof(char*));

    const char *keys[8] = {"1", "2", "3", "4", "5", "6", "7", "8"};
    char **capture_only_files = calloc(8, sizeof(char*));
    int i;
    uint8_t counter=0;
    uint8_t capture_only_counter=0;

    for (i=0; i<8; i++) {
        khiter_t a_iter = kh_get(khStrStr, annotation_files, keys[i]);
        khiter_t c_iter = kh_get(khStrStr, capture_files, keys[i]);
        
        // no keys for both, just skip
        //
        if (a_iter == kh_end(annotation_files) && c_iter == kh_end(capture_files))
            continue;

        // annotation only without capture file
        //
        if (a_iter != kh_end(annotation_files) && c_iter == kh_end(capture_files)) {
            fprintf(stderr, "ERROR: you have entered an annotation file \n%s\n", kh_value(annotation_files, a_iter));
            fprintf(stderr, "without the corresponding capture file\n");
            fprintf(stderr, "Please check your inputs and try again. Thanks!\n");
            exit(EXIT_FAILURE);
        }

        // capture only without annotation file
        //
        if (a_iter == kh_end(annotation_files) && c_iter != kh_end(capture_files)) {
            capture_only_files[capture_only_counter] = 
                (char*) malloc(strlen(kh_value(capture_files, c_iter))+1 * sizeof(char));
            strcpy(capture_only_files[capture_only_counter], kh_value(capture_files, c_iter));
            capture_only_counter++;
        }

        // both have keys
        //
        if (a_iter != kh_end(annotation_files) && c_iter != kh_end(capture_files)) {
            user_inputs->target_files[counter] = 
                (char*) malloc(strlen(kh_value(capture_files, c_iter))+1 * sizeof(char));
            strcpy(user_inputs->target_files[counter], kh_value(capture_files, c_iter));

            user_inputs->user_defined_annotation_files[counter] = 
                (char*) malloc(strlen(kh_value(annotation_files, a_iter))+1 * sizeof(char));
            strcpy(user_inputs->user_defined_annotation_files[counter], kh_value(annotation_files, a_iter));

            counter++;
        }
    }

    for (i=0; i<capture_only_counter; i++) {
        user_inputs->target_files[counter] = (char*) malloc(strlen(capture_only_files[i])+1 * sizeof(char));
        strcpy(user_inputs->target_files[counter], capture_only_files[i]);
        counter++;
    }

    for (i=0; i<capture_only_counter; i++) {
        free(capture_only_files[i]);
    }
    free(capture_only_files);
}

void checkRepeatedCaptureFiles(User_Input *user_inputs) {
    int i, j;
    for (i=0; i< user_inputs->num_of_target_files; i++) {
        for (j=1; j< user_inputs->num_of_target_files; j++) {
            if (i == j)
                continue;

            if (strcmp (user_inputs->target_files[i], user_inputs->target_files[j]) == 0) {
                // same capture file names repeated twice! Error out
                //
                fprintf(stderr, "ERROR: You entered same capture file \n%s\ntwice\n", user_inputs->target_files[i]);
                fprintf(stderr, "Capture file name should be UNIQUE.\n");
                fprintf(stderr, "Please use a different file name instead if this is your intention\n");
                exit(EXIT_FAILURE);
            }
        }
    }
}

// need to write the header line for some of the output files
// type:    1 -> coverage threshold reports;    2 -> capture missing target report; 
//          3 -> Gene Percentage Reports;       4 -> Exon Percentage Reports
//          5 -> Transcript Percentage Reports;
//
void writeHeaderLine(char *file_in, User_Input *user_inputs, int annotation_file_index, uint8_t type) {
    FILE *out_fp = fopen(file_in, "a");

    if (annotation_file_index == 20) {
        // For WGS
        //
        if (user_inputs->wgs_annotation_on) {
            fprintf(out_fp, "##Annotation source: MySQL database\n");
        } else {
            fprintf(out_fp, "##Annotation source: None\n");
        }
    } else {
        // for Capture
        //
        if (user_inputs->num_of_annotation_files <= annotation_file_index) {
            fprintf(out_fp, "##Annotation source: MySQL database\n");
        } else {
            fprintf(out_fp, "##Annotation source: %s\n", user_inputs->user_defined_annotation_files[annotation_file_index]);
        }
    }

    if (type == 1) {
        // for the coverage annotation report (for example: below20x, above10000x coverage reports)
        fprintf(out_fp, "##This file will be produced when the user specifies a low coverage threshold and need detailed annotations for these low coverage regions\n");
        fprintf(out_fp, "##%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Chr", "Start", "End", "Length", "Coverage", "Gene_Symbol", "Synonyms", "Prev_Gene_Symbol", "RefSeq", "CCDS", "VEGA", "miRNA", "Others (SNP, Pseudo-Gene etc.)");
    } else if (type == 2) {
        // for capture missed target file
        fprintf(out_fp, "##This file will be produced when the user specifies the target file (For Capture only)\n");
        fprintf(out_fp, "##It contains the target regions that do not get covered\n");
        fprintf(out_fp, "##%s\t%s\t%s\n", "Chr", "Start", "End");
    } else if (type == 3) {
        // for gene percentage coverage annotation reports
        fprintf(out_fp, "##This file will be produced when the user provides a target bed file (For Capture only)\n");
        fprintf(out_fp, "##NOTE: (Average_)Percentage_Of_Coverage means percentage of bases with coverage >= the user-specified coverage threshold\n");
        fprintf(out_fp, "##%s\t%s\t%s\t%s\n", "Gene_Symbol", "RefSeq_List(Percentage_Of_Coverage)", "Average_Percentage_Of_Coverage", "If_-M_on,_HGMD?");
    } else if (type == 4) {
        // for exon percentage coverage annotation reports
        fprintf(out_fp, "##This file will be produced when the user provides a target bed file (For Capture only)\n");
        fprintf(out_fp, "##NOTE: Percentage_Of_Coverage means percentage of bases with coverage >= the user-specified coverage threshold\n");
        fprintf(out_fp, "##NOTE: Regions_With_Low_Coverage means regions with base coverage below the user-specified coverage threshold\n");
        fprintf(out_fp, "##NOTE: CDS_ID refers to the corresponding coding Exon_ID that users are interested in\n");
        fprintf(out_fp, "##%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Chr", "Gene_Symbol", "RefSeq", "CDS_ID", "Start", "End", "Percentage_Of_Coverage", "Regions_With_Low_Coverage");
    } else if (type == 5) {
        // for transcript percentage coverage reports
        fprintf(out_fp, "##This file will be produced when the user provides a target bed file (For Capture only)\n");
        fprintf(out_fp, "##NOTE: Percentage_Of_Coverage means percentage of bases with coverage >= the user-specified coverage threshold\n");
        fprintf(out_fp, "##NOTE: CDS_Count refers to the corresponding coding exon count that users are interested in\n");
        fprintf(out_fp, "##%s\t%s\t%s\t%s\t%s\t%s\n", "Chr", "Gene_Symbol", "RefSeq", "Length", "CDS_Count", "Percentage_Of_Coverage");
    }

    fclose(out_fp);
}

void outputUserInputOptions(User_Input *user_inputs, khash_t(khStrStr) *capture_files, khash_t(khStrStr) *annotation_files) {
    fprintf(stderr, "The following are the options you have chosen:\n");
    fprintf(stderr, "\tInput bam/cram file: %s\n", user_inputs->bam_file);
    fprintf(stderr, "\tOutput directory is: %s\n", user_inputs->output_dir);
    fprintf(stderr, "\tThe version of official gene annotation is: %s\n", user_inputs->database_version);

    if (user_inputs->n_file)
        fprintf(stderr, "\tThe file that contains all Ns regions is: %s\n", user_inputs->n_file);

    if (user_inputs->chromosome_bed_file)
        fprintf(stderr, "\tThe file that contains the chromosome IDs and regions to be processed: %s\n", user_inputs->chromosome_bed_file);

    if (user_inputs->reference_file)
        fprintf(stderr, "\tThe reference sequence file is: %s\n", user_inputs->reference_file);

    const char *keys[8] = {"1", "2", "3", "4", "5", "6", "7", "8"};
    uint8_t i;
    fprintf(stderr, "\n\tCapture Input Files:\n");
    for (i=0; i<8; i++) {
        khiter_t c_iter = kh_get(khStrStr, capture_files, keys[i]);

        if (c_iter != kh_end(capture_files))
            fprintf(stderr, "\t--t%d capture/target bed file is: %s\n", i+1, kh_value(capture_files, c_iter));
    }
    fprintf(stderr, "\n\tUser Defined Annotation Input Files:\n");

    if (USER_DEFINED_DATABASE) {
        for (i=0; i<8; i++) {
            khiter_t a_iter = kh_get(khStrStr, annotation_files, keys[i]);
            if (a_iter != kh_end(annotation_files))
                fprintf(stderr, "\t--f%d for user provided database is %s.\n", i+1, kh_value(annotation_files, a_iter));
        }
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "\tThe minimum mapping quality is: %d\n", user_inputs->min_map_quality);
    fprintf(stderr, "\tThe minimum base quality is: %d\n", user_inputs->min_base_quality);
    fprintf(stderr, "\tThe number of thread used is: %d\n", user_inputs->num_of_threads);
    fprintf(stderr, "\tThe percentage of reads used for analysis is: %.1f%%\n", user_inputs->percentage*100);
    fprintf(stderr, "\tThe coverage number used for low coverage report is:  < %d\n", user_inputs->low_coverage_to_report);
    fprintf(stderr, "\tThe coverage number used for high coverage report is: > %d\n", user_inputs->high_coverage_to_report);
    fprintf(stderr, "\tThe percentage used for gVCF block grouping is %d%%\n", user_inputs->gVCF_percentage * 100);
    fprintf(stderr, "\tThe buffer size around a target region is %d\n", user_inputs->target_buffer_size);

    fprintf(stderr, "\tThe uniformity data file will be produced\n");
    fprintf(stderr, "\t\tThe uniformity lower bound and upper bound are %d and %d inclusive! \n", user_inputs->lower_bound, user_inputs->upper_bound);

    if (user_inputs->wgs_annotation_on) {
        fprintf(stderr, "\tThe detailed WGS gene annotation is ON\n");
    } else {
        fprintf(stderr, "\tThe detailed WGS gene annotation is OFF\n");
    }

    if (user_inputs->Write_Capture_cov_fasta) {
        fprintf(stderr, "\tThe Capture coverage dump at base level is ON (this will create a cov.fasta file\n");
    } else {
        fprintf(stderr, "\tThe Capture coverage dump is OFF\n");
    }

    if (user_inputs->Write_WGS_cov_fasta) {
        fprintf(stderr, "\tThe WGS coverage dump at base level is ON (this will create a cov.fasta file\n");
    } else {
        fprintf(stderr, "\tThe WGS coverage dump is OFF\n");
    }

    if (user_inputs->Write_WIG) {
        fprintf(stderr, "\tThe WIG file creation is ON\n");
    } else {
        fprintf(stderr, "\tThe WIG file creation is OFF\n");
    }

    if (user_inputs->above_10000_on) {
        fprintf(stderr, "\tThe above_%dx file creation is ON\n", user_inputs->high_coverage_to_report);
    } else {
        fprintf(stderr, "\tThe above_10000x file creation is OFF\n");
    }

    if (user_inputs->wgs_coverage) {
        fprintf(stderr, "\tThe whole genome analysis is ON\n");
    } else {
        fprintf(stderr, "\tThe whole genome analysis is OFF\n");
    }

    if (user_inputs->remove_duplicate) {
        fprintf(stderr, "\tRemove duplicate reads is ON\n");
    } else {
        fprintf(stderr, "\tRemove duplicate reads is OFF\n");
    }

    if (user_inputs->excluding_overlapping_bases) {
        fprintf(stderr, "\tExcluding Overlapping bases is ON\n");
    } else {
        fprintf(stderr, "\tExcluding Overlapping bases is Off\n");
    }

    if (user_inputs->remove_supplementary_alignments) {
        fprintf(stderr, "\tRemove supplementary alignments is ON\n");
    } else {
        fprintf(stderr, "\tRemove supplementary alignments is OFF\n");
    }

    /*if (user_inputs->primary_chromosomes_only) {
        fprintf(stderr, "\tUse primary chromosomes for the Uniformity calculation\n");
    } else {
        fprintf(stderr, "\tUse All chromosomes for the Uniformity calculation\n");
    }*/

    fprintf(stderr, "User Input Options ===> DONE!\n\n");
    //printf("The  is: %d\n", user_inputs->);

}

User_Input * userInputInit() {
    User_Input * user_inputs = calloc(1, sizeof(User_Input));
    if (!user_inputs) {
        fprintf(stderr, "ERROR: Memory allocation failed in line %d!\n", __LINE__);
        exit(EXIT_FAILURE);
    }

    user_inputs->min_map_quality  = 0;
    user_inputs->min_base_quality = 0;
    user_inputs->low_coverage_to_report = 20;
    user_inputs->high_coverage_to_report = 10000;
    user_inputs->lower_bound = 1;
    user_inputs->upper_bound = 150;
    user_inputs->target_buffer_size = 100;
    user_inputs->gVCF_percentage = 10;
    user_inputs->num_of_threads   = 3;
    user_inputs->percentage = 1.0;
    user_inputs->size_of_peak_area = 0;
    user_inputs->num_of_target_files = 0;
    user_inputs->num_of_annotation_files = 0;
    user_inputs->user_set_peak_size_on = false;
    user_inputs->wgs_annotation_on = false;
    user_inputs->above_10000_on = false;
    user_inputs->wgs_coverage = false;
    user_inputs->Write_Capture_cov_fasta = false;
    user_inputs->Write_WGS_cov_fasta = false;
    user_inputs->Write_WIG = false;
    user_inputs->excluding_overlapping_bases = false;
    user_inputs->remove_duplicate = true;
    user_inputs->remove_supplementary_alignments = false;

    user_inputs->n_file = NULL;
    user_inputs->bam_file = NULL;
    user_inputs->output_dir = NULL;
    user_inputs->target_files = NULL;
    user_inputs->reference_file = NULL;
    user_inputs->chromosome_bed_file = NULL;
    user_inputs->annotation_file_basenames = NULL;
    user_inputs->user_defined_annotation_files = NULL;

    // WGS output file
    //
    user_inputs->wgs_cov_file  = NULL;
    user_inputs->wgs_cov_report = NULL;
    user_inputs->wgs_low_cov_file = NULL;
    user_inputs->wgs_high_cov_file = NULL;
    user_inputs->wgs_uniformity_file = NULL;

    // Capture output files
    //
    user_inputs->capture_wig_files = NULL;
    user_inputs->capture_cov_files = NULL;
    user_inputs->capture_cov_bedfiles = NULL;
    user_inputs->capture_cov_reports  = NULL;
    user_inputs->capture_low_cov_files = NULL;
    user_inputs->capture_high_cov_files = NULL;
    user_inputs->capture_all_site_files = NULL;
    user_inputs->low_cov_gene_pct_files = NULL;
    user_inputs->low_cov_exon_pct_files = NULL;
    user_inputs->low_cov_transcript_files = NULL;

    user_inputs->database_version = calloc(10, sizeof(char));
    strcpy(user_inputs->database_version, "hg37");

    user_inputs->non_MC_tag_ON = false;
    
    return user_inputs;
}

void userInputDestroy(User_Input *user_inputs) {

    if (user_inputs->database_version)
        free(user_inputs->database_version);

    if (user_inputs->passwd)
        free(user_inputs->passwd);

    if (user_inputs->user_name)
        free(user_inputs->user_name);

    if (user_inputs->bam_file)
        free(user_inputs->bam_file);

    if (user_inputs->output_dir)
        free(user_inputs->output_dir);

    if (user_inputs->n_file)
        free(user_inputs->n_file);

    // whole genome output files clean-up
    //
    if (user_inputs->wgs_cov_file)
        free(user_inputs->wgs_cov_file);

    if (user_inputs->wgs_cov_report)
        free(user_inputs->wgs_cov_report);

    if (user_inputs->wgs_low_cov_file)
        free(user_inputs->wgs_low_cov_file);

    if (user_inputs->wgs_high_cov_file)
        free(user_inputs->wgs_high_cov_file);

    if (user_inputs->wgs_uniformity_file)
        free(user_inputs->wgs_uniformity_file);

    // Capture (target) output files clean-up
    // cleanCreatedFileArray(user_inputs->num_of_target_files, );
    //
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->target_files);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->target_file_basenames);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->capture_cov_files);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->capture_cov_bedfiles);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->capture_wig_files);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->capture_cov_reports);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->capture_all_site_files);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->capture_low_cov_files);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->capture_high_cov_files);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->low_cov_gene_pct_files);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->low_cov_exon_pct_files);
    cleanCreatedFileArray(user_inputs->num_of_target_files, user_inputs->low_cov_transcript_files);
    cleanCreatedFileArray(user_inputs->num_of_annotation_files, user_inputs->annotation_file_basenames);

    if (user_inputs->chromosome_bed_file)
        free(user_inputs->chromosome_bed_file);

    if (user_inputs->reference_file)
        free(user_inputs->reference_file);

    // For User-Defined Database output files clean-up
    //
    cleanCreatedFileArray(user_inputs->num_of_annotation_files, user_inputs->user_defined_annotation_files);

    if (user_inputs)
        free(user_inputs);
}

void createdFileArray(User_Input *user_inputs, uint8_t f_size, char **file_array, char* description, char* bam_prefix, uint8_t headline_type) {
    if (file_array == NULL) return;

    uint8_t p;
    for (p=0; p<f_size; p++) {
        char* string_to_add = calloc(strlen(user_inputs->target_file_basenames[p]) + strlen(description) + 2, sizeof(char));
        sprintf(string_to_add, ".%s.%s", user_inputs->target_file_basenames[p], description);
        createFileName(user_inputs->output_dir, bam_prefix, &file_array[p], string_to_add, VERSION_);
        writeHeaderLine(file_array[p], user_inputs, p, headline_type);
    }
}

void cleanCreatedFileArray(uint8_t f_size, char **file_array) {
    if (file_array != NULL) {
        uint8_t i;
        for (i=0; i<f_size; i++) {
            if (file_array[i])
                free(file_array[i]);
        }
        free(file_array);
    }
}

void setupMySQLDB(Databases **dbs, User_Input *user_inputs) {
    if (*dbs == NULL) {
        *dbs = calloc(1, sizeof(Databases));
        databaseSetup(*dbs, user_inputs);
    }
}

// Type: 1 for target; 2 for annotation
//
void getBaseFilenameWithoutExtension(User_Input *user_inputs, uint8_t type) {
    // allocate memory
    //
    if (type == 1) {
        user_inputs->target_file_basenames = calloc(user_inputs->num_of_target_files, sizeof(char*));
    } else {
        user_inputs->annotation_file_basenames = calloc(user_inputs->num_of_annotation_files, sizeof(char*));
    }

    uint8_t p, num_of_items;
    type == 1 ? num_of_items = user_inputs->num_of_target_files : user_inputs->num_of_annotation_files;
    for (p=0; p<num_of_items; p++) {
        const char *bn = (type == 1) ? baseFilename(user_inputs->target_files[p]) 
            : baseFilename(user_inputs->user_defined_annotation_files[p]);
        const char* ext = (type == 1) ? getFileExtension(user_inputs->target_files[p]) 
            : getFileExtension(user_inputs->user_defined_annotation_files[p]);

        int ext_len = strlen(ext);
        int bn_len  = strlen(bn);

        int len = bn_len - ext_len;
        char substring[len];
        memcpy(substring, &bn[0], len);
        substring[len] = '\0';

        if (type == 1) { 
            user_inputs->target_file_basenames[p] = calloc(len+1, sizeof(char));
            strcpy(user_inputs->target_file_basenames[p], substring);
        } else {
            user_inputs->annotation_file_basenames[p] = calloc(len+1, sizeof(char));
            strcpy(user_inputs->annotation_file_basenames[p], substring);
        }

        if (bn) free((char*)bn);
        if (ext) free((char*)ext);
    }
}
