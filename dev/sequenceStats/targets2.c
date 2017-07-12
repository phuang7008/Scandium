/*
 * =====================================================================================
 *
 *		Filename:		targets.c
 *
 *		Description:	For the base coverage calculation
 *
 *      Version:		1.0
 *      Created:		02/06/2017 04:45:04 PM
 *      Revision:		none
 *      Compiler:		gcc
 *
 *      Author:			Peiming (Peter) Huang
 *      Company:		Baylor College of Medicine
 *
 * =====================================================================================
 */

#include "targets.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>		// for file access()
#include <time.h>

// It will open the bed-formatted file for the first time. and count the number of lines (items) within it.
//
uint32_t getLineCount(char *bed_file) {
    FILE *fp = fopen(bed_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("target file %s open failed!", bed_file);
        exit(1);
    }

	// Extract characters from file and store in character c
    char c;         // To store a character read from file
    long count = 0; // To store total line number in the target file

    for (c = getc(fp); c != EOF; c = getc(fp)) {
        if (c == '\n')      // Increment count if this character is newline
            count = count + 1;
    }

    fclose(fp);     // Close the file

    return count;
}

// It will open the bed-formatted file and then record all the start and stop positions along with chromosome ids
// for chromosome X and Y are characters, I will use array of chars (ie array of string) to store chromosome ids
//
void loadBedFiles(char * bed_file, Bed_Coords * coords) {
	FILE *fp = fopen(bed_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("target file %s open failed!", bed_file);
        exit(1);
    }

	int num = 0;		// to keep track the tokens from string split using strtok()
	int count = 0;		// index for each target line (item)
	ssize_t read;
	size_t len = 0;
	char *p_token, *line;	// here p_ means a point. so p_token is a point to a token

    while((read = getline(&line, &len, fp)) != -1) {
		//printf("%s\n", line);
		if (line[0] == '\0') continue;		// skip if it is a blank line
		if (line[0] == '#')  continue;		// skip if it is a comment line

		p_token = strtok(line, " \t");
		while (p_token != NULL) {
			if (num == 0) 
				strcpy(coords[count].chr,  p_token);
			else if (num == 1)
				coords[count].start = atoi(p_token);
			else if (num == 2) {
				coords[count].end = atoi(p_token);
				num = 0;
				p_token = NULL;
				break;		// stop 'while' loop as we don't need anything after stop position
			}

			p_token = strtok(NULL, " \t");	// In subsequent calls, strtok expects a null pointer
			num++;
		}
		count++;
	}

	fclose(fp);
	if (line != NULL)
		free(line);
	if (p_token != NULL)
		free(p_token);
}

void processBedFiles(char *bed_file, Bed_Info *bed_info, khash_t(str) *bed_buffer_hash, short type) {
    // load target file again and store the target information (starts, stops and chromosome ids)
    loadBedFiles(bed_file, bed_info->coords);

    // Output stored coordinates for verification
	// int i;
    /*for (i=0; i<bed_info->size; i++) {
        printf("%s\t%ld\t%ld\n", bed_info->coords[i].chr, bed_info->coords[i].start, bed_info->coords[i].end);
    }*/

    // Now we are going to generate target-buffer lookup table for all the loaded targets
    // we will store targets and buffers information based on chromosome ID
    // we will have 22 + X + Y = 24 chromosomes, Here X=23 and Y will be 24
    generateBedBufferLookupTable(bed_info, bed_buffer_hash, type);

	// Check to see if target and buffer is set correctly
	/*khiter_t outer_iter, inner_iter;
	
	for (outer_iter=kh_begin(bed_buffer_hash); outer_iter!=kh_end(bed_buffer_hash); outer_iter++) {
		if (kh_exist(bed_buffer_hash, outer_iter)) {
			// find the key
			printf("Current key is %s\n", kh_key(bed_buffer_hash, outer_iter));
			if (strcmp("1", kh_key(bed_buffer_hash, outer_iter)) == 0) {
				for (inner_iter=kh_begin(kh_value(bed_buffer_hash, outer_iter));
							inner_iter!=kh_end(kh_value(bed_buffer_hash, outer_iter));
								inner_iter++) {
					uint32_t pos = kh_key(kh_value(bed_buffer_hash, outer_iter), inner_iter);
					if ( (pos >= 1219000) && (pos<=1220290) )  {		// for region: 1220086 1220186
						printf("%"PRIu32"\t%c\n",pos, kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter));
					}
				}
			}
		}
	}*/
}

void generateBedBufferLookupTable(Bed_Info * bed_info, khash_t(str) *bed_buffer_hash, short type) {
	int i, j, ret, absent;
	khiter_t outer_iter, inner_iter;
	char *cur_chrom_id  = malloc(15 * sizeof(char));
	strcpy(cur_chrom_id,  "something");

	for (i = 0; i < bed_info->size; i++) {
		//printf("Processing id %d\n", i);

		// check to see if we need to initialize the hash for current chromosome
		if (strcmp(bed_info->coords[i].chr, cur_chrom_id) != 0) {
			khash_t(m32) *tmp_hash = kh_init(m32);
			outer_iter = kh_put(str, bed_buffer_hash, bed_info->coords[i].chr, &absent);
			if (absent)	// insert the key if there is no key
				kh_key(bed_buffer_hash, outer_iter) = strdup(bed_info->coords[i].chr);
			kh_value(bed_buffer_hash, outer_iter) = tmp_hash;

			strcpy(cur_chrom_id,  bed_info->coords[i].chr);
			//printf("hash is added for chromosome: %s\n", coords[i].chr);
		}

		// now get the hash table for current chromosome id
		outer_iter = kh_put(str, bed_buffer_hash, bed_info->coords[i].chr, &ret);
		if (ret)
			fprintf(stderr, "Something went wrong, has key for %s is not added", bed_info->coords[i].chr);

		// type == 1 for target, type == 2 for Ns regions in the reference sequence
		char c_type = type == 1 ? 'T' : 'N';

		// for positions on targets or Ns
		for (j=bed_info->coords[i].start; j<=bed_info->coords[i].end; j++) {
			inner_iter = kh_put(m32, kh_value(bed_buffer_hash, outer_iter), j, &ret);
			kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) = c_type;
		}

		if (type == 1) {
			// for positions on the buffer at the left side
			for (j=bed_info->coords[i].start-BUFFER; j < bed_info->coords[i].start; j++) {
				if (j < 0) continue;
				inner_iter = kh_put(m32, kh_value(bed_buffer_hash, outer_iter), j, &ret);
				kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) = 'B';
			}

			// for positions on the buffer at the right side
			for (j=bed_info->coords[i].end+1; j <= bed_info->coords[i].end+BUFFER; j++ ) {
				inner_iter = kh_put(m32, kh_value(bed_buffer_hash, outer_iter), j, &ret);
				kh_value(kh_value(bed_buffer_hash, outer_iter), inner_iter) = 'B';
			}
		}
	}

	// Need to clean up the memories allocated for seen array
	free(cur_chrom_id);
}

void cleanBedInfo(Bed_Info *bed_info) {
	int i;
	if (bed_info) {
		for(i=0; i<bed_info->size; i++) {
			free(bed_info->coords);
		}

		free(bed_info);
	}
}

void zeroAllNsRegions(char *chrom_id, khash_t(str) *Ns_buffer_hash, Chromosome_Tracking *chrom_tracking) {
	// First, we need to find the index that is used to track current chromosome chrom_id
    uint8_t idx = locateChromosomeIndex(chrom_id, chrom_tracking);
	
    khiter_t outer_iter, inner_iter;
    for (outer_iter=kh_begin(Ns_buffer_hash); outer_iter!=kh_end(Ns_buffer_hash); outer_iter++) {
        if (kh_exist(Ns_buffer_hash, outer_iter)) {
            if (strcmp(chrom_id, kh_key(Ns_buffer_hash, outer_iter)) == 0) {
                // found the chromosome id key in the hash map
                for (inner_iter=kh_begin(kh_value(Ns_buffer_hash, outer_iter));
                        inner_iter!=kh_end(kh_value(Ns_buffer_hash, outer_iter));
                            inner_iter++) {
                    if (kh_exist(kh_value(Ns_buffer_hash, outer_iter), inner_iter)) {
                        uint32_t pos = kh_key(kh_value(Ns_buffer_hash, outer_iter), inner_iter);

                        // set the coverage variable at the pos in chrom_tracking variable to 0
                        chrom_tracking->coverage[idx][pos] = 0;
                    }
                }
            }
        }
    }
}

void getTargetsAndWriteCoverage(char * chrom_id, khash_t(str) *Ns_buffer_hash, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, FILE *cov_fp, FILE *wig_fp, FILE *wgs_fp) {
	// First, we need to find the index that is used to track current chromosome chrom_id
    uint8_t idx = locateChromosomeIndex(chrom_id, chrom_tracking);

	// for the whole genome, we need to use the file that contains regions of all Ns in the reference
	// As they will be not used, so we are going to set the count info in these regions to 0
	if (N_FILE_PROVIDED)
		zeroAllNsRegions(chrom_id, Ns_buffer_hash, chrom_tracking);

	// write to the file that contains whole genome info
    if(user_inputs->wgs_coverage) {

		char header_line[30]; 
		strcpy(header_line, ">chromosome_");
		strcat(header_line, chrom_id);
		fputs(header_line, wgs_fp);

        for(i = 0; i < chrom_tracking->chromosome_lengths[i]; i++) {
            if(i%100==0) fputc('\n', wgs_fp);

			fprintf(wgs_fp, "%d ", chrom_tracking->coverage[idx][i]);
        }
		fputc('\n', wgs_fp);
    }

	// now processing capture targets
	if (TARGET_FILE_PROVIDED) {
	    for(i = 0; i < target_info->size; i++) {
			if ( strcmp(target_info->coords[i].chr, chrom_id) != 0) 
				continue;

	        //stats_info->total_targets++;
			int j, ret;
		    int start = target_info->coords[i].tart;
			int end = target_info->coords.end;
	        int length = end - start + 1;
		    bool collect_target_cov = length > 99 ? true : false ;  // TODO: is it the min length of the exon? Check.

			if (collect_target_cov) {
				for(j = 0; j < PRIMER_SIZE; j++) {
					if ( (start - j) < 0 || (end + j) >= chrom_tracking->chromosome_lengths[idx])
						continue;

					fivePrime[j]  += chrom_tracking->coverage[idx][start-j];
					threePrime[j] += chrom_tracking->coverage[idx][end+j];
				}
			}

			bool target_hit = false;
	        short pc[101], pc2[101];

			char header_line[50];
			strcpy(header_line, ">");
			strcat(header_line, chrom_id);
			fputs(header_line, cov_fp);
			fprintf(cov_fp, " %d %d", start, end);

	        bool space_it = false;
		    if(end - start > 10000) space_it = true;   

			for(j = start; j <= end; j++) {
				if (j > chrom_tracking->chromosome_lengths[idx]) 
					continue;

	            if (space_it && (j-start)%100 == 0) cov_fp->fputc('\n');    // enter a new line after every 100 bases

		        short cov = chrom_tracking->coverage[idx][j];

			    if (cov < 0) {
					cov = 0;
	                printf("Coverage less than 0!!!!!!!\n");
		        }

				addValueToKhashBucket(stats_info->target_cov_histogram, cov, 1);

		        stats_info->total_target_coverage += cov;
			    if (cov > 0) {
					target_hit=true;
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 1, 1);
		        }

				if (cov > 4)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 5, 1);

			    if (cov > 9)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 10, 1);

				if (cov > 14)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 15, 1);

	            if (cov > 19)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 20, 1);

				if (cov > 29)
                    addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 30, 1);

				if (cov > 39)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 40, 1);

	            if (cov > 49)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 50, 1);

				if (cov > 99)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 100, 1);

	            if (cov > 499)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 500, 1);

				if (cov > 999)
					addValueToKhashBucket(stats_info->targeted_base_with_N_coverage, 1000, 1);
                    
				// output to the cov.fasta file
				fprintf(cov_fp, "%d ", cov);

				addValueToKhashBucket(stats_info->target_coverage_for_median, cov, 1);

			    if (collect_target_cov) {
					int pcpos = (int)((double)(i-start)/(double)length*100+0.5);
	                pc[pcpos] += cov;
		            pc2[pcpos]++;
			    }
			}

			// output to the cov.fasta file
			fputc('\n', cov_fp);
    
		    for(int index = 0; index < pc.length; index++) {
				if(pc2[index] != 0)
				{
					int d = (int) (((double)pc[index]/(double)pc2[index])+0.5);
	                pc[index] = (short) d;	
		        }
			}
                
	        for(int i = 0; i < 101; i++)   // TODO: why is this till 101? if its read length shouldn't we get the actual read length?
		    {
			    targetCov[i]+=pc[i];
			}
	        if(targetHit)
		    {
			    hitTargetCount++;
			}
	        else
		    {
			    missTraget.write(targetChrs[j]+"\t"+targetStarts[j]+"\t"+targetStops[j]+"\n");
				boolean hit = false;
	            for(int i = start - BUFFER; i < start && !hit; i++)
		        {
			        if(i < 0) {continue;}
				    if(COVERAGE[i] > 0){
					    hit=true;
					}
				}
	            for(int i = end; i < end+BUFFER && !hit; i++)
		        {
			        if(i >= size) {continue;}
				    if(COVERAGE[i] > 0){
					    hit=true;
					}
				}
	            if(hit){
		            hitTarget_bufferonly_Count++;
			    }
			}
		}
	}
}

void addValueToKhashBucket(khash_t *hash_in, uint32_t pos_key, uint16_t val) {
	int ret;
	khiter_t k_iter = kh_put(m16, stats_info->target_cov_histogram, pos_key, &ret);
    if (ret == 1) {
        kh_value(stats_info->target_cov_histogram, k_iter) = 0;
    } else if (ret == -1) {
        fprintf(stderr, "can't find the key for stats_info->target_cov_histogram at pos %d\n", pos_key);
		exit(1);
    }

	kh_value(stats_info->target_cov_histogram, k_iter) += val;

	return;
}
