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

// It will open the target file for the first time. and count the number of target lines (items) within it.
//
long getTargetCount(char *target_file) {
    FILE *fp = fopen(target_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("target file %s open failed!", target_file);
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

// It will open the target file and then record all the start and stop positions along with chromosome ids
// for chromosome X and Y are characters, I will use array of chars (ie array of string) to store chromosome ids
//
void loadTargets(char * target_file, Target_Coords * targets) {
	FILE *fp = fopen(target_file, "r");
    if (fp == NULL) {       // ensure the target file open correctly!
        printf("target file %s open failed!", target_file);
        exit(1);
    }

	int num = 0;		// to keep track the tokens from string split using strtok()
	int count = 0;		// index for each target line (item)
	ssize_t read;
	size_t len = 0;
	char *p_token, *line;	// here p_ means a point. so p_token is a point to a token

    while((read = getline(&line, &len, fp)) != -1) {
		if (line[0] == '\0') continue;		// skip if it is a blank line
		if (line[0] == '#')  continue;		// skip if it is a comment line

		p_token = strtok(line, " \t");
		while (p_token != NULL) {
			if (num == 0) 
				strcpy(targets[count].chr,  p_token);
			else if (num == 1)
				targets[count].start = atoi(p_token);
			else if (num == 2) {
				targets[count].end = atoi(p_token);
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

/*
void generateTargetBufferLookupTable(Target_Coords * target_coords, long target_line_count, khash_t(32) *target_buffer_hash[]) {
	int i, j, ret, chr_id;
	for (i = 0; i < target_line_count; i++) {
		char *chr = target_coords[i].chr;
		if (strcmp(chr, "X")) {
			chr_id = 23;
		} else if (strcmp(chr, "Y")) {
			chr_id = 24;
		} else {
			chr_id = atoi(chr);
		}
		chr_id--;		// As the array is 0 indexed, need to decrease it by one
		khash_t(32) tmp_hash = target_buffer_hash[chr_id];

		long start, end;
		khiter_t k_iter;

		// for positions on targets
		for (j=start; j<=end; j++) {
			k_iter = kh_put(32, tmp_hash, j, &ret);
			kh_value(tmp_hash, k_iter) = 'T';
		}

		// for positions on the buffer at the left side
		for (j=start-BUFFER; j < start; j++) {
			if (j < 0) continue;
			k_iter = kh_put(32, tmp_hash, j, &ret);
			kh_value(h, k_iter) = 'B';
		}

		// for positions on the buffer at the right side
		for (j=end+1; j < end+BUFFER; j++ ) {
			k_iter = kh_put(32, tmp_hash, j, &ret);
			kh_value(tmp_hash, k_iter) = 'B';
		}
	}
}*/

/*
void getTargetsAndWriteCoverage(char * chromo,  short COVERAGE[], char * covFasta, char * missTraget, char * wgCoverage)
{
    if(wgCoverage != null)
    {
        wgCoverage.write(">"+chromo);
        for(int i = 0; i < COVERAGE.length; i++)
        {
            if(i%100==0) wgCoverage.write("\n");
            wgCoverage.write(COVERAGE[i]+" ");
        }
        wgCoverage.write("\n");
    }
    for(int j = 0; j < targetChrs.length; j++)
    {
        if(!removechr(targetChrs[j]).equals(removechr(chromo))) {continue;}
        //totalTargets++;
        int start = targetStarts[j];
        int end = targetStops[j];
        int length = end - start+1;
        boolean collectTargetCov = length > 99 ;  // TODO: is it the min length of the exon? Check.

        if(supertets)
        {
            System.out.println(targetChrs[j]+" "+start+" "+end);
        }
        if(collectTargetCov)
        {
            for(int i = 0; i < prime_size; i++)
            {
                if((start - i) < 0 || (end+i) >= size){
                    ///System.err.println("The BED Target "+targetChrs[j]+" "+start+" "+end+" is going out of Bound!!!\n");
                    continue;
                }
                fivePrime[i]+=COVERAGE[start-i];
                threePrime[i]+=COVERAGE[end+i];

                ///fivePrime[i]+=COVERAGE[end-i+300];
                //threePrime[i]+=COVERAGE[start+i-300];

            }
        }

        if(supertets)
        {

            for(int i = 0; i < 500; i++)
            {
                if((start-i) < 0) {continue;}
                System.out.print( (start-i)+" ");
            }
            System.out.print("\n");
            for(int i = 0; i < 500; i++)
            {
                if((end+i) >= size) {continue;}
                System.out.print( (end+i)+" ");
            }
            System.out.print("\n");

            supertets= false;
        }

        boolean targetHit = false;
        short[] pc = new short[101];
        short[] pc2 = new short[101];

        covFasta.write(">"+chromo+" "+start+" "+end+"\n");
        boolean spaceit = false;
        if(end - start > 10000) spaceit = true;   
        for(int i = start; i <= end; i++)
        {
            if(i > size) {continue;}
            if(spaceit && (i-start)%100 == 0) covFasta.write("\n");    // enter a new line after every 100 bases
            short cov = COVERAGE[i];
            if(cov < 0)
            {
                cov = 0;
                System.err.println("Coverage less than 0!!!!!!!\n");
            }
            short temp_cov = cov;
            if(temp_cov >= covHisto.length){
                temp_cov = (short) (covHisto.length-1);
            }
            covHisto[temp_cov]++;
            totalTargetCoverage+=cov;
            if(cov > 0)
            {
                targetHit=true;
                basesWithOneHitorMore++;
            }
            if(cov > 9){
                basesWith10HitsorMore++;
			}
            if(cov > 19){
                basesWith20HitsorMore++;
			}
            if(cov > 39){
                basesWith40HitsorMore++;
			}
            if(cov > 49){
                basesWith50HitsorMore++;
			}
            if(cov > 99){
                basesWith100HitsorMore++;
			}
            if(cov > 499){
                basesWith500HitsorMore++;
			}
            if(cov > 999){
                basesWith1000HitsorMore++;
			}
                    
            covFasta.write(cov+" ");

            if(cov < coverage_forMedian.length){
                coverage_forMedian[cov]++;
            }else{
                int[] tmp = new int[coverage_forMedian.length];
                System.arraycopy(coverage_forMedian, 0, tmp, 0, coverage_forMedian.length);
                coverage_forMedian = new int[cov+1];
                System.arraycopy(tmp, 0, coverage_forMedian, 0, tmp.length);
                coverage_forMedian[cov]++;
            }

            if(collectTargetCov)
            {
                int pcpos = (int)((double)(i-start)/(double)length*100+0.5);
                pc[pcpos] += cov;
                pc2[pcpos]++;
            }
        }
        covFasta.write("\n");
    
        for(int index = 0; index < pc.length; index++)
        {
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
}*/

/*
void getTargetAndBufferPositions(char * chromosome_id, int size, char *target_buffer_regions, Target_Coords *target_coords) {
	// remove 'chr' part of the chromosome id 
    chromosome_id = removechr(chromosome_id);
	int i, j;
    for(j = 0; j < target_coords.length; j++) {
		// skip if it is a different chromosome
        if(!strcmp(removechr(target_coords.chr[j]), chromosome_id)) continue;

        int start = target_coords.start[j];
        int end = target_coords.stop[j];

		// hit on the target region
        for(i = start; i <= end; i++) {
            if(i >= size) {
                continue;
            } else {
                target_buffer_regions[i] = 'T';
				TOTAL_TARGETED_BASES++;
            }
        }

		// hit on the left side of the buffer of a target
        for(int i = start - BUFFER; i <start; i++) {
            if (i < 0) {
                continue;
            } else {
                if(target_buffer_regions[i] != 'T'){
                    target_buffer_regions[i] = 'B';
					TOTAL_BUFFER_BASES++;
                }
            }
        }

		// hit on the right side of the buffer of a target
        if (end <= (size-1)) {    //end < (size-1)
            for (int i = end+1; i <= end+BUFFER; i++)
            {
                if (targat_buffer_regions[i] != 'T') {
                    target_buffer_regions[i] = 'B';
					TOTAL_BUFFER_BASES++;
                }
            }
        } else {
            continue;
        }
	}
	return;
}
*/
