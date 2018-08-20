/*
 * =====================================================================================
 *
 *       Filename:  user_defined_annotation.c
 *
 *    Description:  the implementation details for the user_defined_annotation.h
 *
 *        Version:  1.0
 *        Created:  02/12/2018 03:41:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Peiming (Peter) Huang, phuang@bcm.edu
 *        Company:  Baylor College of Medicine, Houston
 *
 * =====================================================================================
 */

#include "user_defined_annotation.h"

void getUserDefinedDatabaseInfo(User_Input *user_inputs, User_Defined_Database_Wrapper *udd_wrapper, khash_t(khStrInt) *cds_lengths, khash_t(khStrInt) *cds_counts, khash_t(khStrInt) *user_defined_targets) {
	// open user-defined-database for read
	//
	FILE *fp = fopen(user_inputs->user_defined_database_file, "r");

	if (fp == NULL) {
		printf("user defined database file %s open failed!", user_inputs->user_defined_database_file);
		exit(EXIT_FAILURE);
	}

	// variables related to User_Defined_Database_Wrapper
	//
	uint64_t total_line_count = 0; 

	// Store all the chromosome info into a KHASH_MAP with string as key and uint32_t as value
	//
	khash_t(khStrInt) *chromosome_info = kh_init(khStrInt);

	// variables used for reading file
	//
	char *gene_symbol = calloc(50,  sizeof(char));
	char *transcript_name   = calloc(70,  sizeof(char));
	char *gene_info   = calloc(100, sizeof(char));
	char *target_line = calloc(100, sizeof(char));
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	char *tokPtr;

	while ((read = getline(&line, &len, fp)) != -1) {
		// printf("%s", line);
		// here we need to handle the case where it is an empty line
		// here we compare char not string (char*), so I need to de-reference it
		//
		if (*line == '\n')
			continue;

		char *savePtr = line;
		int absent = 0;
		khiter_t iter;

		// now need to get cds length and count information
		// 7       87173445        87173591        ABCB1|ENST00000265724|cds_18|gene
		//
		uint32_t j=0, start=0, end=0;

		while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
			if (j==0) {
				strcpy(target_line, tokPtr);
				
				// check to see if the hashkey exists for current chromosome id
				//
				iter = kh_put(khStrInt, chromosome_info, tokPtr, &absent);
				if (absent) {
					kh_key(chromosome_info, iter) = strdup(tokPtr);
					kh_value(chromosome_info, iter) = 0;
				}
				kh_value(chromosome_info, iter)++;	// update the value
			}

			if (j==1) {
				start = (uint32_t) strtol(tokPtr, NULL, 10);
				strcat(target_line, ":");
				strcat(target_line, tokPtr);
			}
	
			if (j==2) {
				end = (uint32_t) strtol(tokPtr, NULL, 10);
				strcat(target_line, ":");
				strcat(target_line, tokPtr);

				// need to handle the case where start = end (convert 1-based annotation to 0-based annotation)
				//
				if (start == end)
					start--;
			}

			if (j==3) {
				strcpy(gene_info, tokPtr);
			}

			j++;
		}

		// Now record the target_line here
		//
		iter = kh_put(khStrInt, user_defined_targets, target_line, &absent);    // get key iterator for target_line
		if (absent) {
			kh_key(user_defined_targets, iter)   = strdup(target_line);
			kh_value(user_defined_targets, iter) = 0;
		}

		kh_value(user_defined_targets, iter)++;

		// now process gene_info to extract gene_symbol and transcript_name
		// ABCB1|ENST00000265724_cds_18
		//
		char *gene_token;
		savePtr = gene_info;
		char *tmp_name=calloc(50, sizeof(char));

		// need to finish the strtok_r() for gene_token
		//
		j=0;
		while ((gene_token = strtok_r(savePtr, "|", &savePtr))) {
			if (j==0)
				strcpy(gene_symbol, gene_token);

			if (j==1)
				strcpy(tmp_name, gene_token);

			j++;
		}

		// Now processing transcript info
		// ENST00000265724_cds_18	or	NM_000531_cds_3
		//
		char *name_token;
		savePtr = tmp_name;
		bool nm_flag=false;

		j=0;
		while ((name_token=strtok_r(savePtr, "_", &savePtr))) {
			if (j==0) {
				strcpy(transcript_name, name_token);

				// some transcript name like the following:
				// NM_005957_cds_6
				//
				if ( (stristr(name_token, "NM") != NULL) || (stristr(name_token, "NR") != NULL) || (stristr(name_token, "XM") != NULL))
					nm_flag=true;
			}

			if (j==1) {
				if (nm_flag) {
					strcat(transcript_name, "_");
					strcat(transcript_name, name_token);
				}
			}
			j++;
		}

		// Now combine both gene_symbol and transcript_name together to form a unique gene ID
		//
		char *gene = calloc(strlen(gene_symbol)+strlen(transcript_name)+10, sizeof(char));
		strcpy(gene, gene_symbol);
		strcat(gene, "-");
		strcat(gene, transcript_name);

		iter = kh_put(khStrInt, cds_lengths, gene, &absent);
		if (absent) {
			kh_key(cds_lengths, iter)   = strdup(gene);
			kh_value(cds_lengths, iter) = 0;
		}
		
		// need to handle the case where start == end (convert 1-based to 0-based)
		// and use it as the key to store total cds length and cds count
		//
		if (start == end)
			start--;
		kh_value(cds_lengths, iter) += end - start;

		iter = kh_put(khStrInt, cds_counts, gene, &absent);
		if (absent) {
			kh_key(cds_counts, iter)   = strdup(gene);
			kh_value(cds_counts, iter) = 0;
		}
		kh_value(cds_counts, iter)++;

		if (gene != NULL) free(gene);
		if (tmp_name  != NULL) free(tmp_name);

		total_line_count++;
	}

	if (line != NULL) free(line);
	if (transcript_name != NULL) free(transcript_name);
	if (gene_info != NULL) free(gene_info);
	if (gene_symbol != NULL) free(gene_symbol);
	if (target_line != NULL) free(target_line);

	udd_wrapper->num_of_lines = total_line_count;

	// Now process hash map to update the udd_wrapper->udd_database (User_Defined_Database)
	//
	uint32_t num_of_chrom = 0;
	khiter_t iter;
	for (iter = kh_begin(chromosome_info); iter != kh_end(chromosome_info); ++iter) {
		if (kh_exist(chromosome_info, iter)) {
			//printf("Key is %s\t", kh_key(chromosome_info, iter));
			//printf("Value is %d\n", kh_value(chromosome_info, iter));
			num_of_chrom++;
		}
	}

	//printf("total number of chromosome is %"PRIu32"\n", num_of_chrom);
	//printf("total number of lines is %"PRIu64"\n", total_line_count);
	
	// Allocate the memory and store the cds count information for each chromosome
	//
	udd_wrapper->num_of_chroms = num_of_chrom;
	udd_wrapper->ud_database_per_chrom = calloc(num_of_chrom, sizeof(User_Defined_Database));

	uint32_t chrom_idx=0;
	char *tmp_cid = calloc(50, sizeof(char));
	for (iter = kh_begin(chromosome_info); iter != kh_end(chromosome_info); ++iter) {
		if (kh_exist(chromosome_info, iter)) {
			udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds = kh_value(chromosome_info, iter);
			udd_wrapper->ud_database_per_chrom[chrom_idx].chrom_id = calloc(strlen(kh_key(chromosome_info, iter))+1, sizeof(char));

			strcpy(udd_wrapper->ud_database_per_chrom[chrom_idx].chrom_id, kh_key(chromosome_info, iter));
			chrom_idx++;
		}
	}

	if (tmp_cid != NULL) free(tmp_cid);

	// clean-up
	//
	cleanKhashStrInt(chromosome_info);

	fclose(fp);
}

void processUserDefinedDatabase(User_Input *user_inputs, Regions_Skip_MySQL *exon_regions, User_Defined_Database_Wrapper *udd_wrapper, Raw_User_Defined_Database * raw_user_defined_database, khash_t(khStrInt) *cds_lengths, khash_t(khStrInt) *cds_counts) {
	// As we already know the total number of line and the cds count for each chromosome
	// we will just need to create the story memory here
	// Note: since the chrom_ids on the User_Defined_Database_Wrapper is not in order, 
	// the exon_regions->chromosome_ids array won't be in order as well
	//
	exon_regions->chrom_list_size = udd_wrapper->num_of_chroms;
	exon_regions->chromosome_ids  = calloc(exon_regions->chrom_list_size, sizeof(char*));

	exon_regions->starts = calloc(exon_regions->chrom_list_size, sizeof(uint32_t*));
	exon_regions->ends   = calloc(exon_regions->chrom_list_size, sizeof(uint32_t*));
	exon_regions->size_r = calloc(exon_regions->chrom_list_size, sizeof(uint32_t*));
	
	exon_regions->prev_search_loc_index   = 0;
	exon_regions->prev_search_chrom_index = 0;

	exon_regions->gene = calloc(exon_regions->chrom_list_size, sizeof(char**));
	exon_regions->Synonymous = calloc(exon_regions->chrom_list_size, sizeof(char**));
	exon_regions->prev_genes = calloc(exon_regions->chrom_list_size, sizeof(char**));
	exon_regions->exon_info = calloc(exon_regions->chrom_list_size, sizeof(char**));

	uint32_t i, j;
	for (i=0; i<exon_regions->chrom_list_size; i++) {
		// assign chromosome id here so that the indices will be the same across different variables
		//
		exon_regions->chromosome_ids[i] = calloc(strlen(udd_wrapper->ud_database_per_chrom[i].chrom_id)+1, sizeof(char));
		strcpy(exon_regions->chromosome_ids[i], udd_wrapper->ud_database_per_chrom[i].chrom_id);

		exon_regions->size_r[i] = udd_wrapper->ud_database_per_chrom[i].num_of_cds;
		exon_regions->starts[i] = calloc(exon_regions->size_r[i], sizeof(uint32_t));
		exon_regions->ends[i]   = calloc(exon_regions->size_r[i], sizeof(uint32_t));
		exon_regions->gene[i]   = calloc(exon_regions->size_r[i], sizeof(char*));
		exon_regions->Synonymous[i] = calloc(exon_regions->size_r[i], sizeof(char*));
		exon_regions->prev_genes[i] = calloc(exon_regions->size_r[i], sizeof(char*));
		exon_regions->exon_info[i]  = calloc(exon_regions->size_r[i], sizeof(char*)); 

		for (j=0; j<exon_regions->size_r[i]; j++) {
			exon_regions->gene[i][j] = NULL;
			exon_regions->Synonymous[i][j] = NULL;
			exon_regions->prev_genes[i][j] = NULL;
			exon_regions->exon_info[i][j]  = NULL;
		}
	}

	// allocate initial memories for raw user defined database
	//
	raw_user_defined_database->chrom_id     = calloc(udd_wrapper->num_of_chroms, sizeof(char*));
	raw_user_defined_database->annotations  = calloc(udd_wrapper->num_of_chroms, sizeof(char**));
	raw_user_defined_database->annotation_size = calloc(udd_wrapper->num_of_chroms, sizeof(uint32_t));
	raw_user_defined_database->num_of_chromosomes  = udd_wrapper->num_of_chroms;

	for (i=0; i<raw_user_defined_database->num_of_chromosomes; i++) {
		// Assign chromosome id here so that the indices will be the same across different variables
		//
		raw_user_defined_database->chrom_id[i] = calloc(strlen(udd_wrapper->ud_database_per_chrom[i].chrom_id)+1, sizeof(char));
		strcpy(raw_user_defined_database->chrom_id[i], udd_wrapper->ud_database_per_chrom[i].chrom_id);

		// allocate memories for the detailed annotations
		//
		raw_user_defined_database->annotations[i] = calloc(udd_wrapper->ud_database_per_chrom[i].num_of_cds, sizeof(char*));

		raw_user_defined_database->annotation_size[i] = 0;
	}

	// now reset the num_of_cds variable inside udd_wrapper->ud_database_per_chrom for the future usage
	//
	for (i=0; i<udd_wrapper->num_of_chroms; i++) {
		udd_wrapper->ud_database_per_chrom[i].num_of_cds = 0;
	}

	// open user-defined-database file for read again and populate the exon_regions
	//
	FILE *fp = fopen(user_inputs->user_defined_database_file, "r");

	if (fp == NULL) {
		printf("user defined database file %s open failed!", user_inputs->user_defined_database_file);
		exit(EXIT_FAILURE);
	}

	char *gene = calloc(100, sizeof(char));				// store tmp gene/transcript name info
	char *buf_cds_length = calloc(50, sizeof(char));	// store tmp cds length to a string
	char *buf_cds_count  = calloc(50, sizeof(char));	// store tmp cds count  to a string
	char *orig_line  = calloc(300, sizeof(char));		// deep copy of a line read as strtok will destroy the string
	char *annotation = calloc(200, sizeof(char));		// store the last part of the annotation for the future usage
	char *tmp_cds_info = calloc(150, sizeof(char));
	char *line = NULL;

	// store unique gene symbol and transcripts
	// number of transcripts will be stored in both cds_lengths and cds_counts
	//
	//khash_t(khStrInt) *list_of_genes = kh_init(khStrInt);
	//khash_t(khStrInt) *list_of_transcripts = kh_init(khStrInt);
	khiter_t iter;

	size_t len = 0;
	ssize_t read;
	
	while ((read = getline(&line, &len, fp)) != -1) {
		// handle empty line here
		//
		if (*line == '\n')
			continue;

		// make a copy first as strtok will destroy the 'line' variable
		// change the last character of the line with '\0' to remove the newline character '\r\n'
		//
		if (line[strlen(line)-1] == '\n')
			line[strlen(line)-1] = '\0';

		if (line[strlen(line)-1] == '\r')
			line[strlen(line)-1] = '\0';

		strcpy(orig_line, line);
		char *savePtr = line;
		char *tokPtr;
		int32_t chrom_idx=-1;

		// line example:
		// 7       87178664        87178834        ABCB1|ENST00000543898_cds_12
		//
		j=0;
		while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
			if (j==0) {
				// need to find the index for current chromosome
				//
				for (i=0; i<exon_regions->chrom_list_size; i++) {
					// need to remove 'chr' from the chromosome id
					//
					if (strcmp(tokPtr, exon_regions->chromosome_ids[i]) == 0) {
						chrom_idx = i;
						break;
					}
				}

				if (chrom_idx == -1) {
					fprintf(stderr, "Something went wrong as the chromosome index shouldn't be -1\n");
					exit(EXIT_FAILURE);
				}
			}

			if ( j == 1) {
				// for start position
				//
				exon_regions->starts[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds] = (uint32_t) strtol(tokPtr, NULL, 10);
			}

			if (j == 2) {
				// for end position
				//
				exon_regions->ends[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds] = (uint32_t) strtol(tokPtr, NULL, 10);

				// need to handle the case where start == end (convert 1-based to 0-based)
				//
				if (exon_regions->starts[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds] == 
						exon_regions->ends[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds]) 
					exon_regions->starts[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds]--;
			}

			if (j == 3) {
				// It seems that I can't do multiple strtok at a time as it will mess up the internal string poiter
				// So I save it first here for the later usage
				//
				strcpy(annotation, tokPtr);
			}
			j++;
		}

		// once we have annotation saved, we will use strtok_r again to do the splitting
		// for annotation: gene|cds_info:
		// ABCB1|ENST00000543898_cds_12
		//
		savePtr = annotation;
		char *token;
		j=0;

		while ((token = strtok_r(savePtr, "|", &savePtr))) {
			if (j==0) {
				exon_regions->gene[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds] = calloc(strlen(token)+1, sizeof(char));
				strcpy(exon_regions->gene[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], token);

				exon_regions->Synonymous[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds] = calloc(3, sizeof(char));
				strcpy(exon_regions->Synonymous[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], ".");

				exon_regions->prev_genes[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds] = calloc(3, sizeof(char));
				strcpy(exon_regions->prev_genes[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], ".");

				// record gene symbol here
				/*
				int absent;
				iter = kh_put(khStrInt, list_of_genes, token, &absent);                                                   
				if (absent) {                                                                                         
					kh_key(list_of_genes, iter)   = strdup(token);                                                        
					kh_value(list_of_genes, iter) = 0;                                                         
				}                                                                                                     
				kh_value(list_of_genes, iter)++;
				*/
				strcpy(gene, token);	// store gene symbol for later usage
			}

			if (j==1)
				strcpy(tmp_cds_info, token);


			if (j > 1) {
				fprintf(stderr, "Something wrong with the user-defined database. Please check the format %s\n", token);
				exit(EXIT_FAILURE);
			}
			j++;
		}

		// allocate memory: here we need to allocate extra space as space-holder for other sources
		// here is an example from the MySQL database annotation
		// NM_000090_exon_21;NR_037401_exon_0	CCDS2297_exon_21	ENST00000304636.7_exon_21   hsa-mir-3606_exon_0
		// Note: for UCSC annotation, I will put it along with the CCDS column
		//
		exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds] = calloc(strlen(tmp_cds_info)+50, sizeof(char));

		// Here I need to find out the source of the annotation
		//
		if ( (stristr(tmp_cds_info, "NM") != NULL) || (stristr(tmp_cds_info, "NR") != NULL) || (stristr(tmp_cds_info, "XM") != NULL)) {
			// refseq
			//
			strcpy(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], tmp_cds_info);
			strcat(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], "\t.\t.\t.\t.");
		} else if ( (stristr(tmp_cds_info, "CCDS") != NULL) || (stristr(tmp_cds_info, "uc") != NULL) ) {
			// CCDS
			//
			strcpy(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], ".\t");
			strcat(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], tmp_cds_info);
			strcat(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], "\t.\t.\t.");
		} else if ( (stristr(tmp_cds_info, "OTT") != NULL) || (stristr(tmp_cds_info, "ENST") != NULL) ) {
			// for VEGA or Gencode
			//
			strcpy(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], ".\t.\t");
			strcat(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], tmp_cds_info);
			strcat(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], "\t.\t.");
		} else if ( (stristr(tmp_cds_info, "hsa") != NULL) || (stristr(tmp_cds_info, "mir") != NULL) ) {
			// for miRNA
			//
			strcpy(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], ".\t.\t.\t");
			strcat(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], tmp_cds_info);
			strcat(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], "\t.");
		} else {
			// everything else goes here
			//
			strcpy(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], ".\t.\t.\t.\t");
			strcat(exon_regions->exon_info[chrom_idx][udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds], tmp_cds_info);
		}

		// The following is for the cds/gene percentage coverage culculation
		//
		char *cds_token;
		savePtr = tmp_cds_info;
		char *transcript_name = calloc(50, sizeof(char));
		j=0;

		while ((cds_token = strtok_r(savePtr, "_", &savePtr))) {
			if (j==0) {
				// get transcript name and add it to the end of gene symbol (would be unique)
				//
				strcat(gene, "-");
				strcat(gene, cds_token);
				strcpy(transcript_name, cds_token);
			}

			if (j==1) {
				// However, some transcript name like the following:
				// NM_005957_cds_6_3_1333
				//
				if ( (stristr(transcript_name, "NM") != NULL) || (stristr(transcript_name, "NR") != NULL) || (stristr(transcript_name, "XM") != NULL)) {
					strcat(gene, "_");
					strcat(gene, cds_token);
				}
			}
			j++;
		}

		if (transcript_name) free(transcript_name);

		// get cds_lengths for current gene/transcript from the cds_lengths and cds_counts khash tables
		//
		iter = kh_get(khStrInt, cds_lengths, gene);
		if (iter == kh_end(cds_lengths)) {  // iter will be equal to kh_end if key not present
			fprintf(stderr, "The gene-transcript %s not present!\n", gene);
		} else {
			sprintf(buf_cds_length, "%"PRIu32, kh_value(cds_lengths, iter));
		}

		// get cds_count information for current gene/transcript
		//
		iter = kh_get(khStrInt, cds_counts, gene);
		if (iter == kh_end(cds_counts)) {
			fprintf(stderr, "The gene-transcript %s not present!\n", gene);
		} else {
			sprintf(buf_cds_count, "%"PRIu32, kh_value(cds_counts, iter));
		}

		// need to store the cds length and count information to the annotation
		//
		strcat(orig_line, "_");
		strcat(orig_line, buf_cds_count);
		strcat(orig_line, "_");
		strcat(orig_line, buf_cds_length);

		// store raw user-defined-database here
		//
		uint32_t anno_idx = raw_user_defined_database->annotation_size[chrom_idx];
		raw_user_defined_database->annotations[chrom_idx][anno_idx] = calloc(strlen(orig_line)+1, sizeof(char));
		strcpy(raw_user_defined_database->annotations[chrom_idx][anno_idx], orig_line);
		raw_user_defined_database->annotation_size[chrom_idx]++;

		// increment the index
		//
		udd_wrapper->ud_database_per_chrom[chrom_idx].num_of_cds++;
	}

	fclose(fp);

	// clean-up
	//
	if (line != NULL) free(line);
	if (gene != NULL) free(gene);
	if (orig_line  != NULL) free(orig_line);
	if (annotation != NULL) free(annotation);
	if (buf_cds_length != NULL) free(buf_cds_length);
	if (buf_cds_count  != NULL) free(buf_cds_count);
	if (tmp_cds_info != NULL) free(tmp_cds_info);
}

// the following is used for the setup of coverage percentage calculation
//
void userDefinedGeneCoverageInit(khash_t(khStrLCG) *user_cds_gene_hash, char *chrom_id, Raw_User_Defined_Database *raw_user_defined_database, khash_t(khStrStrArray) *gene_transcripts) {
	// first need to find out the index used to store the annotation information
	//
	uint32_t i;
   	int32_t chrom_idx=-1;

	for (i=0; i<raw_user_defined_database->num_of_chromosomes; i++) {
		if (strcmp(chrom_id, raw_user_defined_database->chrom_id[i]) == 0) {
			chrom_idx = i;
			break;
		}
	}

	if (chrom_idx == -1)
		return;

	// We need to create a new variable seen_transcript, to avoid add it to the gene_transcripts over and over
	//
	khash_t(khStrInt) *seen_transcript = kh_init(khStrInt);

	// since strtok_r() destroy the string it splits and it stores the pointer for it
	// we will have to finish strtok_r() for one string completely before using another strtok_r()
	// Therefore, we need to use some inter-mediate variables for string split
	//
	char *savePtr, *tokPtr;

	// process one record at a time
	//
	for (i=0; i<raw_user_defined_database->annotation_size[chrom_idx]; i++) {
		char *line = calloc(strlen(raw_user_defined_database->annotations[chrom_idx][i])+5, sizeof(char));
		strcpy(line, raw_user_defined_database->annotations[chrom_idx][i]);

		// Example of a line:
		// 15      75047131        75047429        CYP1A2|NM_000761_cds_6_1_12_1342
		//
		char *gene_annotation = calloc(strlen(raw_user_defined_database->annotations[chrom_idx][i])+2, sizeof(char));
		savePtr = line;
		uint32_t j=0;
		uint32_t cds_target_start=0, cds_target_end=0;

		while ((tokPtr = strtok_r(savePtr, "\t", &savePtr))) {
			if (j==1) {
				cds_target_start = (uint32_t) strtol(tokPtr, NULL, 10);
			}

			if (j==2) {
				cds_target_end = (uint32_t) strtol(tokPtr, NULL, 10);

				// need to handle the case where start = end (convert 1-based annotation to 0-based annotation)
				//
				if (cds_target_start == cds_target_end )
					cds_target_start--;
			}

			if (j == 3) {
				strcpy(gene_annotation, tokPtr);
			}
			j++;
		}

		if (line != NULL) free(line);

		// For debugging
		//
		char *ret = strstr(gene_annotation, "|");
		if (ret == NULL) {
			printf("something is wrong %s and the string is %s in iteration %"PRIu32"\n", 
					gene_annotation, raw_user_defined_database->annotations[chrom_idx][i], i);
		}

		// here is the cds annotation part
		// ANKK1|ENST00000303941_cds_1_8_2290
		//
		char *gene_token=NULL, *cds_annotation=NULL, *gene_symbol=NULL;
		savePtr = gene_annotation;
		j=0;

		while ((gene_token = strtok_r(savePtr, "|", &savePtr))) {
			if (j==0) {
				gene_symbol = calloc(strlen(gene_token)+1, sizeof(char));
				strcpy(gene_symbol, gene_token);
			}

			if (j==1) {
				// get the second part (cds) of the annotation
				// ENST00000303941_cds_1_8_2290
				//
				cds_annotation = calloc(strlen(gene_token)+2, sizeof(char));
				strcpy(cds_annotation, gene_token);
			}
			j++;
		}

		// Now processing the following part
		// ENST00000303941_cds_1_8_2290
		//
		char *cds_token=NULL, *transcript_name=NULL;
		savePtr = cds_annotation;
		bool flag_2nd_field=false;	// if the 2nd field is part of transcript name
		j=0;
		uint16_t exon_id=0, exon_count=0;
		uint32_t cds_length=0;

		while ((cds_token = strtok_r(savePtr, "_", &savePtr))) {
			if (j==0) {
				transcript_name = calloc(strlen(cds_token)+30, sizeof(char));
				strcpy(transcript_name, cds_token);

				// some transcript name like the following:
				// NM_005957_cds_6_3_1333
				//
				if ( (stristr(transcript_name, "NM") != NULL) || (stristr(transcript_name, "NR") != NULL) || (stristr(transcript_name, "XM") != NULL))
					flag_2nd_field=true;
			}

			if (j==1) {
				if (flag_2nd_field) {
					strcat(transcript_name, "_");
					strcat(transcript_name, cds_token);
				}
			}

			if (flag_2nd_field) {
				// j==2 is just 'cds', so skip it
				//
				if (j==3)
					exon_id = (int16_t) strtol(cds_token, NULL, 10);

				if (j==4)
					exon_count = (int16_t) strtol(cds_token, NULL, 10);

				if (j==5)
					cds_length = (int32_t) strtol(cds_token, NULL, 10);
			} else {
				if (j == 2)
					exon_id = (int16_t) strtol(cds_token, NULL, 10);

				if (j == 3)
					exon_count = (int16_t) strtol(cds_token, NULL, 10);
	
				if (j == 4) 
					cds_length = (int32_t) strtol(cds_token, NULL, 10);
			}
			j++;
		}

		// Once we have gene symbol and transcript name, we need to add them into the gene_transcripts
		//
		addToGeneTranscriptKhashTable(gene_symbol, transcript_name, gene_transcripts, seen_transcript);

		// Next we have everything, we need to add them on the hash table, user_cds_gene_hash
		//
		khiter_t iter;
		uint32_t current_key = getHashKey(cds_target_start);
		while (current_key < cds_target_end) {
			char key_str[20];
			sprintf(key_str, "%"PRIu32, current_key);

			// Initialize the bucket with the current key
			//
			lowCoverageGeneHashBucketKeyInit(user_cds_gene_hash, key_str);

			// get the bucket
			//
			iter = kh_get(khStrLCG, user_cds_gene_hash, key_str);
			if (iter == kh_end(user_cds_gene_hash)) {
				// key doesn't exists!
				// This shouldn't happen. So we exit!
				//
				fprintf(stderr, "Memory re-allocation failed at the userDefinedGeneCoverageInit()!\n");
				exit(EXIT_FAILURE);
			}

			// update Gene_Coverage variable ==> user_cds_genes->gene_coverage
			//
			uint32_t idx = kh_value(user_cds_gene_hash, iter)->total_size;
			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].gene_symbol = calloc(strlen(gene_symbol)+1, sizeof(char));
			strcpy(kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].gene_symbol, gene_symbol);

			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].transcript_name = calloc(strlen(transcript_name)+1, sizeof(char));
			strcpy(kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].transcript_name, transcript_name);

			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].cds_start  = cds_target_start;
			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].cds_end    = cds_target_end;
			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].cds_length = cds_length;

			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].cds_target_start = cds_target_start;
			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].cds_target_end   = cds_target_end;
			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].exon_id    = exon_id;
			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].exon_count = exon_count;

			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].targeted = true;
			kh_value(user_cds_gene_hash, iter)->gene_coverage[idx].low_cov_regions = NULL;
			kh_value(user_cds_gene_hash, iter)->total_size++;

			current_key += 1000;
		}

		if (cds_annotation  != NULL) free(cds_annotation);
		if (transcript_name != NULL) free(transcript_name);
		if (gene_annotation != NULL) free(gene_annotation);
		if (gene_symbol     != NULL) free(gene_symbol);
	}

	cleanKhashStrInt(seen_transcript);
}

void recordUserDefinedTargets(khash_t(khStrInt) *user_defined_targets, Bed_Info *user_defined_bed_info) {
	// First, let's find out number of unique targets in the user-defined database
	//
	uint32_t i=0, num_of_targets=0;
	khiter_t iter;
	for (iter = kh_begin(user_defined_targets); iter != kh_end(user_defined_targets); ++iter) {
		if (kh_exist(user_defined_targets, iter)) {
			num_of_targets++;
		}
	}

	user_defined_bed_info->size = num_of_targets;
	user_defined_bed_info->coords = calloc(user_defined_bed_info->size, sizeof(Bed_Coords));

	// now load all the user defined target into the bed_info
	//
	char *line = calloc(100, sizeof(char));
	char *tokPtr, *savePtr;
	for (iter = kh_begin(user_defined_targets); iter != kh_end(user_defined_targets); ++iter) {
		if (kh_exist(user_defined_targets, iter)) {
			// Since strtok_r() will destroy the string value it splits,
			// it would be a good idea to make a copy first and work on the copy one
			//
			strcpy(line, kh_key(user_defined_targets, iter));
			savePtr = line;

			uint32_t j=0;
			while ((tokPtr = strtok_r(savePtr, ":", &savePtr))) {
				if (j==0)
					strcpy(user_defined_bed_info->coords[i].chrom_id,  tokPtr);

				if (j==1)
					user_defined_bed_info->coords[i].start = (uint32_t) strtol(tokPtr, NULL, 10);

				if (j==2) {
					user_defined_bed_info->coords[i].end = (uint32_t) strtol(tokPtr, NULL, 10);

					// need to handle the case where start = end (convert 1-based annotation to 0-based annotation)
					//
					if (user_defined_bed_info->coords[i].start == user_defined_bed_info->coords[i].end)
						user_defined_bed_info->coords[i].start--;
				}
				j++;
			}
			i++;
		}
	}

	if (line != NULL) free(line);
				
}

// for the low coverage gene annotations
//
void writeCoverageForUserDefinedDB(char *chrom_id, Bed_Info *target_info, Chromosome_Tracking *chrom_tracking, User_Input *user_inputs, Stats_Info *stats_info, Regions_Skip_MySQL *exon_regions, Regions_Skip_MySQL *intron_regions) {
	// First, we need to find the index that is used to track current chromosome chrom_id
	//
	int32_t idx = locateChromosomeIndexForChromTracking(chrom_id, chrom_tracking);
	if (idx == -1) return;

	// Need to write information to the following files
	//
	//FILE *missed_target_fp = fopen(user_inputs->missed_targets_file, "a");
	FILE *low_x_fp    = fopen(user_inputs->capture_low_cov_file, "a");                            
	FILE *all_site_fp = fopen(user_inputs->capture_all_site_file, "a");
	
	uint32_t i, j;
	for(i = 0; i < target_info->size; i++) {
		// skip if chromosome id doesn't match
		//
		if ( strcmp(target_info->coords[i].chrom_id, chrom_id) != 0)
			continue;

		stats_info->cov_stats->total_targets++;
		int32_t start = target_info->coords[i].start;
		int32_t end = target_info->coords[i].end;
		int length = end - start + 1;
		bool collect_target_cov = length > 99 ? true : false ;  // TODO: is it the min length of the exon? Check.

		if (collect_target_cov) {
			for(j = 0; j < PRIMER_SIZE; j++) {
				int32_t k = start - j;      // k could go negative, so it is a signed integer
				if ( k < 0 || (end + j) >= chrom_tracking->chromosome_lengths[idx])
					continue;

				stats_info->five_prime[j]  += chrom_tracking->coverage[idx][k];
				stats_info->three_prime[j] += chrom_tracking->coverage[idx][end+j];
			}
		}

		// In the original java code, the following is used for the pc and pc2 definition
		// short pc[101], pc2[101];
		// After discussion with Divya and Qiaoyan, here is what we understand.
		// pc and pc2 is used to setup the 'target_coverage' variable which is defined as an array of 101 in the origain java code
		// here I will use percentage_bin and percentage_count_per_bin
		//
		uint32_t percentage_bin[101], percentage_count_per_bin[101];
		for(j=0; j<101; j++) {
			percentage_bin[j] = 0;
			percentage_count_per_bin[j] = 0;
		}

		bool target_hit = false;

		for(j = 0; j < length; j++) {                                                                     
			// check if it passes the end of the chromosome
			//
			if (j+start >= chrom_tracking->chromosome_lengths[idx])                                       
				continue;

			uint32_t cov = chrom_tracking->coverage[idx][j+start];                                        
			if (!target_hit && (cov > 0))                                                                 
				target_hit = true;
		}

		//if (!target_hit) {
			// need to write to the missed target file
			//
			//fprintf(missed_target_fp, "%s\t%"PRIu32"\t%"PRIu32"\n", chrom_id, start, end);
		//}

		// For All Sites Report. Note: the end position is not included based on the bed format
		// Here we don't have to worry about the NULL pointers passed in as it will be handled 
		//
		produceCaptureAllSitesReport(start, length-1, chrom_tracking, chrom_id, user_inputs, all_site_fp, intron_regions, exon_regions);

		// For low coverage and high coverage Report
		//
		writeLow_HighCoverageReport(start, length, chrom_tracking, chrom_id, user_inputs, low_x_fp, NULL, intron_regions, exon_regions);
	}

	fclose(low_x_fp);
	//fclose(high_x_fp);
	fclose(all_site_fp);
	//fclose(missed_target_fp); 
}

void cleanUserDefinedDatabase(User_Defined_Database_Wrapper *udd_wrapper) {
	uint32_t i;
	for (i=0; i<udd_wrapper->num_of_chroms; i++) {
		if (udd_wrapper->ud_database_per_chrom[i].chrom_id != NULL)
			free(udd_wrapper->ud_database_per_chrom[i].chrom_id);
	}

	if (udd_wrapper->ud_database_per_chrom != NULL)
		free(udd_wrapper->ud_database_per_chrom);

	if (udd_wrapper != NULL)
		free(udd_wrapper);
}

void cleanRawUserDefinedDatabase(Raw_User_Defined_Database * raw_user_defined_database) {
	uint32_t i, j;
	for (i=0; i<raw_user_defined_database->num_of_chromosomes; i++) {
		if (raw_user_defined_database->chrom_id[i] != NULL)
			free(raw_user_defined_database->chrom_id[i]);
	
		for (j=0; j<raw_user_defined_database->annotation_size[i]; j++) {
			if (raw_user_defined_database->annotations[i][j] != NULL)
				free(raw_user_defined_database->annotations[i][j]);
		}

		if (raw_user_defined_database->annotations[i] != NULL)
			free(raw_user_defined_database->annotations[i]);
	}

	if (raw_user_defined_database->annotations != NULL)
		free(raw_user_defined_database->annotations);

	if (raw_user_defined_database->chrom_id != NULL)
		free (raw_user_defined_database->chrom_id);

	if (raw_user_defined_database->annotation_size != NULL)
		free(raw_user_defined_database->annotation_size);

	if (raw_user_defined_database !=NULL )
		free(raw_user_defined_database);
}
