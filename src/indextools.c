/*
 * indextools.c
 *
 *  Created on: Mar 24, 2016
 *      Author: haowen
 */

#include "indextools.h"

uint32_t *occurrence_table;
int *lookup_table;

char *index_file_name;
FILE *index_fp;

FILE *header_fp;

int max_hash_value;

uint32_t hash_table_size;

int initLoadIndex() {
    int tmp;
    index_fp = fopen(index_file_name, "r");
    if (index_fp == NULL) {
        return 0;
    }
    tmp = fread(&(reference.refNum), sizeof(int), 1, index_fp);
    tmp = fread(&window_size, sizeof(int), 1, index_fp);
    tmp = fread(&step_size, sizeof(int), 1, index_fp);
    if (tmp == 0) {
        fprintf(stderr, "Load error while initializing hash table.\n");
        exit(-1);
    }

    max_hash_value = (int)(pow(4, window_size))-1;

    lookup_table = (int*) malloc(sizeof(int) * (max_hash_value+2));
    assert(lookup_table);

    return 1;
}

int loadIndex() {
    double startTime = realtime();
    int tmp = -1;
    for (int refIndex=0;refIndex<reference.refNum;++refIndex) {
        char* refName = reference.names[refIndex];
        tmp = fread(refName, sizeof(char), REF_NAME_LEN_MAX, index_fp);
        //   printf("%s\n",reference.names[refIndex]);
    }
    //load the length of reference bases
    uint32_t *refLookupTable = reference.lookupTable;
    tmp = fread(refLookupTable, sizeof(uint32_t), reference.refNum+1, index_fp);
    //load the int bases of reference
    tmp = fread(reference.bases, sizeof(uint8_t),refLookupTable[reference.refNum], index_fp);

    //load the hash table
    //load the lookup table
    tmp = fread(lookup_table,sizeof(int), max_hash_value + 2, index_fp);
    //load the occurrence table
    tmp = fread(&hash_table_size, sizeof(int), 1, index_fp);

    occurrence_table= (uint32_t*) malloc(sizeof(uint32_t) * hash_table_size);
    assert(occurrence_table);

    tmp = fread(occurrence_table, sizeof(uint32_t), hash_table_size, index_fp);

    fprintf(stderr, "Max hash value: %d, number of reference sequences: %d, number of bases: %d, load in %fs!\n", max_hash_value, reference.refNum, reference.lookupTable[reference.refNum] - reference.lookupTable[0],  realtime() - startTime);
    return tmp;
}

void finalizeLoadIndex() {
    if (occurrence_table != NULL) {
        free(occurrence_table);
        occurrence_table = NULL;
    }

    if (lookup_table != NULL) {
        free(lookup_table);
        lookup_table = NULL;
    }
    fclose(index_fp);
}

int init_saving_index() {
    int tmp = -1;
	index_fp = fopen(index_file_name, "w");
    int refNum = reference.refNum;
    fprintf(stderr, "Number of reference sequences: %d.\n", refNum);
	tmp = fwrite(&refNum, sizeof(int), 1, index_fp);
	tmp = fwrite(&window_size, sizeof(int), 1, index_fp);
	tmp = fwrite(&step_size, sizeof(int), 1, index_fp);
	if (tmp == 0){
		fprintf(stderr, "Write error while initializing hash table.\n");
        exit(-1);
    }
    hash_table_size = 0;
    for(int i = 0; i < reference.refNum; ++i){
        hash_table_size += (reference.lookupTable[i+1] - reference.lookupTable[i] - window_size + 1) / step_size;
    }
    max_hash_value = (int)(pow(4, window_size)) - 1;

    lookup_table = (int*) malloc(sizeof(int) * (max_hash_value+2));
    assert(lookup_table);
    occurrence_table= (uint32_t*) malloc(sizeof(uint32_t) * hash_table_size);
    assert(occurrence_table);

	fprintf(stderr, "The size of hash table is %d, the max hash value is %d.\n", hash_table_size, max_hash_value);
    return 0;
}

__UTILS_HASH_VALUE
int construct_index() {
    double startTime = realtime();
    HashLocationPair* _tmpOccurrenceTable = (HashLocationPair*) malloc(sizeof(HashLocationPair) * hash_table_size);
    //Compute hash value of each element
    for (int j = 0; j < reference.refNum; ++j) {
//#pragma omp parallel for
        for (uint32_t i = reference.lookupTable[j]; i < reference.lookupTable[j+1] - window_size + 1; i += step_size) {
            _tmpOccurrenceTable[i/step_size].location = i;
            int hashVal = hashValue(reference.bases + i, window_size);
            if (hashVal != -1) {
                _tmpOccurrenceTable[i/step_size].hashValue = hashVal;
            } else {
                _tmpOccurrenceTable[i/step_size].hashValue = max_hash_value + 1;
            }
        }
    }
    //Sort
    //tbb::parallel_sort
    qsort(_tmpOccurrenceTable, hash_table_size, sizeof(HashLocationPair), compare_hash_location_pair);
    //build lookup table
    //lookup table gives the start of each list in the occurrence table
    memset(lookup_table, 0, (max_hash_value + 2) * sizeof(int));
    for (uint32_t i = 0; i < hash_table_size; i++) {
        occurrence_table[i] = _tmpOccurrenceTable[i].location;
        int currentHashVal = _tmpOccurrenceTable[i].hashValue + 1;
        if (currentHashVal <= max_hash_value + 1) {
            lookup_table[currentHashVal]++;
        }else{
            break;
        }
    }

    free(_tmpOccurrenceTable);

    int sum = 0;
    for (int i = 0; i <= max_hash_value + 1; i++) {
        sum += lookup_table[i];
        lookup_table[i] = sum;
    }

    hash_table_size = sum;

    fprintf(stderr, "Number of hashed bases: %d\n",lookup_table[max_hash_value]);
    fprintf(stderr, "Generated hash table in %fs.\nMax hash value: %d, hash table size: %u.\n", realtime() - startTime, max_hash_value, hash_table_size);
    return 0;
}

int save_index() {
    int tmp;
	//store the reference
	//store the length of reference name
	//store the name of refenrence
    init_saving_header();
    char hd[1000] = "@HD\tVN:1.4\tSO:unsorted\n";
	tmp = fwrite(hd, sizeof(char), strlen(hd), header_fp);
    for (int i = 0; i < reference.refNum; ++i) {
        char* refName = reference.names[i];
        int refLen = reference.lookupTable[i + 1] - reference.lookupTable[i];
        int length = sprintf(hd, "@SQ\tSN:%s\tLN:%d\n", refName, refLen);
		tmp = fwrite(hd, sizeof(char), length, header_fp);
        tmp = fwrite(refName, sizeof(char), REF_NAME_LEN_MAX, index_fp);
    }
    finalize_saving_header();
    uint32_t *refLookupTable=reference.lookupTable;
    tmp = fwrite(refLookupTable,sizeof(uint32_t),reference.refNum+1,index_fp);
	//store the int bases of reference
	tmp = fwrite(reference.bases, sizeof(uint8_t), refLookupTable[reference.refNum], index_fp);
	//store the hash table
	//store the lookup table
	tmp = fwrite(lookup_table, sizeof(int), max_hash_value + 2, index_fp);
	//store the occurrence table
	tmp = fwrite(&hash_table_size, sizeof(uint32_t), 1, index_fp);
	tmp = fwrite(occurrence_table, sizeof(uint32_t), hash_table_size, index_fp);
	if (tmp == 0){
		fprintf(stderr, "Write error while initializing hash table.\n");
        exit(-1);
    }
    return 0;
}

int finalize_saving_index() {
	fclose(index_fp);
    return 0;
}

void init_saving_header() {
	header_fp = fopen(header_file_name, "w");
	if (header_fp == NULL){
		fprintf(stderr, "Cannot open header file.\n");
        exit(-1);
    }
}

void save_header() {
}

void finalize_saving_header() {
	fclose(header_fp);
}
