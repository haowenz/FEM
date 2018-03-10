/*
 * hashtable.cc
 *
 *  Created on: Mar 24, 2016
 *      Author: howen
 */

#include "indexer.h"

inline int HashTable::hashValue(const uint8_t *seq) {
	int i = 0;
	int val = 0;
	while (i < WINDOW_SIZE) {
		if (seq[i] == (uint8_t) 4) {
			return -1;
		}
		val = (val << 2) | (int) (seq[i]);
		i++;
	}
	return val;
}

inline int HashTable::hashValue(const uint8_t *seq,int windowSize) {
	int i = 0;
	int val = 0;
	while (i < windowSize) {
		if (seq[i] == (uint8_t) 4) {
			return -1;
		}
		val = (val << 2) | (int) (seq[i]);
		i++;
	}
	return val;
}

HashTable::HashTable(Reference *ref, const int windowSize, const int stepSize, const char* fileName) {
	this->ref = ref;

    this->WINDOW_SIZE = windowSize;
	this->stepSize = stepSize;

	maxHashValue = (int) pow(4, WINDOW_SIZE)-1;
	hashTableSize = REF_LEN_MAX/stepSize;
	occurrenceTable = new uint32_t[hashTableSize];
	assert(occurrenceTable);
	lookupTable = new int[maxHashValue+2];
	assert(lookupTable);

	strcpy(this->fileName, fileName);
	strcpy(this->headerFileName, fileName);
	sprintf(headerFileName + strlen(headerFileName), ".header");
	header_fp = NULL;
	fp = NULL;
}

HashTable::~HashTable() {
	clear();
}

void HashTable::clear() {
	if (occurrenceTable != NULL) {
		delete[] occurrenceTable;
		occurrenceTable = NULL;
	}
	if (lookupTable != NULL) {
		delete[] lookupTable;
		lookupTable = NULL;
	}
	

}

void HashTable::config() {
	hashTableSize=0;//(ref->baseLen - WINDOW_SIZE + 1)/stepSize;
    for(int i=0;i<ref->refNum;++i){
        hashTableSize+=(ref->lookupTable[i+1]-ref->lookupTable[i]-WINDOW_SIZE+1)/stepSize;
    }
	printf("hashTableSize=%d.\n",hashTableSize);
}

bool HashTable::generate() {
	double startTime = realtime();
	HashLocationPair* _tmpOccurrenceTable = new HashLocationPair[hashTableSize];
//Compute hash value of each element
    for(int j=0;j<ref->refNum;++j){
#pragma omp parallel for
        for (uint32_t i = ref->lookupTable[j]; i <ref->lookupTable[j+1]-WINDOW_SIZE+1;i+=stepSize) {
            _tmpOccurrenceTable[i/stepSize].location = i;
            int hashVal = hashValue(ref->bases + i, WINDOW_SIZE);
            if (hashVal != -1) {
                _tmpOccurrenceTable[i/stepSize].hashValue = hashVal;
            } else {
                _tmpOccurrenceTable[i/stepSize].hashValue = maxHashValue + 1;
            }
        }
    }
//Sort
	//tbb::parallel_sort
    sort(_tmpOccurrenceTable, _tmpOccurrenceTable + hashTableSize, compare);
//build lookup table
//lookup table gives the start of each list in the occurrence table
	memset(lookupTable, 0, (maxHashValue + 2) * sizeof(int));
	for (uint32_t i = 0; i < hashTableSize; i++) {
		occurrenceTable[i] = _tmpOccurrenceTable[i].location;
		int currentHashVal = _tmpOccurrenceTable[i].hashValue + 1;
		if (currentHashVal <= maxHashValue + 1) {
			lookupTable[currentHashVal]++;
		}else{
			break;
		}
	}

	int sum = 0;
	for (int i = 0; i <= maxHashValue+1; i++) {
		sum += lookupTable[i];
		lookupTable[i] = sum;
	}

	hashTableSize=sum;

	printf("%d\n",lookupTable[maxHashValue]);

	delete[] _tmpOccurrenceTable;

	printf("Generated hash table in %fs.\nMax hash value=%d, hash size=%ud.\n", realtime() - startTime, maxHashValue, hashTableSize);

    return true;
}

void HashTable::initSaving() {
	int tmp;
	fp = fopen(fileName, "w");
    int refNum = ref->refNum;
    printf("refNum=%d.\n",refNum);
	tmp = fwrite(&refNum, sizeof(int), 1, fp);
	tmp = fwrite(&WINDOW_SIZE, sizeof(int), 1, fp);
	tmp = fwrite(&stepSize, sizeof(int), 1, fp);
	if (tmp == 0){
		fprintf(stderr, "Write error while initializing hash table.\n");
        exit(-1);
    }
}

void HashTable::save() {
	int tmp;
	//store the reference
	//store the length of reference name
	//store the name of refenrence
    for(int i=0;i<ref->refNum;++i){
        char* refName = ref->names[i];
        refNames.push(string(refName));
        refLens.push(ref->lookupTable[i+1]-ref->lookupTable[i]);
        tmp = fwrite(refName, sizeof(char), REF_NAME_LEN_MAX, fp);
    }
    uint32_t *refLookupTable=ref->lookupTable;
    tmp = fwrite(refLookupTable,sizeof(uint32_t),ref->refNum+1,fp);
	//store the int bases of reference
	tmp = fwrite(ref->bases, sizeof(uint8_t), refLookupTable[ref->refNum], fp);
	//store the hash table
	//store the lookup table
	tmp = fwrite(lookupTable, sizeof(int), maxHashValue + 2, fp);
	//store the occurrence table
	tmp = fwrite(&hashTableSize, sizeof(uint32_t), 1, fp);
	tmp = fwrite(occurrenceTable, sizeof(uint32_t), hashTableSize, fp);
	if (tmp == 0){
		fprintf(stderr, "Write error while initializing hash table.\n");
        exit(-1);
    }
}
void HashTable::finalizeSaving() {
	fclose(fp);
}

void HashTable::addRef(string name, int length) {
	refNames.push(name);
	refLens.push(length);
}

void HashTable::initSaveHeader() {
	header_fp = fopen(headerFileName, "w");
	if (header_fp == NULL){
		fprintf(stderr, "Cannot open header file.\n");
        exit(-1);
    }
}

void HashTable::saveHeader() {
	char hd[1000] = "@HD\tVN:1.4\tSO:unsorted\n";
	fwrite(hd, sizeof(char), strlen(hd), header_fp);

	if (refNames.empty() && refLens.empty()) {
		return;
	}

	do {
		string tmp = refNames.front();
		const char* refName = tmp.data();
		refNames.pop();
		int refLen = refLens.front();
		refLens.pop();
		int length = sprintf(hd, "@SQ\tSN:%s\tLN:%d\n", refName, refLen);
		fwrite(hd, sizeof(char), length, header_fp);
	} while (!refNames.empty() && !refLens.empty());
}
void HashTable::finalizeSaveHeader() {
	fclose(header_fp);
}

int HashTable::initLoading() {
	int tmp;
	fp = fopen(fileName, "r");
	if (fp == NULL)
		return 0;
	tmp = fread(&(ref->refNum), sizeof(int), 1, fp);
	tmp = fread(&WINDOW_SIZE, sizeof(int), 1, fp);
	tmp = fread(&stepSize, sizeof(int), 1, fp);
	if (tmp == 0) {
		fprintf(stderr, "Load error while initializing hash table.\n");
        exit(-1);
	}
	maxHashValue = (int) pow(4, WINDOW_SIZE)-1;
	return 1;
}

int HashTable::load() {
	double startTime = realtime();
	int tmp;
	//load reference
	//load the length of reference name
	int refNameLength = 0;
	tmp = fread(&refNameLength, sizeof(int), 1, fp);
	if (tmp == 0) {
		return 0;
	}
	//load the name of refenrence
	//char* refName = ref->name;
	//tmp = fread(refName, sizeof(char), refNameLength, fp);
	//refName[refNameLength] = '\0';
	//load the length of reference bases
	int refLength = 0;
	tmp = fread(&refLength, sizeof(int), 1, fp);
	//load the int bases of reference
	tmp = fread(ref->bases, sizeof(uint8_t), refLength, fp);

	//the pointer to intBases should not change
	//ref->nameLen = refNameLength;
	//ref->baseLen = refLength;

	//load the hash table
	//load the lookup table
	tmp = fread(lookupTable, sizeof(int), maxHashValue + 2, fp);
	//load the occurrence table
	tmp = fread(&hashTableSize, sizeof(int), 1, fp);
	//_occurrenceTable = new int[_hashTableSize];
	tmp = fread(occurrenceTable, sizeof(int), hashTableSize, fp);
	printf("Load in %fs!\n", realtime() - startTime);
	return tmp;
}

void HashTable::finalizeLoading() {
	fclose(fp);
}

bool compare(const HashLocationPair &a, const HashLocationPair &b) {
	if (a.hashValue < b.hashValue) {
		return true;
	} else if (a.hashValue == b.hashValue) {
		if (a.location < b.location)
			return true;
		else
			return false;
	} else
		return false;
}

