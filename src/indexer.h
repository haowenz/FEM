/*
 * hashtable.h
 *
 *  Created on: Mar 24, 2016
 *      Author: howen
 */

#ifndef SRC_INDEXER_H_
#define SRC_INDEXER_H_

#include "utils.h"
#include "reftools.h"
#include <math.h>
#include <queue>
#include <string>
#include <algorithm>

using namespace std;
class HashTable {
public:
	HashTable(Reference *ref, const int windowSize,const int stepSize,const char *fileName);
	~HashTable();
	void clear();
	void config();
	bool generate();
	void addRef(string name,int length);

	//for saving and loading
	void initSaving();
	void save();
	void finalizeSaving();
	void initSaveHeader();
	void saveHeader();
	void finalizeSaveHeader();

	int initLoading();
	int load();
	void finalizeLoading();

	inline int hashValue(const uint8_t *seq);
	inline int hashValue(const uint8_t *seq,int windowSize);

	int WINDOW_SIZE;
    int maxHashValue;
	uint32_t hashTableSize;
	int stepSize;

	uint32_t *occurrenceTable;
	int *lookupTable;

	queue<string> refNames;
	queue<int> refLens;

	/*for saving and loading*/
	char fileName[100];
	FILE* fp;
	char headerFileName[100];
	FILE *header_fp;
	Reference *ref;
};

bool compare(const HashLocationPair &a, const HashLocationPair &b);

#endif /* SRC_INDEXER_H_ */
