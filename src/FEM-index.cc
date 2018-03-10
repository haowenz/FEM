/*
 * FEM-index.cc
 *
 *  Created on: Mar 23, 2016
 *      Author: haowen
 */

#include "FEM-index.h"

static int print_usage()
{
    fprintf(stderr, "Usage: FEM index <window_size> <step_size> <ref>\n");
    return 1;
}

int index_main(int argc, char* argv[]) 
{
	if(argc<4)
    {
        return print_usage();
	}

	int WINDOW_SIZE = atoi(argv[1]);
    int stepSize = atoi(argv[2]);
	ref_file_name = argv[3];
    
    /* TODO:add arguments check */ 

	Reference ref;
	initRef(&ref);

	char indexFileName[100];
	strcpy(indexFileName, ref_file_name);
	sprintf(indexFileName + strlen(indexFileName), ".index");
	HashTable *hashTable = new HashTable(&ref,WINDOW_SIZE,stepSize, indexFileName);

	double startTime = 0.0;
	int refLen = 0;

	initRefFile();
    startTime = realtime();
    refLen = getRef(&ref);
    if (refLen > 0) 
    {
        fprintf(stderr, "Reference lenght: %d. Loaded into memory in %fs.\n", refLen,realtime() - startTime);
    }


    startTime = realtime();
    hashTable->initSaving();
    hashTable->config();
    hashTable->generate();
    hashTable->save();
    fprintf(stderr, "Number of bases: %d, number of chromosomes =%d.\n Built index in %fs.\n",ref.lookupTable[ref.refNum],ref.refNum, realtime() - startTime);
    destroyRef(&ref);
    finalizeRefFile();
    hashTable->finalizeSaving();
	/*save header SAM file*/
	hashTable->initSaveHeader();
	hashTable->saveHeader();
	hashTable->finalizeSaveHeader();
	delete hashTable;
    
	return 0;
}
