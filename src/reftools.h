/*
 * reftools.h
 *
 *  Created on: Mar 24, 2016
 *      Author: howen
 */

#ifndef REFTOOLS_H_
#define REFTOOLS_H_

#include "utils.h"
#include "kseq.h"
#include <assert.h>

#define REF_NUM_MAX 128
#define REF_NAME_LEN_MAX 100
#define REF_LEN_MAX 3500000000

typedef struct{
    int refNum;
    uint32_t lookupTable[REF_NUM_MAX];
    char names[REF_NUM_MAX][REF_NAME_LEN_MAX];
    uint8_t *bases;
} Reference;

extern Reference reference;
extern char *ref_file_name;

int countKmer(Reference *ref,uint32_t hashValue,int kmerLength);
void initRefFile();
void finalizeRefFile();
void initRef(Reference *ref);
int getRef(Reference *ref);
void destroyRef(Reference *ref);

#endif /* SRC_REFTOOLS_H_ */
