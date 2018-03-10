/*
 * readtools.h
 *
 *  Created on: Mar 23, 2016
 *      Author: howen
 */

#ifndef READTOOLS_H_
#define READTOOLS_H_

#include "utils.h"
#include "kseq.h"
#include <assert.h>


#define READ_LEN_MAX 300
#define READ_NAME_LEN_MAX 100
#define READBLOCK_NUM_MAX 1000000

extern int read_queue_size_max;

typedef struct {
	int length;
	char name[READ_NAME_LEN_MAX];
	uint8_t bases[READ_LEN_MAX];
	uint8_t rc_bases[READ_LEN_MAX];
	char charBases[READ_LEN_MAX];
} Read;

typedef struct {
	int num;
	int length;
	char *names;
	uint8_t *bases;
	uint8_t *rc_bases;
} ReadBlock;

extern char *read_file_name1;
extern char *read_file_name2;

void initReadBlock(ReadBlock *readBlock);
void initReadFile();
void resetReadFile();
void finalizeReadFile();
int getSingleRead(Read *read);
int getReadBlock(ReadBlock *readBlock);
int popReadQueue(Read *read);
void singleReadQueueThread();
void initReadQueue();
void destroyReadQueue();
void resetReadQueue();

#endif /* READTOOLS_H_ */
