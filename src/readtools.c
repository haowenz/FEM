/*
 * readtools.c
 *
 *  Created on: Mar 23, 2016
 *      Author: howen
 */

#include "readtools.h"
#include "mapper.h"

KSEQ_INIT(int, read)
__UTILS_CHAR_UINT8

FILE *read_seq_fp1;
FILE *read_seq_fp2;
char *read_file_name1;
char *read_file_name2;
kseq_t *read_seq1;
kseq_t *read_seq2;
pthread_mutex_t readMutex1;
pthread_mutex_t readMutex2;
pthread_mutex_t read_queue_mutex;
pthread_cond_t read_queue_pro_cond;
pthread_cond_t read_queue_cos_cond;

Read *read_queue;
int read_queue_size_max = 500000;
int read_queue_front;
int read_queue_rear;
int read_queue_size;
int no_more_read;

void initReadQueue() {
	resetReadQueue();
	read_queue = (Read*) malloc(read_queue_size_max * sizeof(Read));
    assert(read_queue);
	fprintf(stderr, "Initialize read queue successfully!\n");
}

void resetReadQueue() {
	pthread_mutex_init(&read_queue_mutex, NULL);
	pthread_cond_init(&read_queue_pro_cond, NULL);
	pthread_cond_init(&read_queue_cos_cond, NULL);
	no_more_read = 0;
	read_queue_front = 0;
	read_queue_rear = 0;
	read_queue_size = 0;
}

void destroyReadQueue() {
	pthread_mutex_destroy(&read_queue_mutex);
	pthread_cond_destroy(&read_queue_pro_cond);
	pthread_cond_destroy(&read_queue_cos_cond);
	free(read_queue);
	fprintf(stderr, "Destroy read queue successfully!\n");
}

inline void pushReadQueue(Read *read) {
	memcpy(read_queue + read_queue_rear, read, sizeof(Read));
	read_queue_rear = (read_queue_rear + 1) % read_queue_size_max;
	++read_queue_size;
}

int popReadQueue(Read *read) {
	pthread_mutex_lock(&read_queue_mutex);
	if (read_queue_size == 0 && no_more_read == 1) {
        //fprintf(stderr, "Quit.\n");
		pthread_mutex_unlock(&read_queue_mutex);
        pthread_cond_signal(&read_queue_pro_cond);
		return -1;
	}
	while (read_queue_size == 0) {
        pthread_cond_signal(&read_queue_cos_cond);
		pthread_cond_wait(&read_queue_pro_cond, &read_queue_mutex);
		if (read_queue_size == 0 && no_more_read == 1) {
            //fprintf(stderr, "Waked up and quit.\n");
			pthread_mutex_unlock(&read_queue_mutex);
			pthread_cond_signal(&read_queue_pro_cond);
			return -1;
		}
	}
	memcpy(read, read_queue + read_queue_front, sizeof(Read));
	read_queue_front = (read_queue_front + 1) % read_queue_size_max;
	--read_queue_size;
    if (no_more_read == 1) {
        //fprintf(stderr, "Number of left reads: %d.\n", read_queue_size);
        if (read_queue_size > 0) {
            pthread_cond_signal(&read_queue_pro_cond);
        }
    } else {
        pthread_cond_signal(&read_queue_cos_cond);
    }
	pthread_mutex_unlock(&read_queue_mutex);
	return 1;
}

void singleReadQueueThread() {
	Read read;
	while (1) {
		int flag = getSingleRead(&read);
		if (flag == -1) {
			break;
		}
		pthread_mutex_lock(&read_queue_mutex);
		while (read_queue_size == read_queue_size_max) {
			pthread_cond_signal(&read_queue_pro_cond);
			pthread_cond_wait(&read_queue_cos_cond, &read_queue_mutex);
		}
		pushReadQueue(&read);
		pthread_mutex_unlock(&read_queue_mutex);
		pthread_cond_signal(&read_queue_pro_cond);
	}
	pthread_mutex_lock(&read_queue_mutex);
	no_more_read = 1;
    //fprintf(stderr, "No more new read, %d reads left!\n", read_queue_size);
	pthread_mutex_unlock(&read_queue_mutex);
    pthread_cond_signal(&read_queue_pro_cond);
    //pthread_cond_broadcast(&read_queue_pro_cond);
    //fprintf(stderr, "Broadcast!\n");
}

void initReadBlock(ReadBlock *readBlock) {
	readBlock->num = 0;
	readBlock->length = -1;
	readBlock->bases = NULL;
	readBlock->rc_bases = NULL;
	readBlock->names = (char*) malloc(sizeof(char) * READ_NAME_LEN_MAX * READBLOCK_NUM_MAX);
	assert(readBlock->names);
}

void resetReadFile() {
	/*single-end mode*/
	rewind(read_seq_fp1);
	read_seq1 = kseq_init(fileno(read_seq_fp1));
}

void initReadFile() {
	/*single-end mode*/
	readMutex1= PTHREAD_MUTEX_INITIALIZER;

	read_seq_fp1 = fopen(read_file_name1, "r");
	if (read_seq_fp1 == NULL) {
		fprintf(stderr, "Cannnot open read file: %s.\n", read_file_name1);
        exit(-1);
	} else {
		read_seq1 = kseq_init(fileno(read_seq_fp1));
	}
	/*pair-end mode*/
	//eadMutex2= PTHREAD_MUTEX_INITIALIZER;
}

void finalizeReadFile() {
	/*single-end mode*/
	kseq_destroy(read_seq1);
	fclose(read_seq_fp1);
	/*pair-end mode*/

}

int getSingleRead(Read *read) {
	int l = kseq_read(read_seq1);
	while (l == 0) {
		l = kseq_read(read_seq1);
		printf("!");
	}
	if (l > 0) {
		read->length = l;
		memset(read->name, '\0', sizeof(char) * READ_NAME_LEN_MAX);
		memset(read->charBases, '\0', sizeof(char) * READ_LEN_MAX);
		strcpy(read->name, read_seq1->name.s);
		strcpy(read->charBases, read_seq1->seq.s);
		int i;
		for (i = 0; i < l; ++i) {
			char c = read_seq1->seq.s[i];
			uint8_t base = charToUint8(c);
			read->bases[i] = base;
			if (base != 4) {
				read->rc_bases[l - 1 - i] = 3 - base;
			} else {
				read->rc_bases[l - 1 - i] = 4;
			}

		}
		return 1;
	} else if (l == -1) {
		/*end of file*/
		return -1;
	} else {
		/* truncated quality string*/
		fprintf(stderr, "truncated quality string.");
        exit(-1);
	}
	/*never come to here*/
	return 0;
}

int getPairRead(Read *read1, Read *read2) {
	return 0;
}

int getReadBlock(ReadBlock *readBlock) {
	int flag = pthread_mutex_trylock(&readMutex1);
	if (flag != 0) {
		return 0;
	}
	memset(readBlock->names, '\0',
			sizeof(char) * READBLOCK_NUM_MAX * READ_NAME_LEN_MAX);
	int count = 0;
	while (count < READBLOCK_NUM_MAX) {
		int l = kseq_read(read_seq1);
		while (l == 0) {
			l = kseq_read(read_seq1);
		}
		if (l > 0) {
			if (count == 0) {
				if (l != readBlock->length) {
					readBlock->length = l;
					if (readBlock->bases != NULL) {
						_mm_free(readBlock->bases);
					}
					readBlock->bases = (uint8_t*) _mm_malloc(
							sizeof(uint8_t) * l * READBLOCK_NUM_MAX, 16);
					printf("!!!!!!!!!!\n");
					assert(readBlock->bases);
					if (readBlock->rc_bases != NULL) {
						_mm_free(readBlock->rc_bases);
					}
					readBlock->rc_bases = (uint8_t*) _mm_malloc(
							sizeof(uint8_t) * l * READBLOCK_NUM_MAX, 16);
					assert(readBlock->rc_bases);
				}
			} else {
				/*XAM will support variant length later*/
				if (readBlock->length != l) {
					fprintf(stderr,"The length of reads is not equal.\n");
                    exit(-1);
				}

			}
			memcpy(readBlock->names + count * READ_NAME_LEN_MAX,
					read_seq1->name.s, sizeof(uint8_t) * READ_NAME_LEN_MAX);
//			strcpy(readBlock->names + count * READ_NAME_LEN_MAX,
//					read_seq1->name.s);
			int i;
			for (i = 0; i < l; ++i) {
				char c = read_seq1->seq.s[i];
				uint8_t base = charToUint8(c);
				readBlock->bases[count * l + i] = base;
				if (base != 4) {
					readBlock->rc_bases[count * l + l - 1 - i] = 3 - base;
				} else {
					readBlock->rc_bases[count * l + l - 1 - i] = 4;
				}

			}
			count++;
		} else if (l == -1) {
			/*end of file*/
			readBlock->num = count;
			pthread_mutex_unlock(&readMutex1);
			return -1;
		} else {
			pthread_mutex_unlock(&readMutex1);
			/* truncated quality string*/
			fprintf(stderr, "truncated quality string.");
            exit(-1);
		}
	}
	readBlock->num = count;
	pthread_mutex_unlock(&readMutex1);
	return 0;
}
