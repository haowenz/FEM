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

void initialize_ReadBlock(ReadBlock *readBlock);
void initialize_read_file();
void reset_read_file();
void finalize_read_file();
int get_single_read(Read *read);
int get_ReadBlock(ReadBlock *readBlock);
int pop_read_queue(Read *read);
void single_read_queue_thread();
void initialize_read_queue();
void destroy_read_queue();
void reset_read_queue();

#endif /* READTOOLS_H_ */
