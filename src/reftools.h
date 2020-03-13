#ifndef REFTOOLS_H_
#define REFTOOLS_H_

#include "utils.h"
#include "kseq.h"
#include <assert.h>

#define REF_NUM_MAX 128
#define REF_NAME_LEN_MAX 128
#define REF_LEN_MAX 3500000000

typedef struct{
  int refNum;
  uint32_t lookupTable[REF_NUM_MAX];
  char names[REF_NUM_MAX][REF_NAME_LEN_MAX];
  uint8_t *bases;
} Reference;

extern Reference reference;
extern char *ref_file_name;

int count_kmer(Reference *ref, uint32_t hashValue, int kmerLength);
void initialize_ref_file();
void finalize_ref_file();
void initialize_ref(Reference *ref);
int get_ref(Reference *ref);
void destroy_ref(Reference *ref);

#endif /* REFTOOLS_H_ */
