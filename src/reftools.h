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
extern char *reference_file_path;

int count_kmer(Reference *reference, uint32_t hash_value, int kmer_size);
void initialize_ref_file();
void finalize_ref_file();
void initialize_ref(Reference *reference);
int get_ref(Reference *reference);
void destroy_ref(Reference *reference);

#endif /* REFTOOLS_H_ */
