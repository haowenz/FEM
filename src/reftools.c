#include "reftools.h"

KSEQ_INIT(int, read)
__UTILS_CHAR_UINT8

Reference reference;

FILE *ref_seq_fp;
char *ref_file_name;
kseq_t *ref_seq;

void initialize_ref(Reference *ref) {
  ref->refNum = 0;
  ref->lookupTable[0]=0;
  ref->bases = (uint8_t*) _mm_malloc(sizeof(uint8_t) * REF_LEN_MAX, 16);
  assert(ref->bases);
}

void destroy_ref(Reference *ref) {
  if (ref->bases != NULL) {
    _mm_free(ref->bases);
    ref->bases = NULL;
  }
  ref->refNum = 0;
}

void initialize_ref_file() {
  ref_seq_fp = fopen(ref_file_name, "r");
  if (ref_seq_fp == NULL) {
    fprintf(stderr, "Cannnot open read file: %s.\n", ref_file_name);
    exit(-1);
  } else {
    ref_seq = kseq_init(fileno(ref_seq_fp));
  }
}

void finalize_ref_file() {
  kseq_destroy(ref_seq);
  fclose(ref_seq_fp);
}

int get_ref(Reference *ref) {
  int l = kseq_read(ref_seq);
  while (l == 0) {
    l = kseq_read(ref_seq);
  }
  while (l > 0) {
    ref->refNum++;
    ref->lookupTable[ref->refNum] = ref->lookupTable[ref->refNum-1]+l;
    memset(ref->names[ref->refNum-1], '\0', sizeof(char) * REF_NAME_LEN_MAX);
    strcpy(ref->names[ref->refNum-1], ref_seq->name.s);
    uint32_t i;
    for (i = ref->lookupTable[ref->refNum-1]; i < ref->lookupTable[ref->refNum-1]+l; ++i) {
      char c = ref_seq->seq.s[i-ref->lookupTable[ref->refNum-1]];
      ref->bases[i] = charToUint8(c);
    }
    fprintf(stderr, "%s\n", ref->names[ref->refNum-1]);
    l = kseq_read(ref_seq);
  } 
  if (l == -1) {
    /*end of file*/
    return 0;
  } else {
    /* truncated quality string*/
    fprintf(stderr, "truncated quality string.");
    exit(-1);
  }
  /*never come to here*/
  return 0;
}
