#include "reftools.h"

KSEQ_INIT(int, read)
__UTILS_CHAR_UINT8

Reference reference;

FILE *reference_file;
char *reference_file_path;
kseq_t *reference_kseq;

void initialize_ref(Reference *reference) {
  reference->refNum = 0;
  reference->lookupTable[0]=0;
  reference->bases = (uint8_t*) _mm_malloc(sizeof(uint8_t) * REF_LEN_MAX, 16);
  assert(reference->bases);
}

void destroy_ref(Reference *reference) {
  if (reference->bases != NULL) {
    _mm_free(reference->bases);
    reference->bases = NULL;
  }
  reference->refNum = 0;
}

void initialize_ref_file() {
  reference_file = fopen(reference_file_path, "r");
  if (reference_file == NULL) {
    fprintf(stderr, "Cannnot open read file: %s.\n", reference_file_path);
    exit(EXIT_FAILURE);
  } else {
    reference_kseq = kseq_init(fileno(reference_file));
  }
}

void finalize_ref_file() {
  kseq_destroy(reference_kseq);
  fclose(reference_file);
}

int get_ref(Reference *reference) {
  int l = kseq_read(reference_kseq);
  while (l == 0) {
    l = kseq_read(reference_kseq);
  }
  while (l > 0) {
    reference->refNum++;
    reference->lookupTable[reference->refNum] = reference->lookupTable[reference->refNum - 1] + l;
    memset(reference->names[reference->refNum - 1], '\0', sizeof(char) * REF_NAME_LEN_MAX);
    strcpy(reference->names[reference->refNum - 1], reference_kseq->name.s);
    uint32_t i;
    for (i = reference->lookupTable[reference->refNum - 1]; i < reference->lookupTable[reference->refNum - 1] + l; ++i) {
      char c = reference_kseq->seq.s[i-reference->lookupTable[reference->refNum - 1]];
      reference->bases[i] = charToUint8(c);
    }
    fprintf(stderr, "%s\n", reference->names[reference->refNum-1]);
    l = kseq_read(reference_kseq);
  } 
  if (l == -1) {
    /*end of file*/
    return 0;
  } else {
    /* truncated quality string*/
    fprintf(stderr, "truncated quality string.");
    exit(EXIT_FAILURE);
  }
  /*never come to here*/
  return 0;
}
