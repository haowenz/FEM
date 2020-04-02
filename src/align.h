#ifndef ALIGN_H_
#define ALIGN_H_

#include <emmintrin.h>
#include <smmintrin.h>

#include "sequence_batch.h"
#include "utils.h"

#define ALPHABET_SIZE 5
#define NUM_VPU_LANES 8

int verify_candidates(const FEMArgs *fem_args, const SequenceBatch *read_sequence_batch, size_t read_sequence_index, int is_reverse_complement, const SequenceBatch *reference_sequence_batch, const uint64_t *candidates, uint32_t num_candidates, kstring_t *result_string);
int banded_edit_distance(const FEMArgs *fem_args, const char *pattern, const char *text, int read_length, int *mapping_end_position);
void vectorized_banded_edit_distance(const FEMArgs *fem_args, const uint32_t vpu_index, const SequenceBatch *reference_sequence_batch, const char *text, int read_length, const uint64_t *candidates, uint32_t num_candidates, int16_t *mapping_edit_distances, int16_t *mapping_end_positions);

#endif // ALIGN_H_
