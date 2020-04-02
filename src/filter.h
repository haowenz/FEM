#ifndef FILTER_H_
#define FILTER_H_

#include "index.h"
#include "kvec.h"
#include "sequence_batch.h"
#include "utils.h"

uint32_t generate_group_seeding_candidates(const FEMArgs *fem_args, const SequenceBatch *read_sequence_batch, size_t read_index, int is_reverse_complement, const SequenceBatch *reference_sequence_batch, const Index *index, kvec_t_uint64_t *buffer1, kvec_t_uint64_t *buffer2, kvec_t_uint64_t *candidates, uint32_t *num_candidates_without_additonal_qgram_filter);

#endif // FILTER_H_
