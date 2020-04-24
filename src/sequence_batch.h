#ifndef SEQUENCEBATCH_H_
#define SEQUENCEBATCH_H_

#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <zlib.h>

#include "kseq.h"
#include "kvec.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
  uint32_t num_loaded_sequences;
  uint32_t max_num_sequences;
  uint64_t num_bases;
  gzFile sequence_file; 
  kseq_t *sequence_kseq;
  kvec_t(kseq_t*) sequences;
  //kvec_t(kvec_t(char)) negative_sequences;
  kvec_t(kvec_t_char) negative_sequences;
} SequenceBatch;

static inline void swap_sequences_in_sequence_batch(SequenceBatch *a, SequenceBatch *b) {
  uint32_t a_num_loaded_sequences = a->num_loaded_sequences;
  uint32_t a_max_num_sequences = a->max_num_sequences;
  uint64_t a_num_bases = a->num_bases;
  a->num_loaded_sequences = b->num_loaded_sequences;
  a->max_num_sequences = b->max_num_sequences;
  a->num_bases = b->num_bases;
  b->num_loaded_sequences = a_num_loaded_sequences;
  b->max_num_sequences = a_max_num_sequences;
  b->num_bases = a_num_bases;
  kv_swap(kseq_t*, a->sequences, b->sequences);
}

static inline void swap_kstring_t(kstring_t *a, kstring_t *b) {
  kstring_t tmp = *a;
  *a = *b;
  *b = tmp;
}

static inline const char * get_sequence_from_sequence_batch_at(const SequenceBatch *sequence_batch, size_t sequence_index) {
  return kv_A(sequence_batch->sequences, sequence_index)->seq.s;
}

static inline const char * get_negative_sequence_from_sequence_batch_at(const SequenceBatch *sequence_batch, uint32_t sequence_index) {
  return kv_A(sequence_batch->negative_sequences, sequence_index).v.a;
}

static inline uint32_t get_sequence_length_from_sequence_batch_at(const SequenceBatch *sequence_batch, size_t sequence_index) {
  return kv_A(sequence_batch->sequences, sequence_index)->seq.l;
}

static inline char * get_sequence_name_from_sequence_batch_at(const SequenceBatch *sequence_batch, uint32_t sequence_index) {
  return kv_A(sequence_batch->sequences, sequence_index)->name.s;
}

static inline size_t get_sequence_name_length_from_sequence_batch_at(const SequenceBatch *sequence_batch, uint32_t sequence_index) {
  return kv_A(sequence_batch->sequences, sequence_index)->name.l;
}

static inline char * get_sequence_qual_from_sequence_batch_at(const SequenceBatch *sequence_batch, uint32_t sequence_index) {
  return kv_A(sequence_batch->sequences, sequence_index)->qual.s;
}

void initialize_sequence_batch_loading(const char *sequence_file_path, SequenceBatch *sequence_batch);
void finalize_sequence_batch_loading(SequenceBatch *sequence_batch);
void initialize_sequence_batch(SequenceBatch *sequence_batch);
void initialize_sequence_batch_with_max_size(uint32_t max_num_sequences, SequenceBatch *sequence_batch);
void destory_sequence_batch(SequenceBatch *sequence_batch);
void load_batch_of_sequences_into_sequence_batch(SequenceBatch *sequence_batch);
void load_all_sequences_into_sequence_batch(SequenceBatch *sequence_batch);

//bool LoadOneSequenceAndSaveAt(uint32_t sequence_index);

//inline uint32_t GetSequenceNameLengthAt(uint32_t sequence_index) const {
//  return sequence_batch_[sequence_index]->name.l;
//}
//inline uint32_t GetSequenceIdAt(uint32_t sequence_index) const {
//  return sequence_batch_[sequence_index]->id;
//}

//  inline char GetReverseComplementBaseOfSequenceAt(uint32_t sequence_index, uint32_t position) {
//    kseq_t *sequence = sequence_batch_[sequence_index];
//    return Uint8ToChar(((uint8_t)3) ^ (CharToUint8((sequence->seq.s)[sequence->seq.l - position - 1])));
//  }
static inline void prepare_negative_sequence_at(uint32_t sequence_index, SequenceBatch *sequence_batch) {
  kseq_t *sequence = kv_A(sequence_batch->sequences, sequence_index);
  uint32_t sequence_length = sequence->seq.l;
  kv_clear(kv_A(sequence_batch->negative_sequences, sequence_index).v);
  //negative_sequence.reserve(sequence_length);
  for (uint32_t i = 0; i < sequence_length; ++i) {
    kv_push(char, kv_A(sequence_batch->negative_sequences, sequence_index).v, uint8_to_char(((uint8_t)3) ^ (char_to_uint8((sequence->seq.s)[sequence_length - i - 1]))));
  }
}
//inline void TrimSequenceAt(uint32_t sequence_index, int length_after_trim) {
//  kseq_t *sequence = sequence_batch_[sequence_index];
//  negative_sequence_batch_[sequence_index].erase(negative_sequence_batch_[sequence_index].begin(), negative_sequence_batch_[sequence_index].begin() + sequence->seq.l - length_after_trim);
//  sequence->seq.l = length_after_trim;
//}
//inline void SwapSequenceBatch(SequenceBatch &batch) {
//  sequence_batch_.swap(batch.GetSequenceBatch());
//  negative_sequence_batch_.swap(batch.GetNegativeSequenceBatch());
//}


//inline uint64_t GenerateSeedFromSequenceAt(uint32_t sequence_index, uint32_t start_position, uint32_t seed_length) {
//  const char *sequence = GetSequenceAt(sequence_index);
//  uint32_t sequence_length = GetSequenceLengthAt(sequence_index);
//  uint64_t mask = (((uint64_t)1) << (2 * seed_length)) - 1;
//  uint64_t seed = 0;
//  for (uint32_t i = 0; i < seed_length; ++i) {
//    if (start_position + i < sequence_length) {
//      uint8_t current_base = SequenceBatch::CharToUint8(sequence[i + start_position]);
//      if (current_base < 4) { // not an ambiguous base
//        seed = ((seed << 2) | current_base) & mask; // forward k-mer
//      } else {
//        seed = (seed << 2) & mask; // N->A
//      }
//    } else {
//      seed = (seed << 2) & mask; // Pad A
//    }
//  }
//  return seed;
//}

#endif // SEQUENCEBATCH_H_
