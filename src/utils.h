#ifndef UTILS_H_
#define UTILS_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <mm_malloc.h>

#include "ksort.h"
#include "kvec.h"

typedef struct {
  kvec_t(uint64_t) v;
} kvec_t_uint64_t;

typedef struct {
  kvec_t(uint32_t) v;
} kvec_t_uint32_t;

typedef struct {
  kvec_t(uint8_t) v;
} kvec_t_uint8_t;

typedef struct {
  uint8_t edit_distance;
  uint64_t candidate_position;
  int16_t end_position_offset; // end_postion = candiate_position + end_position_offset
  //uint8_t mapq : 6, direction : 1, is_unique : 1;
} Mapping;

#define MappingSortKey(m) ((((uint64_t)(m).edit_distance)<<60)|((m).candidate_position+(m).end_position_offset))
KRADIX_SORT_INIT(mapping, Mapping, MappingSortKey, 8);

//typedef struct {
//  uint32_t read_index;
//  uint8_t edit_distance;
//  uint32_t candidate_position;
//  int16_t end_position_offset; // end_postion = candiate_position + end_position_offset
//  //uint8_t mapq : 6, direction : 1, is_unique : 1;
//} Mapping;

typedef struct {
  kvec_t(Mapping) v;
} kvec_t_Mapping;

typedef struct {
  uint64_t num_reads;
  uint64_t num_mapped_reads;
  uint64_t num_candidates_without_additonal_qgram_filter;
  uint64_t num_candidates;
  uint64_t num_mappings;
} MappingStats;

typedef struct {
  int kmer_size;
  int step_size;
  int error_threshold;
  int num_additional_qgrams;
  int num_threads;
  char seeding_method; // "v" for variable length seeding, "g" for group seeding.
} FEMArgs;

static const uint8_t char_to_uint8_table[256] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
static const char uint8_to_char_table[8] = {'A', 'C', 'G', 'T', 'N', 'N', 'N', 'N'};

static inline uint8_t char_to_uint8(char c) {
  return char_to_uint8_table[(uint8_t)c];
}

static inline char uint8_to_char(uint8_t i) {
  return uint8_to_char_table[i];
}

static inline uint32_t hash_seed_in_sequence(size_t seed_start_position, int seed_length, const char *sequence, size_t sequence_length) {
  uint32_t mask = (((uint32_t)1) << (2 * seed_length)) - 1;
  uint32_t hash_value = 0;
  for (uint32_t i = 0; i < seed_length; ++i) {
    if (seed_start_position + i < sequence_length) {
      uint8_t current_base = char_to_uint8(sequence[i + seed_start_position]);
      if (current_base < 4) { // not an ambiguous base
        hash_value = ((hash_value << 2) | current_base) & mask; // forward k-mer
      } else {
        hash_value = (hash_value << 2) & mask; // N->A
      }
    } else {
      hash_value = (hash_value << 2) & mask; // Pad A
    }
  }
  return hash_value;
}

static inline void hash_all_seeds_in_sequence(size_t seed_start_position, size_t seed_end_position, int seed_length, const char *sequence, size_t sequence_length, int *num_seeds_with_ambiguous_base, uint32_t *seed_hash_values) {
  assert(seed_end_position <= sequence_length);
  *num_seeds_with_ambiguous_base = 0;
  uint32_t mask = (((uint32_t)1) << (2 * seed_length)) - 1;
  uint32_t hash_value = hash_seed_in_sequence(seed_start_position, seed_length, sequence, sequence_length);
  seed_hash_values[0] = hash_value;
  for (size_t i = seed_start_position + 1; i < seed_end_position; ++i) {
    uint8_t current_base = char_to_uint8(sequence[i + seed_length - 1]);
    if (current_base < 4) { // not an ambiguous base
      hash_value = ((hash_value << 2) | current_base) & mask; // forward k-mer
    } else {
      hash_value = (hash_value << 2) & mask; // N->A
      ++(*num_seeds_with_ambiguous_base);
    }
    seed_hash_values[i - seed_start_position] = hash_value;
  }
}

typedef struct {
  uint32_t hash_value;
  uint32_t start_position;
  uint32_t end_position;
  uint32_t num_positions;
} Seed;

static inline int compare_seed(const void *ta, const void *tb) {
  Seed my_ta = *(Seed*) ta;
  Seed my_tb = *(Seed*) tb;
  if (my_ta.num_positions < my_tb.num_positions) {
    return -1;
  } else if (my_ta.num_positions == my_tb.num_positions) {
    return 0;
  } else {
    return 1;
  }
}

static inline double get_cpu_time() {
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static inline double get_real_time() {
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

#endif // UTILS_H_
