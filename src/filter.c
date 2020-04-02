#include "filter.h"

uint32_t generate_optimal_prefix_qgram_for_group_seeding(const FEMArgs *fem_args, const Index *index, int seed_length, int read_length, Seed *seeds, Seed *optimal_seeds) {
  uint32_t num_rows = fem_args->error_threshold + fem_args->num_additional_qgrams + 1 + 1;
  uint32_t num_columns = read_length - (fem_args->error_threshold + fem_args->num_additional_qgrams + 1) * seed_length + 1 + 1; // check if reduce d by one
  uint32_t M[num_rows][num_columns];
  uint32_t D[num_rows][num_columns]; // 3 for stop, 2 for vertical move and 1 for horizontal move
  for (uint32_t i = 1; i < num_rows; ++i) {
    M[i][0] = index->occurrence_table_size; 
    D[i][0] = 3;
  }
  for (uint32_t i = 1; i < num_columns; ++i) {
    M[0][i] = 0;
    D[0][i] = 3;
  }
  for (uint32_t row = 1; row < num_rows; ++row) {
    for (uint32_t column = 1; column < num_columns; ++column) {
      uint32_t position = column + (row - 1) * seed_length - 1;
      uint32_t num_total_occurrences_with_new_seed_at_position = M[row - 1][column] + seeds[position].num_positions;
      if (num_total_occurrences_with_new_seed_at_position < M[row][column - 1]) {
        M[row][column] = num_total_occurrences_with_new_seed_at_position;
        D[row][column] = 2;
      } else {
        M[row][column] = M[row][column - 1];
        D[row][column] = 1;
      }
    }
  }
  // Traceback
  uint32_t seed_row = num_rows - 1;
  uint32_t seed_column = num_columns - 1;
  int num_optimal_seeds = 0;
  while (D[seed_row][seed_column] != 3) {
    if (D[seed_row][seed_column] == 2) {
      optimal_seeds[num_optimal_seeds] = seeds[seed_column + (seed_row - 1) * seed_length - 1];
      ++num_optimal_seeds;
      --seed_row;
    } else if (D[seed_row][seed_column] == 1) {
      --seed_column;
    } 
  }
  return M[num_rows - 1][num_columns - 1];        
}

void merge_kvec_t_uint64_t(const FEMArgs *fem_args, kvec_t_uint64_t *buffer1, kvec_t_uint64_t *buffer2, kvec_t_uint64_t *candidates) {
  size_t buffer1_index = 0;
  size_t buffer2_index = 0;
  while (buffer1_index < kv_size(buffer1->v) || buffer2_index < kv_size(buffer2->v)) {
    if (buffer1_index < kv_size(buffer1->v)) {
      uint64_t buffer1_position = kv_A(buffer1->v, buffer1_index);
      if (buffer2_index < kv_size(buffer2->v)) {
        uint64_t buffer2_position = kv_A(buffer2->v, buffer2_index);
        if (buffer1_position < buffer2_position) {
          if (kv_size(candidates->v) == 0 || buffer1_position > kv_A(candidates->v, kv_size(candidates->v) - 1) + fem_args->error_threshold) {
            kv_push(uint64_t, candidates->v, buffer1_position);
          }
          ++buffer1_index;
        } else {
          if (kv_size(candidates->v) == 0 || buffer2_position > kv_A(candidates->v, kv_size(candidates->v) - 1) + fem_args->error_threshold) {
            kv_push(uint64_t, candidates->v, buffer2_position);
          }
          ++buffer2_index;
        }
      } else {
        if (kv_size(candidates->v) == 0 || buffer1_position > kv_A(candidates->v, kv_size(candidates->v) - 1) + fem_args->error_threshold) {
          kv_push(uint64_t, candidates->v, buffer1_position);
        }
        ++buffer1_index;
      }
    } else {
      uint64_t buffer2_position = kv_A(buffer2->v, buffer2_index);
      if (kv_size(candidates->v) == 0 || buffer2_position > kv_A(candidates->v, kv_size(candidates->v) - 1) + fem_args->error_threshold) {
        kv_push(uint64_t, candidates->v, buffer2_position);
      }
      ++buffer2_index;
    }
  }
}

void merge_candidate_locations(const FEMArgs *fem_args, const Index *index, const Seed *seeds, size_t num_seeds, kvec_t_uint64_t *buffer1, kvec_t_uint64_t *buffer2) {
  for (size_t si = 0; si < num_seeds; ++si) {
    size_t buffer1_index = 0;
    size_t seed_occurrence_index = 0;
    uint64_t *seed_occurrence_list = get_seed_occurrences(index, seeds[si].hash_value);
    while (buffer1_index < kv_size(buffer1->v) || (si != num_seeds - 1 && seed_occurrence_index < seeds[si].num_positions)) { // TODO: for the second case I have to push back one extra
      if (buffer1_index < kv_size(buffer1->v)) {
        uint64_t buffer1_position = kv_A(buffer1->v, buffer1_index);
        if (seed_occurrence_index < seeds[si].num_positions) {
          if ((uint32_t)seed_occurrence_list[seed_occurrence_index] < seeds[si].start_position) {
            ++seed_occurrence_index;
          } else {
            uint64_t seed_position = seed_occurrence_list[seed_occurrence_index] - seeds[si].start_position;
            if (seed_position <= buffer1_position) {
              kv_push(uint64_t, buffer2->v, seed_position);
              ++seed_occurrence_index;
            } else {
              kv_push(uint64_t, buffer2->v, buffer1_position);
              ++buffer1_index;
            }
          }
        } else {
          kv_push(uint64_t, buffer2->v, buffer1_position);
          ++buffer1_index;
        }
      } else {
        if ((uint32_t)seed_occurrence_list[seed_occurrence_index] >= seeds[si].start_position) {
          uint64_t seed_position = seed_occurrence_list[seed_occurrence_index] - seeds[si].start_position;
          kv_push(uint64_t, buffer2->v, seed_position);
        }
        ++seed_occurrence_index;
      }
    }
    kv_swap(uint64_t, buffer1->v, buffer2->v);
    kv_clear(buffer2->v);
  }
}

void additional_qgram_filter(const FEMArgs *fem_args, kvec_t_uint64_t *buffer, kvec_t_uint64_t *candidates) {
  for (size_t ci = 0; ci < kv_size(buffer->v); ++ci) {
    size_t num_candidates_in_range = 1;
    while (ci + num_candidates_in_range < kv_size(buffer->v) && kv_A(buffer->v, ci + num_candidates_in_range) <= kv_A(buffer->v, ci) + fem_args->error_threshold) {
      ++num_candidates_in_range;
      if (num_candidates_in_range > fem_args->num_additional_qgrams) { 
        break;
      }
    }
    if (num_candidates_in_range > fem_args->num_additional_qgrams) { 
      kv_push(uint64_t, candidates->v, kv_A(buffer->v, ci));
    }
  }
}

void remove_out_ranged_candidates(const FEMArgs *fem_args, uint32_t read_length, const SequenceBatch *reference_sequence_batch, kvec_t_uint64_t *buffer, kvec_t_uint64_t *candidates) {
  for (size_t i = 0; i < kv_size(buffer->v); ++i) {
    uint64_t candidate = kv_A(buffer->v, i);
    uint32_t reference_sequence_index = candidate >> 32;
    uint32_t reference_sequence_length = get_sequence_length_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index);
    uint32_t reference_candidate_position = (uint32_t)candidate;
    assert(reference_candidate_position < reference_sequence_length);
    if (reference_candidate_position >= (uint32_t)(fem_args->error_threshold) && reference_candidate_position + read_length + fem_args->error_threshold < reference_sequence_length) {
      kv_push(uint64_t, candidates->v, candidate - fem_args->error_threshold);
    }
  }
}

uint32_t generate_group_seeding_candidates(const FEMArgs *fem_args, const SequenceBatch *read_sequence_batch, size_t read_index, int is_reverse_complement, const SequenceBatch *reference_sequence_batch, const Index *index, kvec_t_uint64_t *buffer1, kvec_t_uint64_t *buffer2, kvec_t_uint64_t *candidates, uint32_t *num_candidates_without_additonal_qgram_filter) {
  //const uint8_t *bases = read->bases;
  //if (is_reverse_complement == 1) {
  //  bases = read->rc_bases;
  //}
  kv_clear(buffer1->v);
  kv_clear(buffer2->v);
  kv_clear(candidates->v);

  uint32_t read_length = get_sequence_length_from_sequence_batch_at(read_sequence_batch, read_index);
  const char *read_sequence = get_sequence_from_sequence_batch_at(read_sequence_batch, read_index);

  // Check if we can select enough seeds in the read
  int seed_length_in_seed_group = fem_args->kmer_size / fem_args->step_size;
  if (fem_args->kmer_size % fem_args->step_size > 0) {
    seed_length_in_seed_group++;
  }
  int num_seeds_in_read = (int)read_length - fem_args->kmer_size + 1;
  assert(num_seeds_in_read > 0);
  int min_num_seeds_in_seed_group = num_seeds_in_read / fem_args->step_size;
  if (fem_args->error_threshold + 1 + fem_args->num_additional_qgrams > min_num_seeds_in_seed_group) {
    // read is too short to be mapped
    return 0;
  }

  // dp for seed selection start
  // Generate seeds
  uint32_t seed_hash_values[num_seeds_in_read];
  //uint32_t seed_frequencies[num_seeds_in_read];
  int num_seeds_with_ambiguous_base = 0;
  hash_all_seeds_in_sequence(0, num_seeds_in_read, fem_args->kmer_size, read_sequence, read_length, &num_seeds_with_ambiguous_base, seed_hash_values);
  if (num_seeds_with_ambiguous_base > fem_args->error_threshold) {
    return 0;
  }
  //for (int si = 0; si < num_seeds_in_read; ++si) {
  //  seed_frequencies[si] = index->lookup_table[seed_hash_values[si] + 1] - lookup_table[seed_hash_values[si]];
  //}

  // Run seeding algorithm in each seed group
  *num_candidates_without_additonal_qgram_filter = 0;
  //uint32_t num_candidates = 0;
  for (int si = 0; si < fem_args->step_size; ++si) {
    // Generate optimal prefix q-gram
    int num_seeds_in_current_seed_group = (read_length - fem_args->kmer_size + 1 - si) / fem_args->step_size;
    Seed seeds_in_current_seed_group[num_seeds_in_current_seed_group];
    Seed optimal_seeds_in_current_seed_group[fem_args->error_threshold + 1 + fem_args->num_additional_qgrams];
    for (int k = 0; k < num_seeds_in_current_seed_group; ++k) {
      int seed_index_in_read = si + k * fem_args->step_size;
      seeds_in_current_seed_group[k].hash_value = seed_hash_values[seed_index_in_read];
      seeds_in_current_seed_group[k].start_position = seed_index_in_read;
      seeds_in_current_seed_group[k].end_position = seed_index_in_read + fem_args->kmer_size;
      seeds_in_current_seed_group[k].num_positions = get_seed_frequency(index, seeds_in_current_seed_group[k].hash_value);
    }
    *num_candidates_without_additonal_qgram_filter += generate_optimal_prefix_qgram_for_group_seeding(fem_args, index, seed_length_in_seed_group, num_seeds_in_current_seed_group, seeds_in_current_seed_group, optimal_seeds_in_current_seed_group);
    // Sort q-grams on their frequency
    qsort(optimal_seeds_in_current_seed_group, fem_args->error_threshold + 1 + fem_args->num_additional_qgrams, sizeof(Seed), compare_seed);
    // Filter seeds with additional q-gram
    kv_clear(buffer1->v);
    kv_clear(buffer2->v);
    merge_candidate_locations(fem_args, index, optimal_seeds_in_current_seed_group, fem_args->error_threshold + 1 + fem_args->num_additional_qgrams, buffer1, buffer2);
    additional_qgram_filter(fem_args, buffer1, buffer2);
    kv_swap(uint64_t, buffer1->v, candidates->v);
    kv_clear(candidates->v);
    merge_kvec_t_uint64_t(fem_args, buffer1, buffer2, candidates);
  }
  for (size_t i = 1; i < kv_size(candidates->v); ++i) {
    if (kv_A(candidates->v, i - 1) == kv_A(candidates->v, i)) {
      fprintf(stderr, "%lu\n", kv_A(candidates->v, i));
    }
  }
  kv_swap(uint64_t, buffer1->v, candidates->v);
  kv_clear(candidates->v);
  remove_out_ranged_candidates(fem_args, read_length, reference_sequence_batch, buffer1, candidates);
  return kv_size(candidates->v);
}
