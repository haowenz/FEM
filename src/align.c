#include "align.h"

int verify_candidates(const FEMArgs *fem_args, const SequenceBatch *read_sequence_batch, size_t read_sequence_index, int is_reverse_complement, const SequenceBatch *reference_sequence_batch, const uint64_t *candidates, uint32_t num_candidates, kstring_t *result_kstring) {
  //const uint8_t *text = read->bases;
  //if (is_reverse_complement == 1) {
  //  text = read->rc_bases;
  //}
  int read_length = get_sequence_length_from_sequence_batch_at(read_sequence_batch, read_sequence_index);
  const char *read_sequence = get_sequence_from_sequence_batch_at(read_sequence_batch, read_sequence_index);

  // Compute how many runs of vectorized code needed
  int num_mappings = 0;
  uint32_t num_vpus = num_candidates / NUM_VPU_LANES;
  uint32_t num_remains = num_candidates % NUM_VPU_LANES;
  //char cigar[read_length + 2 * fem_args->error_threshold];
  //char MD[read_length + 2 * fem_args->error_threshold];
  int16_t mapping_edit_distances[NUM_VPU_LANES];
  int16_t mapping_end_positions[NUM_VPU_LANES]; 

  int64_t previous_mapping_end_position = -read_length;
  for (uint32_t vpu_index = 0; vpu_index < num_vpus; ++vpu_index) {
    for (int li = 0; li < NUM_VPU_LANES; ++li){
      mapping_end_positions[li] = read_length - 1;
    }
    vectorized_banded_edit_distance(fem_args, vpu_index, reference_sequence_batch, read_sequence, read_length, candidates, num_candidates, mapping_edit_distances, mapping_end_positions);
    for (int mi = 0; mi < NUM_VPU_LANES; ++mi) {
      if (mapping_edit_distances[mi] <= fem_args->error_threshold) {
        int64_t mapping_end_position = candidates[vpu_index * NUM_VPU_LANES + mi] + mapping_end_positions[mi] + 1;
        if (previous_mapping_end_position == -read_length) {
          previous_mapping_end_position = mapping_end_position;
        } else if (mapping_end_position == previous_mapping_end_position) {
          continue;
        } else {
          previous_mapping_end_position = mapping_end_position;
        }
        ++num_mappings;
      }
    }
  }
  for (uint32_t ci = 0; ci < num_remains; ++ci) {
    uint64_t candidate = candidates[num_vpus * NUM_VPU_LANES + ci];
    uint32_t reference_sequence_index = candidate >> 32;
    uint32_t reference_candidate_position = (uint32_t)candidate;
    const char *reference_sequence = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index) + (uint32_t)reference_candidate_position;
    int current_mapping_end_position = -read_length;
    int current_mapping_edit_distance = banded_edit_distance(fem_args, reference_sequence, read_sequence, read_length, &current_mapping_end_position);
    if (current_mapping_edit_distance <= fem_args->error_threshold) {
      int64_t mapping_end_position = candidate + current_mapping_end_position + 1;
      if (previous_mapping_end_position == -read_length) {
        previous_mapping_end_position = mapping_end_position;
      } else if (mapping_end_position == previous_mapping_end_position) {
        continue;
      } else {
        previous_mapping_end_position = mapping_end_position;
      }
      ++num_mappings;
    }
  }
  return num_mappings;
}

int banded_edit_distance(const FEMArgs *fem_args, const char *pattern, const char *text, int read_length, int *mapping_end_position) {
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * fem_args->error_threshold; i++) {
    uint8_t base = char_to_uint8(pattern[i]);
    Peq[base] = Peq[base] | (1 << i);
  }
  uint32_t highest_bit_in_band_mask = 1 << (2 * fem_args->error_threshold);
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;
  int num_errors_at_band_start_position = 0;
  for (int i = 0; i < read_length; i++) {
    uint8_t pattern_base = char_to_uint8(pattern[i + 2 * fem_args->error_threshold]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[char_to_uint8(text[i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X & HP;
    VP = HN | ~(X | HP);
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    if (num_errors_at_band_start_position > 3 * fem_args->error_threshold) {
      return fem_args->error_threshold + 1;
    }
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }
  int band_start_position = read_length - 1;
  int min_num_errors = num_errors_at_band_start_position;
  *mapping_end_position = band_start_position;
  for (int i = 0; i < 2 * fem_args->error_threshold; i++) {
    num_errors_at_band_start_position = num_errors_at_band_start_position + ((VP >> i) & (uint32_t) 1);
    num_errors_at_band_start_position = num_errors_at_band_start_position - ((VN >> i) & (uint32_t) 1);
    if (num_errors_at_band_start_position < min_num_errors) {
      min_num_errors = num_errors_at_band_start_position;
      *mapping_end_position = band_start_position + 1 + i;
    }
  }
  return min_num_errors;
}

void vectorized_banded_edit_distance(const FEMArgs *fem_args, const uint32_t vpu_index, const SequenceBatch *reference_sequence_batch, const char *text, int read_length, const uint64_t *candidates, uint32_t num_candidates, int16_t *mapping_edit_distances, int16_t *mapping_end_positions) {
  uint32_t reference_sequence_index0 = candidates[vpu_index * NUM_VPU_LANES + 0] >> 32;
  uint32_t reference_sequence_index1 = candidates[vpu_index * NUM_VPU_LANES + 1] >> 32;
  uint32_t reference_sequence_index2 = candidates[vpu_index * NUM_VPU_LANES + 2] >> 32;
  uint32_t reference_sequence_index3 = candidates[vpu_index * NUM_VPU_LANES + 3] >> 32;
  uint32_t reference_sequence_index4 = candidates[vpu_index * NUM_VPU_LANES + 4] >> 32;
  uint32_t reference_sequence_index5 = candidates[vpu_index * NUM_VPU_LANES + 5] >> 32;
  uint32_t reference_sequence_index6 = candidates[vpu_index * NUM_VPU_LANES + 6] >> 32;
  uint32_t reference_sequence_index7 = candidates[vpu_index * NUM_VPU_LANES + 7] >> 32;
  const char *reference_sequence0 = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index0) + (uint32_t)candidates[vpu_index * NUM_VPU_LANES + 0];
  const char *reference_sequence1 = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index1) + (uint32_t)candidates[vpu_index * NUM_VPU_LANES + 1];
  const char *reference_sequence2 = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index2) + (uint32_t)candidates[vpu_index * NUM_VPU_LANES + 2];
  const char *reference_sequence3 = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index3) + (uint32_t)candidates[vpu_index * NUM_VPU_LANES + 3];
  const char *reference_sequence4 = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index4) + (uint32_t)candidates[vpu_index * NUM_VPU_LANES + 4];
  const char *reference_sequence5 = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index5) + (uint32_t)candidates[vpu_index * NUM_VPU_LANES + 5];
  const char *reference_sequence6 = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index6) + (uint32_t)candidates[vpu_index * NUM_VPU_LANES + 6];
  const char *reference_sequence7 = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index7) + (uint32_t)candidates[vpu_index * NUM_VPU_LANES + 7];
  uint16_t highest_bit_in_band_mask = 1 << (2 * fem_args->error_threshold);
  __m128i highest_bit_in_band_mask_vpu0 = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, highest_bit_in_band_mask);
  __m128i highest_bit_in_band_mask_vpu1 = _mm_set_epi16(0, 0, 0, 0, 0, 0, highest_bit_in_band_mask, 0);
  __m128i highest_bit_in_band_mask_vpu2 = _mm_set_epi16(0, 0, 0, 0, 0, highest_bit_in_band_mask, 0, 0);
  __m128i highest_bit_in_band_mask_vpu3 = _mm_set_epi16(0, 0, 0, 0, highest_bit_in_band_mask, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu4 = _mm_set_epi16(0, 0, 0, highest_bit_in_band_mask, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu5 = _mm_set_epi16(0, 0, highest_bit_in_band_mask, 0, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu6 = _mm_set_epi16(0, highest_bit_in_band_mask, 0, 0, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu7 = _mm_set_epi16(highest_bit_in_band_mask, 0, 0, 0, 0, 0, 0, 0);
  // Init Peq
  __m128i Peq[ALPHABET_SIZE];
  for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
    Peq[ai] = _mm_setzero_si128();
  }
  for (int i = 0; i < 2 * fem_args->error_threshold; i++) {
    uint8_t base0 = char_to_uint8(reference_sequence0[i]);
    uint8_t base1 = char_to_uint8(reference_sequence1[i]);
    uint8_t base2 = char_to_uint8(reference_sequence2[i]);
    uint8_t base3 = char_to_uint8(reference_sequence3[i]);
    uint8_t base4 = char_to_uint8(reference_sequence4[i]);
    uint8_t base5 = char_to_uint8(reference_sequence5[i]);
    uint8_t base6 = char_to_uint8(reference_sequence6[i]);
    uint8_t base7 = char_to_uint8(reference_sequence7[i]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    Peq[base4] = _mm_or_si128(highest_bit_in_band_mask_vpu4, Peq[base4]);
    Peq[base5] = _mm_or_si128(highest_bit_in_band_mask_vpu5, Peq[base5]);
    Peq[base6] = _mm_or_si128(highest_bit_in_band_mask_vpu6, Peq[base6]);
    Peq[base7] = _mm_or_si128(highest_bit_in_band_mask_vpu7, Peq[base7]);
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi16(Peq[ai], 1);
    }
  }

  uint16_t lowest_bit_in_band_mask = 1;
  __m128i lowest_bit_in_band_mask_vpu = _mm_set1_epi16(lowest_bit_in_band_mask);
  __m128i VP = _mm_setzero_si128();
  __m128i VN =  _mm_setzero_si128();
  __m128i X = _mm_setzero_si128();
  __m128i D0 = _mm_setzero_si128();
  __m128i HN = _mm_setzero_si128();
  __m128i HP = _mm_setzero_si128();
  __m128i max_mask_vpu = _mm_set1_epi16(0xffff);
  __m128i num_errors_at_band_start_position_vpu = _mm_setzero_si128();
  __m128i early_stop_threshold_vpu = _mm_set1_epi16(fem_args->error_threshold * 3);
  for (int i = 0; i < read_length; i++) {
    uint8_t base0 = char_to_uint8(reference_sequence0[i + 2 * fem_args->error_threshold]);
    uint8_t base1 = char_to_uint8(reference_sequence1[i + 2 * fem_args->error_threshold]);
    uint8_t base2 = char_to_uint8(reference_sequence2[i + 2 * fem_args->error_threshold]);
    uint8_t base3 = char_to_uint8(reference_sequence3[i + 2 * fem_args->error_threshold]);
    uint8_t base4 = char_to_uint8(reference_sequence4[i + 2 * fem_args->error_threshold]);
    uint8_t base5 = char_to_uint8(reference_sequence5[i + 2 * fem_args->error_threshold]);
    uint8_t base6 = char_to_uint8(reference_sequence6[i + 2 * fem_args->error_threshold]);
    uint8_t base7 = char_to_uint8(reference_sequence7[i + 2 * fem_args->error_threshold]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    Peq[base4] = _mm_or_si128(highest_bit_in_band_mask_vpu4, Peq[base4]);
    Peq[base5] = _mm_or_si128(highest_bit_in_band_mask_vpu5, Peq[base5]);
    Peq[base6] = _mm_or_si128(highest_bit_in_band_mask_vpu6, Peq[base6]);
    Peq[base7] = _mm_or_si128(highest_bit_in_band_mask_vpu7, Peq[base7]);
    X = _mm_or_si128(Peq[char_to_uint8(text[i])], VN);
    D0 = _mm_and_si128(X, VP);
    D0 = _mm_add_epi16(D0, VP);
    D0 = _mm_xor_si128(D0, VP);
    D0 = _mm_or_si128(D0, X);
    HN = _mm_and_si128(VP, D0);
    HP = _mm_or_si128(VP, D0);
    HP = _mm_xor_si128(HP, max_mask_vpu);
    HP = _mm_or_si128(HP, VN);
    X = _mm_srli_epi16(D0, 1);
    VN = _mm_and_si128(X, HP);
    VP = _mm_or_si128(X, HP);
    VP = _mm_xor_si128(VP, max_mask_vpu);
    VP = _mm_or_si128(VP, HN);
    __m128i E = _mm_and_si128(D0, lowest_bit_in_band_mask_vpu);
    E = _mm_xor_si128(E, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu = _mm_add_epi16(num_errors_at_band_start_position_vpu, E);
    __m128i early_stop = _mm_cmpgt_epi16(num_errors_at_band_start_position_vpu, early_stop_threshold_vpu);
    int tmp = _mm_movemask_epi8(early_stop);
    if (tmp == 0xffff) {
      _mm_store_si128((__m128i *)mapping_edit_distances, num_errors_at_band_start_position_vpu);
      return;
    }
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi16(Peq[ai], 1);
    }
  }
  int band_start_position = read_length - 1;
  __m128i min_num_errors_vpu = num_errors_at_band_start_position_vpu;
  for (int i = 0; i < 2 * fem_args->error_threshold; i++) {
    __m128i lowest_bit_in_VP_vpu = _mm_and_si128(VP, lowest_bit_in_band_mask_vpu);
    __m128i lowest_bit_in_VN_vpu = _mm_and_si128(VN, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu = _mm_add_epi16(num_errors_at_band_start_position_vpu, lowest_bit_in_VP_vpu);
    num_errors_at_band_start_position_vpu = _mm_sub_epi16(num_errors_at_band_start_position_vpu, lowest_bit_in_VN_vpu);
    __m128i mapping_end_positions_update_mask_vpu = _mm_cmplt_epi16(num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    int mapping_end_positions_update_mask = _mm_movemask_epi8(mapping_end_positions_update_mask_vpu);
    for (int li = 0; li < NUM_VPU_LANES; ++li) {
      if ((mapping_end_positions_update_mask & 1) == 1) {
        mapping_end_positions[li] = band_start_position + 1 + i;
      }
      mapping_end_positions_update_mask = mapping_end_positions_update_mask >> 2;
    }
    min_num_errors_vpu = _mm_min_epi16(min_num_errors_vpu, num_errors_at_band_start_position_vpu);
    VP = _mm_srli_epi16(VP, 1);
    VN = _mm_srli_epi16(VN, 1);
  }
  _mm_store_si128((__m128i *)mapping_edit_distances, min_num_errors_vpu);
}
