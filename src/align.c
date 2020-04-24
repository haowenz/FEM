#include "align.h"
//#include "kstring.h"

uint32_t verify_candidates(const FEMArgs *fem_args, OutputQueue *output_queue, const SequenceBatch *read_sequence_batch, uint32_t read_sequence_index, int is_reverse_complement, const SequenceBatch *reference_sequence_batch, const uint64_t *candidates, uint32_t num_candidates, kvec_t_Mapping *mappings) {
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
  int16_t mapping_edit_distances[NUM_VPU_LANES];
  int16_t mapping_end_positions[NUM_VPU_LANES]; 
  for (uint32_t vpu_index = 0; vpu_index < num_vpus; ++vpu_index) {
    for (int li = 0; li < NUM_VPU_LANES; ++li){
      mapping_end_positions[li] = read_length - 1;
    }
    vectorized_banded_edit_distance(fem_args, vpu_index, reference_sequence_batch, read_sequence, read_length, candidates, num_candidates, mapping_edit_distances, mapping_end_positions);
    for (int mi = 0; mi < NUM_VPU_LANES; ++mi) {
      if (mapping_edit_distances[mi] <= fem_args->error_threshold) {
        Mapping mapping = {(uint8_t)mapping_edit_distances[mi], candidates[vpu_index * NUM_VPU_LANES + mi], mapping_end_positions[mi]};
        kv_push(Mapping, mappings->v, mapping);
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
      Mapping mapping = {current_mapping_edit_distance, candidate, current_mapping_end_position};
      kv_push(Mapping, mappings->v, mapping);
      ++num_mappings;
    }
  }
  return num_mappings;
}

uint32_t process_mappings(const FEMArgs *fem_args, OutputQueue *output_queue, const SequenceBatch *read_sequence_batch, uint32_t read_sequence_index, int is_reverse_complement, const SequenceBatch *reference_sequence_batch, Mapping *mappings, uint32_t num_mappings, bam1_t **sam_alignment) {
  radix_sort_mapping(mappings, mappings + num_mappings);
  kvec_t_uint32_t cigar_uint32_t;
  kv_init(cigar_uint32_t.v);
  const char *read_sequence = get_sequence_from_sequence_batch_at(read_sequence_batch, read_sequence_index);
  const char *read_qual = get_sequence_qual_from_sequence_batch_at(read_sequence_batch, read_sequence_index);
  int read_length = get_sequence_length_from_sequence_batch_at(read_sequence_batch, read_sequence_index);
  const char *read_name = get_sequence_name_from_sequence_batch_at(read_sequence_batch, read_sequence_index);
  int read_name_length = get_sequence_name_length_from_sequence_batch_at(read_sequence_batch, read_sequence_index);
  num_mappings = 1;
  for (uint32_t mi = 0; mi < num_mappings; ++mi) {
    uint8_t edit_distance = mappings[mi].edit_distance;
    uint64_t candidate_position = mappings[mi].candidate_position;
    uint32_t reference_sequence_index = candidate_position >> 32;
    const char *reference_sequence = get_sequence_from_sequence_batch_at(reference_sequence_batch, reference_sequence_index) + (uint32_t)candidate_position;
    kv_clear(cigar_uint32_t.v);
    int mapping_start_position = generate_alignment(fem_args, reference_sequence, read_sequence, read_length, mappings[mi].edit_distance, mappings[mi].end_position_offset, NULL, &cigar_uint32_t);
    mapping_start_position += (uint32_t)candidate_position;
    uint8_t mapping_quality = 255;
    uint16_t flag = 0;
    if (mi > 0) {
      flag |= BAM_FSECONDARY;
      generate_bam1_t(edit_distance, mapping_start_position, reference_sequence_index, mapping_quality, flag, read_name, read_name_length, cigar_uint32_t.v.a, kv_size(cigar_uint32_t.v), read_sequence, read_qual, 0, *sam_alignment);
    } else {
      generate_bam1_t(edit_distance, mapping_start_position, reference_sequence_index, mapping_quality, flag, read_name, read_name_length, cigar_uint32_t.v.a, kv_size(cigar_uint32_t.v), read_sequence, read_qual, read_length, *sam_alignment);
    }
    push_output_queue(sam_alignment, output_queue);
  }
  kv_destroy(cigar_uint32_t.v);
  return num_mappings;
}

// Banded Myers bit-parallel algorithm
/* @param fem_args             FEM parameters 
   @param pattern              Reference sequence starting from (candidate position - error threshold) 
   @param text                 Read sequence
   @param text_length          Read length
   @param mapping_end_position 0-based mapping end position (inclusive) 
   @return Edit distance of the mapping.
   */
int banded_edit_distance(const FEMArgs *fem_args, const char *pattern, const char *text, int text_length, int *mapping_end_position) {
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
  for (int i = 0; i < text_length; i++) {
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
  int band_start_position = text_length - 1;
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

int generate_alignment(const FEMArgs *fem_args, const char *pattern, const char *text, int read_length, int mapping_edit_distance, int mapping_end_position, kstring_t *cigar, kvec_t_uint32_t *cigar_uint32_t) {
  // Note that we do a semi-global alignemnt, that is, errors at two ends of ref are not penalized and read is aligned globally
  // Also note that cigar operations are on ref 
  // M/I/S/=/X operations shall equal the length of SEQ
  // When using ED, "ID" or "DI" won't happen

  int mapping_start_position = mapping_end_position - read_length + 1;
  assert(mapping_start_position >= 0);
  int num_errors = 0;
  // Check if there are only substitutions
  for (int i = 0; i < read_length; i++) {
    if (text[i] != pattern[mapping_start_position + i]) {
      ++num_errors;
    }
  }
  if (num_errors == mapping_edit_distance) {
    //ksprintf(cigar, "%d%c", read_length, 'M');
    uint32_t cigar_uint = read_length << 4;
    kv_push(uint32_t, cigar_uint32_t->v, cigar_uint | 0);
    return mapping_start_position;
  }

  // Alignment traceback
  uint32_t D0s[read_length];
  uint32_t HPs[read_length];
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
  //int num_errors_at_band_start_position = 0;
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
    D0s[i] = D0;
    HPs[i] = HP;
    //num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    //if (num_errors_at_band_start_position > 3 * fem_args->error_threshold) {
    //  return fem_args->error_threshold + 1;
    //}
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }

  int pattern_bit_position = mapping_end_position - read_length + 1; // position of ending bit in bit vector 
  int text_position = read_length - 1; // start from the read end
  char pre_operation = 'S';
  int pre_num_operations = 1;
  num_errors = 0; // # errors in alignment (including soft clip)
  if (((D0s[text_position] >> pattern_bit_position) & lowest_bit_in_band_mask) && (pattern[mapping_end_position] == text[text_position])) { // match
    --text_position;
    --mapping_end_position;
    pre_operation = 'M';
    pre_num_operations = 1;
  } else if (!((D0s[text_position] >> pattern_bit_position) & lowest_bit_in_band_mask)) { // mismatch
    assert(pattern[mapping_end_position] != text[text_position]);
    --text_position;
    --mapping_end_position;
    ++num_errors;
    pre_operation = 'S';
    //pre_operation = 'M';
    pre_num_operations = 1;
  } else if (((D0s[text_position] >> pattern_bit_position) & lowest_bit_in_band_mask) && ((HPs[text_position] >> pattern_bit_position) & lowest_bit_in_band_mask)) { // insertion
    --text_position;
    ++pattern_bit_position;
    ++num_errors;
    pre_operation = 'S';
    //pre_operation = 'I';
    pre_num_operations = 1;
  } else { // deletion
    assert(1 == 0);
  }

  int cigar_operation_index = 0;
  char cigar_operations[read_length];
  int num_cigar_operations[read_length];
  while (text_position >= 0) {
    if (num_errors == mapping_edit_distance) {
      break;
    }
    if (((D0s[text_position] >> pattern_bit_position) & lowest_bit_in_band_mask) && (pattern[mapping_end_position] == text[text_position])) { // match: consume one base from both target and query
      --text_position;
      --mapping_end_position;
      //if (pre_operation == 'S') {
      //  pre_operation = 'M';
      //  ++pre_num_operations;
      //} else 
      if (pre_operation != 'M') {
        cigar_operations[cigar_operation_index] = pre_operation;
        num_cigar_operations[cigar_operation_index] = pre_num_operations;
        ++cigar_operation_index;
        pre_operation = 'M';
        pre_num_operations = 1;
      } else {
        ++pre_num_operations;
      }
    } else if (!((D0s[text_position] >> pattern_bit_position) & lowest_bit_in_band_mask)) { // mismatch
      assert(pattern[mapping_end_position] != text[text_position]);
      --text_position;
      --mapping_end_position;
      ++num_errors;
      if (pre_operation == 'S') {
        ++pre_num_operations;
      } else if (pre_operation != 'M') {
        cigar_operations[cigar_operation_index] = pre_operation;
        num_cigar_operations[cigar_operation_index] = pre_num_operations;
        ++cigar_operation_index;
        pre_operation = 'M';
        pre_num_operations = 1;
      } else {
        ++pre_num_operations;
      }
    } else if (((D0s[text_position] >> pattern_bit_position) & lowest_bit_in_band_mask) && (HPs[text_position] >> pattern_bit_position) & lowest_bit_in_band_mask) { // Insertion: consume one base in query
      --text_position;
      ++pattern_bit_position;
      ++num_errors;
      if (pre_operation == 'S') {
        ++pre_num_operations;
      } else if (pre_operation != 'I') {
        cigar_operations[cigar_operation_index] = pre_operation;
        num_cigar_operations[cigar_operation_index] = pre_num_operations;
        ++cigar_operation_index;
        pre_operation = 'I';
        pre_num_operations = 1;
      } else {
        ++pre_num_operations;
      }
      ++mapping_start_position;
    } else { // Deletion: consume one base in ref
      --pattern_bit_position;
      --mapping_end_position;
      ++num_errors;
      if (pre_operation != 'D') {
        cigar_operations[cigar_operation_index] = pre_operation;
        num_cigar_operations[cigar_operation_index] = pre_num_operations;
        ++cigar_operation_index;
        pre_operation = 'D';
        pre_num_operations = 1;
      } else {
        ++pre_num_operations;
      }
      --mapping_start_position;
    }
  }
  // After all errors are consumed, the rest must be matches
  //cigar_operations[cigar_operation_index] = pre_operation;
  //num_cigar_operations[cigar_operation_index] = pre_num_operations;
  //++cigar_operation_index;
  if (text_position >= 0) {
    if (pre_operation != 'M') {
      cigar_operations[cigar_operation_index] = pre_operation;
      num_cigar_operations[cigar_operation_index] = pre_num_operations;
      ++cigar_operation_index;
      cigar_operations[cigar_operation_index] = 'M';
      num_cigar_operations[cigar_operation_index] = text_position + 1;
    } else {
      cigar_operations[cigar_operation_index] = 'M';
      num_cigar_operations[cigar_operation_index] = pre_num_operations + text_position + 1;
    }
  } else {
    cigar_operations[cigar_operation_index] = pre_operation;
    num_cigar_operations[cigar_operation_index] = pre_num_operations;
  }

  //CIGAR data is stored as in the BAM format, i.e. (op_len << 4) | op
  //where op_len is the length in bases and op is a value between 0 and 8
  //representing one of the operations "MIDNSHP=X" (M = 0; X = 8)
  //int size_SM = 0;
  int cigar_operation_index_end = 0;
  if (cigar_operations[0] == 'S') {
    num_cigar_operations[1] += num_cigar_operations[0]; 
    cigar_operation_index_end = 1;
  }
  for (int i = cigar_operation_index; i >= cigar_operation_index_end; i--) {
    //if (Route_Char_Whole[j] == 'M' || Route_Char_Whole[j] == 'S') {
    //  size_SM = 0;
    //  while (j > 0 && (Route_Char_Whole[j] == 'M' || Route_Char_Whole[j] == 'S')) {
    //    size_SM = size_SM + Route_Size_Whole[j];
    //    j--;
    //  }
    //  j++;
    //  ksprintf(cigar, "%d%c", size_SM, 'M');
    //} else {
    //ksprintf(cigar, "%d%c", num_cigar_operations[i], cigar_operations[i]);
    //fprintf(stderr, "%d%c", num_cigar_operations[i], cigar_operations[i]);
    uint32_t cigar_uint = num_cigar_operations[i] << 4;
    if (cigar_operations[i] == 'M') {
      kv_push(uint32_t, cigar_uint32_t->v, cigar_uint | 0);
    } else if (cigar_operations[i] == 'I') {
      kv_push(uint32_t, cigar_uint32_t->v, cigar_uint | 1);
    } else if (cigar_operations[i] == 'D') {
      kv_push(uint32_t, cigar_uint32_t->v, cigar_uint | 2);
    } else if (cigar_operations[i] == 'N') {
      assert(1 == 0);
    } else if (cigar_operations[i] == 'S') {
      assert(1 == 0);
    } else if (cigar_operations[i] == 'H') {
      assert(1 == 0);
    } else if (cigar_operations[i] == 'P') {
      assert(1 == 0);
    } else if (cigar_operations[i] == '=') {
      assert(1 == 0);
    } else if (cigar_operations[i] == 'X') {
      assert(1 == 0);
    } else {
      //fprintf(stderr, "%d %d %c %d\n", cigar_operation_index, i, cigar_operations[i], num_cigar_operations[i]);
      assert(1 == 0);
    }
    //}
  }
  //fprintf(stderr, "\n");
  return mapping_start_position;
}

void generate_bam1_t(uint8_t edit_distance, uint32_t mapping_start_position, int32_t reference_sequence_index, uint8_t mapping_quality, uint16_t flag, const char *query_name, uint16_t query_name_length, uint32_t *cigar, uint32_t num_cigar_operations, const char *query, const char *query_qual, int32_t query_length, bam1_t *sam_alignment) {
  /*! @typedef
   *  @abstract Structure for core alignment information.
   *  @field  pos     0-based leftmost coordinate
   *  @field  tid     chromosome ID, defined by sam_hdr_t
   *  @field  bin     bin calculated by bam_reg2bin()
   *  @field  qual    mapping quality
   *  @field  l_extranul length of extra NULs between qname & cigar (for alignment)
   *  @field  flag    bitwise flag
   *  @field  l_qname length of the query name
   *  @field  n_cigar number of CIGAR operations
   *  @field  l_qseq  length of the query sequence (read)
   *  @field  mtid    chromosome ID of next read in template, defined by sam_hdr_t
   *  @field  mpos    0-based leftmost coordinate of next read in template
   *  @field  isize   observed template length ("insert size")
   */
  sam_alignment->core.pos = mapping_start_position;
  sam_alignment->core.tid = reference_sequence_index;
  sam_alignment->core.qual = mapping_quality;
  sam_alignment->core.l_extranul = 0; 
  if ((query_name_length + 1) % 4 != 0) {
    sam_alignment->core.l_extranul = 4 - ((query_name_length + 1) % 4);
  }
  sam_alignment->core.flag = flag;
  sam_alignment->core.l_qname = query_name_length + 1 + sam_alignment->core.l_extranul;
  sam_alignment->core.n_cigar = num_cigar_operations;
  sam_alignment->core.l_qseq = query_length;
  sam_alignment->core.mtid = -1;
  sam_alignment->core.mpos = -1;
  
  /*! @typedef
   @abstract Structure for one alignment.
   @field  core       core information about the alignment
   @field  id
   @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
   @field  l_data     current length of bam1_t::data
   @field  m_data     maximum length of bam1_t::data
   @field  mempolicy  memory handling policy, see bam_set_mempolicy()
   @discussion Notes:
   1. The data blob should be accessed using bam_get_qname, bam_get_cigar,
   bam_get_seq, bam_get_qual and bam_get_aux macros.  These returns pointers
   to the start of each type of data.
   2. qname is terminated by one to four NULs, so that the following
   cigar data is 32-bit aligned; core.l_qname includes these trailing NULs,
   while core.l_extranul counts the excess NULs (so 0 <= l_extranul <= 3).
   3. Cigar data is encoded 4 bytes per CIGAR operation.
   See the bam_cigar_* macros for manipulation.
   4. seq is nibble-encoded according to bam_nt16_table.
   See the bam_seqi macro for retrieving individual bases.
   5. Per base qualilties are stored in the Phred scale with no +33 offset.
   Ie as per the BAM specification and not the SAM ASCII printable method.
   */
  /*
   @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
   8 for T and 15 for N. Two bases are packed in one byte with the base
   at the higher 4 bits having smaller coordinate on the read. It is
   recommended to use bam_seqi() macro to get the base.
  */
  // First calculate the length of data
  sam_alignment->l_data = sam_alignment->core.l_qname + (sam_alignment->core.n_cigar << 2) + ((sam_alignment->core.l_qseq + 1) >> 1) + sam_alignment->core.l_qseq;
  if (sam_alignment->l_data > sam_alignment->m_data) {
    sam_alignment->m_data = sam_alignment->l_data;
    free(sam_alignment->data);
    sam_alignment->data = (uint8_t*)calloc(sam_alignment->m_data, sizeof(uint8_t));
  }
  // copy qname
  memcpy(bam_get_qname(sam_alignment), query_name, query_name_length * sizeof(char));
  // add NULs after qname and before cigar, let me add one nul at the moment and see if okay.
  bam_get_qname(sam_alignment)[query_name_length] = '\0';
  // copy cigar
  memcpy(bam_get_cigar(sam_alignment), cigar, sam_alignment->core.n_cigar * sizeof(uint32_t));
  // set seq
  uint8_t *seq = bam_get_seq(sam_alignment);
  for (size_t i = 0; i < query_length; ++i) {
    bam_set_seqi(seq, i, seq_nt16_table[(uint8_t)query[i]]);
  }
  // copy seq qual
  uint8_t *seq_qual = bam_get_qual(sam_alignment);
  memcpy(seq_qual, query_qual, query_length * sizeof(char));
  // remove +33 offset
  for (int i = 0; i < query_length; ++i) {
    seq_qual[i] -= 33;
  }
  bam_aux_update_int(sam_alignment, "NM", edit_distance);
}
