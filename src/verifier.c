#include "verifier.h"

inline void permute(__m128i *patterns) {
  __m128i tp[V_CPU];
  tp[0] = _mm_unpacklo_epi8(patterns[0], patterns[1]);
  tp[1] = _mm_unpackhi_epi8(patterns[0], patterns[1]);
  tp[2] = _mm_unpacklo_epi8(patterns[2], patterns[3]);
  tp[3] = _mm_unpackhi_epi8(patterns[2], patterns[3]);
  tp[4] = _mm_unpacklo_epi8(patterns[4], patterns[5]);
  tp[5] = _mm_unpackhi_epi8(patterns[4], patterns[5]);
  tp[6] = _mm_unpacklo_epi8(patterns[6], patterns[7]);
  tp[7] = _mm_unpackhi_epi8(patterns[6], patterns[7]);
  patterns[0] = _mm_unpacklo_epi16(tp[0], tp[2]);
  patterns[1] = _mm_unpacklo_epi16(tp[1], tp[3]);
  patterns[2] = _mm_unpackhi_epi16(tp[0], tp[2]);
  patterns[3] = _mm_unpackhi_epi16(tp[1], tp[3]);
  patterns[4] = _mm_unpacklo_epi16(tp[4], tp[6]);
  patterns[5] = _mm_unpacklo_epi16(tp[5], tp[7]);
  patterns[6] = _mm_unpackhi_epi16(tp[4], tp[6]);
  patterns[7] = _mm_unpackhi_epi16(tp[5], tp[7]);
  tp[0] = _mm_unpacklo_epi32(patterns[0], patterns[4]);
  tp[1] = _mm_unpackhi_epi32(patterns[0], patterns[4]);
  tp[2] = _mm_unpacklo_epi32(patterns[2], patterns[6]);
  tp[3] = _mm_unpackhi_epi32(patterns[2], patterns[6]);
  tp[4] = _mm_unpacklo_epi32(patterns[1], patterns[5]);
  tp[5] = _mm_unpackhi_epi32(patterns[1], patterns[5]);
  tp[6] = _mm_unpacklo_epi32(patterns[3], patterns[7]);
  tp[7] = _mm_unpackhi_epi32(patterns[3], patterns[7]);
  for (int pi = 0; pi < V_CPU; ++pi) {
    _mm_storeu_si128(patterns + pi, tp[pi]);
  }
}

int verify_candidates(const Read *read, const int is_reverse_complement, const uint32_t *candidates, const uint32_t num_candidates, int *no_left_result, uint32_t *result_string_length, char **result_string) {
  const uint8_t *text = read->bases;
  if (is_reverse_complement == 1) {
    text = read->rc_bases;
  }
  int num_mappings = 0;

  uint32_t num_registers = num_candidates / V_CPU;
  uint32_t remains = num_candidates % V_CPU;
  if (remains != 0) {
    ++num_registers;
  }

  char cigar[READ_LEN_MAX * 2];
  char MD[READ_LEN_MAX * 2];
  int16_t locations[V_CPU]; 
  int16_t errs[V_CPU];
  uint32_t last_location = 0xffffffff;
  for (uint32_t register_index = 0; register_index < num_registers; ++register_index) {
    vectorized_banded_edit_distance(register_index, text, read->length, candidates, num_candidates, errs, locations);
    for (int mi = 0; mi < V_CPU; ++mi) {
      if (register_index * V_CPU + mi >= num_candidates) {
        break;
      }
      if ((int) errs[mi] <= error_threshold) {
        uint8_t *pattern = reference.bases + candidates[register_index * V_CPU + mi];
        locations[mi] = (int16_t) generate_alignment(pattern, text, (int) locations[mi], cigar, read->length);
        uint32_t location = candidates[register_index * V_CPU + mi] + (uint32_t) locations[mi] + 1;
        if (last_location == 0xffffffff) {
          last_location = location;
        } else if (location == last_location) {
          continue;
        } else {
          last_location = location;
        }
        generate_MD_tag(pattern, text, (int) locations[mi], cigar, MD);

        int secondary = 0;
        if (num_mappings > 0) {
          secondary = 1;
        }
        for(int reference_index = 0; reference_index < reference.refNum; ++reference_index){
          if (location < reference.lookupTable[reference_index + 1]) {
            append_result_string(reference_index, read, is_reverse_complement, secondary, location, cigar, (int) errs[mi], MD, result_string_length, result_string);
            break;
          }
        }
        *no_left_result = 0;
        ++num_mappings;
      }
    }
  }

  if ((*no_left_result) == 0 && push_output_queue(*result_string)) {
    *no_left_result = 1;
    *result_string_length = 4096;
    *result_string = NULL;
    *result_string = (char*) calloc((*result_string_length), sizeof(char));
  }

  return num_mappings;
}

int banded_edit_distance(const uint8_t* pattern, const uint8_t *text, int *location, const int read_length) {
  int reference_length = read_length + 2 * error_threshold;
  int band_down = 2 * error_threshold;
  int band_length = 2 * error_threshold + 1;

  uint32_t Peq[ALPHABET_SIZE];

  for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
    Peq[ai] = (uint32_t) 0;
  }

  uint32_t tmp = (uint32_t) 1;
  for (int i = 0; i < band_length; i++) {
    Peq[pattern[i]] = Peq[pattern[i]] | tmp;
    tmp = tmp << 1;
  }

  uint32_t Mask = (uint32_t) 1 << (band_length - 1);
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;

  int err = 0;

  uint32_t err_mask = (uint32_t) 1;
  int i_bd = band_down;
  int last_high = band_length - read_length + reference_length - band_down - 1;
  for (int i = 0; i < read_length; i++) {
    X = Peq[text[i]] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X & HP;
    VP = HN | ~(X | HP);
    if (!(D0 & err_mask)) {
      ++err;
      if ((err - last_high) > error_threshold)
        return error_threshold + 1;
    }

    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] >>= 1;
    }

    ++i_bd;
    Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;
  }

  int site = reference_length - last_high - 1;

  *location = -1;
  int error = error_threshold + 1;
  if ((err <= error_threshold) && (err < error)) {
    error = err;
    *location = site;
  }

  for (int i = 0; i < last_high; i++) {
    err = err + ((VP >> i) & (uint32_t) 1);
    err = err - ((VN >> i) & (uint32_t) 1);
    if ((err <= error_threshold) && (err < error)) {
      error = err;
      *location = site + i + 1;
    }
  }

  return error;
}

//#pragma GCC push_options
//#pragma GCC optimize ("unroll-loops")
void vectorized_banded_edit_distance(const uint32_t register_index, const uint8_t *text, const int read_length, const uint32_t *candidates, const uint32_t num_candidates, int16_t *errors, int16_t *locations) {
  int reference_length = read_length + 2 * error_threshold;
  int band_length = 2 * error_threshold + 1;
  int bandLenMask = 0;
  //1 on least significant bits
  for (int i = 0; i < band_length; ++i) {
    bandLenMask <<= 1;
    ++bandLenMask;
  }
  __m128i Peq[ALPHABET_SIZE];
  for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
    Peq[ai] = _mm_set1_epi16(0);
  }
  __m128i patterns[V_CPU];
  //#pragma unroll (V_CPU)
  //for (int pi = 0; pi < V_CPU; ++pi) {
  //	if (register_index * V_CPU + pi >= num_candidates) {
  //			patterns[pi] = _mm_set1_epi8(4);
  //		} else {
  //			patterns[pi] =
  //					_mm_loadu_si128(
  //							(__m128i *) (refs[refIndex].bases
  //									+ candidates[register_index * V_CPU + pi]));
  //		}
  //		__m128i result;
  //		int tmp = 0;
  //#pragma unroll (ALPHABET_SIZE-1)
  //		for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
  //			__m128i letter = _mm_set1_epi8((uint8_t) ai);
  //			result = _mm_cmpeq_epi8(letter, patterns[pi]);
  //			tmp = _mm_movemask_epi8(result);
  //			tmp &= bandLenMask;
  //			Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, pi);
  //		}
  //	}

  if (register_index * V_CPU + 0 >= num_candidates) {
    patterns[0] = _mm_set1_epi8(4);
  } else {    
    patterns[0] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + 0]));
  }
  __m128i result;
  int tmp = 0;
  for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
    __m128i letter = _mm_set1_epi8((uint8_t) ai);
    result = _mm_cmpeq_epi8(letter, patterns[0]);
    tmp = _mm_movemask_epi8(result);
    tmp &= bandLenMask;
    Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 0);
  }

  if (register_index * V_CPU + 1 >= num_candidates) {
    patterns[1] = _mm_set1_epi8(4);
  } else {                                   
    patterns[1] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + 1]));
  }                                                    
  for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
    __m128i letter = _mm_set1_epi8((uint8_t) ai);
    result = _mm_cmpeq_epi8(letter, patterns[1]);
    tmp = _mm_movemask_epi8(result);            
    tmp &= bandLenMask;                                      
    Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 1);
  }

  if (register_index * V_CPU + 2 >= num_candidates) {
    patterns[2] = _mm_set1_epi8(4);
  } else {                                   
    patterns[2] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + 2]));
  }                                                    
  for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
    __m128i letter = _mm_set1_epi8((uint8_t) ai);
    result = _mm_cmpeq_epi8(letter, patterns[2]);
    tmp = _mm_movemask_epi8(result);            
    tmp &= bandLenMask;                                      
    Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 2);
  }

  if (register_index * V_CPU + 3 >= num_candidates) {
    patterns[3] = _mm_set1_epi8(4);
  } else {                                   
    patterns[3] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + 3]));
  }                                                    
  for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
    __m128i letter = _mm_set1_epi8((uint8_t) ai);
    result = _mm_cmpeq_epi8(letter, patterns[3]);
    tmp = _mm_movemask_epi8(result);            
    tmp &= bandLenMask;                                      
    Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 3);
  }

  if (register_index * V_CPU + 4 >= num_candidates) {
    patterns[4] = _mm_set1_epi8(4);
  } else {                                   
    patterns[4] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + 4]));
  }                                                    
  for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
    __m128i letter = _mm_set1_epi8((uint8_t) ai);
    result = _mm_cmpeq_epi8(letter, patterns[4]);
    tmp = _mm_movemask_epi8(result);            
    tmp &= bandLenMask;                                      
    Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 4);
  }

  if (register_index * V_CPU + 5 >= num_candidates) {
    patterns[5] = _mm_set1_epi8(4);
  } else {                                   
    patterns[5] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + 5]));
  }                                                    
  for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
    __m128i letter = _mm_set1_epi8((uint8_t) ai);
    result = _mm_cmpeq_epi8(letter, patterns[5]);
    tmp = _mm_movemask_epi8(result);            
    tmp &= bandLenMask;                                      
    Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 5);
  }

  if (register_index * V_CPU + 6 >= num_candidates) {
    patterns[6] = _mm_set1_epi8(4);
  } else {                                   
    patterns[6] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + 6]));
  }                                                    
  for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
    __m128i letter = _mm_set1_epi8((uint8_t) ai);
    result = _mm_cmpeq_epi8(letter, patterns[6]);
    tmp = _mm_movemask_epi8(result);            
    tmp &= bandLenMask;                                      
    Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 6);
  }


  if (register_index * V_CPU + 7 >= num_candidates) {
    patterns[7] = _mm_set1_epi8(4);
  } else {                                   
    patterns[7] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + 7]));
  }                                                    
  for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
    __m128i letter = _mm_set1_epi8((uint8_t) ai);
    result = _mm_cmpeq_epi8(letter, patterns[7]);
    tmp = _mm_movemask_epi8(result);            
    tmp &= bandLenMask;                                      
    Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 7);
  }

  //a0,a1,a2...a7,a8,a9...->a0,b0,c0...a1,b1,c1...
  permute(patterns);
  uint16_t Mask = (uint16_t) 1 << (band_length - 1);
  __m128i PeqMask = _mm_set1_epi16(Mask);
  __m128i VP = _mm_set1_epi16(0);
  __m128i VN = _mm_set1_epi16(0);
  __m128i X = _mm_set1_epi16(0);
  __m128i D0 = _mm_set1_epi16(0);
  __m128i HN = _mm_set1_epi16(0);
  __m128i HP = _mm_set1_epi16(0);
  __m128i maxMask = _mm_set1_epi16(0xffff);
  __m128i errMask = _mm_set1_epi16(1);
  __m128i errs = _mm_set1_epi16(0);
  int i_bd = 2 * error_threshold;
  int last_high = 2 * error_threshold;
  __m128i threshold = _mm_set1_epi16(error_threshold * 3);
  for (int i = 0; i < read_length; i++) {
    X = _mm_or_si128(Peq[text[i]], VN);
    D0 = _mm_and_si128(X, VP);
    D0 = _mm_add_epi16(D0, VP);
    D0 = _mm_xor_si128(D0, VP);
    D0 = _mm_or_si128(D0, X);
    HN = _mm_and_si128(VP, D0);
    HP = _mm_or_si128(VP, D0);
    HP = _mm_xor_si128(HP, maxMask);
    HP = _mm_or_si128(HP, VN);
    X = _mm_srli_epi16(D0, 1);
    VN = _mm_and_si128(X, HP);
    VP = _mm_or_si128(X, HP);
    VP = _mm_xor_si128(VP, maxMask);
    VP = _mm_or_si128(VP, HN);
    __m128i E = _mm_and_si128(D0, errMask);
    E = _mm_xor_si128(E, errMask);
    errs = _mm_add_epi16(errs, E);
    __m128i earlyEnd = _mm_cmpgt_epi16(errs, threshold);
    int tmp = _mm_movemask_epi8(earlyEnd);
    if (tmp == 0xffff) {
      _mm_store_si128((__m128i *) errors, errs);
      return;
    }
    ++i_bd;
    if (i_bd % 16 == 0) {
      //#pragma unroll (V_CPU)
      for (int pi = 0; pi < V_CPU; ++pi) {
        if (register_index * V_CPU + pi >= num_candidates) {
          patterns[pi] = _mm_set1_epi8(4);
        } else {
          patterns[pi] = _mm_loadu_si128( (__m128i *) (reference.bases + candidates[register_index * V_CPU + pi] + i_bd));
        }
      }
      permute(patterns);
    }
    int pi = (i_bd % 16) / 2;
    if ((i_bd % 16) % 2 != 0) {
      patterns[pi] = _mm_unpackhi_epi64(patterns[pi], patterns[pi]);
    }
    //#pragma unroll (ALPHABET_SIZE-1)
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
      Peq[ai] = _mm_srli_epi16(Peq[ai], 1);
      __m128i letter = _mm_set1_epi16((int16_t) ai);
      __m128i tmpPat = _mm_cvtepu8_epi16(patterns[pi]);
      __m128i mask = _mm_cmpeq_epi16(letter, tmpPat);
      mask = _mm_and_si128(mask, PeqMask);
      Peq[ai] = _mm_or_si128(Peq[ai], mask);
    }
  }
  int site = reference_length - last_high - 1;
  __m128i tmpLocs = _mm_set1_epi16(site);
  __m128i tmpErrs = _mm_set1_epi16(error_threshold + 1);
  tmpErrs = _mm_min_epu16(tmpErrs, errs);
  threshold = _mm_set1_epi16(error_threshold + 1);
  for (int i = 0; i < last_high; i++) {
    __m128i tmpVP = _mm_and_si128(VP, errMask);
    __m128i tmpVN = _mm_and_si128(VN, errMask);
    errs = _mm_add_epi16(errs, tmpVP);
    errs = _mm_sub_epi16(errs, tmpVN);
    __m128i tmpMask1 = _mm_cmplt_epi16(errs, threshold);
    __m128i tmpMask2 = _mm_cmplt_epi16(errs, tmpErrs);
    //__m128i tmpMask3 = _mm_cmpeq_epi16(errs,tmpErrs);
    //tmpMask2 =  _mm_or_si128(tmpMask2, tmpMask3);
    tmpErrs = _mm_min_epu16(tmpErrs, errs);
    tmpMask1 = _mm_and_si128(tmpMask1, tmpMask2);
    int tmp = _mm_movemask_epi8(tmpMask1);
    //#pragma unroll (V_CPU)
    //for (int pi = 0; pi < V_CPU; ++pi) {
    //  if (tmp & (1 << (pi * 2))) {
    //    tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, pi);
    //  }
    //}
    if (tmp & (1 << (0 * 2))) {
      tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, 0);
    }
    if (tmp & (1 << (1 * 2))) {
      tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, 1);
    }
    if (tmp & (1 << (2 * 2))) {
      tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, 2);
    }
    if (tmp & (1 << (3 * 2))) {
      tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, 3);
    }
    if (tmp & (1 << (4 * 2))) {
      tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, 4);
    }
    if (tmp & (1 << (5 * 2))) {
      tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, 5);
    }
    if (tmp & (1 << (6 * 2))) {
      tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, 6);
    }
    if (tmp & (1 << (7 * 2))) {
      tmpLocs = _mm_insert_epi16(tmpLocs, site + i + 1, 7);
    }

    VP = _mm_srli_epi16(VP, 1);
    VN = _mm_srli_epi16(VN, 1);
  }
  _mm_store_si128((__m128i *) errors, tmpErrs);
  _mm_store_si128((__m128i *) locations, tmpLocs);
}
//#pragma GCC pop_options

int generate_alignment(const uint8_t *pattern, const uint8_t *text, int match_site, char *cigar, const int read_length) {
  int reference_length = read_length + 2 * error_threshold;
  int band_down = 2 * error_threshold;
  int band_length = 2 * error_threshold + 1;

  cigar[0] = 0;

  int start_location = match_site - read_length + 1;
  int tmp_err = 0;

  for (int i = 0; i < read_length; i++)
    if (text[i] != pattern[i + start_location])
      tmp_err++;

  if (tmp_err == error_threshold) {
    sprintf(cigar + strlen(cigar), "%d%c", read_length, 'M');
    return start_location;
  }

  uint32_t D0_arry_64[READ_LEN_MAX];
  uint32_t HP_arry_64[READ_LEN_MAX];
  int Route_Size_Whole[READ_LEN_MAX];
  char Route_Char_Whole[READ_LEN_MAX];

  uint32_t Peq[ALPHABET_SIZE];

  for (int ai = 0; ai < ALPHABET_SIZE - 1; ai++) {
    Peq[ai] = (uint32_t) 0;
  }

  uint32_t tmp = (uint32_t) 1;
  for (int i = 0; i < band_length; i++) {
    Peq[pattern[i]] = Peq[pattern[i]] | tmp;
    tmp = tmp << 1;
  }

  uint32_t Mask = (uint32_t) 1 << (band_length - 1);
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t HN = 0;

  int j = 0;
  int i_bd = band_down;
  int last_high = band_length - read_length + reference_length - band_down - 1;

  uint32_t tmp_D0 = 0;
  uint32_t tmp_HP = 0;
  int i = 0;
  for (i = 0; i < read_length; i++) {
    X = Peq[text[i]] | VN;
    tmp_D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & tmp_D0;
    tmp_HP = VN | ~(VP | tmp_D0);
    X = tmp_D0 >> 1;
    VN = X & tmp_HP;
    VP = HN | ~(X | tmp_HP);
    D0_arry_64[i] = tmp_D0;
    HP_arry_64[i] = tmp_HP;

    for (int ai = 0; ai < ALPHABET_SIZE - 1; ai++) {
      Peq[ai] >>= 1;
    }

    ++i_bd;
    if ((i_bd) < reference_length)
      Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;
  }

  int site = reference_length - last_high - 1;
  int search_site = match_site - site;
  int pre_size = 1;
  char pre_char = 'N';
  uint32_t Mask_1 = (uint32_t) 1;
  i = read_length - 1;
  int sum_err = 0;

  j = 1;
  if (((D0_arry_64[i] >> search_site) & Mask_1)
      && (pattern[match_site] == text[i])) {
    i--;
    pre_size = 1;
    pre_char = 'M';
    match_site--;
  } else if (!((D0_arry_64[i] >> search_site) & Mask_1)
      && (pattern[match_site] != text[i])) {
    i--;
    pre_size = 1;
    pre_char = 'S';
    match_site--;
    sum_err++;
  } else if ((HP_arry_64[i] >> search_site) & Mask_1) {
    i--;
    search_site++;
    pre_size = 1;

    pre_char = 'S';
    sum_err++;
    start_location++;
  } else {
    search_site--;
    pre_size = 1;
    pre_char = 'D';

    match_site--;
    sum_err++;
    start_location--;
  }

  while (i >= 0) {
    if (sum_err == error_threshold)
      break;

    if (((D0_arry_64[i] >> search_site) & Mask_1)
        && (pattern[match_site] == text[i])) {
      i--;
      match_site--;

      if (pre_char != 'M') {
        Route_Size_Whole[j] = pre_size;
        Route_Char_Whole[j++] = pre_char;
        pre_size = 1;
        pre_char = 'M';
      } else
        pre_size++;
    } else if (!((D0_arry_64[i] >> search_site) & Mask_1)
        && (pattern[match_site] != text[i])) {
      i--;
      match_site--;
      sum_err++;

      if (pre_char != 'S') {
        Route_Size_Whole[j] = pre_size;
        Route_Char_Whole[j++] = pre_char;
        pre_size = 1;
        pre_char = 'S';
      } else
        pre_size++;
    } else if ((HP_arry_64[i] >> search_site) & Mask_1) {
      i--;
      search_site++;
      sum_err++;
      if (pre_char != 'I') {
        Route_Size_Whole[j] = pre_size;
        Route_Char_Whole[j++] = pre_char;
        pre_size = 1;
        pre_char = 'I';
      } else
        pre_size++;
      start_location++;
    } else {
      search_site--;
      match_site--;
      sum_err++;
      if (pre_char != 'D') {
        Route_Size_Whole[j] = pre_size;
        Route_Char_Whole[j++] = pre_char;
        pre_size = 1;
        pre_char = 'D';
      } else
        pre_size++;
      start_location--;
    }
  }

  Route_Size_Whole[j] = pre_size;
  Route_Char_Whole[j++] = pre_char;
  if (i >= 0) {
    Route_Size_Whole[j] = i + 1;
    Route_Char_Whole[j++] = 'M';
  }

  int size_SM = 0;
  for (j = j - 1; j > 0; j--) {
    if (Route_Char_Whole[j] == 'M' || Route_Char_Whole[j] == 'S') {
      size_SM = 0;
      while (j > 0
          && (Route_Char_Whole[j] == 'M' || Route_Char_Whole[j] == 'S')) {
        size_SM = size_SM + Route_Size_Whole[j];
        j--;
      }
      j++;
      sprintf(cigar + strlen(cigar), "%d%c", size_SM, 'M');
    } else
      sprintf(cigar + strlen(cigar), "%d%c", Route_Size_Whole[j], Route_Char_Whole[j]);
  }

  return start_location;
}

void generate_MD_tag(const uint8_t *pattern, const uint8_t *text, int start_location, char *cigar, char *MD) {
  MD[0] = '\0';
  int cigar_length = strlen(cigar);
  char num_operations_string[READ_LEN_MAX];
  int num_operation_digit_index = 0;
  int num_matches = 0;
  const uint8_t *read = text;
  const uint8_t *reference = pattern + start_location;
  int read_index = 0;
  int reference_index = 0;

  for (int ci = 0; ci < cigar_length; ci++) {
    char c = cigar[ci];
    if (c >= '0' && c <= '9') {
      num_operations_string[num_operation_digit_index++] = c;
    } else {
      num_operations_string[num_operation_digit_index] = '\0';
      int num_operations = atoi(num_operations_string);
      if (c == 'M') {
        for (int opi = 0; opi < num_operations; opi++) {
          if (reference[reference_index] != read[read_index]) {
            //a mismatch
            if (num_matches != 0) {
              sprintf(MD + strlen(MD), "%d", num_matches);
              num_matches = 0;
            }
            switch (reference[reference_index]) {
              case 0:
                sprintf(MD + strlen(MD), "%c", 'A');
                break;
              case 1:
                sprintf(MD + strlen(MD), "%c", 'C');
                break;
              case 2:
                sprintf(MD + strlen(MD), "%c", 'G');
                break;
              case 3:
                sprintf(MD + strlen(MD), "%c", 'T');
                break;
            }
          } else {
            ++num_matches;
          }
          ++reference_index;
          ++read_index;
        }
      } else if (c == 'I') {
        //matchNum += num_operations;
        read_index += num_operations;
      } else if (c == 'D') {
        if (num_matches != 0) {
          sprintf(MD + strlen(MD), "%d", num_matches);
          num_matches = 0;
        }
        sprintf(MD + strlen(MD), "%c", '^');
        for (int opi = 0; opi < num_operations; opi++) {
          switch (reference[reference_index]) {
            case 0:
              sprintf(MD + strlen(MD), "%c", 'A');
              break;
            case 1:
              sprintf(MD + strlen(MD), "%c", 'C');
              break;
            case 2:
              sprintf(MD + strlen(MD), "%c", 'G');
              break;
            case 3:
              sprintf(MD + strlen(MD), "%c", 'T');
              break;
          }
          reference_index++;
        }
      }
      num_operation_digit_index = 0;
    }
  }
  if (num_matches != 0) {
    sprintf(MD + strlen(MD), "%d", num_matches);
    num_matches = 0;
  }
}

void append_result_string(const int reference_index, const Read *read, const int is_reverse_complement, const int secondary, const int location, const char *cigar, const int err, const char *MD, uint32_t *result_string_length, char **result_string) {
  char result[1000];
  result[0] = '\0';
  strcpy(result, read->name);
  int flag = is_reverse_complement * READ_REVERSE_MAPPED + secondary * NOT_PRIMARY_ALIGN;
  sprintf(result + strlen(result), "\t%d\t", flag);
  strcpy(result + strlen(result), reference.names[reference_index]);
  sprintf(result + strlen(result), "\t%d\t%d\t", location - reference.lookupTable[reference_index], 255);
  strcpy(result + strlen(result), cigar);
  strcpy(result + strlen(result), "\t*\t0\t0\t");
  if (secondary == 0) {
    strcpy(result + strlen(result), read->charBases);
  } else {
    strcpy(result + strlen(result), "*");
  }
  sprintf(result + strlen(result), "\t*\tNM:i:%d\tMD:Z:", err);
  strcpy(result + strlen(result), MD);
  strcpy(result + strlen(result), "\n");
  if ((strlen(*result_string) + strlen(result)) >= *result_string_length) {
    (*result_string_length) *= 2;
    char *temp = (char*) realloc(*result_string, sizeof(char) * (*result_string_length));
    assert(temp);
    *result_string = temp;
  }
  strcpy((*result_string) + strlen((*result_string)), result);
}

//uint8_t reverseComplementInt8(uint8_t base) {
//  uint8_t ret;
//  if (base != (uint8_t) 4) {
//    ret = (uint8_t) 3 - base;
//  } else {
//    ret = (uint8_t) 4;
//  }
//  return ret;
//}
