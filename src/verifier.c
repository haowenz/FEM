/*
 * verifier.c
 *
 *  Created on: Mar 25, 2016
 *      Author: howen
 */

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

int CPUVBM(const Read *read, const int rc, const uint32_t *candidates, const uint32_t candiNum, int *noLeftResult, uint32_t *resultLen, char** resultStr) {
	const uint8_t *text = read->bases;
	if (rc == 1) {
		text = read->rc_bases;
	}
	int mappingNum = 0;

	uint32_t regNum = candiNum / V_CPU;
	uint32_t remains = candiNum % V_CPU;
	if (remains != 0) {
		++regNum;
	}

	char cigar[READ_LEN_MAX * 2];
	char MD[READ_LEN_MAX * 2];
	int16_t locations[V_CPU]; 
	int16_t errs[V_CPU];
	uint32_t lastLocation = 0xffffffff;
	for (uint32_t regIndex = 0; regIndex < regNum; ++regIndex) {
		CPUVMBMKernel(regIndex, text, read->length, candidates, candiNum, errs, locations);
		for (int mi = 0; mi < V_CPU; ++mi) {
			if (regIndex * V_CPU + mi >= candiNum) {
				break;
			}
			if ((int) errs[mi] <= edit_distance) {
				uint8_t *pattern = reference.bases + candidates[regIndex * V_CPU + mi];
				locations[mi] = (int16_t) CPUAlignKernel(pattern, text, (int) locations[mi], cigar, read->length);
				uint32_t location = candidates[regIndex * V_CPU + mi] + (uint32_t) locations[mi] + 1;
				if (lastLocation == 0xffffffff) {
					lastLocation = location;
				} else if (location == lastLocation) {
					continue;
				} else {
					lastLocation = location;
				}
				CPUNMTag(pattern, text, (int) locations[mi], cigar, MD);

				int secondary = 0;
				if (mappingNum > 0) {
					secondary = 1;
				}
                for(int refIndex=0;refIndex<reference.refNum;++refIndex){
                    if(location<reference.lookupTable[refIndex+1]){
                        genOneResult(refIndex,read, rc, secondary, location, cigar, (int) errs[mi], MD, resultLen, resultStr);
                        break;
                    }
                }
                *noLeftResult = 0;
                ++mappingNum;
			}
		}
	}

	if ((*noLeftResult) == 0 && pushOutputQueue(*resultStr)) {
		*noLeftResult = 1;
		*resultLen = 4096;
		*resultStr = NULL;
		*resultStr = (char*) calloc((*resultLen), sizeof(char));
	}

	return mappingNum;
}

int CPUSBMKernel(const uint8_t* pattern, const uint8_t *text, int *location,
		const int readLen) {
	int refLen = readLen + 2 * edit_distance;
	int band_down = 2 * edit_distance;
	int band_length = 2 * edit_distance + 1;

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
	int last_high = band_length - readLen + refLen - band_down - 1;
	for (int i = 0; i < readLen; i++) {

		X = Peq[text[i]] | VN;
		D0 = ((VP + (X & VP)) ^ VP) | X;
		HN = VP & D0;
		HP = VN | ~(VP | D0);

		X = D0 >> 1;
		VN = X & HP;
		VP = HN | ~(X | HP);
		if (!(D0 & err_mask)) {
			++err;
			if ((err - last_high) > edit_distance)
				return edit_distance + 1;
		}

		for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
			Peq[ai] >>= 1;
		}

		++i_bd;
		Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;
	}

	int site = refLen - last_high - 1;

	*location = -1;
	int error = edit_distance + 1;
	if ((err <= edit_distance) && (err < error)) {
		error = err;
		*location = site;
	}

	for (int i = 0; i < last_high; i++) {
		err = err + ((VP >> i) & (uint32_t) 1);
		err = err - ((VN >> i) & (uint32_t) 1);

		if ((err <= edit_distance) && (err < error)) {
			error = err;
			*location = site + i + 1;
		}
	}

	return error;
}

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
inline void CPUVMBMKernel(const uint32_t regIndex, const uint8_t* text, const int readLen, const uint32_t *candidates, const uint32_t candiNum,int16_t *errors, int16_t *locations) {
	int refLen = readLen + 2 * edit_distance;
	int band_length = 2 * edit_distance + 1;
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
	//	if (regIndex * V_CPU + pi >= candiNum) {
//			patterns[pi] = _mm_set1_epi8(4);
//		} else {
//			patterns[pi] =
//					_mm_loadu_si128(
//							(__m128i *) (refs[refIndex].bases
//									+ candidates[regIndex * V_CPU + pi]));
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

    if (regIndex * V_CPU + 0 >= candiNum) {
        patterns[0] = _mm_set1_epi8(4);
    } else {    
        patterns[0] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + 0]));
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

    if (regIndex * V_CPU + 1 >= candiNum) {
        patterns[1] = _mm_set1_epi8(4);
    } else {                                   
        patterns[1] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + 1]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[1]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 1);
    }


    if (regIndex * V_CPU + 2 >= candiNum) {
        patterns[2] = _mm_set1_epi8(4);
    } else {                                   
        patterns[2] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + 2]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[2]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 2);
    }


    if (regIndex * V_CPU + 3 >= candiNum) {
        patterns[3] = _mm_set1_epi8(4);
    } else {                                   
        patterns[3] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + 3]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[3]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 3);
    }


    if (regIndex * V_CPU + 4 >= candiNum) {
        patterns[4] = _mm_set1_epi8(4);
    } else {                                   
        patterns[4] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + 4]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[4]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 4);
    }


    if (regIndex * V_CPU + 5 >= candiNum) {
        patterns[5] = _mm_set1_epi8(4);
    } else {                                   
        patterns[5] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + 5]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[5]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 5);
    }


    if (regIndex * V_CPU + 6 >= candiNum) {
        patterns[6] = _mm_set1_epi8(4);
    } else {                                   
        patterns[6] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + 6]));
    }                                                    
    for (int ai = 0; ai < ALPHABET_SIZE - 1; ++ai) {
        __m128i letter = _mm_set1_epi8((uint8_t) ai);
        result = _mm_cmpeq_epi8(letter, patterns[6]);
        tmp = _mm_movemask_epi8(result);            
        tmp &= bandLenMask;                                      
        Peq[ai] = _mm_insert_epi16(Peq[ai], tmp, 6);
    }


    if (regIndex * V_CPU + 7 >= candiNum) {
        patterns[7] = _mm_set1_epi8(4);
    } else {                                   
        patterns[7] =  _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + 7]));
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
	int i_bd = 2 * edit_distance;
	int last_high = 2 * edit_distance;
	__m128i threshold = _mm_set1_epi16(edit_distance * 3);
	for (int i = 0; i < readLen; i++) {
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
				if (regIndex * V_CPU + pi >= candiNum) {
					patterns[pi] = _mm_set1_epi8(4);
				} else {
					patterns[pi] = _mm_loadu_si128( (__m128i *) (reference.bases + candidates[regIndex * V_CPU + pi] + i_bd));
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
	int site = refLen - last_high - 1;
	__m128i tmpLocs = _mm_set1_epi16(site);
	__m128i tmpErrs = _mm_set1_epi16(edit_distance + 1);
	tmpErrs = _mm_min_epu16(tmpErrs, errs);
	threshold = _mm_set1_epi16(edit_distance + 1);
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
#pragma GCC pop_options

inline int CPUAlignKernel(const uint8_t* pattern, const uint8_t* text, int match_site, char* cigar, const int readLen) {
	int refLen = readLen + 2 * edit_distance;
	int band_down = 2 * edit_distance;
	int band_length = 2 * edit_distance + 1;

	cigar[0] = 0;

	int start_location = match_site - readLen + 1;
	int tmp_err = 0;

	for (int i = 0; i < readLen; i++)
		if (text[i] != pattern[i + start_location])
			tmp_err++;

	if (tmp_err == edit_distance) {
		sprintf(cigar + strlen(cigar), "%d%c", readLen, 'M');
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
	int last_high = band_length - readLen + refLen - band_down - 1;

	uint32_t tmp_D0 = 0;
	uint32_t tmp_HP = 0;
	int i = 0;
	for (i = 0; i < readLen; i++) {
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
		if ((i_bd) < refLen)
			Peq[pattern[i_bd]] = Peq[pattern[i_bd]] | Mask;
	}

	int site = refLen - last_high - 1;
	int search_site = match_site - site;
	int pre_size = 1;
	char pre_char = 'N';
	uint32_t Mask_1 = (uint32_t) 1;
	i = readLen - 1;
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
		if (sum_err == edit_distance)
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
			sprintf(cigar + strlen(cigar), "%d%c", Route_Size_Whole[j],
					Route_Char_Whole[j]);
	}

	return start_location;

}

inline void CPUNMTag(const uint8_t* pattern, const uint8_t* text, int startLocation, char* cigar, char* NM) {
	NM[0] = '\0';
	int cigarLen = strlen(cigar);
	char num[READ_LEN_MAX];
	int numIndex = 0;
	int matchNum = 0;
	const uint8_t *read = text;
	const uint8_t *reference = pattern + startLocation;
	int readIndex = 0;
	int refIndex = 0;

	for (int ci = 0; ci < cigarLen; ci++) {
		char c = cigar[ci];
		if (c >= '0' && c <= '9') {
			num[numIndex++] = c;
		} else {
			num[numIndex] = '\0';
			int opNum = atoi(num);
			if (c == 'M') {
				for (int opi = 0; opi < opNum; opi++) {
					if (reference[refIndex] != read[readIndex]) {
						//a mismatch
						if (matchNum != 0) {
							sprintf(NM + strlen(NM), "%d", matchNum);
							matchNum = 0;
						}
						switch (reference[refIndex]) {
						case 0:
							sprintf(NM + strlen(NM), "%c", 'A');
							break;
						case 1:
							sprintf(NM + strlen(NM), "%c", 'C');
							break;
						case 2:
							sprintf(NM + strlen(NM), "%c", 'G');
							break;
						case 3:
							sprintf(NM + strlen(NM), "%c", 'T');
							break;
						}
					} else {
						++matchNum;
					}
					++refIndex;
					++readIndex;
				}
			} else if (c == 'I') {
				//matchNum += opNum;
				readIndex += opNum;
			} else if (c == 'D') {
				if (matchNum != 0) {
					sprintf(NM + strlen(NM), "%d", matchNum);
					matchNum = 0;
				}
				sprintf(NM + strlen(NM), "%c", '^');
				for (int opi = 0; opi < opNum; opi++) {
					switch (reference[refIndex]) {
					case 0:
						sprintf(NM + strlen(NM), "%c", 'A');
						break;
					case 1:
						sprintf(NM + strlen(NM), "%c", 'C');
						break;
					case 2:
						sprintf(NM + strlen(NM), "%c", 'G');
						break;
					case 3:
						sprintf(NM + strlen(NM), "%c", 'T');
						break;
					}
					refIndex++;
				}
			}
			numIndex = 0;
		}
	}
	if (matchNum != 0) {
		sprintf(NM + strlen(NM), "%d", matchNum);
		matchNum = 0;
	}
}

inline void genOneResult(const int refIndex,const Read *read, const int rc, const int secondary,
		const int location, const char* cigar, const int err, const char* MD,
		uint32_t *resultStrLen, char** resultStr) {
	char result[1000];
	result[0] = '\0';
	strcpy(result, read->name);
	int flag = rc * READ_REVERSE_MAPPED + secondary * NOT_PRIMARY_ALIGN;
	sprintf(result + strlen(result), "\t%d\t", flag);
	strcpy(result + strlen(result), reference.names[refIndex]);
	sprintf(result + strlen(result), "\t%d\t%d\t", location-reference.lookupTable[refIndex], 255);
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
	if ((strlen(*resultStr) + strlen(result)) >= *resultStrLen) {
		(*resultStrLen) *= 2;
		char *temp = (char*) realloc(*resultStr, sizeof(char) * (*resultStrLen));
		assert(temp);
		*resultStr = temp;
	}
	strcpy((*resultStr) + strlen((*resultStr)), result);
}

uint8_t reverseComplementInt8(uint8_t base) {
	uint8_t ret;
	if (base != (uint8_t) 4) {
		ret = (uint8_t) 3 - base;
	} else {
		ret = (uint8_t) 4;
	}
	return ret;
}
