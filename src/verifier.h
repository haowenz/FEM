/*
 * verifier.h
 *
 *  Created on: Mar 25, 2016
 *      Author: howen
 */

#ifndef VERIFIER_H_
#define VERIFIER_H_

#include "utils.h"
#include "readtools.h"
#include "reftools.h"
#include "filter.h"
#include "mapper.h"
#include "outputer.h"
#include <emmintrin.h>
#include <smmintrin.h>

/*we consider alphabet set as {A, T, C, G, N}*/
#define ALPHABET_SIZE 5
/*v_c and v_p are 8 and 16 respectively*/
#define V_CPU 8
#define V_PHI 16

inline void permute(__m128i *patterns);
int CPUVBM(const Read *read, const int rc, const uint32_t *candidates, const uint32_t candiNum,int *noLeftResult, uint32_t *resultStrLen, char** resultStr);
int CPUSBMKernel(const uint8_t* pattern, const uint8_t *text, int *location, const int readLen);
inline void CPUVMBMKernel(const uint32_t regIndex, const uint8_t* text, const int readLen, const uint32_t *candidates, const uint32_t candiNum, int16_t *errors, int16_t *locations);
inline int CPUAlignKernel(const uint8_t* pattern, const uint8_t* text, int match_site, char* cigar, const int readLen); 
inline void CPUNMTag(const uint8_t* pattern, const uint8_t* text, int startLocation, char* cigar, char* NM);
inline void genOneResult(const int refIndex,const Read *read, const int rc, const int secondary, const int location, const char* cigar, const int err,const char* MD, uint32_t *resulStrtLen, char** resultStr);
uint8_t reverseComplementInt8(uint8_t base);

#endif /* VERIFIER_H_ */
