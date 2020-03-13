#ifndef UTILS_H_
#define UTILS_H_

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

#define V_CPU_WIDE 4

#define __UTILS_CHAR_UINT8 	\
  static inline uint8_t charToUint8(const char c) {\
    switch (c) {\
      case 'A':\
               return 0;\
      case 'C':\
               return 1;\
      case 'G':\
               return 2;\
      case 'T':\
               return 3;\
      default: {\
                 return 4;\
               }\
    }\
  }

#define __UTILS_HASH_VALUE 	 \
  static inline int hashValue(const uint8_t *seq,int windowSize) {\
    int i = 0;\
    int val = 0;\
    while (i < windowSize) {\
      if (seq[i] == (uint8_t) 4) {\
        return -1;\
      }\
      val = (val << 2) | (int) (seq[i]);\
      i++;\
    }\
    return val;\
  }

typedef struct {
  int hashValue;
  uint32_t location;
} HashLocationPair;

typedef struct {
  uint64_t cost;
  int offset; 
  int ptr; 
} entry;

typedef struct {
  uint32_t a;
  int b;
} TwoTuple;

typedef struct {
  int site;
  int locationNum;
} SubSeed;

typedef struct {
  int hashValue;
  uint32_t site;
  uint32_t locationsNum;
  uint32_t *locations;
} gCandidate;

typedef struct {
  int hashValue;
  uint32_t startSite;
  uint32_t endSite;
  uint32_t locationsNum;
  uint32_t *locations;
} vlCandidate;

double get_real_time();
double get_cpu_time();
int compare_sub_seed(const void *ta, const void *tb); 
int compare_gCandidate(const void *ta, const void *tb);
int compare_hash_location_pair(const void *a, const void *b); 

#endif /* UTILS_H_ */
