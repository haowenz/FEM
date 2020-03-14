#ifndef VERIFIER_H_
#define VERIFIER_H_

#include <emmintrin.h>
#include <smmintrin.h>

#include "utils.h"
#include "readtools.h"
#include "reftools.h"
#include "filter.h"
#include "mapper.h"
#include "outputer.h"

/*we consider alphabet set as {A, T, C, G, N}*/
#define ALPHABET_SIZE 5
/*v_c and v_p are 8 and 16 respectively*/
#define V_CPU 8
#define V_PHI 16

void permute(__m128i *patterns);
int verify_candidates(const Read *read, const int is_reverse_complement, const uint32_t *candidates, const uint32_t num_candidates, int *no_left_result, uint32_t *result_string_length, char** result_string);
int banded_edit_distance(const uint8_t *pattern, const uint8_t *text, int *location, const int read_length);
void vectorized_banded_edit_distance(const uint32_t register_index, const uint8_t *text, const int read_length, const uint32_t *candidates, const uint32_t num_candidates, int16_t *errors, int16_t *locations);
int generate_alignment(const uint8_t *pattern, const uint8_t *text, int match_site, char *cigar, const int read_length); 
void generate_MD_tag(const uint8_t *pattern, const uint8_t *text, int start_location, char *cigar, char *MD);
void append_result_string(const int reference_index, const Read *read, const int is_reverse_complement, const int secondary, const int location, const char *cigar, const int err, const char *MD, uint32_t *result_string_length, char **result_string);
//uint8_t reverseComplementInt8(uint8_t base);

#endif /* VERIFIER_H_ */
