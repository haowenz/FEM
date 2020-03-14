#ifndef FILTER_H_
#define FILTER_H_

#include "indextools.h"
#include "utils.h"
#include "mapper.h"
#include "readtools.h"

uint32_t generate_group_seeding_candidates(const Read *read, const int is_reverse_complement, uint32_t **candidates, uint32_t **swap_candidates, uint32_t *max_num_candidates, TwoTuple **candis, TwoTuple **temp_candis, uint32_t *num_tuples, uint32_t *num_candidates_without_additonal_qgram_filter);

uint32_t generate_variable_length_seeding_candidates(const Read *read, const int is_reverse_complement, uint32_t **candidates, uint32_t **swap_candidates, uint32_t *max_num_candidates, TwoTuple **candis, TwoTuple **temp_candis, uint32_t *num_tuples, uint32_t *num_candidates_without_additonal_qgram_filter);
#endif /* FILTER_H_ */
