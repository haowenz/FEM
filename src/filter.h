#ifndef FILTER_H_
#define FILTER_H_

#include "indextools.h"
#include "utils.h"
#include "mapper.h"
#include "readtools.h"

uint32_t generate_group_seeding_candidates(const Read *read, const int rc, uint32_t **candidates,uint32_t **swapCandidates,uint32_t *maxCandiNum,TwoTuple **candis, TwoTuple **tempCandis, uint32_t *tupleNum,uint32_t *candiNumWithoutAddFilter);

uint32_t generate_variable_length_seeding_candidates(const Read *read, const int rc, uint32_t **candidates,uint32_t **swapCandidates,uint32_t *maxCandiNum,TwoTuple **candis, TwoTuple **tempCandis, uint32_t *tupleNum,uint32_t *candiNumWithoutAddFilter);
#endif /* FILTER_H_ */
