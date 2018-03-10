/*
 * filter.h
 *
 *  Created on: Mar 24, 2016
 *      Author: howen
 */

#ifndef FILTER_H_
#define FILTER_H_

#include "indextools.h"
#include "utils.h"
#include "mapper.h"
#include "readtools.h"

uint32_t groupSeedingCandidateGenerator(const Read *read, const int rc, uint32_t **candidates,uint32_t **swapCandidates,uint32_t *maxCandiNum,TwoTuple **candis, TwoTuple **tempCandis, uint32_t *tupleNum,uint32_t *candiNumWithoutAddFilter);

uint32_t variableLengthSeedingCandidateGenerator(const Read *read, const int rc, uint32_t **candidates,uint32_t **swapCandidates,uint32_t *maxCandiNum,TwoTuple **candis, TwoTuple **tempCandis, uint32_t *tupleNum,uint32_t *candiNumWithoutAddFilter);
#endif /* FILTER_H_ */
