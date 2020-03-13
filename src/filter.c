#include "filter.h"

__UTILS_HASH_VALUE

void generate_optimal_prefix_qgram_for_group_seeding_new(uint32_t readLen,uint32_t gramLen,gCandidate *optCandidates, gCandidate *swapOptCandidates, uint32_t *totalCandidatesNum) {
  uint32_t nGrams = readLen - gramLen + 1;
  //if(nGrams > invertedLists.size()){ cerr << "Assertion Failed: not enough inverted lists" << endl; exit(1); }
  int nSeeds = edit_distance + 2;

  entry matrix[nSeeds][nGrams - nSeeds*gramLen + gramLen];

  uint32_t offset;
  uint32_t min_size, bsize, esize;
  for(int r = 0; r < nSeeds; r++){
    int start = (nSeeds - r - 1)*gramLen;
    int end = nGrams - r*gramLen;
    for(int c = start; c < end; c++){
      int idx = c - start;

      offset = 0;
      min_size = bsize = optCandidates[c].locationsNum;//invertedLists[c].second->size();
      uint64_t min_cost = INT_MAX; 

      for(int i = 0; i < end - start - idx; i++){
        esize = optCandidates[c+i].locationsNum;//invertedLists[c+i].second->size();

        //cost estimation
        uint64_t cost = (i == 0 ? bsize : bsize + esize);
        //if(type == SELECT_INDEL || type == SELECT_HEURISTIC_INDEL){
        if(esize < min_size) min_size = esize;
        uint32_t factor = ((int)gramLen < i ? 2*gramLen : gramLen+i);
        factor *= 10;
        cost += min_size*(readLen + edit_distance)*edit_distance/factor;
        //}
        //else if(type == SELECT_SUBST) cost += 0;

        if(r != 0) cost += matrix[r-1][matrix[r-1][i + idx].ptr].cost;

        if(cost < min_cost){
          min_cost = cost;
          offset = i;
        }
      }

      matrix[r][idx].cost = min_cost;
      matrix[r][idx].offset = offset;
    }

    matrix[r][end-start-1].ptr = end-start-1;
    for(int i = end - start - 2; i >= 0; i--){
      if(matrix[r][i].cost > matrix[r][matrix[r][i+1].ptr].cost)
        matrix[r][i].ptr = matrix[r][i+1].ptr;
      else matrix[r][i].ptr = i;
    }
  }

  // calculate prefix lists and save it into prefixLists
  //prefixLists.clear();

  int s = 0;
  for(int i = 0; i < nSeeds; i++){
    int r = nSeeds - i - 1;

    s = matrix[r][s].ptr;
    offset = matrix[r][s].offset;

    swapOptCandidates[i].site = optCandidates[i*gramLen + s + offset].site;
    swapOptCandidates[i].hashValue = optCandidates[i*gramLen + s + offset].hashValue;
    swapOptCandidates[i].locationsNum = optCandidates[i*gramLen + s + offset].locationsNum;
    *totalCandidatesNum += swapOptCandidates[i].locationsNum;
    // prefixLists.push_back(ListIterWrapper(i*gramLen + s, i*gramLen + s + offset, invertedLists));

    s += offset;
  }

  // sort prefix lists in ascending order of their sizes
  //for (uint32_t i = 0; i < prefixLists.size(); i++) 
  //    for (int j = i + 1; j < nSeeds; j++) 
  //        if (prefixLists[i].size > prefixLists[j].size){
  //            ListIterWrapper tmp = prefixLists[i];
  //            prefixLists[i] = prefixLists[j];
  //            prefixLists[j] = tmp;
  //        }
}



void generate_optimal_prefix_qgram_for_variable_length_seeding_new(uint32_t readLen,uint32_t gramLen,vlCandidate *optCandidates, vlCandidate *swapOptCandidates, uint32_t *totalCandidatesNum){
  uint32_t nGrams = readLen - gramLen + 1;
  //if(nGrams > invertedLists.size()){ cerr << "Assertion Failed: not enough inverted lists" << endl; exit(1); }
  int nSeeds = edit_distance + 2;

  entry matrix[nSeeds][nGrams - nSeeds*gramLen + gramLen];

  uint32_t offset;
  uint32_t min_size, bsize, esize;
  for(int r = 0; r < nSeeds; r++){
    int start = (nSeeds - r - 1)*gramLen;
    int end = nGrams - r*gramLen;
    for(int c = start; c < end; c++){
      int idx = c - start;

      offset = 0;
      min_size = bsize = optCandidates[c].locationsNum;//invertedLists[c].second->size();
      uint64_t min_cost = INT_MAX; 

      for(int i = 0; i < end - start - idx; i++){
        esize = optCandidates[c+i].locationsNum;//invertedLists[c+i].second->size();

        //cost estimation
        uint64_t cost = (i == 0 ? bsize : bsize + esize);
        //if(type == SELECT_INDEL || type == SELECT_HEURISTIC_INDEL){
        if(esize < min_size) min_size = esize;
        uint32_t factor = ((int)gramLen < i ? 2*gramLen : gramLen+i);
        factor *= 10;
        cost += min_size*(readLen + edit_distance)*edit_distance/factor;
        //}
        //else if(type == SELECT_SUBST) cost += 0;

        if(r != 0) cost += matrix[r-1][matrix[r-1][i + idx].ptr].cost;

        if(cost < min_cost){
          min_cost = cost;
          offset = i;
        }
      }

      matrix[r][idx].cost = min_cost;
      matrix[r][idx].offset = offset;
    }

    matrix[r][end-start-1].ptr = end-start-1;
    for(int i = end - start - 2; i >= 0; i--){
      if(matrix[r][i].cost > matrix[r][matrix[r][i+1].ptr].cost)
        matrix[r][i].ptr = matrix[r][i+1].ptr;
      else matrix[r][i].ptr = i;
    }
  }

  // calculate prefix lists and save it into prefixLists
  //prefixLists.clear();

  int s = 0;
  for(int i = 0; i < nSeeds; i++){
    int r = nSeeds - i - 1;

    s = matrix[r][s].ptr;
    offset = matrix[r][s].offset;

    swapOptCandidates[i].startSite = optCandidates[i*gramLen + s + offset].startSite;
    swapOptCandidates[i].endSite = optCandidates[i*gramLen + s + offset].endSite;
    swapOptCandidates[i].hashValue = optCandidates[i*gramLen + s + offset].hashValue;
    swapOptCandidates[i].locationsNum = optCandidates[i*gramLen + s + offset].locationsNum;
    *totalCandidatesNum += swapOptCandidates[i].locationsNum;
    // prefixLists.push_back(ListIterWrapper(i*gramLen + s, i*gramLen + s + offset, invertedLists));

    s += offset;
  }

  // sort prefix lists in ascending order of their sizes
  //for (uint32_t i = 0; i < prefixLists.size(); i++) 
  //    for (int j = i + 1; j < nSeeds; j++) 
  //        if (prefixLists[i].size > prefixLists[j].size){
  //            ListIterWrapper tmp = prefixLists[i];
  //            prefixLists[i] = prefixLists[j];
  //            prefixLists[j] = tmp;
  //        }
}



void generate_optimal_prefix_qgram_for_group_seeding(int kmerSize,int readLen, gCandidate *kmerCandidates, gCandidate *swapKmerCandidates, uint32_t *totalCandiNum) {
  int row = edit_distance + 1 + additional_gram_num + 1;
  int col = readLen - (edit_distance + 1 + additional_gram_num) * kmerSize + 1 + 1;
  int M[row][col];
  int direction[row][col];
  for (int i = 1; i < row; ++i) {
    M[i][0] = hash_table_size;
    direction[i][0] = 3;
  }
  for (int i = 1; i < col; ++i) {
    M[0][i] = 0;
    direction[0][i] = 3;
  }
  for (int ri = 1; ri < row; ++ri) {
    for (int ci = 1; ci < col; ++ci) {
      int positionIndex = ci + (ri - 1) * kmerSize - 1;
      int cmp = M[ri - 1][ci] + kmerCandidates[positionIndex].locationsNum;
      if (cmp < M[ri][ci - 1]) {
        M[ri][ci] = cmp;
        direction[ri][ci] = 2;
      } else {
        M[ri][ci] = M[ri][ci - 1];
        direction[ri][ci] = 1;
      }
    }
  }
  *totalCandiNum = M[row - 1][col - 1];        
  int count = 0;
  int i = row - 1;
  int j = col - 1;
  while (1) {
    if (direction[i][j] == 2) {
      swapKmerCandidates[count].hashValue = kmerCandidates[j + (i - 1) * kmerSize - 1].hashValue;
      swapKmerCandidates[count].site = kmerCandidates[j + (i - 1) * kmerSize - 1].site;
      swapKmerCandidates[count].locationsNum = kmerCandidates[j + (i - 1) * kmerSize - 1].locationsNum;
      ++count;
      --i;
    } else if (direction[i][j] == 1) {
      --j;
    } else {
      break;
    }
  }
}


void generate_optimal_prefix_qgram_for_variable_length_seeding(int kmerSize,int readLen, vlCandidate *kmerCandidates, vlCandidate *swapKmerCandidates, uint32_t *totalCandiNum) {
  int row = edit_distance + 1 + additional_gram_num + 1;
  int col = readLen - (edit_distance + 1 + additional_gram_num) * kmerSize + 1 + 1;
  int M[row][col];
  int direction[row][col];
  for (int i = 1; i < row; ++i) {
    M[i][0] = hash_table_size;
    direction[i][0] = 3;
  }
  for (int i = 1; i < col; ++i) {
    M[0][i] = 0;
    direction[0][i] = 3;
  }
  for (int ri = 1; ri < row; ++ri) {
    for (int ci = 1; ci < col; ++ci) {
      int positionIndex = ci + (ri - 1) * kmerSize - 1;
      int cmp = M[ri - 1][ci] + kmerCandidates[positionIndex].locationsNum;
      if (cmp < M[ri][ci - 1]) {
        M[ri][ci] = cmp;
        direction[ri][ci] = 2;
      } else {
        M[ri][ci] = M[ri][ci - 1];
        direction[ri][ci] = 1;
      }
    }
  }
  *totalCandiNum = M[row - 1][col - 1];        
  int count = 0;
  int i = row - 1;
  int j = col - 1;
  while (1) {
    if (direction[i][j] == 2) {
      swapKmerCandidates[edit_distance + additional_gram_num - count].startSite = kmerCandidates[j + (i - 1) * kmerSize - 1].startSite;
      swapKmerCandidates[edit_distance + additional_gram_num - count].endSite = kmerCandidates[j + (i - 1) * kmerSize - 1].endSite;
      swapKmerCandidates[edit_distance + additional_gram_num - count].locationsNum = kmerCandidates[j + (i - 1) * kmerSize - 1].locationsNum;
      ++count;
      --i;
    } else if (direction[i][j] == 1) {
      --j;
    } else {
      break;
    }
  }
}



uint32_t generate_group_seeding_candidates(const Read *read, const int reverseComplemented, uint32_t **candidatesList,uint32_t **swapCandidatesList, uint32_t *candidatesNumMax, TwoTuple **candis, TwoTuple **tempCandis, uint32_t *tupleNum,uint32_t *candiNumWithoutAddFilter) {
  const uint8_t *bases = read->bases;
  if (reverseComplemented == 1) {
    bases = read->rc_bases;
  }

  //dp for seed selection start
  int scanSize = read->length - window_size + 1;
  uint32_t tempCandidateNums[scanSize];
  int tempHashValues[scanSize];

  int mask = max_hash_value;
  int Nflag = 1;
  int hashVal = -1;
  int NNum = 0;
  for (int ri = 0; ri < scanSize; ++ri) {
    if (Nflag == 1) {
      hashVal = hashValue(bases + ri, window_size);
    } else {
      int nextBase = (int) bases[ri + window_size - 1];
      if (nextBase == 4) {
        ++NNum;
        if (NNum > edit_distance) {
          return 0;
        }
        hashVal = -1;
      }else{
        hashVal = ((hashVal << 2) & mask) + nextBase;
      }
    }
    tempHashValues[ri] = hashVal;
    //deal with N in read
    if (hashVal != -1) {
      Nflag = 0;
      tempCandidateNums[ri] = lookup_table[hashVal+1] - lookup_table[hashVal];
    } else {
      Nflag = 1;
      tempCandidateNums[ri] = 0;
    }
  }

  int coveredNum = window_size/step_size;
  if(window_size%step_size>0){
    coveredNum++;
  }

  int minTotalSubSeedNum = (read->length - window_size + 1 - step_size)/step_size;
  if (coveredNum * (edit_distance + 1 +additional_gram_num) > minTotalSubSeedNum) {
    return 0;
  }

  uint32_t totalCandidatesCount = 0;
  *candiNumWithoutAddFilter=0;
  uint32_t count = 0;
  for(int si=0;si<step_size;++si){
    int totalSubSeedNum = (read->length - window_size + 1 - si)/step_size;
    gCandidate optCandidates[totalSubSeedNum];
    gCandidate swapOptCandidates[edit_distance + 1 + additional_gram_num];
    for(int k=0;k<totalSubSeedNum;++k){
      optCandidates[k].site = k*step_size + si; 
      optCandidates[k].hashValue = tempHashValues[si+k*step_size];
      optCandidates[k].locationsNum = tempCandidateNums[si+k*step_size];
    }
    uint32_t totalPositionNum = 0;
    generate_optimal_prefix_qgram_for_group_seeding(coveredNum, totalSubSeedNum, optCandidates, swapOptCandidates,&totalPositionNum);

    //generate_optimal_prefix_qgram_for_group_seeding_new(totalSubSeedNum,coveredNum,optCandidates,swapOptCandidates,&totalPositionNum);

    totalCandidatesCount += totalPositionNum;

    *candiNumWithoutAddFilter=totalCandidatesCount;
    qsort(swapOptCandidates, edit_distance + 1 + additional_gram_num, sizeof(gCandidate), compare_gCandidate);
    if (totalPositionNum > *tupleNum) {
      *tupleNum = totalPositionNum + (V_CPU_WIDE - (totalPositionNum % V_CPU_WIDE));
      TwoTuple* t_candis = (TwoTuple*) realloc(*candis, sizeof(TwoTuple) * (*tupleNum));
      assert(t_candis);
      *candis = t_candis;
      TwoTuple* t_tempCandis = (TwoTuple*) realloc(*tempCandis, sizeof(TwoTuple) * (*tupleNum));
      assert(t_tempCandis);
      *tempCandis = t_tempCandis;
    }
    uint32_t locationNum = swapOptCandidates[0].locationsNum;
    uint32_t *locations =  occurrence_table + lookup_table[swapOptCandidates[0].hashValue];
    for (uint32_t pi = 0; pi < locationNum; ++pi) {
      TwoTuple tmpPosCount;
      tmpPosCount.a = locations[pi] - swapOptCandidates[0].site;
      tmpPosCount.b = 1;
      (*candis)[pi] = tmpPosCount;
    }
    uint32_t candiSize = locationNum;
    for (int ki = 1; ki < edit_distance + 1 + additional_gram_num; ++ki) {
      locations = occurrence_table + lookup_table[swapOptCandidates[ki].hashValue];
      locationNum = swapOptCandidates[ki].locationsNum;

      uint32_t tempIndex = 0;
      uint32_t ci = 0;
      uint32_t j = 0;
      int additionalPrefix = ki - edit_distance;

      while (ci < candiSize) {
        uint32_t location = locations[j] - swapOptCandidates[ki].site;

        int detel = (*candis)[ci].a - location;
        if (j < locationNum && detel >edit_distance) {
          /*if this is dealing with the last one,
           ** then there is no need to do this.*/
          if (additionalPrefix <= 0) {
            TwoTuple tmpPosCount;
            tmpPosCount.a = location;
            tmpPosCount.b = 1;
            (*tempCandis)[tempIndex++] = tmpPosCount;
          }
          ++j;
        } else if (j == locationNum || (-detel) > edit_distance) {
          if (additionalPrefix <= 0 || (*candis)[ci].b > additional_gram_num) {
            (*tempCandis)[tempIndex++] = (*candis)[ci];
          }
          ++ci;
        } else {
          if (detel == 0) {
            TwoTuple tmpPosCount;
            tmpPosCount.a = location;
            tmpPosCount.b = (*candis)[ci].b + 1;
            (*tempCandis)[tempIndex++] = tmpPosCount;
            ++ci;
            ++j;
          } else if (detel > 0) {
            TwoTuple tmpPosCount;
            tmpPosCount.a = location;
            tmpPosCount.b = 2;
            (*tempCandis)[tempIndex++] = tmpPosCount;
            ++j;
          } else {
            TwoTuple tmpPosCount;
            tmpPosCount.a = (*candis)[ci].a;
            tmpPosCount.b = (*candis)[ci].b + 1;
            (*tempCandis)[tempIndex++] = tmpPosCount;
            ++ci;
          }
        }
      }
      /* if this is dealing with the last one,
       ** then there is no need to do this.*/
      while (additionalPrefix <= 0 && j < locationNum) {
        uint32_t location = locations[j] - swapOptCandidates[ki].site;
        TwoTuple tmpPosCount;
        tmpPosCount.a = location;
        tmpPosCount.b = 1;
        (*tempCandis)[tempIndex++] = tmpPosCount;
        ++j;
      }
      TwoTuple *tmpTuple = *candis;
      *candis = *tempCandis;
      *tempCandis = tmpTuple;
      candiSize = tempIndex;
    }

    /* we have to initialize maxCandiNum with a constant before we use this function.
     ** I think maybe 4096 is OK.*/
    for (uint32_t i = 0; i < candiSize; ++i) {
      if ((*candis)[i].b > additional_gram_num) {
        if ((*candis)[i].a >= (uint32_t)edit_distance) {
          uint32_t location = (*candis)[i].a - edit_distance;
          (*candidatesList)[count++] = location;
          if (count >= *candidatesNumMax) {
            *candidatesNumMax= count * 2;
            uint32_t *t_candidates = (uint32_t*) realloc(*candidatesList, sizeof(uint32_t) * count * 2);
            assert(t_candidates);
            *candidatesList = t_candidates;
          }
        }
      }
    }

  }
  return count;
}

uint32_t generate_variable_length_seeding_candidates(const Read *read, const int reverseComplemented, uint32_t **candidatesList,uint32_t **swapCandidatesList, uint32_t *candidatesNumMax, TwoTuple **candis, TwoTuple **tempCandis, uint32_t *tupleNum,uint32_t *candiNumWithoutAddFilter) {
  const uint8_t *bases = read->bases;
  if (reverseComplemented == 1) {
    bases = read->rc_bases;
  }

  //dp for seed selection start
  int seedLengthMin = window_size + step_size - 1;
  int scanSize = read->length - window_size + 1;

  if (seedLengthMin * (edit_distance + 1 + additional_gram_num) > read->length) {
    return 0;
  }

  uint32_t tempCandidateNums[scanSize];
  int tempHashValues[scanSize];

  int mask = max_hash_value;
  int Nflag = 1;
  int hashVal = -1;
  int NNum = 0;
  for (int ri = 0; ri < scanSize; ++ri) {
    if (Nflag == 1) {
      hashVal = hashValue(bases + ri, window_size);
    } else {
      int nextBase = (int) bases[ri + window_size - 1];
      if (nextBase == 4) {
        ++NNum;
        if (NNum > edit_distance) {
          return 0;
        }
        hashVal = -1;
      }else{
        hashVal = ((hashVal << 2) & mask) + nextBase;
      }
    }
    tempHashValues[ri] = hashVal;
    //deal with N in read
    if (hashVal != -1) {
      Nflag = 0;
      tempCandidateNums[ri] = lookup_table[hashVal+1] - lookup_table[hashVal];
    } else {
      Nflag = 1;
      tempCandidateNums[ri] = 0;
    }
  }
  vlCandidate optCandidates[scanSize];
  vlCandidate swapOptCandidates[edit_distance+additional_gram_num+1];

  uint32_t estimatedLocationNum[scanSize];
  memset(estimatedLocationNum,0,sizeof(int)*scanSize);
  for(int i=0;i<read->length-window_size+1-step_size+1;++i){
    estimatedLocationNum[i]=tempCandidateNums[i];
    for(int j=1;j<step_size;++j){
      estimatedLocationNum[i]+=tempCandidateNums[i+j];
    }
    optCandidates[i].startSite = i;
    optCandidates[i].endSite = i + seedLengthMin - 1;
    optCandidates[i].locationsNum = estimatedLocationNum[i];
  }

  uint32_t totalPositionNum = 0;

  //generate_optimal_prefix_qgram_for_variable_length_seeding_new(read->length,seedLengthMin,optCandidates, swapOptCandidates, &totalPositionNum);

  generate_optimal_prefix_qgram_for_variable_length_seeding(seedLengthMin, read->length, optCandidates, swapOptCandidates,&totalPositionNum);
  if(totalPositionNum==0){
    return 0;
  }


  swapOptCandidates[0].startSite = 0;
  swapOptCandidates[edit_distance + additional_gram_num].endSite = read->length - 1;
  for (int ki = 1; ki < edit_distance + 1 + additional_gram_num; ++ki){
    if(swapOptCandidates[ki].locationsNum>=swapOptCandidates[ki-1].locationsNum){
      swapOptCandidates[ki].startSite = swapOptCandidates[ki - 1].endSite + 1;
    }else{
      swapOptCandidates[ki - 1].endSite = swapOptCandidates[ki].startSite - 1;
    }
  }


  uint32_t totalCandidatesCount = 0;
  *candiNumWithoutAddFilter=0;
  int totalKmerNum=(edit_distance+1+additional_gram_num)*(step_size);
  SubSeed minSubSeeds[totalKmerNum];
  int seedNum = edit_distance+1+additional_gram_num;
  for(int ei=0; ei<edit_distance+1+additional_gram_num; ++ei){
    for(int si=0;si<step_size;++si){
      minSubSeeds[ei*step_size+si].site=swapOptCandidates[ei].startSite+si;
      minSubSeeds[ei*step_size+si].locationNum=tempCandidateNums[minSubSeeds[ei*step_size+si].site];
    }
    int seedLength = swapOptCandidates[ei].endSite - swapOptCandidates[ei].startSite + 1;
    int endIndex = step_size-1;
    for(int sj=step_size;sj<seedLength-window_size+1;++sj){
      int index = sj%step_size;
      if(tempCandidateNums[swapOptCandidates[ei].startSite + sj]<(uint32_t)minSubSeeds[ei*step_size+index].locationNum){
        minSubSeeds[ei*step_size+index].site=swapOptCandidates[ei].startSite + sj;
        minSubSeeds[ei*step_size+index].locationNum=tempCandidateNums[minSubSeeds[ei*step_size+index].site];
        endIndex=sj;
      }
    }

    if(step_size>1){
      if(endIndex>step_size -1 && ei<edit_distance+additional_gram_num){
        swapOptCandidates[ei+1].startSite = swapOptCandidates[ei].startSite + endIndex+ window_size;
      }
    }
    for(int sk=0;sk<step_size;++sk){
      totalCandidatesCount+=minSubSeeds[ei*step_size+sk].locationNum;
    }
  }

  *candiNumWithoutAddFilter=totalCandidatesCount;
  if (totalCandidatesCount > *tupleNum) {
    *tupleNum = totalCandidatesCount + (V_CPU_WIDE - (totalCandidatesCount % V_CPU_WIDE));
    TwoTuple* t_candis = (TwoTuple*) realloc(*candis, sizeof(TwoTuple) * (*tupleNum));
    assert(t_candis);
    *candis = t_candis;
    TwoTuple* t_tempCandis = (TwoTuple*) realloc(*tempCandis, sizeof(TwoTuple) * (*tupleNum));
    assert(t_tempCandis);
    *tempCandis = t_tempCandis;
  }

  uint32_t candiSize = 0;
  qsort(minSubSeeds, (edit_distance+additional_gram_num+1)*step_size, sizeof(SubSeed), compare_sub_seed);
  uint32_t locationNum = minSubSeeds[0].locationNum;
  uint32_t *locations=occurrence_table + lookup_table[tempHashValues[minSubSeeds[0].site]];
  for(uint32_t i=0;i<locationNum;++i){
    uint32_t location = locations[i] - minSubSeeds[0].site;
    TwoTuple tmpPosCount;
    tmpPosCount.a = location;
    tmpPosCount.b = 1;
    (*candis)[candiSize++] = tmpPosCount;
  }

  for(int ei = 0; ei<seedNum; ++ei){
    for(int si=0;si<step_size;++si){
      if(ei==0&&si==0){
        continue;
      }

      locationNum = minSubSeeds[ei*step_size+si].locationNum;
      locations =  occurrence_table + lookup_table[tempHashValues[minSubSeeds[ei*step_size+si].site]];

      uint32_t tempIndex = 0;
      uint32_t ci = 0;
      uint32_t j = 0;
      int additionalPrefix =  ei * step_size+si- ((edit_distance+1+additional_gram_num) * step_size );

      while (ci < candiSize) {
        uint32_t location = locations[j] - minSubSeeds[ei*step_size+si].site;
        //__builtin_prefetch((const void*) (locations + ((ci + 2))), 0, 0);
        int detel = (*candis)[ci].a - location;
        if (j < locationNum && detel >edit_distance) {
          /*if this is dealing with the last one,
           ** then there is no need to do this.*/
          if (additionalPrefix <= 0) {
            TwoTuple tmpPosCount;
            tmpPosCount.a = location;
            tmpPosCount.b = 1;
            (*tempCandis)[tempIndex++] = tmpPosCount;
          }
          ++j;
        } else if (j == locationNum || (-detel) > edit_distance) {
          if (additionalPrefix <= 0 || (*candis)[ci].b > additional_gram_num) {
            (*tempCandis)[tempIndex++] = (*candis)[ci];
          }
          ++ci;
        } else {
          if (detel == 0) {
            TwoTuple tmpPosCount;
            tmpPosCount.a = location;
            tmpPosCount.b = (*candis)[ci].b + 1;
            (*tempCandis)[tempIndex++] = tmpPosCount;
            ++ci;
            ++j;
          } else if (detel > 0) {
            TwoTuple tmpPosCount;
            tmpPosCount.a = location;
            tmpPosCount.b = 2;
            (*tempCandis)[tempIndex++] = tmpPosCount;
            ++j;
          } else {
            TwoTuple tmpPosCount;
            tmpPosCount.a = (*candis)[ci].a;
            tmpPosCount.b = (*candis)[ci].b + 1;
            (*tempCandis)[tempIndex++] = tmpPosCount;
            ++ci;
          }
        }
      }
      /*if this is dealing with the last one,
       ** then there is no need to do this.*/
      while (additionalPrefix <= 0 && j < locationNum) {
        uint32_t location = locations[j] - minSubSeeds[ei*step_size+si].site;
        TwoTuple tmpPosCount;
        tmpPosCount.a = location;
        tmpPosCount.b = 1;
        (*tempCandis)[tempIndex++] = tmpPosCount;
        ++j;
      }
      TwoTuple *tmpTuple = *candis;
      *candis = *tempCandis;
      *tempCandis = tmpTuple;
      candiSize = tempIndex;
    }
  }

  /*we have to initialize maxCandiNum with a constant before we use this function.
   ** I think maybe 4096 is OK.*/
  uint32_t count = 0;
  for (uint32_t i = 0; i < candiSize; ++i) {
    if ((*candis)[i].b > additional_gram_num) {
      if ((*candis)[i].a >= (uint32_t)edit_distance&& (*candis)[i].a<reference.lookupTable[reference.refNum] ) {
        uint32_t location = (*candis)[i].a - edit_distance;
        (*candidatesList)[count++] = location;
        if (count >= *candidatesNumMax) {
          *candidatesNumMax= count * 2;
          uint32_t *t_candidates = (uint32_t*) realloc(*candidatesList, sizeof(uint32_t) * count * 2);
          assert(t_candidates);
          *candidatesList = t_candidates;
        }
      }
    }
  }
  return count;
}
