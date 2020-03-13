#include "mapper.h"

int window_size;
int step_size;
int edit_distance;
int additional_gram_num;
int cpu_thread_num;
int finished_thread_num;
int cpu_thread_num_per_group;

uint32_t candidate_num_without_add_filter[256];
uint32_t mapped_read_num[256];
uint32_t candidate_num[256];
uint32_t mapping_num[256];
uint32_t read_count[256];

uint32_t (*novelCandidateGenerator)(const Read *read, const int rc, uint32_t **candidates,uint32_t **swapCandidates,uint32_t *maxCandiNum,TwoTuple **candis, TwoTuple **tempCandis, uint32_t *tupleNum,uint32_t *candiNumWithoutAddFilter);


void initialize_mapper() {
  initialize_ref(&reference);
  initialize_read_file();
  initialize_read_queue();
  initialize_loading_index();
  initialize_output_file();
  initialize_output_queue();
  output_header_from_index_file();
}

void finalize_mapper() {
  destroy_ref(&reference);
  destroy_read_queue();
  destroy_output_queue();
  finalize_read_file();
  finalize_loading_index();
  finalize_output_file();
}

int CPU_map(int threadID) {
  uint32_t mappingNum = 0;
  uint32_t localCandiNum = 0;
  uint32_t totalCandidateNum = 0;
  uint32_t candiNumWithoutAddFilter = 0;
  uint32_t totalCandidateNumWithoutAddFilter = 0;
  int noLeftResult = 0;
  uint32_t readCount = 0;
  uint32_t mappedReadCount = 0;
  uint32_t maxCandiNum = 4096;
  uint32_t tupleNum = 4096;
  uint32_t resultStrLen = 4096;
  uint32_t *candidates = (uint32_t*) calloc(maxCandiNum, sizeof(uint32_t));
  uint32_t *swapCandidates = (uint32_t*) calloc(maxCandiNum, sizeof(uint32_t));
  TwoTuple *candis = (TwoTuple*) calloc(tupleNum, sizeof(TwoTuple));
  TwoTuple *tempCandis = (TwoTuple*) calloc(tupleNum, sizeof(TwoTuple));
  char *resultStr = (char*) calloc(resultStrLen, sizeof(char));
  Read read;

  while (pop_read_queue(&read) != -1) {
    readCount++;
    int tmpMappedReadNum = 0;
    localCandiNum = novelCandidateGenerator(&read, 0, &candidates, &swapCandidates, &maxCandiNum, &candis, &tempCandis, &tupleNum, &candiNumWithoutAddFilter);
    totalCandidateNumWithoutAddFilter += candiNumWithoutAddFilter;

    if (localCandiNum > 0) {
      totalCandidateNum += localCandiNum;
      tmpMappedReadNum += verify_candidates(&read, 0, candidates, localCandiNum, &noLeftResult, &resultStrLen, &resultStr);
    }

    localCandiNum = novelCandidateGenerator(&read, 1, &candidates,&swapCandidates, &maxCandiNum,&candis,&tempCandis,&tupleNum ,&candiNumWithoutAddFilter);
    totalCandidateNumWithoutAddFilter += candiNumWithoutAddFilter;

    if (localCandiNum > 0) {
      totalCandidateNum += localCandiNum;
      tmpMappedReadNum += verify_candidates(&read, 1, candidates, localCandiNum, &noLeftResult, &resultStrLen, &resultStr);
    }

    if(tmpMappedReadNum>0) {
      mappingNum += tmpMappedReadNum; 
      mappedReadCount++;
    }
  }

  //fprintf(stderr, "Thread %d is pushing local results to output queue...\n", threadID);

  while (noLeftResult != 1) {
    if (push_output_queue(resultStr)) {
      noLeftResult = 1;
      //free(resultStr);
      resultStr = NULL;
    }
  }

  //fprintf(stderr, "Thread %d pushed local results to output queue successfully.\n", threadID);

  free(candidates);
  free(swapCandidates);
  free(candis);
  free(tempCandis);
  candidate_num_without_add_filter[threadID] = totalCandidateNumWithoutAddFilter;
  candidate_num[threadID] = totalCandidateNum;
  mapping_num[threadID] = mappingNum;
  read_count[threadID] = readCount;
  mapped_read_num[threadID] = mappedReadCount;
  fprintf(stderr, "Thread %d completed.\n", threadID);

  return mappingNum;
}


void *start_CPU_thread(void *arg) {
  int threadID = *((int*) arg);
  CPU_map(threadID);
  return NULL;
}

void *start_single_read_queue_thread(void *arg) {
  single_read_queue_thread();
  return NULL;
}

void *start_output_queue_thread(void *arg) {
  output_queue_thread();
  return NULL;
}
