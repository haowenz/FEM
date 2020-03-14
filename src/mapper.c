#include "mapper.h"

int kmer_size;
int step_size;
int error_threshold;
int num_additional_qgrams;
int num_threads;
int num_finished_threads;
int num_threads_per_group;

uint32_t num_candidates_without_filter[256];
uint32_t num_mapped_reads[256];
uint32_t num_candidates[256];
uint32_t num_mappings[256];
uint32_t num_reads[256];

uint32_t (*generate_candidates)(const Read *read, const int is_reverse_complement, uint32_t **candidates, uint32_t **swap_candidates, uint32_t *max_num_candidates, TwoTuple **candis, TwoTuple **temp_candis, uint32_t *num_tuples, uint32_t *num_candidates_without_additonal_qgram_filter);

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

int CPU_map(int thread_id) {
  uint32_t mappingNum = 0;
  uint32_t tmp_num_candidates = 0;
  uint32_t total_num_candidates = 0;
  uint32_t tmp_num_candidates_without_additional_qgram_filter = 0;
  uint32_t total_num_candidates_without_additional_qgram_filter = 0;
  int no_left_result = 0;
  uint32_t tmp_num_reads = 0;
  uint32_t tmp_num_mapped_reads = 0;
  uint32_t max_num_candidates = 4096;
  uint32_t num_tuples = 4096;
  uint32_t result_string_length = 4096;
  uint32_t *candidates = (uint32_t*) calloc(max_num_candidates, sizeof(uint32_t));
  uint32_t *swap_candidates = (uint32_t*) calloc(max_num_candidates, sizeof(uint32_t));
  TwoTuple *candis = (TwoTuple*) calloc(num_tuples, sizeof(TwoTuple));
  TwoTuple *tmp_candis = (TwoTuple*) calloc(num_tuples, sizeof(TwoTuple));
  char *result_string = (char*) calloc(result_string_length, sizeof(char));
  Read read;

  while (pop_read_queue(&read) != -1) {
    tmp_num_reads++;
    int tmpMappedReadNum = 0;
    tmp_num_candidates = generate_candidates(&read, 0, &candidates, &swap_candidates, &max_num_candidates, &candis, &tmp_candis, &num_tuples, &tmp_num_candidates_without_additional_qgram_filter);
    total_num_candidates_without_additional_qgram_filter += tmp_num_candidates_without_additional_qgram_filter;

    if (tmp_num_candidates > 0) {
      total_num_candidates += tmp_num_candidates;
      tmpMappedReadNum += verify_candidates(&read, 0, candidates, tmp_num_candidates, &no_left_result, &result_string_length, &result_string);
    }

    tmp_num_candidates = generate_candidates(&read, 1, &candidates, &swap_candidates, &max_num_candidates, &candis, &tmp_candis, &num_tuples, &tmp_num_candidates_without_additional_qgram_filter);
    total_num_candidates_without_additional_qgram_filter += tmp_num_candidates_without_additional_qgram_filter;

    if (tmp_num_candidates > 0) {
      total_num_candidates += tmp_num_candidates;
      tmpMappedReadNum += verify_candidates(&read, 1, candidates, tmp_num_candidates, &no_left_result, &result_string_length, &result_string);
    }

    if (tmpMappedReadNum>0) {
      mappingNum += tmpMappedReadNum; 
      tmp_num_mapped_reads++;
    }
  }

  //fprintf(stderr, "Thread %d is pushing local results to output queue...\n", thread_id);

  while (no_left_result != 1) {
    if (push_output_queue(result_string)) {
      no_left_result = 1;
      //free(result_string);
      result_string = NULL;
    }
  }

  //fprintf(stderr, "Thread %d pushed local results to output queue successfully.\n", thread_id);

  free(candidates);
  free(swap_candidates);
  free(candis);
  free(tmp_candis);
  num_candidates_without_filter[thread_id] = total_num_candidates_without_additional_qgram_filter;
  num_candidates[thread_id] = total_num_candidates;
  num_mappings[thread_id] = mappingNum;
  num_reads[thread_id] = tmp_num_reads;
  num_mapped_reads[thread_id] = tmp_num_mapped_reads;
  fprintf(stderr, "Thread %d completed.\n", thread_id);

  return mappingNum;
}


void *start_CPU_thread(void *arg) {
  int thread_id = *((int*) arg);
  CPU_map(thread_id);
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
