#ifndef MAPPER_H_
#define MAPPER_H_

#include "reftools.h"
#include "readtools.h"
#include "utils.h"
#include "indextools.h"
#include "filter.h"
#include "verifier.h"

extern int kmer_size;
extern int step_size;
extern int error_threshold;
extern int num_threads;
extern int num_additional_qgrams;
extern int num_finished_threads;

extern uint32_t num_candidates_without_filter[256];
extern uint32_t num_mapped_reads[256];
extern uint32_t num_candidates[256];
extern uint32_t num_mappings[256];
extern uint32_t num_reads[256];
extern uint32_t (*generate_candidates)(const Read *read, const int is_reverse_complement, uint32_t **candidates, uint32_t **swap_candidates, uint32_t *max_num_candidates, TwoTuple **candis, TwoTuple **temp_candis, uint32_t *num_tuples, uint32_t *num_candidates_without_additonal_qgram_filter);

void initialize_mapper();
void finalize_mapper();
int CPU_map(int threadID);
void *start_CPU_thread(void *arg);
void *start_single_read_queue_thread(void *arg);
void *start_output_queue_thread(void *arg);

#endif /* MAPPER_H_ */
