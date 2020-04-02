#ifndef MAP_H_
#define MAP_H_

#include "align.h"
#include "filter.h"
#include "index.h"
#include "input_queue.h"
#include "output_queue.h"
#include "sequence_batch.h"
#include "utils.h"

//typedef struct {
//  int kmer_size;
//  int step_size;
//  int error_threshold;
//  int num_additional_qgrams;
//  int num_threads;
//  char seeding_method; // "v" for variable length seeding, "g" for group seeding.
//} FEMArgs;
//
//void initialize_fem_args(FEMArgs *fem_args) {
//  fem_args.kmer_size = 12;
//  fem_args.step_size = 3;
//  fem_args.error_threshold = 2;
//  fem_args.num_additional_qgrams = 1;
//  fem_args.num_threads = 1;
//  fem_args.seeding_method = 'g'; // "v" for variable length seeding, "g" for group seeding.
//}
//
//void check_fem_args(FEMArgs *fem_args) {
//  if (fem_args.error_threshold < 0 || fem_args.error_threshold > 7) {
//    fprintf(stderr, "%s\n", "Wrong error threshold.");
//    exit(EXIT_FAILURE);
//  } else if (num_threads <= 0) {
//    fprintf(stderr, "%s\n", "Wrong number of threads.");
//    exit(EXIT_FAILURE);
//  } else if (num_additional_qgrams < 0 || num_additional_qgrams > 2) {
//    fprintf(stderr, "%s\n", "Wrong number of additional q-grams.");
//    exit(EXIT_FAILURE);
//  }
//}

typedef struct {
  uint64_t num_reads;
  uint64_t num_mapped_reads;
  uint64_t num_candidates_without_additonal_qgram_filter;
  uint64_t num_candidates;
  uint64_t num_mappings;
} MappingStats;

typedef struct {
  int thread_id;
  uint32_t max_read_batch_size;
  FEMArgs *fem_args;
  SequenceBatch *reference_sequence_batch;
  Index *index;
  InputQueue *input_queue;
  OutputQueue *output_queue;
  MappingStats mapping_stats;
} MappingArgs;

//void initialize_mapping_args(MappingArgs *mapping_args, int thread_id, const FEMArgs *fem_args, const SequenceBatch *reference_sequence_batch, const Index *index, const InputQueue *input_queue, const OutputQueue *output_queue) {
//  mapping_args->thread_id = thread_id;
//  mapping_args->fem_args = fem_args;
//  mapping_args->reference_sequence_batch = reference_sequence_batch;
//  mapping_args->index = index;
//  mapping_args->input_queue = input_queue;
//  mapping_args->output_queue = output_queue;
//  mapping_args->mapping_stats = {0, 0, 0, 0, 0};
//}

void *single_end_read_mapping_thread(void *mapping_args);

#endif // MAP_H_
