#ifndef OUTPUTQUEUE_H_
#define OUTPUTQUEUE_H_

#include <pthread.h>

#include "htslib/sam.h"
#include "kstring.h"
#include "sequence_batch.h"
#include "utils.h"

typedef struct {
  size_t size;
  size_t front;
  size_t rear;
  size_t max_size;
  int num_mapping_threads;
  int num_finished_mapping_threads;
  pthread_mutex_t queue_mutex;
  pthread_cond_t pro_cond;
  pthread_cond_t con_cond;
  //FILE *output_file;
  samFile *output_sam_file;
  sam_hdr_t *sam_header;
  bam1_t **sam_alignments;
  //kstring_t *output_kstrings;
} OutputQueue;

void initialize_output_queue(const char *output_file_path, const SequenceBatch *sequence_batch, int num_mapping_threads, size_t max_queue_size, OutputQueue *output_queue);
void destroy_output_queue(OutputQueue *output_queue) ;
void push_output_queue(bam1_t **sam_alignment, OutputQueue *output_queue);
void *output_queue_thread(void *output_queue);
void output_sam_header(const char* output_file_path, const SequenceBatch *reference_sequence_batch, OutputQueue *output_queue);
void swap_bam1_t(bam1_t **a, bam1_t **b);

#endif // OUTPUTQUEUE_H_
