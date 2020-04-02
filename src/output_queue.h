#ifndef OUTPUTQUEUE_H_
#define OUTPUTQUEUE_H_

#include <pthread.h>

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
  FILE *output_file;
  kstring_t *output_kstrings;
} OutputQueue;

void initialize_output_queue(const char *output_file_path, int num_mapping_threads, size_t max_queue_size, OutputQueue *output_queue);
void destroy_output_queue(OutputQueue *output_queue) ;
void push_output_queue(kstring_t *result, OutputQueue *output_queue);
void *output_queue_thread(void *output_queue);

#endif // OUTPUTQUEUE_H_
