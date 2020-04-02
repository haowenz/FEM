#ifndef INPUTQUEUE_H_
#define INPUTQUEUE_H_

#include <assert.h>

#include "sequence_batch.h"
#include "utils.h"

typedef struct {
  size_t max_size;
  size_t front;
  size_t rear;
  size_t size;
  int no_more_read;
  pthread_mutex_t queue_mutex;
  pthread_cond_t pro_cond;
  pthread_cond_t con_cond;
  SequenceBatch read1_batch_for_loading;
  SequenceBatch read2_batch_for_loading;
  SequenceBatch *batches;
} InputQueue;

void initialize_input_queue(const char *read1_file_path, uint32_t max_batch_size, size_t max_queue_size, InputQueue *intput_queue);
void destroy_input_queue(InputQueue *input_queue);

void initialize_input_queue_read_file(InputQueue *input_queue);
void finalize_input_queue_read_file(InputQueue *input_queue);

void *input_queue_thread(void *input_queue);
void pop_input_queue(SequenceBatch *read1_batch, InputQueue *input_queue);

#endif // INPUTQUEUE_H_
