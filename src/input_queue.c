#include "input_queue.h"

void initialize_input_queue(const char *read1_file_path, uint32_t max_batch_size, size_t max_queue_size, InputQueue *input_queue) {
  input_queue->max_size = max_queue_size;
  input_queue->front = 0;
  input_queue->rear = 0;
  input_queue->size = 0;
  input_queue->no_more_read = 0;
  pthread_mutex_init(&(input_queue->queue_mutex), NULL);
  pthread_cond_init(&(input_queue->pro_cond), NULL);
  pthread_cond_init(&(input_queue->con_cond), NULL);
  initialize_sequence_batch_with_max_size(max_batch_size, &(input_queue->read1_batch_for_loading));
  initialize_sequence_batch_loading(read1_file_path, &(input_queue->read1_batch_for_loading));
  input_queue->batches = (SequenceBatch*)malloc(input_queue->max_size * sizeof(SequenceBatch));
  for (size_t i = 0; i < input_queue->max_size; ++i) {
    initialize_sequence_batch_with_max_size(max_batch_size, &(input_queue->batches[i]));
  }
  fprintf(stderr, "Initialize input queue successfully!\n");
}

void destroy_input_queue(InputQueue *input_queue) {
  pthread_mutex_destroy(&(input_queue->queue_mutex));
  pthread_cond_destroy(&(input_queue->pro_cond));
  pthread_cond_destroy(&(input_queue->con_cond));
  finalize_sequence_batch_loading(&(input_queue->read1_batch_for_loading));
  destory_sequence_batch(&(input_queue->read1_batch_for_loading));
  for (size_t i = 0; i < input_queue->max_size; ++i) {
    destory_sequence_batch(&(input_queue->batches[i]));
  }
  free(input_queue->batches);
  fprintf(stderr, "Destroy input queue successfully!\n");
}

void pop_input_queue(SequenceBatch *read1_batch, InputQueue *input_queue) {
  pthread_mutex_lock(&(input_queue->queue_mutex));
  while (input_queue->size == 0) {
    if (input_queue->no_more_read == 1) {
      pthread_mutex_unlock(&(input_queue->queue_mutex));
      read1_batch->num_loaded_sequences = 0;
      return;
    }
    // Check if queue is empty, if yes, wait for input thread
    pthread_cond_wait(&(input_queue->pro_cond), &(input_queue->queue_mutex));
  }
  // Take away one batch
  swap_sequences_in_sequence_batch(read1_batch, &(input_queue->batches[input_queue->front]));
  input_queue->front = (input_queue->front + 1) % input_queue->max_size;
  --(input_queue->size);
  pthread_cond_signal(&(input_queue->con_cond));
  pthread_mutex_unlock(&(input_queue->queue_mutex));
}

void *input_queue_thread(void *input_queue_v) {
  InputQueue *input_queue = (InputQueue*)input_queue_v;
  while (1) {
    // Load a batch
    load_batch_of_sequences_into_sequence_batch(&(input_queue->read1_batch_for_loading));
    if (input_queue->read1_batch_for_loading.num_loaded_sequences == 0) {
      // If no more reads, set no_more_read to 1
      pthread_mutex_lock(&(input_queue->queue_mutex));
      input_queue->no_more_read = 1;
      // Wake up all mapping threads so that they can proceed and end
      pthread_cond_broadcast(&(input_queue->pro_cond));
      pthread_mutex_unlock(&(input_queue->queue_mutex));
      return NULL;
    }
    pthread_mutex_lock(&(input_queue->queue_mutex));
    while (input_queue->size == input_queue->max_size) {
      // If the queue is full. wait for mapping thread to take batch
      pthread_cond_wait(&(input_queue->con_cond), &(input_queue->queue_mutex));
    }
    // Push a batch into the queue
    swap_sequences_in_sequence_batch(&(input_queue->read1_batch_for_loading), &(input_queue->batches[input_queue->rear]));
    input_queue->rear = (input_queue->rear + 1) % input_queue->max_size;
    ++(input_queue->size);
    pthread_cond_signal(&(input_queue->pro_cond));
    pthread_mutex_unlock(&(input_queue->queue_mutex));
  }
}
