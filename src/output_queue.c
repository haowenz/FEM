#include "output_queue.h"

void initialize_output_queue(const char *output_file_path, int num_mapping_threads, size_t max_queue_size, OutputQueue *output_queue) {
  output_queue->size = 0;
  output_queue->max_size = max_queue_size;
  output_queue->front = 0;
  output_queue->rear = 0;
  output_queue->num_mapping_threads = num_mapping_threads;
  output_queue->num_finished_mapping_threads = 0;
  pthread_mutex_init(&(output_queue->queue_mutex), NULL);
  pthread_cond_init(&(output_queue->pro_cond), NULL);
  pthread_cond_init(&(output_queue->con_cond), NULL);
  output_queue->output_kstrings = (kstring_t*)malloc(output_queue->max_size * sizeof(kstring_t));
  for (size_t i = 0; i < output_queue->max_size; ++i) {
    (output_queue->output_kstrings)[i].l = 0;
    (output_queue->output_kstrings)[i].m = 0;
    (output_queue->output_kstrings)[i].s = NULL;
  }
  output_queue->output_file = fopen(output_file_path, "w");
  if (output_queue->output_file == NULL) {
    fprintf(stderr, "Cannot open output file.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "Initialize output queue successfully!\n");
}

void destroy_output_queue(OutputQueue *output_queue) {
  pthread_mutex_destroy(&(output_queue->queue_mutex));
  pthread_cond_destroy(&(output_queue->pro_cond));
  pthread_cond_destroy(&(output_queue->con_cond));
  if (output_queue != NULL) {
    for (size_t i = 0; i < output_queue->max_size; ++i) {
      if (output_queue->output_kstrings[i].s != NULL) {
        free(output_queue->output_kstrings[i].s);
      }
    }
    free(output_queue->output_kstrings);
  }
  fclose(output_queue->output_file);
  fprintf(stderr, "Destroy output queue successfully!\n");
}

void push_output_queue(kstring_t *result, OutputQueue *output_queue) {
  pthread_mutex_lock(&(output_queue->queue_mutex));
  // Wait if the queue is full
  while (output_queue->size == output_queue->max_size) {
    pthread_cond_wait(&(output_queue->con_cond), &(output_queue->queue_mutex));
  }
  // Push one result into the queue
  swap_kstring_t(&(output_queue->output_kstrings[output_queue->rear]), result);
  output_queue->rear = (output_queue->rear + 1) % output_queue->max_size;
  ++(output_queue->size);
  pthread_cond_signal(&(output_queue->pro_cond));
  pthread_mutex_unlock(&(output_queue->queue_mutex));
}

void *output_queue_thread(void *output_queue_v) {
  OutputQueue *output_queue = (OutputQueue*)output_queue_v;
  while (1) {
    //char *result_string;
    pthread_mutex_lock(&(output_queue->queue_mutex));
    while (output_queue->size == 0) {
      // If all mapping threads finished, then the output queue thread can stop
      if (output_queue->num_finished_mapping_threads == output_queue->num_mapping_threads) {
        pthread_mutex_unlock(&(output_queue->queue_mutex));
        return NULL;
      }
      // Note that when the mapping threads finished, they must signal output_queue_pro_cond, so that the output queue thread won't get stuck here forever
      pthread_cond_wait(&(output_queue->pro_cond), &(output_queue->queue_mutex));
    }
    // Take one result from queue
    kstring_t output_string = {0, 0, NULL};
    swap_kstring_t(&output_string, &(output_queue->output_kstrings[output_queue->front]));
    //output_queue[output_queue_front] = NULL;
    output_queue->front = (output_queue->front + 1) % output_queue->max_size;
    --(output_queue->size);
    pthread_cond_signal(&(output_queue->con_cond));
    pthread_mutex_unlock(&(output_queue->queue_mutex));
    // Save the result into the file
  }
}

//void initialize_output_file(const char *output_file_path, OutputQueue *output_queue) {
//  output_queue->output_file = fopen(output_file_path, "w");
//  if (output_queue->output_file == NULL) {
//    fprintf(stderr, "Cannot open output file.\n");
//    exit(EXIT_FAILURE);
//  }
//}
//
//void finalize_output_file(OutputQueue *output_queue) {
//  fclose(output_queue->output_file);
//}
//void output_header_from_index_file() {
//  int buffer_size = 4096;
//  char buffer[buffer_size];
//  header_file = fopen(header_file_path, "r");
//  int num = 0;
//  do {
//    num = fread(buffer, sizeof(char), buffer_size, header_file);
//    fwrite(buffer, sizeof(char), num, output_file);
//  } while (num > 0);
//  fclose(header_file);
//}
