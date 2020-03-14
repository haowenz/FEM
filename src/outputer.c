#include "outputer.h"

FILE *output_file;
FILE *header_file;
char *output_file_path;
char *header_file_path;
pthread_mutex_t output_mutex;
pthread_mutex_t output_queue_mutex;
pthread_cond_t output_queue_pro_cond;
pthread_cond_t output_queue_cos_cond;
pthread_cond_t not_empty_cond;

int output_queue_size;
int output_queue_size_max = 100000;

int output_queue_front;
int output_queue_rear;
/*output queue*/
char** output_queue;

void initialize_output_queue() {
  pthread_mutex_init(&output_queue_mutex, NULL);
  pthread_cond_init(&output_queue_pro_cond, NULL);
  output_queue_size = 0;
  output_queue = (char**) malloc(output_queue_size_max * sizeof(char*));
}

void destroy_output_queue() {
  pthread_mutex_destroy(&output_queue_mutex);
  pthread_cond_destroy(&output_queue_pro_cond);
  free(output_queue);
}

int push_output_queue(char* result) {
  int flag = pthread_mutex_trylock(&output_queue_mutex);
  if (flag != 0) {
    return 0;
  }

  if (output_queue_size >= output_queue_size_max) {
    pthread_mutex_unlock(&output_queue_mutex);
    pthread_cond_signal(&output_queue_pro_cond);
    return 0;
  }

  output_queue[output_queue_rear] = result;
  output_queue_rear = (output_queue_rear + 1) % output_queue_size_max;
  ++output_queue_size;
  pthread_mutex_unlock(&output_queue_mutex);
  pthread_cond_signal(&output_queue_pro_cond);

  return 1;
}

void output_queue_thread() {
  while (num_finished_threads < num_threads) {
    char *result_string;
    pthread_mutex_lock(&output_queue_mutex);
    while (output_queue_size == 0) {
      pthread_cond_wait(&output_queue_pro_cond, &output_queue_mutex);
      if (num_finished_threads == num_threads) {
        pthread_mutex_unlock(&output_queue_mutex);
        break;
      }
    }

    result_string = output_queue[output_queue_front];
    output_queue[output_queue_front] = NULL;
    output_queue_front = (output_queue_front + 1) % output_queue_size_max;
    --output_queue_size;
    pthread_mutex_unlock(&output_queue_mutex);

    if (result_string != NULL) {
      fwrite(result_string, sizeof(char), strlen(result_string), output_file);
      free(result_string);
      result_string = NULL;
    }
  }

  while (output_queue_size > 0) {
    fwrite(output_queue[output_queue_front], sizeof(char), strlen(output_queue[output_queue_front]), output_file);
    free(output_queue[output_queue_front]);
    output_queue[output_queue_front] = NULL;
    output_queue_front = (output_queue_front + 1) % output_queue_size_max;
    --output_queue_size;
    //printf("queue_size=%d.\n", output_queue_size);
  }
}

void initialize_output_file() {
  //outputMutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_init(&output_mutex, NULL);
  output_file = fopen(output_file_path, "w");
  if (output_file == NULL) {
    fprintf(stderr, "Cannot open output file.\n");
    exit(EXIT_FAILURE);
  }
}

void output_header_from_index_file() {
  int buffer_size = 4096;
  char buffer[buffer_size];
  header_file = fopen(header_file_path, "r");
  int num = 0;
  do {
    num = fread(buffer, sizeof(char), buffer_size, header_file);
    fwrite(buffer, sizeof(char), num, output_file);
  } while (num > 0);
  fclose(header_file);
}

int output_string(const char* result) {
  int out = pthread_mutex_trylock(&output_mutex);
  if (out == 0) {
    fwrite(result, sizeof(char), strlen(result), output_file);
    pthread_mutex_unlock(&output_mutex);
    return 1;
  } else {
    return 0;
  }
}

void finalize_output_file() {
  fclose(output_file);
}
