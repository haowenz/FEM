#ifndef OUTPUTER_H_
#define OUTPUTER_H_

#include "utils.h"
#include "readtools.h"
#include "reftools.h"
#include "mapper.h"
#include "pthread.h"

#include <sys/queue.h>

extern char *result_file_name;
extern char *header_file_name;
extern pthread_cond_t output_queue_pro_cond;

enum SAMFlags {
  READ_UNMAPPED = 4,
  READ_REVERSE_MAPPED = 16,
  READ_PAIR_UNMAPPED1 = 77,
  READ_PAIR_UNMAPPED2 = 141,
  READ_PAIR_FORWARD_MAPPED1 = 99,
  READ_PAIR_REVERSE_MAPPED1 = 83,
  READ_PAIR_FORWARD_MAPPED2 = 147,
  READ_PAIR_REVERSE_MAPPED2 = 163,
  READ_PAIR_FIRST_FORWARD_MAPPED = 0,
  READ_PAIR_FIRST_REVERSE_MAPPED = 16,
  READ_PAIR_FIRST_UNMAPPED = 4,
  READ_PAIR_SECOND_FORWARD_MAPPED = 0,
  READ_PAIR_SECOND_REVERSE_MAPPED = 16,
  READ_PAIR_SECOND_UNMAPPED = 4,

  NOT_PRIMARY_ALIGN = 256
};

typedef struct {
  char resultStr[1024];
} OutputResult;

void initialize_output_queue();
void destroy_output_queue() ;
int push_output_queue(char* result);
void initialize_output_file();
void output_header_from_index_file();
int output_string(const char* result);
void finalize_output_file();
void output_queue_thread();

#endif /* OUTPUTER_H_ */
