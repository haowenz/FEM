#ifndef INDEX_H_
#define INDEX_H_

#include "sequence_batch.h"
#include "utils.h"

typedef struct {
  int kmer_size;
  int step_size;
  FILE *index_file;
  uint32_t *lookup_table;
  size_t occurrence_table_size;
  uint64_t *occurrence_table;
} Index;

void initialize_index(Index *index);
void destroy_index(Index *index);
void construct_index(const SequenceBatch *sequence_batch, Index *index);
void load_index(const char *index_file_path, Index *index);
void save_index(const char *index_file_path, Index *index); 

static inline uint32_t get_seed_frequency(const Index *index, uint32_t hash_value) {
  return index->lookup_table[hash_value + 1] - index->lookup_table[hash_value];
}

static inline uint64_t* get_seed_occurrences(const Index *index, uint32_t hash_value) {
  return &(index->occurrence_table[index->lookup_table[hash_value]]);
}
//extern char *index_file_path;
//extern int max_hash_value;
//extern uint32_t hash_table_size;
//extern uint32_t *occurrence_table;
//extern int *lookup_table;


//void initialize_loading_index(Index *index);
//void finalize_loading_index(Index *index);

//void initialize_saving_header();
//void finalize_saving_header();

//void initialize_saving_index(Index *index); 
//void finalize_saving_index(Index *index); 
#endif // INDEX_H_
