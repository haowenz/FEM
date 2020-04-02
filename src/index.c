#include "index.h"

#include <assert.h>

#include "ksort.h"

void initialize_index(Index *index) {
  index->index_file = NULL;
  index->lookup_table = NULL;
  index->occurrence_table = NULL;
}

void destroy_index(Index *index) {
  if (index->lookup_table != NULL) {
    free(index->lookup_table);
    index->lookup_table = NULL;
  }
  if (index->occurrence_table != NULL) {
    free(index->occurrence_table);
    index->occurrence_table = NULL;
  }
}

typedef struct {
  uint32_t hash_value;
  uint64_t location;
} HashTableEntry;

#define HashTableSortKey(m) ((m).hash_value)
KRADIX_SORT_INIT(hash_table, HashTableEntry, HashTableSortKey, 4);

#define OccurrenceTableSortKey(m) (m)
KRADIX_SORT_INIT(occurrence_table, uint64_t, OccurrenceTableSortKey, 8);

//static int compare_hash_table_entry(const void *a, const void *b) {
//  if (((HashTableEntry*)a)->hash_value < ((HashTableEntry*)b)->hash_value) {
//    return -1;
//  } else if (((HashTableEntry*)a)->hash_value == ((HashTableEntry*)b)->hash_value) {
//    return 0;
//  } else {
//    return 1;
//  }
//}
//
//int compare_uint64_t(const void * a, const void * b) {
//  uint64_t x = *(uint64_t*)a;
//  uint64_t y = *(uint64_t*)b;
//  if (x < y) {
//    return -1;
//  } else if (x == y) {
//    return 0;
//  } else {
//    return 1;
//  }
//}

void construct_index(const SequenceBatch *sequence_batch, Index *index) {
  double real_start_time = get_real_time();
  HashTableEntry* tmp_hash_table = (HashTableEntry*) malloc(sizeof(HashTableEntry) * (sequence_batch->num_bases / index->step_size + 1));
  size_t num_tmp_hash_table_entries = 0;
  //Compute hash value of each element
  for (uint32_t sequence_index = 0; sequence_index < sequence_batch->num_loaded_sequences; ++sequence_index) {
    uint32_t sequence_length = get_sequence_length_from_sequence_batch_at(sequence_batch, sequence_index);
    const char *sequence = get_sequence_from_sequence_batch_at(sequence_batch, sequence_index);
    for (uint32_t sequence_position = 0; sequence_position + index->kmer_size - 1 < sequence_length; sequence_position += index->step_size) {
      tmp_hash_table[num_tmp_hash_table_entries].hash_value = hash_seed_in_sequence(sequence_position, index->kmer_size, sequence, sequence_length);
      tmp_hash_table[num_tmp_hash_table_entries].location = ((uint64_t)sequence_index) << 32 | (uint32_t)sequence_position;
      ++num_tmp_hash_table_entries;
    }
  }
  fprintf(stderr, "Collected %ld seeds.\n", num_tmp_hash_table_entries);
  // Sort hash table
  //qsort(tmp_hash_table, num_tmp_hash_table_entries, sizeof(HashTableEntry), compare_hash_table_entry);
  radix_sort_hash_table(tmp_hash_table, tmp_hash_table + num_tmp_hash_table_entries);
  fprintf(stderr, "Sorted all the seeds.\n");
  // Build occ table and lookup table
  size_t lookup_table_size = (1 << (2 * index->kmer_size)) + 1;
  index->lookup_table = (uint32_t*)calloc(lookup_table_size, sizeof(uint32_t));
  assert(index->lookup_table);
  index->occurrence_table_size = num_tmp_hash_table_entries;
  index->occurrence_table = (uint64_t*)malloc(sizeof(uint64_t) * index->occurrence_table_size);
  assert(index->occurrence_table);
  for (uint32_t i = 0; i < num_tmp_hash_table_entries; i++) {
    index->lookup_table[tmp_hash_table[i].hash_value + 1]++;
    index->occurrence_table[i] = tmp_hash_table[i].location;
  }
  free(tmp_hash_table);
  size_t occurrence_table_size = 0;
  for (size_t i = 1; i < lookup_table_size; i++) {
    occurrence_table_size += index->lookup_table[i];
    //qsort(index->occurrence_table + index->lookup_table[i - 1], index->lookup_table[i], sizeof(uint64_t), compare_uint64_t);
    index->lookup_table[i] = occurrence_table_size;
    radix_sort_occurrence_table(index->occurrence_table + index->lookup_table[i - 1], index->occurrence_table + index->lookup_table[i]);
  }
  assert(occurrence_table_size == num_tmp_hash_table_entries);
  fprintf(stderr, "Lookup table size: %ld, occurrence table size: %ld.\n", lookup_table_size, occurrence_table_size);
  fprintf(stderr, "Built index in %fs.\n", get_real_time() - real_start_time);
}

void load_index(const char *index_file_path, Index *index) {
  double real_start_time = get_real_time();
  index->index_file = fopen(index_file_path, "rb");
  if (index->index_file == NULL) {
    fprintf(stderr, "Failed to open index file %s\n", index_file_path);
    exit(EXIT_FAILURE);
  }
  size_t num_read_elements = fread(&(index->kmer_size), sizeof(int), 1, index->index_file); //TODO: check fread return value
  if (num_read_elements != 1) {
    fprintf(stderr, "Load error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
  num_read_elements = fread(&(index->step_size), sizeof(int), 1, index->index_file);
  if (num_read_elements != 1) {
    fprintf(stderr, "Load error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
  size_t lookup_table_size = (1 << (2 * index->kmer_size)) + 1;
  index->lookup_table = (uint32_t*) malloc(sizeof(uint32_t) * lookup_table_size);
  assert(index->lookup_table);
  //load the hash table
  //load the lookup table
  num_read_elements = fread(index->lookup_table, sizeof(uint32_t), lookup_table_size, index->index_file);
  //load the occurrence table
  num_read_elements = fread(&(index->occurrence_table_size), sizeof(size_t), 1, index->index_file);
  assert(num_read_elements == 1);
  index->occurrence_table = (uint64_t*) malloc(sizeof(uint64_t) * index->occurrence_table_size);
  assert(index->occurrence_table);
  num_read_elements = fread(index->occurrence_table, sizeof(uint64_t), index->occurrence_table_size, index->index_file);
  fclose(index->index_file);
  fprintf(stderr, "Loaded index in %fs!\n", get_real_time() - real_start_time);
}

void save_index(const char *index_file_path, Index *index) {
  index->index_file = fopen(index_file_path, "wb");
  size_t num_written_elements = 0;
  num_written_elements = fwrite(&(index->kmer_size), sizeof(int), 1, index->index_file);
  num_written_elements = fwrite(&(index->step_size), sizeof(int), 1, index->index_file);
  if (num_written_elements == 0) {
    fprintf(stderr, "Write error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
  //store the reference
  //store the length of reference name
  //store the name of refenrence
  //initialize_saving_header();
  //char hd[1000] = "@HD\tVN:1.4\tSO:unsorted\n";
  //num_written_elements = fwrite(hd, sizeof(char), strlen(hd), header_file);
  //for (int i = 0; i < reference.refNum; ++i) {
  //  char* reference_sequence_name = reference.names[i];
  //  int reference_sequence_length = reference.lookupTable[i + 1] - reference.lookupTable[i];
  //  int length = sprintf(hd, "@SQ\tSN:%s\tLN:%d\n", reference_sequence_name, reference_sequence_length);
  //  num_written_elements = fwrite(hd, sizeof(char), length, header_file);
  //  num_written_elements = fwrite(reference_sequence_name, sizeof(char), REF_NAME_LEN_MAX, index_file);
  //}
  //finalize_saving_header();
  //store the hash table
  //store the lookup table
  size_t lookup_table_size = (1 << (2 * index->kmer_size)) + 1;
  num_written_elements = fwrite(index->lookup_table, sizeof(uint32_t), lookup_table_size, index->index_file);
  //store the occurrence table
  num_written_elements = fwrite(&(index->occurrence_table_size), sizeof(size_t), 1, index->index_file);
  num_written_elements = fwrite(index->occurrence_table, sizeof(uint64_t), index->occurrence_table_size, index->index_file);
  if (num_written_elements == 0){
    fprintf(stderr, "Write error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
  fclose(index->index_file);
}
