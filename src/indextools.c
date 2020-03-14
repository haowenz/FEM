#include "indextools.h"

uint32_t *occurrence_table;
int *lookup_table;

char *index_file_path;
FILE *index_file;

FILE *header_file;

int max_hash_value;
uint32_t hash_table_size;

void initialize_loading_index() {
  index_file = fopen(index_file_path, "r");
  if (index_file == NULL) {
    fprintf(stderr, "%s %s\n", "Failed to open index file", index_file_path);
    exit(EXIT_FAILURE);
  }
  size_t num_read_elements = fread(&(reference.refNum), sizeof(int), 1, index_file);
  if (num_read_elements != 1) {
    fprintf(stderr, "Load error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
  num_read_elements = fread(&kmer_size, sizeof(int), 1, index_file);
  if (num_read_elements != 1) {
    fprintf(stderr, "Load error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
  num_read_elements = fread(&step_size, sizeof(int), 1, index_file);
  if (num_read_elements != 1) {
    fprintf(stderr, "Load error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
  max_hash_value = (int)(pow(4, kmer_size)) - 1;
  lookup_table = (int*) malloc(sizeof(int) * (max_hash_value + 2));
  assert(lookup_table);
}

void load_index() {
  double real_start_time = get_real_time();
  size_t num_read_elements = 0; //TODO: check fread return value
  for (int reference_sequence_index = 0; reference_sequence_index < reference.refNum; ++reference_sequence_index) {
    char* reference_sequence_name = reference.names[reference_sequence_index];
    num_read_elements = fread(reference_sequence_name, sizeof(char), REF_NAME_LEN_MAX, index_file);
  }
  //load the length of reference bases
  uint32_t *reference_lookup_table = reference.lookupTable;
  num_read_elements = fread(reference_lookup_table, sizeof(uint32_t), reference.refNum + 1, index_file);
  //load the int bases of reference
  num_read_elements = fread(reference.bases, sizeof(uint8_t), reference_lookup_table[reference.refNum], index_file);

  //load the hash table
  //load the lookup table
  num_read_elements = fread(lookup_table, sizeof(int), max_hash_value + 2, index_file);
  //load the occurrence table
  num_read_elements = fread(&hash_table_size, sizeof(int), 1, index_file);
  assert(num_read_elements == 1);

  occurrence_table = (uint32_t*) malloc(sizeof(uint32_t) * hash_table_size);
  assert(occurrence_table);

  num_read_elements = fread(occurrence_table, sizeof(uint32_t), hash_table_size, index_file);

  fprintf(stderr, "Max hash value: %d, number of reference sequences: %d, number of bases: %d, load in %fs!\n", max_hash_value, reference.refNum, reference.lookupTable[reference.refNum] - reference.lookupTable[0],  get_real_time() - real_start_time);
}

void finalize_loading_index() {
  if (occurrence_table != NULL) {
    free(occurrence_table);
    occurrence_table = NULL;
  }
  if (lookup_table != NULL) {
    free(lookup_table);
    lookup_table = NULL;
  }
  fclose(index_file);
}

void initialize_saving_index() {
  size_t num_written_elements = 0;
  index_file = fopen(index_file_path, "w");
  int refNum = reference.refNum;
  fprintf(stderr, "Number of reference sequences: %d.\n", refNum);
  num_written_elements = fwrite(&refNum, sizeof(int), 1, index_file);
  num_written_elements = fwrite(&kmer_size, sizeof(int), 1, index_file);
  num_written_elements = fwrite(&step_size, sizeof(int), 1, index_file);
  if (num_written_elements == 0) {
    fprintf(stderr, "Write error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
  hash_table_size = 0;
  for (int i = 0; i < reference.refNum; ++i) {
    hash_table_size += (reference.lookupTable[i+1] - reference.lookupTable[i] - kmer_size + 1) / step_size;
  }
  max_hash_value = (int)(pow(4, kmer_size)) - 1;

  lookup_table = (int*) malloc(sizeof(int) * (max_hash_value+2));
  assert(lookup_table);
  occurrence_table= (uint32_t*) malloc(sizeof(uint32_t) * hash_table_size);
  assert(occurrence_table);

  fprintf(stderr, "The size of hash table is %d, the max hash value is %d.\n", hash_table_size, max_hash_value);
}

__UTILS_HASH_VALUE
void construct_index() {
  double real_start_time = get_real_time();
  HashLocationPair* tmp_occurrence_table = (HashLocationPair*) malloc(sizeof(HashLocationPair) * hash_table_size);
  //Compute hash value of each element
  for (int j = 0; j < reference.refNum; ++j) {
    //#pragma omp parallel for
    for (uint32_t i = reference.lookupTable[j]; i < reference.lookupTable[j + 1] - kmer_size + 1; i += step_size) {
      tmp_occurrence_table[i / step_size].location = i;
      int hashVal = hashValue(reference.bases + i, kmer_size);
      if (hashVal != -1) {
        tmp_occurrence_table[i / step_size].hashValue = hashVal;
      } else {
        tmp_occurrence_table[i / step_size].hashValue = max_hash_value + 1;
      }
    }
  }
  //Sort
  //tbb::parallel_sort
  qsort(tmp_occurrence_table, hash_table_size, sizeof(HashLocationPair), compare_hash_location_pair);
  //build lookup table
  //lookup table gives the start of each list in the occurrence table
  memset(lookup_table, 0, (max_hash_value + 2) * sizeof(int));
  for (uint32_t i = 0; i < hash_table_size; i++) {
    occurrence_table[i] = tmp_occurrence_table[i].location;
    int currentHashVal = tmp_occurrence_table[i].hashValue + 1;
    if (currentHashVal <= max_hash_value + 1) {
      lookup_table[currentHashVal]++;
    }else{
      break;
    }
  }
  free(tmp_occurrence_table);

  int sum = 0;
  for (int i = 0; i <= max_hash_value + 1; i++) {
    sum += lookup_table[i];
    lookup_table[i] = sum;
  }

  hash_table_size = sum;

  fprintf(stderr, "Number of hashed bases: %d\n", lookup_table[max_hash_value]);
  fprintf(stderr, "Generated hash table in %fs.\nMax hash value: %d, hash table size: %u.\n", get_real_time() - real_start_time, max_hash_value, hash_table_size);
}

void save_index() {
  size_t num_written_elements = 0;
  //store the reference
  //store the length of reference name
  //store the name of refenrence
  initialize_saving_header();
  char hd[1000] = "@HD\tVN:1.4\tSO:unsorted\n";
  num_written_elements = fwrite(hd, sizeof(char), strlen(hd), header_file);
  for (int i = 0; i < reference.refNum; ++i) {
    char* reference_sequence_name = reference.names[i];
    int reference_sequence_length = reference.lookupTable[i + 1] - reference.lookupTable[i];
    int length = sprintf(hd, "@SQ\tSN:%s\tLN:%d\n", reference_sequence_name, reference_sequence_length);
    num_written_elements = fwrite(hd, sizeof(char), length, header_file);
    num_written_elements = fwrite(reference_sequence_name, sizeof(char), REF_NAME_LEN_MAX, index_file);
  }
  finalize_saving_header();
  uint32_t *reference_lookup_table = reference.lookupTable;
  num_written_elements = fwrite(reference_lookup_table, sizeof(uint32_t), reference.refNum + 1, index_file);
  //store the int bases of reference
  num_written_elements = fwrite(reference.bases, sizeof(uint8_t), reference_lookup_table[reference.refNum], index_file);
  //store the hash table
  //store the lookup table
  num_written_elements = fwrite(lookup_table, sizeof(int), max_hash_value + 2, index_file);
  //store the occurrence table
  num_written_elements = fwrite(&hash_table_size, sizeof(uint32_t), 1, index_file);
  num_written_elements = fwrite(occurrence_table, sizeof(uint32_t), hash_table_size, index_file);
  if (num_written_elements == 0){
    fprintf(stderr, "Write error while initializing hash table.\n");
    exit(EXIT_FAILURE);
  }
}

void finalize_saving_index() {
  fclose(index_file);
}

void initialize_saving_header() {
  header_file = fopen(header_file_path, "w");
  if (header_file == NULL){
    fprintf(stderr, "Cannot open header file.\n");
    exit(EXIT_FAILURE);
  }
}

void save_header() {
}

void finalize_saving_header() {
  fclose(header_file);
}
