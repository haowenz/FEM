#include "FEM_index.h"

#include "index.h"
#include "sequence_batch.h"
#include "utils.h"

static inline void print_usage() {
  fprintf(stderr, "Usage: FEM index <window_size> <step_size> <reference> <output>\n");
}

int index_main(int argc, char* argv[]) {
  if (argc < 5) {
    fprintf(stderr, "%s\n", "Too few args!");
    print_usage();
    exit(EXIT_FAILURE);
  }

  int kmer_size = atoi(argv[1]);
  int step_size = atoi(argv[2]);
  const char *reference_file_path = argv[3];
  const char *index_file_path = argv[4];
  fprintf(stderr, "k: %d, step size: %d, reference: %s, output: %s\n", kmer_size, step_size, reference_file_path, index_file_path);
  // TODO: check arguments

  SequenceBatch reference_sequence_batch;
  initialize_sequence_batch(&reference_sequence_batch);
  initialize_sequence_batch_loading(reference_file_path, &reference_sequence_batch);
  load_all_sequences_into_sequence_batch(&reference_sequence_batch);
  Index index;
  initialize_index(&index);
  index.kmer_size = kmer_size;
  index.step_size = step_size;
  construct_index(&reference_sequence_batch, &index);
  save_index(index_file_path, &index);
  destroy_index(&index);
  finalize_sequence_batch_loading(&reference_sequence_batch);
  destory_sequence_batch(&reference_sequence_batch);
  return 0;
}
