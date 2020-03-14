#include "FEM-index.h"

static inline void print_usage() {
  fprintf(stderr, "Usage: FEM index <window_size> <step_size> <ref>\n");
}

int index_main(int argc, char* argv[]) {
  if (argc < 4) {
    fprintf(stderr, "%s\n", "Too few args!");
    print_usage();
    exit(EXIT_FAILURE);
  }

  kmer_size = atoi(argv[1]);
  step_size = atoi(argv[2]);
  reference_file_path = argv[3];
  /* TODO:add arguments check */ 

  initialize_ref(&reference);

  char tmp_index_file_path[128]; // TODO: remove hard-coded length
  tmp_index_file_path[0] = '\0';
  strcpy(tmp_index_file_path, reference_file_path);
  sprintf(tmp_index_file_path + strlen(tmp_index_file_path), ".fem");
  index_file_path = tmp_index_file_path;
  char tmp_header_file_path[128];
  tmp_header_file_path[0] = '\0';
  strcpy(tmp_header_file_path, index_file_path);
  sprintf(tmp_header_file_path + strlen(tmp_header_file_path), ".header");
  header_file_path = tmp_header_file_path;

  initialize_ref_file();
  double real_start_time = get_real_time();
  get_ref(&reference);
  fprintf(stderr, "Reference length: %u. Loaded into memory in %fs.\n", reference.lookupTable[reference.refNum], get_real_time() - real_start_time);

  real_start_time = get_real_time();
  initialize_saving_index();
  construct_index();
  save_index();
  fprintf(stderr, "Number of bases: %u, number of chromosomes: %d.\n Built index in %fs.\n", reference.lookupTable[reference.refNum], reference.refNum, get_real_time() - real_start_time);
  destroy_ref(&reference);
  finalize_ref_file();
  finalize_saving_index();
  return 0;
}
