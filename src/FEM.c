#include <stdio.h>
#include <string.h>

#include "FEM_index.h"
#include "FEM_map.h"
#include "utils.h"

#define FEM_VERSION "0.1"

static inline void print_usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: FEM (Fast and Efficient short read Mapper)\n");
  fprintf(stderr, "Version: %s\n", FEM_VERSION);
  fprintf(stderr, "Contact: Haowen Zhang <hwzhang@gatech.edu>\n\n");
  fprintf(stderr, "Usage:   FEM <command> [options]\n\n");
  fprintf(stderr, "Command: index   build index for reference\n");
  fprintf(stderr, "         map     map reads\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: To use FEM, you need to first index the genome with `FEM index'.\n\n");
}


int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "%s\n", "Too few arguements.");
    print_usage();
    exit(EXIT_FAILURE);
  }

  int return_value = 0;
  double real_start_time = get_real_time();
  double cpu_start_time = get_cpu_time();
  if (strcmp(argv[1], "index") == 0) {
    return_value = index_main(argc - 1, argv + 1);
  } else if (strcmp(argv[1], "map") == 0) {
    return_value = map_main(argc - 1, argv + 1);
  } else {
    fprintf(stderr, "[%s] unrecognized command '%s'\n", __func__, argv[1]);
    exit(EXIT_FAILURE);
  }

  if (return_value == 0) {
    fprintf(stderr, "[%s] Version: %s\n", __func__, FEM_VERSION);
    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
      fprintf(stderr, " %s", argv[i]);
    }
    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, get_real_time() - real_start_time, get_cpu_time() - cpu_start_time);
  }
  return return_value;
}
