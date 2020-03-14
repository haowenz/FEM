#include "FEM-align.h"

static inline void print_usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   FEM align [options] \n\n");
  fprintf(stderr, "Options:");
  fprintf(stderr, "\n");
  fprintf(stderr, "         -e        INT    error threshold \n");
  fprintf(stderr, "         -t        INT    number of threads \n");
  fprintf(stderr, "         -f        STR    seeding algorithm: \"g\" for group seeding and \"vl\" for variable-length seeding \n");
  fprintf(stderr, "         -a               select one additional q-gram (only for test)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input/output: ");
  fprintf(stderr, "\n");
  fprintf(stderr, "         --ref     STR    Input reference file\n");
  fprintf(stderr, "         --read    STR    Input read file\n");
  fprintf(stderr, "         -o        STR    Output SAM file \n");
  fprintf(stderr, "\n");
}

static inline int check_args() {
  if (error_threshold < 0 || error_threshold > 7) {
    fprintf(stderr, "%s\n", "Wrong error threshold.");
    return 0;
  } else if (num_threads <= 0) {
    fprintf(stderr, "%s\n", "Wrong number of threads.");
    return 0;
  } else if (num_additional_qgrams < 0 || num_additional_qgrams > 2) {
    fprintf(stderr, "%s\n", "The number of additional q-gram is too large.");
    return 0;
  } else if (!reference_file_path || !strcmp(reference_file_path,"")) {
    fprintf(stderr, "%s\n", "Reference file path is required.");
    return 0;
  } else if (!read1_file_path || !strcmp(read1_file_path,"")) {
    fprintf(stderr, "%s\n", "Read file path is required.");
    return 0;
  } else if (!output_file_path || !strcmp(output_file_path,"")) {
    fprintf(stderr, "%s\n", "Output file path is required.");
    return 0;
  } else {
    fprintf(stderr, "e: %d, t:%d, a: %d, ref: %s, read: %s, output: %s\n", error_threshold, num_threads, num_additional_qgrams, reference_file_path, read1_file_path, output_file_path);
    return 1;
  }
}

int align_main(int argc, char *argv[]) {
  // Parse args
  const char *short_opt = "ahf:e:t:o:r:i:";
  struct option long_opt[] = 
  {
    {"help", no_argument, NULL, 'h'},
    {"ref",required_argument, NULL, 'r'},
    {"read",required_argument, NULL,'i'}
  };
  int c, option_index;
  num_additional_qgrams = 0;
  num_threads = 1;
  while((c = getopt_long(argc, argv, short_opt, long_opt, &option_index)) >= 0) {
    switch(c) {
      case 'r':
        fprintf(stderr, "name: %s, ref: %s\n", long_opt[option_index].name, optarg);
        reference_file_path = optarg;
        break;
      case 'i':
        fprintf(stderr, "name: %s, read: %s\n", long_opt[option_index].name, optarg);
        read1_file_path = optarg;
        break;
      case 'e':
        error_threshold = atoi(optarg);
        break;
      case 't':
        num_threads = atoi(optarg);
        break;
      case 'a':
        num_additional_qgrams = 1;
        break;
      case 'f':
        if(strcmp(optarg, "vl") == 0) {
          generate_candidates = generate_variable_length_seeding_candidates;
        } else if(strcmp(optarg,"g") == 0) {
          generate_candidates = generate_group_seeding_candidates;
        } else {
          fprintf(stderr, "%s\n", "Wrong name of seeding algorithm!");
          print_usage();
          exit(EXIT_FAILURE);
        }
        break;
      case 'o':
        output_file_path = optarg;
        break;
      default:
        print_usage();
        exit(EXIT_SUCCESS);
    }
  }

  // Check args
  if (!check_args()) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  // Start mapping
  char temp_index_file_path[256];
  strcpy(temp_index_file_path, reference_file_path);
  sprintf(temp_index_file_path + strlen(temp_index_file_path), ".fem");
  index_file_path = temp_index_file_path;
  char temp_header_file_path[256];
  strcpy(temp_header_file_path, temp_index_file_path);
  sprintf(temp_header_file_path + strlen(temp_header_file_path), ".header");
  header_file_path = temp_header_file_path;
  pthread_t thread_handles[num_threads];
  pthread_t read_queue_thread_handle;
  pthread_t output_queue_thread_handle;
  initialize_mapper();
  load_index();

  double startTime = get_real_time();
  fprintf(stderr, "kmer size: %d, step size: %d.\n", kmer_size, step_size);
  num_finished_threads = 0;
  int thread_id_array[num_threads];
  for (int i = 0; i < num_threads; ++i) {
    num_candidates[i] = 0;
    num_candidates_without_filter[i] = 0;
    num_mappings[i] = 0;
    num_reads[i] = 0;
    num_mapped_reads[i] = 0;
    thread_id_array[i] = i;
  }

  int pthread_err = 0;
  pthread_err = pthread_create(&read_queue_thread_handle, NULL, start_single_read_queue_thread, NULL);
  if (pthread_err == 0) {
    fprintf(stderr, "Created read queue successfully.\n");
  }
  pthread_err = pthread_create(&output_queue_thread_handle, NULL, start_output_queue_thread, NULL);
  if (pthread_err == 0) {
    fprintf(stderr, "Created output queue successfully.\n");
  }
  for (int i = 0; i < num_threads; ++i) {
    pthread_err = pthread_create(thread_handles + i, NULL, start_CPU_thread, thread_id_array + i);
  }
  for (int i = 0; i < num_threads; ++i) {
    pthread_err = pthread_join(thread_handles[i], NULL);
    ++num_finished_threads;
  }
  pthread_err = pthread_cond_signal(&output_queue_pro_cond);
  assert(pthread_err == 0);
  pthread_err = pthread_join(output_queue_thread_handle, NULL);
  if (pthread_err == 0) {
    fprintf(stderr, "Output queue thread joint successfully.\n");
  }
  pthread_err = pthread_join(read_queue_thread_handle, NULL);
  if (pthread_err == 0) {
    fprintf(stderr, "Read queue thread joint successfully.\n");
  }

  // Generate mapping statistics
  uint32_t total_num_reads = 0;
  uint32_t total_num_mapped_reads =0;
  uint32_t total_num_candidates_without_filter = 0;
  uint32_t total_num_candidates = 0;
  uint32_t total_num_mappings = 0;
  for (int i = 0; i < num_threads; ++i) {
    total_num_reads += num_reads[i];
    total_num_mapped_reads += num_mapped_reads[i];
    total_num_candidates_without_filter += num_candidates_without_filter[i];
    total_num_candidates += num_candidates[i];
    total_num_mappings += num_mappings[i];
  }
  fprintf(stderr, "The number of read: %u.\n", total_num_reads);
  fprintf(stderr, "The number of mapped read: %u.\n", total_num_mapped_reads);
  fprintf(stderr, "The number of candidate before additional q-gram filter: %u.\n", total_num_candidates_without_filter);
  fprintf(stderr, "The number of candidate: %u.\n", total_num_candidates);
  fprintf(stderr, "The number of mapping: %u.\n", total_num_mappings);
  fprintf(stderr, "Time: %fs.\n", get_real_time() - startTime);

  finalize_mapper();
  return 0;
}
