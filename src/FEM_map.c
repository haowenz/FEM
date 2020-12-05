#include "FEM_map.h"

#include <assert.h>
#include <getopt.h>

#include "map.h"
#include "input_queue.h"
#include "output_queue.h"

static inline void print_usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:  FEM map [options] \n\n");
  fprintf(stderr, "Options:");
  fprintf(stderr, "\n");
  fprintf(stderr, "        -e       INT  error threshold \n");
  fprintf(stderr, "        -t       INT  number of threads \n");
  fprintf(stderr, "        -f       STR  seeding algorithm: \"g\" for group seeding and \"v\" for variable-length seeding \n");
  fprintf(stderr, "        -a       INT  # additional q-grams (only for test)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input/output: ");
  fprintf(stderr, "\n");
  fprintf(stderr, "        --ref    STR  Input reference file\n");
  fprintf(stderr, "        --index  STR  Input index file\n");
  fprintf(stderr, "        --read1  STR  Input read1 file\n");
  fprintf(stderr, "        -o       STR  Output SAM file \n");
  fprintf(stderr, "\n");
}

static inline int check_args(const FEMArgs *fem_args, const char *reference_file_path, const char *index_file_path, const char *read1_file_path, const char *output_file_path) {
  if (fem_args->error_threshold < 0 || fem_args->error_threshold > 7) {
    fprintf(stderr, "%s\n", "Wrong error threshold.");
    return 0;
  } 
  if (fem_args->num_threads <= 0) {
    fprintf(stderr, "%s\n", "Wrong number of threads.");
    return 0;
  } 
  if (fem_args->num_additional_qgrams < 0 || fem_args->num_additional_qgrams > 2) {
    fprintf(stderr, "%s\n", "Wrong number of additional q-grams.");
    return 0;
  }
  if (reference_file_path == NULL) {
    fprintf(stderr, "%s\n", "Reference file path is required.");
    return 0;
  } 
  if (read1_file_path == NULL) {
    fprintf(stderr, "%s\n", "Read file path is required.");
    return 0;
  } 
  if (output_file_path == NULL) {
    fprintf(stderr, "%s\n", "Output file path is required.");
    return 0;
  }
  return 1;
}

int map_main(int argc, char *argv[]) {
  {
  int a = 2, b = 3, c = 5, d = 9, e = 17, f = 33, g = 139;
  fprintf(stderr, "%d %d %d %d %d %d %d\n", kv_roundup32(a), kv_roundup32(b), kv_roundup32(c), kv_roundup32(d), kv_roundup32(e), kv_roundup32(f), kv_roundup32(g));
  }
  char *reference_file_path = NULL;
  char *index_file_path = NULL;
  char *read1_file_path = NULL;
  char *output_file_path = NULL;
  FEMArgs fem_args;
  fem_args.kmer_size = 12;
  fem_args.step_size = 3;
  fem_args.error_threshold = 2;
  fem_args.num_additional_qgrams = 1;
  fem_args.num_threads = 1;
  fem_args.seeding_method = 'g'; // "v" for variable length seeding, "g" for group seeding.

  //initialize_fem_args(&fem_args);
  // Parse args
  const char *short_opt = "ha:f:e:t:o:r:i:b:";
  struct option long_opt[] = 
  {
    {"help", no_argument, NULL, 'h'},
    {"ref", required_argument, NULL, 'r'},
    {"index", required_argument, NULL, 'i'},
    {"read1", required_argument, NULL,'b'}
  };
  int c, option_index;
  while((c = getopt_long(argc, argv, short_opt, long_opt, &option_index)) >= 0) {
    switch (c) {
      case 'r':
        reference_file_path = optarg;
        fprintf(stderr, "name: %s, ref: %s\n", long_opt[option_index].name, reference_file_path);
        break;
      case 'i':
        index_file_path = optarg;
        fprintf(stderr, "name: %s, index: %s\n", long_opt[option_index].name, index_file_path);
        break;
      case 'b':
        read1_file_path = optarg;
        fprintf(stderr, "name: %s, read1: %s\n", long_opt[option_index].name, read1_file_path);
        break;
      case 'e':
        fem_args.error_threshold = atoi(optarg);
        break;
      case 't':
        fem_args.num_threads = atoi(optarg);
        break;
      case 'a':
        fem_args.num_additional_qgrams = atoi(optarg);
        break;
      case 'f':
        if (strcmp(optarg, "v") == 0) {
          //generate_candidates = generate_variable_length_seeding_candidates;
        } else if(strcmp(optarg,"g") == 0) {
          //generate_candidates = generate_group_seeding_candidates;
        } else {
          fprintf(stderr, "%s\n", "Wrong name of seeding algorithm!");
          print_usage();
          exit(EXIT_FAILURE);
        }
        break;
      case 'o':
        output_file_path = optarg;
        fprintf(stderr, "output: %s\n", output_file_path);
        break;
      default:
        print_usage();
        exit(EXIT_SUCCESS);
    }
  }

  // Check args
  if (!check_args(&fem_args, reference_file_path, index_file_path, read1_file_path, output_file_path)) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  // Load reference
  SequenceBatch reference_sequence_batch;
  initialize_sequence_batch(&reference_sequence_batch);
  initialize_sequence_batch_loading(reference_file_path, &reference_sequence_batch);
  load_all_sequences_into_sequence_batch(&reference_sequence_batch);

  // Load index
  Index index;
  load_index(index_file_path, &index);

  pthread_t mapping_thread_handles[fem_args.num_threads];
  pthread_t input_queue_thread_handle;
  pthread_t output_queue_thread_handle;

  InputQueue input_queue;
  uint32_t input_queue_max_size = 100;
  uint32_t read_batch_max_size = 10000; 
  initialize_input_queue(read1_file_path, read_batch_max_size, input_queue_max_size, &input_queue);
  OutputQueue output_queue;
  uint32_t output_queue_max_size = 100000;
  initialize_output_queue(output_file_path, &reference_sequence_batch, fem_args.num_threads, output_queue_max_size, &output_queue);
  MappingArgs mapping_args[fem_args.num_threads];
  for (int i = 0; i < fem_args.num_threads; ++i) {
    mapping_args[i].thread_id = i;
    mapping_args[i].max_read_batch_size = read_batch_max_size;
    mapping_args[i].fem_args = &fem_args;
    mapping_args[i].reference_sequence_batch = &reference_sequence_batch;
    mapping_args[i].index = &index;
    mapping_args[i].input_queue = &input_queue;
    mapping_args[i].output_queue = &output_queue;
    mapping_args[i].mapping_stats.num_reads = 0;
    mapping_args[i].mapping_stats.num_mapped_reads = 0;
    mapping_args[i].mapping_stats.num_candidates_without_additonal_qgram_filter = 0;
    mapping_args[i].mapping_stats.num_candidates = 0;
    mapping_args[i].mapping_stats.num_mappings = 0;
  }

  double startTime = get_real_time();
  int pthread_err = 0;
  pthread_err = pthread_create(&input_queue_thread_handle, NULL, input_queue_thread, &input_queue);
  if (pthread_err == 0) {
    fprintf(stderr, "Created input queue queue successfully.\n");
  }
  pthread_err = pthread_create(&output_queue_thread_handle, NULL, output_queue_thread, &output_queue);
  if (pthread_err == 0) {
    fprintf(stderr, "Created output queue successfully.\n");
  }
  for (int i = 0; i < fem_args.num_threads; ++i) {
    pthread_err = pthread_create(mapping_thread_handles + i, NULL, single_end_read_mapping_thread, &(mapping_args[i]));
    assert(pthread_err == 0);
  }
  for (int i = 0; i < fem_args.num_threads; ++i) {
    pthread_err = pthread_join(mapping_thread_handles[i], NULL);
    assert(pthread_err == 0);
  }
  assert(pthread_err == 0);
  pthread_err = pthread_join(output_queue_thread_handle, NULL);
  if (pthread_err == 0) {
    fprintf(stderr, "Output queue thread joint successfully.\n");
  }
  pthread_err = pthread_join(input_queue_thread_handle, NULL);
  if (pthread_err == 0) {
    fprintf(stderr, "Input queue thread joint successfully.\n");
  }

  // Generate mapping statistics
  uint64_t num_reads = 0;
  uint64_t num_mapped_reads = 0;;
  uint64_t num_candidates_without_additonal_qgram_filter = 0;
  uint64_t num_candidates = 0;
  uint64_t num_mappings = 0;
  for (int i = 0; i < fem_args.num_threads; ++i) {
    num_reads += mapping_args[i].mapping_stats.num_reads;
    num_mapped_reads += mapping_args[i].mapping_stats.num_mapped_reads;
    num_candidates_without_additonal_qgram_filter += mapping_args[i].mapping_stats.num_candidates_without_additonal_qgram_filter;
    num_candidates += mapping_args[i].mapping_stats.num_candidates;
    num_mappings += mapping_args[i].mapping_stats.num_mappings;
  }

  fprintf(stderr, "The number of read: %"PRIu64"\n", num_reads);
  fprintf(stderr, "The number of mapped read: %"PRIu64"\n", num_mapped_reads);
  fprintf(stderr, "The number of candidate before additional q-gram filter: %"PRIu64"\n", num_candidates_without_additonal_qgram_filter);
  fprintf(stderr, "The number of candidate: %"PRIu64"\n", num_candidates);
  fprintf(stderr, "The number of mapping: %"PRIu64"\n", num_mappings);
  fprintf(stderr, "Time: %fs\n", get_real_time() - startTime);

  destroy_output_queue(&output_queue);
  destroy_input_queue(&input_queue);
  destroy_index(&index);
  finalize_sequence_batch_loading(&reference_sequence_batch);
  destory_sequence_batch(&reference_sequence_batch);
  return 0;
}
