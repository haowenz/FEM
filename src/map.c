#include "map.h"

void *single_end_read_mapping_thread(void *mapping_args_v) {
  MappingArgs *mapping_args = (MappingArgs*)mapping_args_v;
  kvec_t_uint64_t candidates;
  kv_init(candidates.v);
  kvec_t_uint64_t buffer1;
  kv_init(buffer1.v);
  kvec_t_uint64_t buffer2;
  kv_init(buffer2.v);
  SequenceBatch read_batch;
  initialize_sequence_batch_with_max_size(mapping_args->max_read_batch_size, &read_batch);
  kstring_t result_kstring = {0, 0, NULL};
  while (1) {
    // Get read batch
    pop_input_queue(&read_batch, mapping_args->input_queue);
    // If no more reads, end the mapping loop
    if (read_batch.num_loaded_sequences == 0) {
      break;
    }
    double real_start_time = get_real_time();
    mapping_args->mapping_stats.num_reads += read_batch.num_loaded_sequences;
    // Generate candidates
    for (uint32_t read_index = 0; read_index < read_batch.num_loaded_sequences; ++read_index) {
      uint32_t num_candidates_without_additonal_qgram_filter = 0;
      uint32_t num_candidates = generate_group_seeding_candidates(mapping_args->fem_args, &read_batch, read_index, 0, mapping_args->reference_sequence_batch, mapping_args->index, &buffer1, &buffer2, &candidates, &num_candidates_without_additonal_qgram_filter);
      mapping_args->mapping_stats.num_candidates_without_additonal_qgram_filter += num_candidates_without_additonal_qgram_filter;
      mapping_args->mapping_stats.num_candidates += num_candidates;
      if (num_candidates > 0) {
        // Verify candidates
        uint32_t num_mappings = verify_candidates(mapping_args->fem_args, &read_batch, read_index, 0, mapping_args->reference_sequence_batch, candidates.v.a, num_candidates, &result_kstring);
        mapping_args->mapping_stats.num_mappings += num_mappings;
        if (num_mappings > 0) {
          ++(mapping_args->mapping_stats.num_mapped_reads);
          // Output mappings
          push_output_queue(&result_kstring, mapping_args->output_queue);
        }
      }
    }
    fprintf(stderr, "Mapped read batch in %fs.\n", get_real_time() - real_start_time);
  }
  pthread_mutex_lock(&(mapping_args->output_queue->queue_mutex));
  ++(mapping_args->output_queue->num_finished_mapping_threads);
  pthread_cond_signal(&(mapping_args->output_queue->pro_cond));
  pthread_mutex_unlock(&(mapping_args->output_queue->queue_mutex));
  destory_sequence_batch(&read_batch);
  assert(result_kstring.s == NULL);
  kv_destroy(candidates.v);
  kv_destroy(buffer1.v);
  kv_destroy(buffer2.v);
  fprintf(stderr, "Thread %d completed.\n", mapping_args->thread_id);
  return NULL;
}
