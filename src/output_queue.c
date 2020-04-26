#include "output_queue.h"

void initialize_output_queue(const char *output_file_path, const SequenceBatch *reference_sequence_batch, int num_mapping_threads, size_t max_queue_size, OutputQueue *output_queue) {
  output_queue->size = 0;
  output_queue->max_size = max_queue_size;
  output_queue->front = 0;
  output_queue->rear = 0;
  output_queue->num_mapping_threads = num_mapping_threads;
  output_queue->num_finished_mapping_threads = 0;
  pthread_mutex_init(&(output_queue->queue_mutex), NULL);
  pthread_cond_init(&(output_queue->pro_cond), NULL);
  pthread_cond_init(&(output_queue->con_cond), NULL);
  output_queue->sam_alignment_kvec = (kvec_t_bam1_t_ptr*)malloc(output_queue->max_size * sizeof(kvec_t_bam1_t_ptr));
  for (size_t i = 0; i < output_queue->max_size; ++i) {
    kv_init(output_queue->sam_alignment_kvec[i].v);
  }
  output_queue->output_sam_file = sam_open_format(output_file_path, "w", NULL);
  output_queue->sam_header = sam_hdr_init();
  output_sam_header(output_file_path, reference_sequence_batch, output_queue);
  fprintf(stderr, "Initialize output queue successfully!\n");
}

void destroy_output_queue(OutputQueue *output_queue) {
  pthread_mutex_destroy(&(output_queue->queue_mutex));
  pthread_cond_destroy(&(output_queue->pro_cond));
  pthread_cond_destroy(&(output_queue->con_cond));
  // free alignments
  for (size_t i = 0; i < output_queue->max_size; ++i) {
    for (size_t si = 0; si < kv_size(output_queue->sam_alignment_kvec[i].v); ++si) {
      bam_destroy1(kv_A(output_queue->sam_alignment_kvec[i].v, si));
    }
    kv_destroy(output_queue->sam_alignment_kvec[i].v);
  }
  free(output_queue->sam_alignment_kvec);
  // free header
  //free(output_queue->sam_header->target_len);
  //free(output_queue->sam_header->target_name);
  for (uint32_t ri = 0; ri < output_queue->sam_header->n_targets; ++ri) {
    output_queue->sam_header->target_name[ri] = NULL;
  }
  sam_hdr_destroy(output_queue->sam_header);
  sam_close(output_queue->output_sam_file);
  fprintf(stderr, "Destroy output queue successfully!\n");
}

void push_output_queue(kvec_t_bam1_t_ptr *sam_alignment_kvec, OutputQueue *output_queue) {
  pthread_mutex_lock(&(output_queue->queue_mutex));
  // Wait if the queue is full
  while (output_queue->size == output_queue->max_size) {
    pthread_cond_wait(&(output_queue->con_cond), &(output_queue->queue_mutex));
  }
  // Push one result into the queue
  kv_swap(bam1_t*, sam_alignment_kvec->v, output_queue->sam_alignment_kvec[output_queue->rear].v);
  output_queue->rear = (output_queue->rear + 1) % output_queue->max_size;
  ++(output_queue->size);
  pthread_cond_signal(&(output_queue->pro_cond));
  pthread_mutex_unlock(&(output_queue->queue_mutex));
}

void *output_queue_thread(void *output_queue_v) {
  OutputQueue *output_queue = (OutputQueue*)output_queue_v;
  kvec_t_bam1_t_ptr sam_alignment_kvec;
  kv_init(sam_alignment_kvec.v);
  while (1) {
    pthread_mutex_lock(&(output_queue->queue_mutex));
    while (output_queue->size == 0) {
      // If all mapping threads finished, then the output queue thread can stop
      if (output_queue->num_finished_mapping_threads == output_queue->num_mapping_threads) {
        pthread_mutex_unlock(&(output_queue->queue_mutex));
        return NULL;
      }
      // Note that when the mapping threads finished, they must signal output_queue_pro_cond, so that the output queue thread won't get stuck here forever
      pthread_cond_wait(&(output_queue->pro_cond), &(output_queue->queue_mutex));
    }
    // Take one result from queue
    kv_swap(bam1_t*, sam_alignment_kvec.v, output_queue->sam_alignment_kvec[output_queue->front].v);
    output_queue->front = (output_queue->front + 1) % output_queue->max_size;
    --(output_queue->size);
    pthread_cond_signal(&(output_queue->con_cond));
    pthread_mutex_unlock(&(output_queue->queue_mutex));
    // Save the result into the file
    for (size_t si = 0; si < kv_size(sam_alignment_kvec.v); ++si) {
      int htslib_err = sam_write1(output_queue->output_sam_file, output_queue->sam_header, kv_A(sam_alignment_kvec.v, si));
      assert(htslib_err >= 0);
    }
  }
  for (size_t si = 0; si < kv_size(sam_alignment_kvec.v); ++si) {
    bam_destroy1(kv_A(sam_alignment_kvec.v, si));
  }
  kv_destroy(sam_alignment_kvec.v);
}

void output_sam_header(const char *output_file_path, const SequenceBatch *reference_sequence_batch, OutputQueue *output_queue) {
  //sam_hdr_t sam_header;
  output_queue->sam_header->n_targets = reference_sequence_batch->num_loaded_sequences;
  output_queue->sam_header->target_len = (uint32_t*) malloc(output_queue->sam_header->n_targets * sizeof(uint32_t));
  output_queue->sam_header->target_name = (char**) malloc(output_queue->sam_header->n_targets * sizeof(char*));
  kstring_t header_kstr = {0, 0, NULL};
  // TODO(Haowen): add PG later
  //ksprintf(&header_kstr, "@PG\tID:FEM\tPN:FEM\tVN:%s\tCL:%s", FEM_PACKAGE, argv[0]);
  //for (int i = 1; i < argc; ++i) {
  //  ksprintf(&header_kstr, " %s", argv[i]);
  //}
  for (uint32_t ri = 0; ri < output_queue->sam_header->n_targets; ++ri) {
    output_queue->sam_header->target_len[ri] = get_sequence_length_from_sequence_batch_at(reference_sequence_batch, ri);
    output_queue->sam_header->target_name[ri] = get_sequence_name_from_sequence_batch_at(reference_sequence_batch, ri);
    ksprintf(&header_kstr, "@SQ\tSN:%s\tLN:%d\n", output_queue->sam_header->target_name[ri], output_queue->sam_header->target_len[ri]);
  }
  output_queue->sam_header->l_text = ks_len(&header_kstr);
  output_queue->sam_header->text = ks_str(&header_kstr);
  output_queue->sam_header->sdict = NULL;
  output_queue->sam_header->hrecs = NULL;
  output_queue->sam_header->ref_count = 1;
  int htslib_err = sam_hdr_write(output_queue->output_sam_file, output_queue->sam_header);
  assert(htslib_err == 0);
}

void swap_bam1_t(bam1_t **a, bam1_t **b) {
  bam1_t *tmp = *a;
  *a = *b;
  *b = tmp;
}
