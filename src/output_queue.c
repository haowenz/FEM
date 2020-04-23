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
  output_queue->sam_alignments = (bam1_t**)malloc(output_queue->max_size * sizeof(bam1_t*));
  //output_queue->output_kstrings = (kstring_t*)malloc(output_queue->max_size * sizeof(kstring_t));
  for (size_t i = 0; i < output_queue->max_size; ++i) {
    //output_queue->sam_alignments[i].l_data = 0;
    //output_queue->sam_alignments[i].m_data = 0;
    //output_queue->sam_alignments[i].data = NULL;
    //bam_set_mempolicy(&(output_queue->sam_alignments[i]), BAM_USER_OWNS_STRUCT|BAM_USER_OWNS_DATA);
    output_queue->sam_alignments[i] = bam_init1();
    assert(output_queue->sam_alignments[i]);
    //bam_set_mempolicy(output_queue->sam_alignments[i], BAM_USER_OWNS_STRUCT|BAM_USER_OWNS_DATA);
    //(output_queue->output_kstrings)[i].l = 0;
    //(output_queue->output_kstrings)[i].m = 0;
    //(output_queue->output_kstrings)[i].s = NULL;
  }
  output_queue->output_sam_file = sam_open_format(output_file_path, "w", NULL);
  output_queue->sam_header = sam_hdr_init();
  output_sam_header(output_file_path, reference_sequence_batch, output_queue);
  //output_queue->output_file = fopen(output_file_path, "w");
  //if (output_queue->output_file == NULL) {
  //  fprintf(stderr, "Cannot open output file.\n");
  //  exit(EXIT_FAILURE);
  //}
  fprintf(stderr, "Initialize output queue successfully!\n");
}

void destroy_output_queue(OutputQueue *output_queue) {
  pthread_mutex_destroy(&(output_queue->queue_mutex));
  pthread_cond_destroy(&(output_queue->pro_cond));
  pthread_cond_destroy(&(output_queue->con_cond));
  // free alignments
  if (output_queue != NULL) {
    for (size_t i = 0; i < output_queue->max_size; ++i) {
      //if (output_queue->output_kstrings[i].s != NULL) {
      //  free(output_queue->output_kstrings[i].s);
      //}
      //if (output_queue->sam_alignments[i].data != NULL) {
      //  free(output_queue->sam_alignments[i].data);
      //}
      bam_destroy1(output_queue->sam_alignments[i]);
    }
    free(output_queue->sam_alignments);
    //free(output_queue->output_kstrings);
  }
  // free header
  //free(output_queue->sam_header->target_len);
  //free(output_queue->sam_header->target_name);
  for (uint32_t ri = 0; ri < output_queue->sam_header->n_targets; ++ri) {
    output_queue->sam_header->target_name[ri] = NULL;
  }
  sam_hdr_destroy(output_queue->sam_header);
  sam_close(output_queue->output_sam_file);
  //fclose(output_queue->output_file);
  fprintf(stderr, "Destroy output queue successfully!\n");
}

void push_output_queue(bam1_t **sam_alignment, OutputQueue *output_queue) {
  pthread_mutex_lock(&(output_queue->queue_mutex));
  // Wait if the queue is full
  while (output_queue->size == output_queue->max_size) {
    pthread_cond_wait(&(output_queue->con_cond), &(output_queue->queue_mutex));
  }
  // Push one result into the queue
  //swap_kstring_t(&(output_queue->output_kstrings[output_queue->rear]), result);
  swap_bam1_t(sam_alignment, &(output_queue->sam_alignments[output_queue->rear]));
  output_queue->rear = (output_queue->rear + 1) % output_queue->max_size;
  ++(output_queue->size);
  pthread_cond_signal(&(output_queue->pro_cond));
  pthread_mutex_unlock(&(output_queue->queue_mutex));
}

void *output_queue_thread(void *output_queue_v) {
  OutputQueue *output_queue = (OutputQueue*)output_queue_v;
  while (1) {
    //char *result_string;
    bam1_t *sam_alignment = bam_init1();
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
    //kstring_t output_string = {0, 0, NULL};
    //swap_kstring_t(&output_string, &(output_queue->output_kstrings[output_queue->front]));
    swap_bam1_t(&sam_alignment, &(output_queue->sam_alignments[output_queue->front]));
    output_queue->front = (output_queue->front + 1) % output_queue->max_size;
    --(output_queue->size);
    pthread_cond_signal(&(output_queue->con_cond));
    pthread_mutex_unlock(&(output_queue->queue_mutex));
    // Save the result into the file
    int htslib_err = sam_write1(output_queue->output_sam_file, output_queue->sam_header, sam_alignment);
    assert(htslib_err >= 0);
    //fprintf(stderr, "Output one sam\n");
    bam_destroy1(sam_alignment);
  }
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
