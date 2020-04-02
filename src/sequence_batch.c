#include "sequence_batch.h"

#include "utils.h"

void initialize_sequence_batch(SequenceBatch *sequence_batch) { 
  kv_init(sequence_batch->sequences);
}

void initialize_sequence_batch_with_max_size(uint32_t max_num_sequences, SequenceBatch *sequence_batch) { 
  sequence_batch->max_num_sequences = max_num_sequences;
  kv_init(sequence_batch->sequences);
  kv_resize(kseq_t*, sequence_batch->sequences, max_num_sequences);
  for (uint32_t i = 0; i < max_num_sequences; ++i) {
    kv_push(kseq_t*, sequence_batch->sequences, (kseq_t*)calloc(1, sizeof(kseq_t)));
  }
}

void destory_sequence_batch(SequenceBatch *sequence_batch) {
  kv_destroy(sequence_batch->sequences);
}

void initialize_sequence_batch_loading(const char *sequence_file_path, SequenceBatch *sequence_batch) {
  sequence_batch->sequence_file = gzopen(sequence_file_path, "r");
  if (sequence_batch->sequence_file == NULL) {
    fprintf(stderr, "Cannot find sequence file!");
    exit(EXIT_FAILURE);
  }
  sequence_batch->sequence_kseq = kseq_init(sequence_batch->sequence_file);
}

void finalize_sequence_batch_loading(SequenceBatch *sequence_batch) {
  kseq_destroy(sequence_batch->sequence_kseq);
  gzclose(sequence_batch->sequence_file);
}

void load_batch_of_sequences_into_sequence_batch(SequenceBatch *sequence_batch) {
  double real_start_time = get_real_time();
  sequence_batch->num_loaded_sequences = 0;
  sequence_batch->num_bases = 0;
  for (uint32_t sequence_index = 0; sequence_index < sequence_batch->max_num_sequences; ++sequence_index) { 
    int length = kseq_read(sequence_batch->sequence_kseq);
    while (length == 0) { // Skip the sequences of length 0
      length = kseq_read(sequence_batch->sequence_kseq);
    }
    if (length > 0) {
      kseq_t *sequence = kv_A(sequence_batch->sequences, sequence_index);
      swap_kstring_t(&(sequence_batch->sequence_kseq->seq), &(sequence->seq));
      swap_kstring_t(&(sequence_batch->sequence_kseq->name), &(sequence->name));
      swap_kstring_t(&(sequence_batch->sequence_kseq->comment), &(sequence->comment));
      if (sequence_batch->sequence_kseq->qual.l != 0) { // fastq file
        swap_kstring_t(&(sequence_batch->sequence_kseq->qual), &(sequence->qual));
      }
      //sequence->id = sequence_batch->num_loaded_sequences_;
      ++(sequence_batch->num_loaded_sequences);
      sequence_batch->num_bases += length;
    } else {
      if (length != -1) {
        fprintf(stderr, "Didn't reach the end of sequence file, which might be corrupted!");
        exit(EXIT_FAILURE);
      }
      // make sure to reach the end of file rather than meet an error
      break;
    }
  }
  if (sequence_batch->num_loaded_sequences != 0) {
    fprintf(stderr, "Number of sequences: %d\n", sequence_batch->num_loaded_sequences);
    fprintf(stderr, "Number of bases: %ld\n", sequence_batch->num_bases);
    fprintf(stderr, "Loaded sequence batch successfully in %fs\n", get_real_time() - real_start_time);
  } else {
    fprintf(stderr, "No more sequences.\n");
  }
}

void load_all_sequences_into_sequence_batch(SequenceBatch *sequence_batch) {
  double real_start_time = get_real_time();
  sequence_batch->num_loaded_sequences = 0;
  sequence_batch->num_bases = 0;
  int length = kseq_read(sequence_batch->sequence_kseq);
  while (length >= 0) { 
    if (length == 0) { // Skip the sequences of length 0
      continue;
    } else if (length > 0) {
      kseq_t *sequence = (kseq_t*)calloc(1, sizeof(kseq_t));
      kv_push(kseq_t*, sequence_batch->sequences, sequence);
      swap_kstring_t(&(sequence_batch->sequence_kseq->seq), &(sequence->seq));
      swap_kstring_t(&(sequence_batch->sequence_kseq->name), &(sequence->name));
      swap_kstring_t(&(sequence_batch->sequence_kseq->comment), &(sequence->comment));
      if (sequence_batch->sequence_kseq->qual.l != 0) { // fastq file
        swap_kstring_t(&(sequence_batch->sequence_kseq->qual), &(sequence->qual));
      }
      //sequence->id = sequence_batch->num_loaded_sequences;
      ++(sequence_batch->num_loaded_sequences);
      sequence_batch->num_bases += length;
    } else {
      if (length != -1) {
        fprintf(stderr, "Didn't reach the end of sequence file, which might be corrupted!");
        exit(EXIT_FAILURE);
      }
      // make sure to reach the end of file rather than meet an error
      break;
    }
    length = kseq_read(sequence_batch->sequence_kseq);
  }
  fprintf(stderr, "Number of sequences: %d\n", sequence_batch->num_loaded_sequences);
  fprintf(stderr, "Number of bases: %ld\n", sequence_batch->num_bases);
  fprintf(stderr, "Loaded all sequences successfully in %fs\n", get_real_time() - real_start_time);
}

//bool LoadOneSequenceAndSaveAt(uint32_t sequence_index) {
//  //double real_start_time = Chromap::GetRealTime();
//  bool no_more_sequence = false;
//  int length = kseq_read(sequence_kseq_);
//  while (length == 0) { // Skip the sequences of length 0
//    length = kseq_read(sequence_kseq_);
//  }
//  if (length > 0) {
//    kseq_t *sequence = sequence_batch_[sequence_index];
//    std::swap(sequence_kseq_->seq, sequence->seq);
//    std::swap(sequence_kseq_->name, sequence->name);
//    std::swap(sequence_kseq_->comment, sequence->comment);
//    sequence->id = num_loaded_sequences_;
//    ++num_loaded_sequences_;
//    if (sequence_kseq_->qual.l != 0) { // fastq file
//      std::swap(sequence_kseq_->qual, sequence->qual);
//    } 
//  } else {
//    if (length != -1) {
//      Chromap<>::ExitWithMessage("Didn't reach the end of sequence file, which might be corrupted!");
//    }
//    // make sure to reach the end of file rather than meet an error
//    no_more_sequence = true;
//  }
//  return no_more_sequence;
//}
