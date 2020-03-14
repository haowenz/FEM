#include "readtools.h"
#include "mapper.h"

KSEQ_INIT(int, read)
__UTILS_CHAR_UINT8

FILE *read1_file;
FILE *read2_file;
char *read1_file_path;
char *read2_file_path;
kseq_t *read1_kseq;
kseq_t *read2_kseq;
pthread_mutex_t read1_mutex;
pthread_mutex_t read2_mutex;
pthread_mutex_t read_queue_mutex;
pthread_cond_t read_queue_pro_cond;
pthread_cond_t read_queue_cos_cond;

Read *read_queue;
int read_queue_size_max = 500000;
int read_queue_front;
int read_queue_rear;
int read_queue_size;
int no_more_read;

void initialize_read_queue() {
  reset_read_queue();
  read_queue = (Read*) malloc(read_queue_size_max * sizeof(Read));
  assert(read_queue);
  fprintf(stderr, "Initialize read queue successfully!\n");
}

void reset_read_queue() {
  pthread_mutex_init(&read_queue_mutex, NULL);
  pthread_cond_init(&read_queue_pro_cond, NULL);
  pthread_cond_init(&read_queue_cos_cond, NULL);
  no_more_read = 0;
  read_queue_front = 0;
  read_queue_rear = 0;
  read_queue_size = 0;
}

void destroy_read_queue() {
  pthread_mutex_destroy(&read_queue_mutex);
  pthread_cond_destroy(&read_queue_pro_cond);
  pthread_cond_destroy(&read_queue_cos_cond);
  free(read_queue);
  fprintf(stderr, "Destroy read queue successfully!\n");
}

inline void push_read_queue(Read *read) {
  memcpy(read_queue + read_queue_rear, read, sizeof(Read));
  read_queue_rear = (read_queue_rear + 1) % read_queue_size_max;
  ++read_queue_size;
}

int pop_read_queue(Read *read) {
  pthread_mutex_lock(&read_queue_mutex);
  if (read_queue_size == 0 && no_more_read == 1) {
    //fprintf(stderr, "Quit.\n");
    pthread_mutex_unlock(&read_queue_mutex);
    pthread_cond_signal(&read_queue_pro_cond);
    return -1;
  }
  while (read_queue_size == 0) {
    pthread_cond_signal(&read_queue_cos_cond);
    pthread_cond_wait(&read_queue_pro_cond, &read_queue_mutex);
    if (read_queue_size == 0 && no_more_read == 1) {
      //fprintf(stderr, "Waked up and quit.\n");
      pthread_mutex_unlock(&read_queue_mutex);
      pthread_cond_signal(&read_queue_pro_cond);
      return -1;
    }
  }
  memcpy(read, read_queue + read_queue_front, sizeof(Read));
  read_queue_front = (read_queue_front + 1) % read_queue_size_max;
  --read_queue_size;
  if (no_more_read == 1) {
    //fprintf(stderr, "Number of left reads: %d.\n", read_queue_size);
    if (read_queue_size > 0) {
      pthread_cond_signal(&read_queue_pro_cond);
    }
  } else {
    pthread_cond_signal(&read_queue_cos_cond);
  }
  pthread_mutex_unlock(&read_queue_mutex);
  return 1;
}

void single_read_queue_thread() {
  Read read;
  while (1) {
    int flag = get_single_read(&read);
    if (flag == -1) {
      break;
    }
    pthread_mutex_lock(&read_queue_mutex);
    while (read_queue_size == read_queue_size_max) {
      pthread_cond_signal(&read_queue_pro_cond);
      pthread_cond_wait(&read_queue_cos_cond, &read_queue_mutex);
    }
    push_read_queue(&read);
    pthread_mutex_unlock(&read_queue_mutex);
    pthread_cond_signal(&read_queue_pro_cond);
  }
  pthread_mutex_lock(&read_queue_mutex);
  no_more_read = 1;
  //fprintf(stderr, "No more new read, %d reads left!\n", read_queue_size);
  pthread_mutex_unlock(&read_queue_mutex);
  pthread_cond_signal(&read_queue_pro_cond);
  //pthread_cond_broadcast(&read_queue_pro_cond);
  //fprintf(stderr, "Broadcast!\n");
}

void initialize_ReadBlock(ReadBlock *read_block) {
  read_block->num = 0;
  read_block->length = -1;
  read_block->bases = NULL;
  read_block->rc_bases = NULL;
  read_block->names = (char*) malloc(sizeof(char) * READ_NAME_LEN_MAX * READBLOCK_NUM_MAX);
  assert(read_block->names);
}

void reset_read_file() {
  /*single-end mode*/
  rewind(read1_file);
  read1_kseq = kseq_init(fileno(read1_file));
}

void initialize_read_file() {
  /*single-end mode*/
  //readMutex1= PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_init(&read1_mutex, NULL);

  read1_file = fopen(read1_file_path, "r");
  if (read1_file == NULL) {
    fprintf(stderr, "Cannnot open read file: %s.\n", read1_file_path);
    exit(EXIT_FAILURE);
  } else {
    read1_kseq = kseq_init(fileno(read1_file));
  }
  /*pair-end mode*/
  //eadMutex2= PTHREAD_MUTEX_INITIALIZER;
}

void finalize_read_file() {
  /*single-end mode*/
  kseq_destroy(read1_kseq);
  fclose(read1_file);
  /*pair-end mode*/
}

int get_single_read(Read *read) {
  int l = kseq_read(read1_kseq);
  while (l == 0) {
    l = kseq_read(read1_kseq);
    printf("!");
  }
  if (l > 0) {
    read->length = l;
    memset(read->name, '\0', sizeof(char) * READ_NAME_LEN_MAX);
    memset(read->charBases, '\0', sizeof(char) * READ_LEN_MAX);
    strcpy(read->name, read1_kseq->name.s);
    strcpy(read->charBases, read1_kseq->seq.s);
    int i;
    for (i = 0; i < l; ++i) {
      char c = read1_kseq->seq.s[i];
      uint8_t base = charToUint8(c);
      read->bases[i] = base;
      if (base != 4) {
        read->rc_bases[l - 1 - i] = 3 - base;
      } else {
        read->rc_bases[l - 1 - i] = 4;
      }

    }
    return 1;
  } else if (l == -1) {
    /*end of file*/
    return -1;
  } else {
    /* truncated quality string*/
    fprintf(stderr, "truncated quality string.");
    exit(-1);
  }
  /*never come to here*/
  return 0;
}

int get_ReadBlock(ReadBlock *read_block) {
  int flag = pthread_mutex_trylock(&read1_mutex);
  if (flag != 0) {
    return 0;
  }
  memset(read_block->names, '\0', sizeof(char) * READBLOCK_NUM_MAX * READ_NAME_LEN_MAX);
  int count = 0;
  while (count < READBLOCK_NUM_MAX) {
    int l = kseq_read(read1_kseq);
    while (l == 0) {
      l = kseq_read(read1_kseq);
    }
    if (l > 0) {
      if (count == 0) {
        if (l != read_block->length) {
          read_block->length = l;
          if (read_block->bases != NULL) {
            _mm_free(read_block->bases);
          }
          read_block->bases = (uint8_t*) _mm_malloc(sizeof(uint8_t) * l * READBLOCK_NUM_MAX, 16);
          printf("!!!!!!!!!!\n");
          assert(read_block->bases);
          if (read_block->rc_bases != NULL) {
            _mm_free(read_block->rc_bases);
          }
          read_block->rc_bases = (uint8_t*) _mm_malloc(sizeof(uint8_t) * l * READBLOCK_NUM_MAX, 16);
          assert(read_block->rc_bases);
        }
      } else {
        /*XAM will support variant length later*/
        if (read_block->length != l) {
          fprintf(stderr,"The length of reads is not equal.\n");
          exit(-1);
        }
      }
      memcpy(read_block->names + count * READ_NAME_LEN_MAX, read1_kseq->name.s, sizeof(uint8_t) * READ_NAME_LEN_MAX);
      //			strcpy(read_block->names + count * READ_NAME_LEN_MAX,
      //					read_seq1->name.s);
      int i;
      for (i = 0; i < l; ++i) {
        char c = read1_kseq->seq.s[i];
        uint8_t base = charToUint8(c);
        read_block->bases[count * l + i] = base;
        if (base != 4) {
          read_block->rc_bases[count * l + l - 1 - i] = 3 - base;
        } else {
          read_block->rc_bases[count * l + l - 1 - i] = 4;
        }
      }
      count++;
    } else if (l == -1) {
      /*end of file*/
      read_block->num = count;
      pthread_mutex_unlock(&read1_mutex);
      return -1;
    } else {
      pthread_mutex_unlock(&read1_mutex);
      /* truncated quality string*/
      fprintf(stderr, "truncated quality string.");
      exit(EXIT_FAILURE);
    }
  }
  read_block->num = count;
  pthread_mutex_unlock(&read1_mutex);
  return 0;
}
