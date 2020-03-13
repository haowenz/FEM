#ifndef INDEXTOOLS_H_
#define INDEXTOOLS_H_

#include "utils.h"
#include "mapper.h"

extern char *index_file_name;
extern int max_hash_value;
extern uint32_t hash_table_size;
extern uint32_t *occurrence_table;
extern int *lookup_table;

int initialize_loading_index();
int load_index();
void finalize_loading_index();
void initialize_saving_header();
void finalize_saving_header();

int initialize_saving_index(); 
int construct_index();
int save_index(); 
int finalize_saving_index(); 
#endif /* INDEXTOOLS_H_ */
