#ifndef INDEXTOOLS_H_
#define INDEXTOOLS_H_

#include "utils.h"
#include "mapper.h"

extern char *index_file_path;
extern int max_hash_value;
extern uint32_t hash_table_size;
extern uint32_t *occurrence_table;
extern int *lookup_table;

void initialize_loading_index();
void load_index();
void finalize_loading_index();

void initialize_saving_header();
void finalize_saving_header();

void initialize_saving_index(); 
void construct_index();
void save_index(); 
void finalize_saving_index(); 
#endif /* INDEXTOOLS_H_ */
