/*
 * indextools.h
 *
 *  Created on: Mar 24, 2016
 *      Author: howen
 */

#ifndef INDEXTOOLS_H_
#define INDEXTOOLS_H_

#include "utils.h"
#include "mapper.h"

extern char *index_file_name;
extern int max_hash_value;
extern uint32_t hash_table_size;
extern uint32_t *occurrence_table;
extern int *lookup_table;

int initLoadIndex();
int loadIndex();
void finalizeLoadIndex();
void init_saving_header();
void finalize_saving_header();

int init_saving_index(); 
int construct_index();
int save_index(); 
int finalize_saving_index(); 
#endif /* INDEXTOOLS_H_ */
