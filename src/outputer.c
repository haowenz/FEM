/*
 * outputer.c
 *
 *  Created on: Mar 25, 2016
 *      Author: howen
 */

#include "outputer.h"

FILE *result_file_fp;
FILE *header_file_fp;
char *result_file_name;
char *header_file_name;
pthread_mutex_t outputMutex;
pthread_mutex_t output_queue_mutex;
pthread_cond_t output_queue_pro_cond;
pthread_cond_t output_queue_cos_cond;
pthread_cond_t not_empty_cond;

int output_queue_size;
int output_queue_size_max = 100000;

int output_queue_front;
int output_queue_rear;
/*output queue*/
char** output_queue;

void initOutputQueue() {
	pthread_mutex_init(&output_queue_mutex, NULL);
	pthread_cond_init(&output_queue_pro_cond, NULL);
	output_queue_size = 0;
	output_queue = (char**) malloc(output_queue_size_max * sizeof(char*));
}

void destroyOutputQueue() {
	pthread_mutex_destroy(&output_queue_mutex);
	pthread_cond_destroy(&output_queue_pro_cond);
	free(output_queue);
}

int pushOutputQueue(char* result) {
	int flag = pthread_mutex_trylock(&output_queue_mutex);
	if (flag != 0) {
		return 0;
	}

	if (output_queue_size >= output_queue_size_max) {
		pthread_mutex_unlock(&output_queue_mutex);
		pthread_cond_signal(&output_queue_pro_cond);
		return 0;
	}

	output_queue[output_queue_rear] = result;
	output_queue_rear = (output_queue_rear + 1) % output_queue_size_max;
	++output_queue_size;
	pthread_mutex_unlock(&output_queue_mutex);
	pthread_cond_signal(&output_queue_pro_cond);

	return 1;
}

inline int popOutputQueue(char* result) {

	return 1;
}

void outputQueueThread() {
	while (finished_thread_num < cpu_thread_num) {
		char *resultStr;
		pthread_mutex_lock(&output_queue_mutex);
		while (output_queue_size == 0) {
			pthread_cond_wait(&output_queue_pro_cond, &output_queue_mutex);
			if (finished_thread_num == cpu_thread_num) {
                pthread_mutex_unlock(&output_queue_mutex);
				break;
			}
		}

		resultStr = output_queue[output_queue_front];
		output_queue[output_queue_front] = NULL;
		output_queue_front = (output_queue_front + 1) % output_queue_size_max;
		--output_queue_size;
		pthread_mutex_unlock(&output_queue_mutex);

		if (resultStr != NULL) {
			fwrite(resultStr, sizeof(char), strlen(resultStr), result_file_fp);
			free(resultStr);
            resultStr = NULL;
		}
	}

	while (output_queue_size > 0) {
		fwrite(output_queue[output_queue_front], sizeof(char), strlen(output_queue[output_queue_front]), result_file_fp);
		free(output_queue[output_queue_front]);
		output_queue[output_queue_front] = NULL;
		output_queue_front = (output_queue_front + 1) % output_queue_size_max;
		--output_queue_size;
		//printf("queue_size=%d.\n", output_queue_size);
	}
}

void clearOutputQueue() {

}

void initOutput() {
	outputMutex = PTHREAD_MUTEX_INITIALIZER;
	result_file_fp = fopen(result_file_name, "w");
	if (result_file_fp == NULL) {
		fprintf(stderr, "Cannot open output file.\n");
        exit(-1);
	}
}

void outputHeaderFromIndexFile() {
	int bufferSize = 4096;
	char buffer[bufferSize];
	header_file_fp = fopen(header_file_name, "r");
	int num = 0;
	do {
		num = fread(buffer, sizeof(char), bufferSize, header_file_fp);
		fwrite(buffer, sizeof(char), num, result_file_fp);
	} while (num > 0);
	fclose(header_file_fp);
}

int outputString(const char* result) {
	int out = pthread_mutex_trylock(&outputMutex);
	if (out == 0) {
		fwrite(result, sizeof(char), strlen(result), result_file_fp);
		pthread_mutex_unlock(&outputMutex);
		return 1;
	} else {
		return 0;
	}
}

void finalizeOutput() {
	fclose(result_file_fp);
}
