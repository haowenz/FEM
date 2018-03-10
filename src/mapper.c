/*
 * mapper.c
 *
 *  Created on: Mar 24, 2016
 *      Author: howen
 */
#include "mapper.h"

int window_size;
int step_size;
int edit_distance;
int additional_gram_num;
int cpu_thread_num;
int finished_thread_num;
int cpu_thread_num_per_group;

uint32_t candidate_num_without_add_filter[256];
uint32_t mapped_read_num[256];
uint32_t candidate_num[256];
uint32_t mapping_num[256];
uint32_t read_count[256];

uint32_t (*novelCandidateGenerator)(const Read *read, const int rc, uint32_t **candidates,uint32_t **swapCandidates,uint32_t *maxCandiNum,TwoTuple **candis, TwoTuple **tempCandis, uint32_t *tupleNum,uint32_t *candiNumWithoutAddFilter);


void initMapper() 
{
    initRef(&reference);
    initReadFile();
    initReadQueue();
    initLoadIndex();
    initOutput();
    initOutputQueue();
    outputHeaderFromIndexFile();
}

void finalizeMapper() 
{
    destroyRef(&reference);
    destroyReadQueue();
    destroyOutputQueue();
    finalizeReadFile();
    finalizeLoadIndex();
    finalizeOutput();
}

int CPUMap(int threadID) 
{
    uint32_t mappingNum = 0;
    uint32_t localCandiNum = 0;
    uint32_t totalCandidateNum = 0;
    uint32_t candiNumWithoutAddFilter = 0;
    uint32_t totalCandidateNumWithoutAddFilter = 0;
    int noLeftResult = 0;
    uint32_t readCount = 0;
    uint32_t mappedReadCount = 0;
    uint32_t maxCandiNum = 4096;
    uint32_t tupleNum = 4096;
    uint32_t resultStrLen = 4096;
    uint32_t *candidates = (uint32_t*) calloc(maxCandiNum, sizeof(uint32_t));
    uint32_t *swapCandidates = (uint32_t*) calloc(maxCandiNum, sizeof(uint32_t));
    TwoTuple *candis = (TwoTuple*) calloc(tupleNum, sizeof(TwoTuple));
    TwoTuple *tempCandis = (TwoTuple*) calloc(tupleNum, sizeof(TwoTuple));
    char *resultStr = (char*) calloc(resultStrLen, sizeof(char));
    Read read;
    
    while (popReadQueue(&read) != -1) 
    {
        readCount++;
        int tmpMappedReadNum = 0;
        localCandiNum = novelCandidateGenerator(&read, 0, &candidates, &swapCandidates, &maxCandiNum, &candis, &tempCandis, &tupleNum, &candiNumWithoutAddFilter);
        totalCandidateNumWithoutAddFilter += candiNumWithoutAddFilter;

        if (localCandiNum > 0) 
        {
            totalCandidateNum += localCandiNum;
            tmpMappedReadNum += CPUVBM(&read, 0, candidates, localCandiNum, &noLeftResult, &resultStrLen, &resultStr);
        }

        localCandiNum = novelCandidateGenerator(&read, 1, &candidates,&swapCandidates, &maxCandiNum,&candis,&tempCandis,&tupleNum ,&candiNumWithoutAddFilter);
        totalCandidateNumWithoutAddFilter += candiNumWithoutAddFilter;

        if (localCandiNum > 0) 
        {
            totalCandidateNum += localCandiNum;
            tmpMappedReadNum += CPUVBM(&read, 1, candidates, localCandiNum, &noLeftResult, &resultStrLen, &resultStr);
        }

        if(tmpMappedReadNum>0)
        {
            mappingNum += tmpMappedReadNum; 
            mappedReadCount++;
        }
    }

    //fprintf(stderr, "Thread %d is pushing local results to output queue...\n", threadID);

    while (noLeftResult != 1) 
    {
        if (pushOutputQueue(resultStr)) 
        {
            noLeftResult = 1;
            //free(resultStr);
            resultStr = NULL;
        }
    }

    //fprintf(stderr, "Thread %d pushed local results to output queue successfully.\n", threadID);

    free(candidates);
    free(swapCandidates);
    free(candis);
    free(tempCandis);
    candidate_num_without_add_filter[threadID] = totalCandidateNumWithoutAddFilter;
    candidate_num[threadID] = totalCandidateNum;
    mapping_num[threadID] = mappingNum;
    read_count[threadID] = readCount;
    mapped_read_num[threadID] = mappedReadCount;
    fprintf(stderr, "Thread %d completed.\n", threadID);

    return mappingNum;
}


void *startCPUThread(void *arg) 
{
    int threadID = *((int*) arg);
    CPUMap(threadID);
    return NULL;
}

void *startSingleReadQueueThread(void *arg) 
{
    singleReadQueueThread();
    return NULL;
}

void *startOutputQueueThread(void *arg) 
{
    outputQueueThread();
    return NULL;
}
