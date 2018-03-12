#include "FEM-align.h"

static int print_usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   FEM align [options] \n\n");
    fprintf(stderr, "Options:");
    fprintf(stderr, "\n");
    fprintf(stderr, "         -e        INT    error threshold \n");
    fprintf(stderr, "         -t        INT    number of threads \n");
    fprintf(stderr, "         -f        STR    seeding algorithm: \"g\" for group seeding and \"vl\" for variable-length seeding \n");
    fprintf(stderr, "         -a               select one additional q-gram (only for test)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input/output: ");
    fprintf(stderr, "\n");
    fprintf(stderr, "         --ref     STR    Input reference file\n");
    fprintf(stderr, "         --read    STR    Input read file\n");
    fprintf(stderr, "         -o        STR    Output SAM file \n");
    fprintf(stderr, "\n");

    return 1;
}

static int check_args() {
    if(edit_distance<0 || edit_distance>7) {
        printf("Wrong edit distance.\n");
        return 0;
    } else if(cpu_thread_num<=0) {
        printf("Wrong number of threads.\n");
        return 0;
    } else if(additional_gram_num<0||additional_gram_num>2) {
        printf("The number of additional q-gram is too large");
        return 0;
    } else if(!ref_file_name||!strcmp(ref_file_name,"")) {
        printf("Reference file path is required.\n");
        return 0;
    } else if(!read_file_name1||!strcmp(read_file_name1,"")) {
        printf("Read file path is required.\n");
        return 0;
    } else if(!result_file_name||!strcmp(result_file_name,"")) {
        printf("Output file path is required.\n");
        return 0;
    } else {
        printf("e: %d, t:%d, a: %d, ref: %s, read: %s, output: %s\n ", edit_distance, cpu_thread_num, additional_gram_num, ref_file_name, read_file_name1, result_file_name);
        return 1;
    }
}


int align_main(int argc, char *argv[]) {
    const char *short_opt = "ahf:e:t:o:r:i:";
    struct option long_opt[] = 
    {
        {"help", no_argument, NULL, 'h'},
        {"ref",required_argument, NULL, 'r'},
        {"read",required_argument, NULL,'i'}
    };
    int c, option_index;
    additional_gram_num = 0;
    cpu_thread_num = 1;
    while((c = getopt_long(argc, argv, short_opt, long_opt, &option_index))>=0) {
        switch(c) {
            case 'r':
                printf("name: %s, ref: %s\n",long_opt[option_index].name,optarg);
                ref_file_name =optarg;
                break;
            case 'i':
                printf("name: %s, read: %s\n",long_opt[option_index].name, optarg);
                read_file_name1= optarg;
                break;
            case 'e':
                edit_distance = atoi(optarg);
                break;
            case 't':
                cpu_thread_num = atoi(optarg);
                break;
            case 'a':
                additional_gram_num = 1;
                break;
            case 'f':
                if(strcmp(optarg, "vl")==0)
                    novelCandidateGenerator = variableLengthSeedingCandidateGenerator;
                else if(strcmp(optarg,"g")==0)
                    novelCandidateGenerator = groupSeedingCandidateGenerator;
                else {
                    fprintf(stderr, "Wrong name of seeding algorithm!");
                    return print_usage();
                }
                break;
            case 'o':
                result_file_name = optarg;
                break;
            default:
                return print_usage();
        }
    }

    if (!check_args()) {
        return print_usage();
    }

    char tempIndexFileName[256];
    strcpy(tempIndexFileName, ref_file_name);
    sprintf(tempIndexFileName + strlen(tempIndexFileName), ".fem");
    index_file_name =tempIndexFileName;
    char tempHeaderFileName[256];
    strcpy(tempHeaderFileName, tempIndexFileName);
    sprintf(tempHeaderFileName + strlen(tempHeaderFileName), ".header");
    header_file_name = tempHeaderFileName;
    pthread_t CPUTaskHandle[cpu_thread_num];
    pthread_t readQueueHandle;
    pthread_t outputQueueHandle;
    initMapper();
    loadIndex();

    double startTime = realtime();
    fprintf(stderr, "Window_size: %d.\n", window_size);
    finished_thread_num = 0;
    int thread_id_array[cpu_thread_num];
    for (int i = 0; i < cpu_thread_num; ++i) {
        candidate_num[i] = 0;
        candidate_num_without_add_filter[i] = 0;
        mapping_num[i] = 0;
        read_count[i] = 0;
        mapped_read_num[i] = 0;
        thread_id_array[i] = i;
    }

    int err = 0;
    err = pthread_create(&readQueueHandle, NULL, startSingleReadQueueThread, NULL);
    if (err == 0) {
        fprintf(stderr, "Created read queue successfully.\n");
    }
    err = pthread_create(&outputQueueHandle, NULL, startOutputQueueThread, NULL);
    if (err == 0) {
        fprintf(stderr, "Created output queue successfully.\n");
    }

    for (int i = 0; i < cpu_thread_num; ++i) {
        err = pthread_create(CPUTaskHandle + i, NULL, startCPUThread, thread_id_array + i);
    }

    for (int i = 0; i < cpu_thread_num; ++i) {
        err = pthread_join(CPUTaskHandle[i], NULL);
        ++finished_thread_num;
    }

    pthread_cond_signal(&output_queue_pro_cond);
    err = pthread_join(outputQueueHandle, NULL);
    if (err == 0) {
        fprintf(stderr, "Output queue thread joint successfully.\n");
    }


    err = pthread_join(readQueueHandle, NULL);
    if (err == 0) {
        fprintf(stderr, "Read queue thread joint successfully.\n");
    }



    uint32_t totalCandidateNum = 0;
    uint32_t totalMappingNum = 0;
    uint32_t totalReadNum=0;
    uint32_t totalCanddidateNumWithoutFilter = 0;
    uint32_t totalMappedReadNum =0;
    for (int i = 0; i < cpu_thread_num; ++i) {
        totalCandidateNum += candidate_num[i];
        totalMappingNum += mapping_num[i];
        totalCanddidateNumWithoutFilter += candidate_num_without_add_filter[i];
        totalMappedReadNum += mapped_read_num[i];
        totalReadNum += read_count[i];
    }

    fprintf(stderr, "The number of read: %u.\n", totalReadNum);
    fprintf(stderr, "The number of mapped read: %u.\n", totalMappedReadNum);
    fprintf(stderr, "The number of candidate before additional q-gram filter: %u.\n", totalCanddidateNumWithoutFilter);
    fprintf(stderr, "The number of candidate: %u.\n", totalCandidateNum);
    fprintf(stderr, "The number of mapping: %u.\n", totalMappingNum);
    fprintf(stderr, "Time: %fs.\n", realtime() - startTime);

    finalizeMapper();

    return 0;
}
