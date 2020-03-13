/*
 * utils.c
 *
 *  Created on: Mar 23, 2016
 *      Author: howen
 */

#include "utils.h"

double get_cpu_time() {
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double get_real_time() {
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

int compare_sub_seed(const void *ta, const void *tb) {
    SubSeed my_ta = *(SubSeed *) ta;
    SubSeed my_tb = *(SubSeed *) tb;
    if (my_ta.locationNum > my_tb.locationNum) {
        return 1;
    } else if (my_ta.locationNum == my_tb.locationNum) {
        return 0;
    } else {
        return -1;
    }
}

int compare_gCandidate(const void *ta, const void *tb) {
    gCandidate my_ta = *(gCandidate *) ta;
    gCandidate my_tb = *(gCandidate *) tb;
    if (my_ta.locationsNum > my_tb.locationsNum) {
        return 1;
    } else if (my_ta.locationsNum == my_tb.locationsNum) {
        return 0;
    } else {
        return -1;
    }
}

int compare_hash_location_pair(const void *a, const void *b) {
    HashLocationPair ta = *((HashLocationPair*)a);
    HashLocationPair tb = *((HashLocationPair*)b);
	if (ta.hashValue < tb.hashValue) {
		return -1;
	} else if (ta.hashValue == tb.hashValue) {
		if (ta.location < tb.location)
			return -1;
		else
			return 1;
	} else
		return 1;
}
