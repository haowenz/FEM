#include <stdio.h>
#include <string.h>
#include "FEM-align.h"
#include "FEM-index.h"
#include "utils.h"

#define PACKAGE_VERSION "1.1"

static int print_usage() {
    fprintf(stderr, "\n");
	fprintf(stderr, "Program: FEM (Fast and Efficient short read Mapper)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Haowen Zhang <hwzhang@gatech.edu>\n\n");
	fprintf(stderr, "Usage:   FEM <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
	fprintf(stderr, "         align         FEM algorithm\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Note: To use FEM, you need to first index the genome with `FEM index'.\n\n");
    return 1;
}


int main(int argc, char *argv[]) {
    if(argc<2) {
        return print_usage();
    }

    int ret;
    double t_real = realtime();
    if(strcmp(argv[1],"index")==0) {
        ret = index_main(argc-1,argv+1);
    } else if(strcmp(argv[1],"align")==0) {
        ret = align_main(argc-1,argv+1);
    } else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
        return 1;
    }

    if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (int i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}

    return ret;
}
