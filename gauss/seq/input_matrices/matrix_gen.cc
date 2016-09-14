#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) 
{
    if (argc != 2) {
	fprintf(stderr, "usage: %s <matrixorder>\n", argv[0]);
	exit(-1);
    }

    int nsize = atoi(argv[1]);
    
    fprintf(stdout, "%d %d %d\n", nsize, nsize, nsize*nsize);

    for (int col=1; col<=nsize; col++) {
	for (int row=1; row<=nsize; row++) {
	    double value = (row < col) ? 2*row : 2*col;
	    fprintf(stdout, "%d %d %f\n", row, col, value);
	}
    }

    fprintf(stdout, "0 0 0\n");

    return 0;
}
