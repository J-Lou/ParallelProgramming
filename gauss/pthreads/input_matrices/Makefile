all: matrix_gen

CC = g++
CCFLAGS = -g -Wall -O3
LDFLAGS = -g -Wall -O3

matrix_gen: matrix_gen.o
	$(CC) $(LDFLAGS) -o matrix_gen matrix_gen.o

matrix_gen.o: matrix_gen.cc
	$(CC) $(CCFLAGS) -c matrix_gen.cc -o matrix_gen.o

clean:
	-rm *.o *~ matrix_gen
