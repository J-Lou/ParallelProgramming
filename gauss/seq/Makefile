all: gauss

CC = gcc
CCFLAGS = -g -Wall -O3
LDFLAGS = -g -Wall -O3

gauss: gauss.o
	$(CC) $(LDFLAGS) -o gauss gauss.o

gauss.o: gauss.c
	$(CC) $(CCFLAGS) -c gauss.c -o gauss.o

clean:
	-rm *.o *~ gauss
