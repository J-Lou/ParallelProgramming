# makefile for digits recognizer

.SUFFIXES: .o .c

CC = gcc
CFLAGS = -g -O3
LFLAGS = -lm 

CTAGS =     ctags
BIN = sor
CFILES =   sor.c
OBJECTS =  sor.o
HFILES = 
OTHERSOURCES =  
SOURCES =   $(HFILES) $(CFILES) $(OTHERSOURCES)

.c.o:  
	$(CC) -c $(CFLAGS) $*.c

$(BIN): $(OBJECTS)
	$(CC)  -o $(BIN) $(OBJECTS) $(LFLAGS)

tags: $(CFILES) $(HFILES) $(EXTRACTCFILES) $(RECOGCFILES)
# for vi:
	$(CTAGS) -d -t -T -w $(CFILES) $(HFILES) 
	sort -o tags tags

clean:
	-rm $(BIN) $(OBJECTS) 
