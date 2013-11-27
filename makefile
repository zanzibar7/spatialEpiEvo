CC := g++-4.3
CINCLUDE:=
#HD:=$(HOME)/include/4.3
#CLIB := -lcurses -lm -lgsl -lgslcblas $(HD)/option.o $(HD)/mymath.o 
#CFLAG := -g -I$(HD)
CLIB := -lcurses -lm -lgsl -lgslcblas 
CFLAG := -O4 
OBJS :=

all: epi 

epi: main.o
	$(CC) $(CFLAG) $(CLIB) main.cc -o $@

main.o: main.cc headers.h
	$(CC) $(CFLAG) -c $<

ark: 
	tar -czvf epidemic_c_code.tar.gz main.cc headers.h makefile
	scp epidemic_c_code.tar.gz math:~/ &> /dev/null

clean:
	rm -f -v *.o

