CC := g++
CINCLUDE:=
CLIB :=  -L/opt/local/lib -lcurses -lm -lgsl -lgslcblas 
CFLAG := -g -I/opt/local/include 
OBJS :=

all: epi 

epi: main.o option.o
	$(CC) $(CFLAG) $(CLIB) $^ -o $@

main.o: main.cc headers.h
	$(CC) $(CFLAG) -c $<

option.o: option.cc 
	$(CC) $(CFLAG) -c $<

ark: 
	tar -czvf epidemic_c_code.tar.gz main.cc headers.h option.h option.cc makefile
	scp epidemic_c_code.tar.gz math:~/ &> /dev/null

clean:
	rm -f -v epi *.o

