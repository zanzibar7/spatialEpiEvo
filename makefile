#CC:=g++
CC:=  clang
CPP:=  clang++
CINCLUDE:=
CLIB:= -L /lib/x86_64-linux-gnu -lncurses -lgsl -lgslcblas
CFLAG := -O2 -Os -Wall -Wextra -pedantic -I/opt/local/include 
OBJS :=

all: epi 

epi: main.o option.o
	$(CPP) $^ $(CFLAG) $(CLIB) -o $@

#epi -x120 -y120 -b.01 -g30 -t100000 -R18134238 >! data1

main.o: main.cc headers.h simulation.cc
	$(CPP) $(CFLAG) -c $<

option.o: option.cc 
	$(CPP) $(CFLAG) -c $<

ark: 
	tar -czvf epidemic_c_code.tar.gz main.cc headers.h option.h option.cc makefile
	scp epidemic_c_code.tar.gz math:~/ &> /dev/null

clean:
	rm -f -v epi *.o

# g++ main.cc -L /lib/x86_64-linux-gnu -lncurses -lgsl -lgslcblas option.o 
