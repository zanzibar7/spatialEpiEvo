CPP:=clang++
CINCLUDE:=
CLIB:= -L /lib/x86_64-linux-gnu -lncurses -lgsl -lgslcblas
CFLAG := -O2 -Os -Wall -Wextra -pedantic -I/opt/local/include 
OBJS :=

all: spatialEpiEvo 

spatialEpiEvo: main.o option.o
	$(CPP) $^ $(CFLAG) $(CLIB) -o $@

main.o: main.cc headers.h simulation.cc
	$(CPP) $(CFLAG) -c $<

option.o: option.cc 
	$(CPP) $(CFLAG) -c $<

clean:
	rm -f -v spatialEpiEvo *.o
