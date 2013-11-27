#import <unistd.h>
#import <iostream>
#import <gsl/gsl_errno.h>
#import <gsl/gsl_randist.h>
#import "curses.h"

extern gsl_rng* random_number_generator;
extern int output_mode;

void Assert(short test, const char* s) {
	if ( !test ) {
		fprintf(stderr, "Error: %s\n", s );
		//exit(EXIT_FAILURE);
		abort();
	}
}


