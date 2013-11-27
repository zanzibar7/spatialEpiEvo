#include <unistd.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include "curses.h"
#include "option.h"

extern gsl_rng* random_number_generator;
extern int output_mode;

void Assert(short test, const char* s) {
	if ( !test ) {
		fprintf(stderr, "Error: %s\n", s );
		//exit(EXIT_FAILURE);
		abort();
	}
}


