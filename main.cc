#include "headers.h"


gsl_rng* random_number_generator;

bool WITH_GRAPHICS = false;

#include "simulation.cc"

int main(int argc, const char** argv) {
	const char* opt_message[] = {
		"Program: epi, February 12, 2010",
		"Started:  July 11, 2009",
		"Purpose: evolution of a spatially structured endemic disease",
		"Usage: epi <options>",
		"Options:",
		"\t-v\tdisplays this message.",
		NULL
	}; Options opt(opt_message); //opt.exit_on_verbose=false;

	int nx = 40;
	opt.AddOption(new OptionInt('x',&nx, (const char*)
		"\t-x#\tnumber of columns in lattice(int)"));

	int ny = 40;
	opt.AddOption(new OptionInt('y',&ny, (const char*)
		"\t-y#\tnumber of rows in lattice(int)"));

	int sweeps = 1000;
	opt.AddOption(new OptionInt('t',&sweeps, (const char*)
		"\t-t#\tnumber of time sweeps(int)"));

	double tauR = 100.;
	opt.AddOption( new OptionDouble('a',&tauR, (const char*)
		"\t-a#\tduration of immunity immunity(double)"));

	double mu = .01; // mutation size
	opt.AddOption( new OptionDouble('m',&mu, (const char*)
		"\t-m#\tmutation width(double)"));

	double first_beta  = 0.01; // initial transmission
	opt.AddOption( new OptionDouble('b',&first_beta, (const char*)
		"\t-b#\tinitial transmission probability(double)"));

	double first_tauI = 30.; // initial recovery
	opt.AddOption( new OptionDouble('g',&first_tauI, (const char*)
		"\t-g#\tinitial infection period(double)"));

	int random_number_seed = (int) time(NULL) + ( (int) getpid() << 14  ) ;
	opt.AddOption( new OptionInt('R',&random_number_seed, (const char*)
		"\t-R#\tset the pseudorandom number seed value(int)"));

	int num_trials = 1;
	opt.AddOption( new OptionInt('N',&num_trials, (const char*)
		"\t-N#\tset the number of repeats(int)"));

	opt.AddOption( new OptionBool('W',&WITH_GRAPHICS, (const char*)
		"\t-W\tturns on curses graphics"));

	opt.CheckOptions(argc, argv);

	gsl_rng_env_setup();
	//gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );
	//gsl_rng * rng = gsl_rng_alloc( gsl_rng_random64_glibc2 );
	random_number_generator = gsl_rng_alloc( gsl_rng_default );
	gsl_rng_set(random_number_generator, random_number_seed );
	printf("Random seed %d\n",random_number_seed);

	for (int i = 0; i < num_trials; ++i ) {
		simulation( nx, ny, sweeps, mu, first_beta, first_tauI );
		printf("\n\n# Trial %d\n",i);
	}
	exit(0);
}

