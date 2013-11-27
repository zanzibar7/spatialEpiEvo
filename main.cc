#include "headers.h"

enum state { SUSCE, INFEC, IMMUN };

gsl_rng* random_number_generator;

bool WITH_GRAPHICS = false;

int simulation(
	const int nx,
	const int ny,
	const int sweeps,
	const double alpha,
	const double mu,
	const double first_beta,
	const double first_gamma
	) {

	const int nxny = nx*ny;
	const int t_max = sweeps*nxny;

	int num_INFEC = 0;
	int num_IMMUN = 0;
	int num_SUSCE = nxny;

	state*  __restrict__ pop   = new state[nxny];
	double* __restrict__ beta  = new double[nxny];
	double* __restrict__ gamma = new double[nxny];
	int neighbors[8];

	pop[0] = INFEC;
	beta[0] = first_beta;
	gamma[0] = first_gamma;
	--num_SUSCE;
	++num_INFEC;

	// start windowing
	int refresh_period;
	if (WITH_GRAPHICS) {
		initscr();
		clear();
		refresh_period = 10;
	} else {
		refresh_period = nxny;
	}

	for ( int t = 0; t < t_max; ++t ) {
		int i = gsl_rng_uniform_int(random_number_generator,nxny);
		switch (pop[i]) {
			case SUSCE:
				break;
			case IMMUN:
				if ( gsl_ran_bernoulli(random_number_generator, alpha) ) {
					pop[i] = SUSCE;
					--num_IMMUN;
					++num_SUSCE;
				}
				break;
			case INFEC:
				neighbors[0] = i+nx+1;
				neighbors[1] = i+nx;
				neighbors[2] = i+nx-1;
				neighbors[3] = i+1;
				neighbors[4] = i-1;
				neighbors[5] = i-nx+1;
				neighbors[6] = i-nx;
				neighbors[7] = i-nx-1;
				for (int j = 0; j < 8; ++j ) {
					int h = neighbors[j];
					if ( h < 0 or h >= nxny ) { h = (h + nxny ) % nxny ; }
					if ( SUSCE == pop[h] ) {
						if (gsl_ran_bernoulli(random_number_generator,beta[i])){
							pop[h] = INFEC;

							beta[h] = beta[i]*(1+gsl_ran_gaussian(random_number_generator, mu));
							if (beta[h] < 0.) { beta[h] = 0.; }
							else if (beta[h] > 1.) { beta[h] = 1.; }

							gamma[h] = gamma[i]*(1+gsl_ran_gaussian(random_number_generator, mu));
							if (gamma[h] < 0.) { gamma[h] = 0.; }
							else if (gamma[h] > 1.) { gamma[h] = 1.; }

							--num_SUSCE;
							++num_INFEC;
						}
					}
				}
				if ( gsl_ran_bernoulli(random_number_generator, gamma[i]) ) {
					pop[i] = IMMUN;
					--num_INFEC;
					++num_IMMUN;
				}
				break;
			default:
				Assert(false, "invalid state");
	
		}
		if ( 0 == num_INFEC ) {
			//printf( "%d\t%d\t%d\t%d\n", t, num_SUSCE, num_INFEC, num_IMMUN );
			break;
		}

		if ( 0 == t%refresh_period ) {
			double sum_gamma = 0.;
			double sum_beta = 0.;
			for (int i = 0; i < nxny; ++i ) {
				switch (pop[i]) {
					case INFEC:
						sum_gamma += gamma[i];
						sum_beta += beta[i];
				}
			}
			if (WITH_GRAPHICS) {
				for (int i = 0; i < nxny; ++i ) {
					move( i/nx, i%nx );
					switch (pop[i]) {
						case IMMUN:
							addch(' ');
							break;
						case INFEC:
							addch('*');
							break;
						default:
							break;
							addch(' ');
					}
				}
				mvprintw(ny+1,0," beta: %.4f",sum_beta/double(num_INFEC));
				mvprintw(ny+2,0,"gamma: %.4f",sum_gamma/double(num_INFEC));
				refresh();
			} else {
				printf("%d\t%d\t%d\t%d\t%.4f\t%.4f\n",
					t/nxny,
					num_SUSCE,
					num_INFEC,
					num_IMMUN,
					sum_beta/double(num_INFEC),
					sum_gamma/double(num_INFEC));
			}
		}

		
	}
	if (WITH_GRAPHICS) { endwin(); }
}

int main(int argc, const char** argv) {
	const char* opt_message[] = {
		"Program: epi, July 11, 2009",
		"Started:  July 11, 2009",
		"Purpose: evolution of a spatially structured endemic disease",
		"Usage: epi <options>",
		"Options:",
		"\t-v\tdisplays this message.",
		NULL
	}; Options opt(opt_message); //opt.exit_on_verbose=false;

	int nx = 80;
	opt.AddOption(new OptionInt('x',&nx, (const char*)
		"\t-x#\tnumber of columns in lattice(int)"));

	int ny = 40;
	opt.AddOption(new OptionInt('y',&ny, (const char*)
		"\t-y#\tnumber of rows in lattice(int)"));

	int sweeps = 1000;
	opt.AddOption(new OptionInt('t',&sweeps, (const char*)
		"\t-t#\tnumber of time sweeps(int)"));

	double alpha = .02;
	opt.AddOption( new OptionDouble('a',&alpha, (const char*)
		"\t-a#\tprobability of lossing immunity(double)"));

	double mu = .02; // mutation size
	opt.AddOption( new OptionDouble('m',&mu, (const char*)
		"\t-m#\tmutation width(double)"));

	double first_beta  = 0.5; // initial transmission
	opt.AddOption( new OptionDouble('b',&first_beta, (const char*)
		"\t-b#\tinitial transmission probability(double)"));

	double first_gamma = 0.3; // initial recovery
	opt.AddOption( new OptionDouble('g',&first_gamma, (const char*)
		"\t-g#\tinitial recovery probability(double)"));

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

	for (int i = 0; i < num_trials; ++i ) {
		simulation( nx, ny, sweeps, alpha, mu, first_beta, first_gamma );
	}
	exit(0);
}

