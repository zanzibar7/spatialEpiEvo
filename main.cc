#include "headers.h"

enum state { SUSCE, INFEC, IMMUN };

gsl_rng* random_number_generator;

bool WITH_GRAPHICS = false;

int simulation(
	const int nx,
	const int ny,
	const int t_max,
	const double mu,
	const double first_beta,
	double first_tauI
	) {
	// start windowing
	int refresh_period;
	if (WITH_GRAPHICS) {
		initscr();
		start_color(); 
		clear();
		refresh_period = 1;
   		init_pair (1, COLOR_WHITE, COLOR_BLACK);
   		init_pair (2, COLOR_RED, COLOR_BLACK);
   		init_pair (3, COLOR_BLUE, COLOR_BLACK);

	} else {
		refresh_period = 100;
	}

	const int tauR = 100;

	const int nxny = nx*ny;

	double* __restrict__ beta  = new double[nxny];
	double* __restrict__ tauI  = new double[nxny];
	int* __restrict__ tic_INFEC = new int[nxny];
	int* __restrict__ tic_IMMUN = new int[nxny];

	state**  popf = new state*[2];
	popf[0]   = new state[nxny];
	popf[1]   = new state[nxny];
	state* pop    = popf[0];
	state* newpop = popf[1];

	int i;
	for ( i = 0; i < nxny; ++i ) {
		popf[0][i] = SUSCE;
		popf[1][i] = SUSCE;
		tic_INFEC[i] = 0;
		tic_IMMUN[i] = 0;
	}
	int neighbors[8];

	int num_INFEC = 0;
	int num_IMMUN = 0;
	int num_SUSCE = nxny;
	double sum_tauI = 0.;
	double sum_beta = 0.;

	i = nxny/2+nx/2; pop[i] = INFEC; beta[i] = first_beta; tauI[i] = first_tauI;
	--num_SUSCE; ++num_INFEC;

	i = 2; pop[i] = IMMUN; tic_IMMUN[i] = tauR/2;
	--num_SUSCE; ++num_IMMUN;
	i = 3; pop[i] = IMMUN; tic_IMMUN[i] = tauR/2;
	--num_SUSCE; ++num_IMMUN;

	for ( int t = 0; t < t_max; ++t ) {
		pop = popf[t%2];
		newpop = popf[(t+1)%2];
		Assert( pop != newpop, "flip collision");
		
		for ( int i = 0; i < nxny; ++i) {
			switch (pop[i]) {
				case IMMUN:
					tic_IMMUN[i]--;
					if ( 0 >= tic_IMMUN[i] ) {
						newpop[i] = SUSCE;
						--num_IMMUN;
						++num_SUSCE;
					} else { newpop[i] = IMMUN; }
				//printf("%d\t%d\t%d\t%d\t%.4f\t%.4f\n", t, num_SUSCE, num_INFEC, num_IMMUN, sum_beta/double(num_INFEC), sum_tauI/double(num_INFEC));
					break;
				case INFEC:
					tic_INFEC[i]--;
					if ( 0 >= tic_INFEC[i] ) {
						newpop[i] = IMMUN;
						tic_IMMUN[i] = tauR;
						--num_INFEC;
						++num_IMMUN;
						//printf("%dh\t%d\t%d\t%d\t%.4f\t%.4f\n", t, num_SUSCE, num_INFEC, num_IMMUN, sum_beta/double(num_INFEC), sum_tauI/double(num_INFEC));
					} else { newpop[i] = INFEC; }
					break;
				case SUSCE:
					newpop[i] = SUSCE;
					{
						neighbors[0] = i+nx+1;
						neighbors[1] = i+nx;
						neighbors[2] = i+nx-1;
						neighbors[3] = i+1;
						neighbors[4] = i-1;
						neighbors[5] = i-nx+1;
						neighbors[6] = i-nx;
						neighbors[7] = i-nx-1;
					}
					{
						double tic = 1.0; // one time step
						for (int j = 0; j < 8; ++j ) { // BUG: not randomized
							int h = neighbors[j];
							if ( h < 0 or h >= nxny ) { h = (h + nxny ) % nxny ; }

							if (INFEC!=pop[h]) continue;

							double toc = gsl_ran_exponential(random_number_generator,1./beta[h]);
							if (toc > tic) {
								continue;
							} else {
								tic = toc;
								newpop[i] = INFEC;
								beta[i] = beta[h]*(1+gsl_ran_gaussian(random_number_generator, mu));
								if (beta[i] < 0.) { beta[i] = 0.; }
								if (beta[i] > 4.) { beta[i] = 4.; }

								tauI[i] = tauI[h]*(1+gsl_ran_gaussian(random_number_generator, mu));
								if (tauI[i] < 0.) { tauI[i] = 0.; }
								if (tauI[i] > 400) { tauI[i] = 400; }
							}
						}
						if ( INFEC == newpop[i] ) {
							--num_SUSCE;
							++num_INFEC;
							tic_INFEC[i] = int(tauI[i]);
							break;
						}
					}
				//printf("%d\t%d\t%d\t%d\t%.4f\t%.4f\n", t, num_SUSCE, num_INFEC, num_IMMUN, sum_beta/double(num_INFEC), sum_tauI/double(num_INFEC));
					break;
				default:
					Assert(false, "invalid state");
			}
		}

		/////////////////////////////////////////////////////

		if ( 0 == num_INFEC ) {
			//printf( "%d\t%d\t%d\t%d\n", t, num_SUSCE, num_INFEC, num_IMMUN );
			break;
		}

		if ( 0 == t%refresh_period ) {
			sum_tauI = 0.;
			sum_beta = 0.;
			num_INFEC = 0;
			num_IMMUN = 0;
			num_SUSCE = 0;
			for ( int i = 0; i < nxny; ++i) {
				switch (newpop[i]) {
					case IMMUN: ++num_IMMUN; break;
					case SUSCE: ++num_SUSCE; break;
					case INFEC: ++num_INFEC;
						sum_tauI += tauI[i];
						sum_beta += beta[i];
						break;
					default: Assert(false,"982");
				}
			}
			if (WITH_GRAPHICS) {
				for (int i = 0; i < nxny; ++i ) {
					move( i/nx+1, i%nx );
					switch (newpop[i]) {
						case IMMUN:
							addch(' ');
							break;
						case INFEC:
							attron (COLOR_PAIR (2));
							addch('*');
							break;
						default:
							attron (COLOR_PAIR (3));
							addch('.');
							break;
					}
				}
				attron (COLOR_PAIR (1));
				mvprintw(0,0," beta: %.4f",sum_beta/double(num_INFEC));
				mvprintw(0,20,"tauI: %.4f",sum_tauI/double(num_INFEC));
				mvprintw(0,40,"time: %d",t);
				refresh();
				usleep(1000);

			} else {
				printf("%d\t%d\t%d\t%d\t%.8f\t%.8f\n",
					t,
					num_SUSCE,
					num_INFEC,
					num_IMMUN,
					sum_beta/double(num_INFEC),
					sum_tauI/double(num_INFEC));
			}
		}

		Assert( num_SUSCE + num_INFEC + num_IMMUN == nxny , "conservation" );
		Assert( num_SUSCE >=0 , "positive SUSCE" );
		Assert( num_INFEC >=0 , "positive INFEC" );
		Assert( num_IMMUN >=0 , "positive IMMUN" );
	}
	if (WITH_GRAPHICS) { endwin(); }
}

int main(int argc, const char** argv) {
	const char* opt_message[] = {
		"Program: epi, November 21, 2009",
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

