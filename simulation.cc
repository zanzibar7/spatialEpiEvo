/*
2010-02-06
WARNING: there seems to be a bug in this
code, making it difficult to identify the
threshold values of the transmission rate
and recovery rate parameters beta and tauI
*/

enum state { SUSCE, INFEC, IMMUN };

double min(double x,double y){ return (y<x) ? y : x;}

int simulation(
	const int nx,
	const int ny,
	const int t_max,
	const double mu, // mutation rate
	const double first_beta,
	const double first_tauI,
	const double tauR,
	const int icStyle, // initial condition style
	const int sleeptime  // time to sleep between updates
	) {
	// start windowing
	int refresh_period = -1;
	if (WITH_GRAPHICS) {
		refresh_period = 10;
		initscr();
		start_color(); 
		clear();
   		init_pair (1, COLOR_WHITE, COLOR_BLACK);
   		init_pair (2, COLOR_RED, COLOR_BLACK);
   		init_pair (3, COLOR_BLUE, COLOR_BLACK);
	}
	/*
	Currently use a fixed time-step of 1 for update
	sweeps.  All other time scales must be relative 
	to this, and longer timescales lead to higher
	resolution.
	*/


	// configuration

	const double timestep_length = .0025;

	const double beta_lowerbound = 0.;
	const double beta_upperbound = 10.;
	const double tauI_lowerbound = 0.;
	const double tauI_upperbound = 4.;

	const int nxny = nx*ny;

	const int num_neighbors=8;
	int neighbors[num_neighbors];
	{
		neighbors[0] = +nx+1;
		neighbors[1] = +nx;
		neighbors[2] = +nx-1;
		neighbors[3] = +1;
		neighbors[4] = -1;
		neighbors[5] = -nx+1;
		neighbors[6] = -nx;
		neighbors[7] = -nx-1;
	}


	// ---- declare data structures -------------
	/*
	Model state is represented by these 5 arrays.
	Each site has genotype (beta, tauI).
	The main state is controlled by pop.  tic_INFEC
	and tic_IMMUN control auxillary state variables.

	Indexing is based on a torus without boundaries or jumps.
	Any lattice for could be used.
	*/
	double* __restrict__ beta  = new double[nxny];
	double* __restrict__ tauI  = new double[nxny];
	double* __restrict__ tic_INFEC = new double[nxny];
	double* __restrict__ tic_IMMUN = new double[nxny];

	/* For fast updating, we keep 2 storage arrays
	   and switch between them to do the updates
	*/
	state**  popf = new state*[2];
	popf[0]   = new state[nxny];
	popf[1]   = new state[nxny];
	state* pop    = popf[0];
	state* newpop = popf[1];


	// ----   initialization --------------------

	int num_INFEC = 0;
	int num_IMMUN = 0;
	int num_SUSCE = nxny;
	double sum_tauI = 0.;
	double sum_beta = 0.;

	if (not WITH_GRAPHICS) {
		refresh_period = 1./timestep_length/8.;
	}

	{
		int i;
		for ( i = 0; i < nxny; ++i ) {
			popf[0][i] = SUSCE;
			popf[1][i] = SUSCE;
			tic_INFEC[i] = 0;
			tic_IMMUN[i] = 0;
		}

		switch ( icStyle ) {
			case 1:
				// start with 1 infected cell
				// in the middle of the lattice.
				i = nxny/2+nx/2;
				newpop[i] = INFEC;
				--num_SUSCE;
				++num_INFEC;
				beta[i] = first_beta;
				tauI[i] = first_tauI;
				tic_INFEC[i] = tauI[i];

				break;

			case 2:
				// start with a line of infection
				// and a 1-sided wall of resistance
				for ( i=nx/2; i<nxny; i+=nx) {
					for (int j =1; j < nx/2; ++j ) { 
						newpop[i-j] = IMMUN;
						--num_SUSCE;
						++num_IMMUN;
						tic_IMMUN[i-j] = tauR;
					}
					newpop[i] = INFEC;
					--num_SUSCE;
					++num_INFEC;
					beta[i] = first_beta;
					tauI[i] = first_tauI;
					tic_INFEC[i] = tauI[i];
				}
				break;
			default:
				Assert(false,"initial condition failure");
		}
	}


	// main loop
	double t_time = 0.;
	for ( int t_sweep = 1; t_sweep < t_max and num_INFEC > 0; ++t_sweep ) {

		// update outputs
		if ( 0 == t_sweep%refresh_period ) {
			sum_tauI = 0.;
			sum_beta = 0.;
			//num_INFEC = 0;
			//num_IMMUN = 0;
			//num_SUSCE = 0;
			for ( int i = 0; i < nxny; ++i) {
				switch (pop[i]) {
					case IMMUN:
						//++num_IMMUN;
						break;
					case SUSCE:
						//++num_SUSCE;
						break;
					case INFEC:
						//++num_INFEC;
						sum_tauI += tauI[i];
						sum_beta += beta[i];
						break;
					default: Assert(false,"982");
				}
			}
			if (WITH_GRAPHICS) {
				for (int i = 0; i < nxny; ++i ) {
					move( i/nx+1, i%nx );
					switch (pop[i]) {
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
				mvprintw(0,0," beta: %.3f ",sum_beta/double(num_INFEC));
				mvprintw(0,16,"tauI: %.3f ",sum_tauI/double(num_INFEC));
				mvprintw(0,30,"R: %.2f ",
					num_neighbors*sum_beta/double(num_INFEC)*sum_tauI/double(num_INFEC) );
				mvprintw(0,40,"time: %.0f",t_time);
				refresh();
				usleep(sleeptime);

			} else {
				printf("%f\t%d\t%d\t%d\t%.8f\t%.8f\t%d\n",
					t_time,
					num_SUSCE,
					num_INFEC,
					num_IMMUN,
					sum_beta/double(num_INFEC),
					sum_tauI/double(num_INFEC),
					t_sweep);
			}
		}

		// test sanity
		Assert( num_SUSCE + num_INFEC + num_IMMUN == nxny, "conservation" );
		Assert( num_SUSCE >=0, "SUSCE must be positive" );
		Assert( num_INFEC >=0, "INFEC must be positive" );
		Assert( num_IMMUN >=0, "IMMUN must be positive" );

		// update and blit
		t_time += timestep_length;
		pop = popf[t_sweep%2];
		newpop = popf[(t_sweep+1)%2];
		Assert( pop != newpop, "flip collision");

		// for each lattice node
		for ( int i = 0; i < nxny; ++i) {
			double tic = timestep_length; // one time step
			switch (pop[i]) {
				case IMMUN:
					tic_IMMUN[i] -= tic;
					if ( 0 >= tic_IMMUN[i] ) {
						newpop[i] = SUSCE;
						--num_IMMUN;
						++num_SUSCE;
					} else {
						newpop[i] = IMMUN;
					}
				//printf("%d\t%d\t%d\t%d\t%.4f\t%.4f\n", t, num_SUSCE, num_INFEC, num_IMMUN, sum_beta/double(num_INFEC), sum_tauI/double(num_INFEC));
					break;
				case INFEC:
					tic_INFEC[i] -= tic;
					if ( 0 >= tic_INFEC[i] ) {
						newpop[i] = IMMUN;
						--num_INFEC;
						++num_IMMUN;
						tic_IMMUN[i] = tauR;
						//printf("%dh\t%d\t%d\t%d\t%.4f\t%.4f\n", t, num_SUSCE, num_INFEC, num_IMMUN, sum_beta/double(num_INFEC), sum_tauI/double(num_INFEC));
					} else {
						newpop[i] = INFEC;
					}
					break;
				case SUSCE:
					newpop[i] = SUSCE;
					{
						// Check if a neighbor infects you.
						// Use random numbers to determine which
						// neighbor infects you first.

						for (int j = 0; j < num_neighbors; ++j ) { // for each neighbor
							int h = i+neighbors[j];
							if ( h < 0 or h >= nxny ) { h = (h + nxny ) % nxny ; }

							if (INFEC!=pop[h]) continue;

							double toc = gsl_ran_exponential(random_number_generator,1./beta[h]);
							if (toc < tic) {
								tic = toc;
								newpop[i] = INFEC;

								// assign bounded mutant beta
								beta[i] = beta[h]*(1+gsl_ran_gaussian(random_number_generator, mu));
								if (beta[i] < beta_lowerbound) { beta[i] = beta_lowerbound; }
								else
								if (beta[i] > beta_upperbound) { beta[i] = beta_upperbound; }

								// assign bounded mutant tauI
								tauI[i] = tauI[h]*(1+gsl_ran_gaussian(random_number_generator, mu));
								if (tauI[i] < tauI_lowerbound) { tauI[i] = tauI_lowerbound; }
								else
								if (tauI[i] > tauI_upperbound) { tauI[i] = tauI_upperbound; }

								// continue with the rest of the for-loop in case
								// a different neighbor infected you first.
							}
						}
						if ( INFEC == newpop[i] ) {
							// if the site is newly infected, do the rest of the state update
							--num_SUSCE;
							++num_INFEC;
							tic_INFEC[i] = tauI[i];
						}
					}
				//printf("%d\t%d\t%d\t%d\t%.4f\t%.4f\n", t_sweep, num_SUSCE, num_INFEC, num_IMMUN, sum_beta/double(num_INFEC), sum_tauI/double(num_INFEC));
					break;
				default:
					Assert(false, "invalid state");
			}
		}
	}
	if (WITH_GRAPHICS) { usleep(100000); endwin(); }
}

