enum state { SUSCE, INFEC, IMMUN };

// 8*8 = 64 long.  Each block of 8 is a random permutation of integers 0-7.
// Used to introduce a little disorder to transmission process
const int32_t rdzer[] = {
	7, 6, 1, 3, 2, 5, 0, 4,
	5, 0, 6, 1, 4, 7, 3, 2,
	1, 7, 6, 2, 5, 0, 3, 4,
	4, 3, 7, 6, 5, 2, 0, 1,
	0, 4, 2, 1, 6, 5, 3, 7,
	3, 2, 4, 1, 5, 6, 7, 0,
	2, 5, 6, 7, 4, 0, 3, 1,
	6, 1, 2, 7, 5, 3, 4, 0,
};

int simulation(
	const int32_t nx,
	const int32_t ny,
	const uint64_t t_max,
	const double mu,
	const double first_beta,
	double first_tauI
	) {
	// start windowing
	uint32_t refresh_period;
	if (WITH_GRAPHICS) {
		initscr();
		start_color(); 
		refresh_period = 1;
   		init_pair (1, COLOR_WHITE, COLOR_BLACK);
   		init_pair (2, COLOR_MAGENTA, COLOR_BLACK);
   		init_pair (3, COLOR_GREEN, COLOR_BLACK);
		timeout(0);
	} else {
		refresh_period = 100;
	}

	const uint32_t tauR = 255;
	const double beta_max = 3.;
	const double beta_min = 0.;
	const double tauI_max = 100.;
	const double tauI_min = 0.;

	const int32_t nxny = nx*ny;

	uint32_t sleeptime=10000;

	double* __restrict__ beta  = new double[nxny];
	double* __restrict__ tauI  = new double[nxny];
	int16_t* __restrict__ tic_INFEC = new int16_t[nxny];
	int16_t* __restrict__ tic_IMMUN = new int16_t[nxny];

	state**  popf = new state*[2];
	popf[0]   = new state[nxny];
	popf[1]   = new state[nxny];
	state* pop    = popf[0];
	state* newpop = popf[1];

	int32_t neighbors[8] = {nx+1, nx, nx-1, 1, -1, 0-nx+1, 0-nx, 0-nx-1};

	int32_t num_INFEC = 0;
	int32_t num_IMMUN = 0;
	int32_t num_SUSCE = nxny;
	double sum_tauI = 0.;
	double sum_beta = 0.;

	// initial conditions
	for (int32_t i = 0; i < nxny; ++i ) {
		popf[0][i] = SUSCE;
		popf[1][i] = SUSCE;
		tic_INFEC[i] = 0;
		tic_IMMUN[i] = 0;
	}
	{ uint32_t i;
		i = nxny/2+nx/2;
		pop[i] = INFEC;
		beta[i] = first_beta;
		tauI[i] = first_tauI;
		--num_SUSCE;
		++num_INFEC;

		i = 2; pop[i] = IMMUN; tic_IMMUN[i] = tauR/2;
		--num_SUSCE; ++num_IMMUN;
		i = 3; pop[i] = IMMUN; tic_IMMUN[i] = tauR/2;
		--num_SUSCE; ++num_IMMUN;
	}


	for ( uint64_t t = 0; t < t_max; ++t ) {
		pop = popf[t%2];
		newpop = popf[(t+1)%2];
		Assert( pop != newpop, "flip collision");

		for ( int32_t i = 0; i < nxny; ++i) {
			/*
			// Added holes to the lattice to simulate reduction of 
			// erratic heartbeat by introduction of heterogeneity.
			// 2010-04-17 TCR
			if ( 2 > i % 4 and 2 > (i/nx)%4 ) {
				pop[i] = IMMUN;
				tic_IMMUN[i] = 1;
				continue;
			}
			*/
			switch (pop[i]) {
				case IMMUN:
					tic_IMMUN[i]--;
					if ( 0 >= tic_IMMUN[i] ) {
						newpop[i] = SUSCE;
						--num_IMMUN;
						++num_SUSCE;
					} else { newpop[i] = IMMUN; }
					break;
				case INFEC:
					tic_INFEC[i]--;
					if ( 0 >= tic_INFEC[i] ) {
						newpop[i] = IMMUN;
						tic_IMMUN[i] = tauR;
						--num_INFEC;
						++num_IMMUN;
					} else { newpop[i] = INFEC; }
					break;
				case SUSCE:
					newpop[i] = SUSCE;
					{
						double tic = 1.0; // one time step
						for (uint8_t j = 0; j < 8; ++j ) { 
							int32_t h = i+neighbors[rdzer[(((i+t)%8)*8)+j]]; // Pseudo-randomization of neighbor update order
							if ( h < 0 or h >= nxny ) { h = (h + nxny ) % nxny ; }

							if (INFEC!=pop[h]) continue;

							double toc = gsl_ran_exponential(random_number_generator,1./beta[h]);
							if (toc < tic) {
								tic = toc;
								newpop[i] = INFEC;
								beta[i] = beta[h]*(1+gsl_ran_gaussian(random_number_generator, mu));
								if (beta[i] < beta_min) { beta[i] = beta_min; }
								else if (beta[i] > beta_max) { beta[i] = beta_max; }

								tauI[i] = tauI[h]*(1+gsl_ran_gaussian(random_number_generator, mu));
								if (tauI[i] < tauI_min) { tauI[i] = tauI_min; }
								if (tauI[i] > tauI_max) { tauI[i] = tauI_max; }
							}
						}
						if ( INFEC == newpop[i] ) {
							--num_SUSCE;
							++num_INFEC;
							tic_INFEC[i] = uint32_t(tauI[i]);
							break;
						}
					}
					break;
				default:
					Assert(false, "invalid state");
			}
		}

		/////////////////////////////////////////////////////

		if ( 0 == num_INFEC ) { break; }

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
							addch('*');
							break;
					}
				}
				attron (COLOR_PAIR (1));
				mvprintw(0,0," beta: %.3f",sum_beta/double(num_INFEC));
				mvprintw(0,16,"tauI: %.3f",sum_tauI/double(num_INFEC)/double(tauR));
				mvprintw(0,30,"R: %.2f",8*(sum_tauI/double(num_INFEC))*(sum_beta/double(num_INFEC)));
				mvprintw(0,43,"time: %d",t);
				move(0,0);
				refresh();
				if ( sleeptime ) { usleep(sleeptime); }
				{ char c = getch(); 
					switch (c) {
						case ' ':
							for (char cc; (cc=getchar()); ) {
								usleep(10);
								if (' ' == cc ) { break; }
							}
							break;
						case '0': sleeptime = 0; break;
						case '1': sleeptime =100; break;
						case '2': sleeptime =1000; break;
						case '3': sleeptime =10000; break;
						case '4': sleeptime =100000; break;
						case '5': sleeptime =1000000; break;
						case '6': sleeptime =10000000; break;
						case 'Q' :  //End program at control-d or Q
						case '': { endwin(); return 1; break; }
						default:
							break;
					}
				}
			} else {
				printf("%lu\t%d\t%d\t%d\t%.8f\t%.8f\n",
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
	return 0;
}

