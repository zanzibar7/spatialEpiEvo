/*
Added holes to the lattice to simulate reduction of 
erratic heartbeat by introduction of heterogeneity.
2010-04-17 TCR
*/
enum state { SUSCE, INFEC, IMMUN };

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
		//clear();
		refresh_period = 1;
   		init_pair (1, COLOR_WHITE, COLOR_BLACK);
   		init_pair (2, COLOR_MAGENTA, COLOR_BLACK);
   		init_pair (3, COLOR_GREEN, COLOR_BLACK);
		timeout(0);
	} else {
		refresh_period = 100;
	}

	const int tauR = 100;

	const int nxny = nx*ny;

	int sleeptime=10000;

	double* __restrict__ beta  = new double[nxny];
	double* __restrict__ tauI  = new double[nxny];
	int* __restrict__ tic_INFEC = new int[nxny];
	int* __restrict__ tic_IMMUN = new int[nxny];

	state**  popf = new state*[2];
	popf[0]   = new state[nxny];
	popf[1]   = new state[nxny];
	state* pop    = popf[0];
	state* newpop = popf[1];

	for (int i = 0; i < nxny; ++i ) {
		popf[0][i] = SUSCE;
		popf[1][i] = SUSCE;
		tic_INFEC[i] = 0;
		tic_IMMUN[i] = 0;
	}
	int neighbors[8];
	{
		neighbors[0] = 0+nx+1;
		neighbors[1] = 0+nx;
		neighbors[2] = 0+nx-1;
		neighbors[3] = 0+1;
		neighbors[4] = 0-1;
		neighbors[5] = 0-nx+1;
		neighbors[6] = 0-nx;
		neighbors[7] = 0-nx-1;
	}

	int num_INFEC = 0;
	int num_IMMUN = 0;
	int num_SUSCE = nxny;
	double sum_tauI = 0.;
	double sum_beta = 0.;

	{ int i;
	i = nxny/2+nx/2; pop[i] = INFEC; beta[i] = first_beta; tauI[i] = first_tauI;
	--num_SUSCE; ++num_INFEC;

	i = 2; pop[i] = IMMUN; tic_IMMUN[i] = tauR/2;
	--num_SUSCE; ++num_IMMUN;
	i = 3; pop[i] = IMMUN; tic_IMMUN[i] = tauR/2;
	--num_SUSCE; ++num_IMMUN;
	}


	for ( int t = 0; t < t_max; ++t ) {
		pop = popf[t%2];
		newpop = popf[(t+1)%2];
		Assert( pop != newpop, "flip collision");

		for ( int i = 0; i < nxny; ++i) {
			if ( 2 > i % 4 and 2 > (i/nx)%4 ) {
				pop[i] = IMMUN;
				tic_IMMUN[i] = 1;
				continue;
			}
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
						for (int j = 0; j < 8; ++j ) { // BUG: not randomized
							int h = i+neighbors[j];
							if ( h < 0 or h >= nxny ) { h = (h + nxny ) % nxny ; }

							if (INFEC!=pop[h]) continue;

							double toc = gsl_ran_exponential(random_number_generator,1./beta[h]);
							if (toc < tic) {
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
				mvprintw(0,0," beta: %.3f",sum_beta/double(num_INFEC)*double(tauR));
				mvprintw(0,16,"tauI: %.3f",sum_tauI/double(num_INFEC)/double(tauR));
				mvprintw(0,30,"R: %.3f",8*sum_tauI/double(num_INFEC)*sum_beta/double(num_INFEC));
				mvprintw(0,40,"time: %d",t);
				move(0,0);
				refresh();
				if ( sleeptime ) { usleep(sleeptime); }
				//nodelay(stdscr,1);
				{ char c = getch(); 
					switch (c) {
						case ' ':
							for (char cc; cc=getchar(); ) {
								usleep(10);
								if (' ' == cc ) { break; }
							}
							break;
						case '0': sleeptime = 0; break;
						case '1': sleeptime =1000; break;
						case '2': sleeptime =10000; break;
						case '3': sleeptime =100000; break;
						case '4': sleeptime =1000000; break;
						case 'Q' :  //End program at control-d or Q
						case '': { endwin(); return 1; break; }
						default:
							break;
					}
				}
			// for (int j =0;!j;){usleep(1);j=getch(); printw("%d ",j); if(j>0) j=1; else j=0; }
			// endwin(); printf("\nhitkb end\n"); return 0;


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

