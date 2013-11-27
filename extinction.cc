int extinction(
	const int nx,
	const int ny,
	const int t_max,
	const double mu,
	const double first_beta,
	double first_tauI,
	const int icStyle,
	const int sleeptime
	) {

	const double tauR = 1.;

	const int timeRes = 24.;
	const int refresh_period = 10;

	const double beta_lowerbound = 0.;
	const double beta_upperbound = 10.;
	const double tauI_lowerbound = 0.;
	const double tauI_upperbound = 4.;

	const int nxny = nx*ny;

	/* declare data structures */
	double* __restrict__ beta  = new double[nxny];
	double* __restrict__ tauI  = new double[nxny];
	double* __restrict__ tic_INFEC = new double[nxny];
	double* __restrict__ tic_IMMUN = new double[nxny];

	state**  popf = new state*[2];
	popf[0]   = new state[nxny];
	popf[1]   = new state[nxny];
	state* pop    = popf[0];
	state* newpop = popf[1];

	/* initialization */

	double tic_baseline = 1./(first_beta*timeRes);

	int neighbors[8];
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

	int num_INFEC = 0;
	int num_IMMUN = 0;
	int num_SUSCE = nxny;
	double sum_tauI = 0.;
	double sum_beta = 0.;

	{ int i;
	for ( i = 0; i < nxny; ++i ) {
		popf[0][i] = SUSCE;
		popf[1][i] = SUSCE;
		tic_INFEC[i] = 0;
		tic_IMMUN[i] = 0;
	}

	switch ( 1 ) {
		case 1:
			i = nxny/2+nx/2; pop[i] = INFEC;
			--num_SUSCE; ++num_INFEC;
			beta[i] = first_beta;
			tauI[i] = first_tauI;

			break;

		case 2:
			for ( i=nx/2; i<nxny; i+=nx) {
				for (int j =1; j < nx/2; ++j ) { 
					pop[i-j] = IMMUN;
					--num_SUSCE; ++num_IMMUN;
					tic_IMMUN[i-j] = tauR;
				}
				pop[i] = INFEC;
				--num_SUSCE; ++num_INFEC;
				beta[i] = first_beta;
				tauI[i] = first_tauI;
				tic_INFEC[i] = tauI[i];
			}
			break;
	} }


	/* main loop */
	double time = 0.;
	for ( int t = 0; t < t_max; ++t ) {
		pop = popf[t%2];
		newpop = popf[(t+1)%2];
		Assert( pop != newpop, "flip collision");

		time += tic_baseline;
		
		for ( int i = 0; i < nxny; ++i) {
			double tic = tic_baseline; // one time step
			switch (pop[i]) {
				case IMMUN:
					tic_IMMUN[i] -= tic;
					if ( 0 >= tic_IMMUN[i] ) {
						newpop[i] = SUSCE;
						--num_IMMUN;
						++num_SUSCE;
					} else { newpop[i] = IMMUN; }
					break;
				case INFEC:
					tic_INFEC[i] -= tic;
					if ( 0 >= tic_INFEC[i] ) {
						newpop[i] = IMMUN;
						--num_INFEC;
						++num_IMMUN;
						tic_IMMUN[i] = tauR;
					} else { newpop[i] = INFEC; }
					break;
				case SUSCE:
					newpop[i] = SUSCE;
					{
						// Check if a neighbor infects you.
						// Use random numbers to determine which
						// neighbor infects you first.

						for (int j = 0; j < 8; ++j ) { // for each neighbor
							int h = i+neighbors[j];
							if ( h < 0 or h >= nxny ) { h = (h + nxny ) % nxny ; }

							if (INFEC!=pop[h]) continue;

							double toc = gsl_ran_exponential(random_number_generator,1./beta[h]);
							if (toc > tic) {
								continue;
							} else {
								tic = toc;
								newpop[i] = INFEC;

								beta[i] = beta[h]*(1+gsl_ran_gaussian(random_number_generator, mu));
								if (beta[i] < beta_lowerbound) { beta[i] = beta_lowerbound; }
								else
								if (beta[i] > beta_upperbound) { beta[i] = beta_upperbound; }

								tauI[i] = tauI[h]*(1+gsl_ran_gaussian(random_number_generator, mu));
								if (tauI[i] < tauI_lowerbound) { tauI[i] = tauI_lowerbound; }
								else
								if (tauI[i] > tauI_upperbound) { tauI[i] = tauI_upperbound; }
							}
						}
						if ( INFEC == newpop[i] ) { // if the site is newly infected
							--num_SUSCE;
							++num_INFEC;
							tic_INFEC[i] = tauI[i];
						}
					}
					break;
				default:
					Assert(false, "invalid state");
			}
		}

		/////////////////////////////////////////////////////

		if ( 0 == num_INFEC ) {
			printf("%f\n", time);
			delete beta; delete tauI; delete tic_INFEC; delete tic_IMMUN; delete popf[0]; delete popf[1]; delete popf;
			return 1;
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
				}
			}
			tic_baseline = 1./(sum_beta/double(num_INFEC)*timeRes);
		}

		Assert( num_SUSCE + num_INFEC + num_IMMUN == nxny , "conservation" );
		Assert( num_SUSCE >=0 , "positive SUSCE" );
		Assert( num_INFEC >=0 , "positive INFEC" );
		Assert( num_IMMUN >=0 , "positive IMMUN" );
	}

	printf("%f\n", time);
	delete beta; delete tauI; delete tic_INFEC; delete tic_IMMUN; delete popf[0]; delete popf[1]; delete popf;
	return 1;
}

