#include "headers.h"

const int nx = 80;
const int ny = 40;
const int nxny = nx*ny;

enum state { SUSCE, INFEC, IMMUN };



gsl_rng* random_number_generator;

int main(int argc, const char** argv) {
	const double beta = .8; // large values are more likely
	const double gamma = .3;
	const double alpha = .02;
	const int t_max = 1000*nxny;

	/*
	initscr();
	clear();
	*/

	int random_number_seed = (int) time(NULL) + ( (int) getpid() << 14  ) ;

	gsl_rng_env_setup();
	//gsl_rng * rng = gsl_rng_alloc( gsl_rng_taus );
	//gsl_rng * rng = gsl_rng_alloc( gsl_rng_random64_glibc2 );
	random_number_generator = gsl_rng_alloc( gsl_rng_default );
	gsl_rng_set(random_number_generator, random_number_seed );
	
	int num_INFEC = 0;
	int num_IMMUN = 0;
	int num_SUSCE = nxny;
	int neighbors[8];

	state* __restrict__ pop = new state[nxny];

	pop[0] = INFEC;
	--num_SUSCE;
	++num_INFEC;
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
					if ( h < 0 or h >= nxny ) {
						h = (h + nxny ) % nxny ;
					}
					if ( SUSCE == pop[h] ) {
						if ( gsl_ran_bernoulli(random_number_generator, beta) ) {
							pop[h] = INFEC;
							--num_SUSCE;
							++num_INFEC;
						}
					}
				}
				if ( gsl_ran_bernoulli(random_number_generator, gamma) ) {
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
		if ( 0 == t%nxny ) {
			printf( "%d\t%d\t%d\t%d\n", t/nxny, num_SUSCE, num_INFEC, num_IMMUN );
		}
		/*
		if ( 0 == t%10 ) {
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
						addch(' ');
				}
			}
			refresh();
		}
		*/
		
	}
}


/*

discrete-time (days)


*/
