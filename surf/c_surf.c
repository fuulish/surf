#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define sqr(x) ((x)*(x))
#define DIM 3
#define PI 3.141592653589793

double c_coarse_grained_density( double *mepos, double *pos, double *zeta, long int natoms,
                                   double *pbc, long int periodic, long int calc_grad, double *grad );

void calc_distance(double * dst, double * x, double * pos, double * pbc);

double c_coarse_grained_density( double *mepos, double *pos, double *zeta, long int natoms,
                                   double *pbc, long int periodic, long int calc_grad, double *grad )
{
    int a;
    // double prefactor, dummy, cutshft;

    double distance;
    double density = 0.0;

    if ( grad != NULL ) {
      grad[0] = 0.;
      grad[1] = 0.;
      grad[2] = 0.;
    }

// #ifdef OPENMP
// #pragma omp parallel for default(none) \
//     private(a,distance) shared(atoms,pbc,periodic,surface,natoms,mask,zeta,mepos,density,grad) // \
//         // schedule(guided, surface.n[2])
//     // schedule(dynamic)
// #endif
//     //this natoms here is already the one accounting for number of atoms in mask only

    for ( a=0; a<natoms; a++ ) {

        double sqzeta = sqr(zeta[a]);
        double trplzt = 3*zeta[a];

        double dst[DIM];
        calc_distance(dst, mepos, &(pos[a*DIM]), pbc);

        distance = sqrt( sqr(dst[0]) + sqr(dst[1]) + sqr(dst[2]) );

        if ( distance > trplzt )
            continue;

        double mttsqzeta = -2. * sqzeta;
        double dummy = 2. * PI * sqzeta;

        double prefactor = 1. / dummy / (sqrt(dummy));
        double cutshft = prefactor * exp( sqr( trplzt ) / (mttsqzeta));

        double tmpdens = prefactor * exp( sqr( distance ) / (mttsqzeta)) - cutshft;
        // the below formula can be simplified, check here
#pragma omp atomic update
        // density += scale * ( prefactor * exp( sqr( distance ) / (mttsqzeta)) - cutshft );
        density += tmpdens;

        if ( calc_grad ) {
          /* calculate gradient and save result in grad array
           * can make use of tmpdens and distance vector */

          double fact = -2. / mttsqzeta;

          grad[0] += dst[0] * fact * tmpdens;
          grad[1] += dst[1] * fact * tmpdens;
          grad[2] += dst[2] * fact * tmpdens;
        }
    }

    return density;

}

void calc_distance(double * dst, double * x, double * pos, double * pbc)
{
  dst[0] = x[0] - pos[0];
  dst[1] = x[1] - pos[1];
  dst[2] = x[2] - pos[2];

  dst[0] -= pbc[0] * round( dst[0] / pbc[0] );
  dst[1] -= pbc[1] * round( dst[1] / pbc[1] );
  dst[2] -= pbc[2] * round( dst[2] / pbc[2] );
}
