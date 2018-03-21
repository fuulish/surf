#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #ifdef HAVE_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
int multivariate_f (const gsl_vector * x, void *params, gsl_vector * f);
double c_opt_distance_to_surface_gsl( double *init_guess, double *mepos, double *pos,
                                      double *zeta, long int natoms, double surfcut, double *pbc,
                                      long int periodic, double *bnds, double xtol, double ctol );
void print_state (size_t iter, gsl_multiroot_fsolver * s);
// #endif

typedef struct {
  double *ref;
  double *pbc;
} my_func_data_type;

typedef struct {
  double *pos;
  double *zeta;
  double surfcut;
  double * pbc;
  long int periodic;
  long int natoms;
  double *ref;
} my_constraint_data_type;

#define sqr(x) ((x)*(x))
#define DIM 3
#define PI 3.141592653589793

double c_coarse_grained_density( double *mepos, double *pos, double *zeta, long int natoms,
                                   double *pbc, long int periodic, long int calc_grad, double *grad );

void calc_distance(double * dst, double * x, double * pos, double * pbc);
double vector_norm( double * vec, int len );

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
#pragma omp parallel for default(none) \
    private(a,distance) shared(pos,pbc,periodic,natoms,zeta,mepos,density,grad,calc_grad) \
        // schedule(guided, surface.n[2])
    // schedule(dynamic)
// #endif
    //this natoms here is already the one accounting for number of atoms in mask only

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

#pragma omp critical
{

        // density += scale * ( prefactor * exp( sqr( distance ) / (mttsqzeta)) - cutshft );

// #pragma omp atomic update
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

double vector_norm( double *vec, int len )
{
  double nrm = 0.;

  for( int i=0; i<len; ++i )
    nrm += sqr(vec[i]);

  return sqrt(nrm);
}

// #ifdef HAVE_GSL
int multivariate_f (const gsl_vector * x, void *params, gsl_vector * f)
{

  my_constraint_data_type *d = (my_constraint_data_type *) params;

  double grad[DIM];
  double spos[DIM];

  spos[0] = gsl_vector_get (x, 0);
  spos[1] = gsl_vector_get (x, 1);
  spos[2] = gsl_vector_get (x, 2);

  double lambda = gsl_vector_get( x, 3 );

  double density = c_coarse_grained_density( spos, d->pos, d->zeta, d->natoms, d->pbc, d->periodic, 1, grad );

  double dx[DIM];
  calc_distance( dx, d->ref, spos, d->pbc );
  double dst = vector_norm(dx, DIM);
  double dinv = 1. / dst;

  const double y0 = grad[0] * lambda + dinv * dx[0];
  const double y1 = grad[1] * lambda + dinv * dx[1];
  const double y2 = grad[2] * lambda + dinv * dx[2];
  const double y3 = density - d->surfcut;

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  gsl_vector_set (f, 2, y2);
  gsl_vector_set (f, 3, y3);

  return GSL_SUCCESS;
}

// #define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-06

double c_opt_distance_to_surface_gsl( double *init_guess, double *mepos, double *pos,
                                      double *zeta, long int natoms, double surfcut, double *pbc,
                                      long int periodic, double *bnds, double xtol, double ctol )
{
  my_constraint_data_type cons_data;
  cons_data.pos = pos ;
  cons_data.zeta = zeta;
  cons_data.surfcut = surfcut;
  cons_data.pbc = pbc;
  cons_data.periodic = periodic;
  cons_data.ref = mepos;
  cons_data.natoms = natoms;

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t iter = 0;

  const size_t n = 4;

  gsl_multiroot_function f = {&multivariate_f, n, &cons_data};

  double x_init[4] = {init_guess[0], init_guess[1], init_guess[2], 0.};
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  gsl_vector_set (x, 2, x_init[2]);
  gsl_vector_set (x, 3, x_init[3]);

  // T = gsl_multiroot_fsolver_hybrid;
  T = gsl_multiroot_fsolver_hybrids;
  // T = gsl_multiroot_fsolver_dnewton;
  // T = gsl_multiroot_fsolver_broyden;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);

  // print_state (iter, s);

  do {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      if (status)   /* check if solver is stuck */
        break;

      status = gsl_multiroot_test_residual (s->f, xtol);

    } while (status == GSL_CONTINUE && iter < 1000);

  if( iter == 1000 )
    printf("Exceeded maximum number of iterations!\n");
  // printf ("status = %s\n", gsl_strerror (status));

  init_guess[0] = gsl_vector_get( s->x, 0 );
  init_guess[1] = gsl_vector_get( s->x, 1 );
  init_guess[2] = gsl_vector_get( s->x, 2 );

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  double dst[DIM];
  calc_distance( dst, init_guess, mepos, pbc );

  // FUDO| in principle could do sign-flip here as well
  // real density = get_coarse_grained_density( init_guess, cons_data.mask, cons_data.atoms, cons_data.zeta, cons_data.pbc, cons_data.periodic, NULL );
  // if ( ( density - surfcut ) > xtol )
  //   printf("density value: %f surfcut: %f iterations: %i distance: %g\n", density, surfcut, iter, dst);

  return vector_norm( dst, DIM );

}

void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3lu x = % .3f % .3f % .3f % .3f"
          "f(x) = % .3e % .3e % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2),
          gsl_vector_get (s->f, 3));
}
// #endif
