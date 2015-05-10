AC_INIT([surf], [0.1], [uhlig.frank@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC([gcc icc])
AM_PROG_CC_C_O
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])

AC_SEARCH_LIBS([sqrt], [m], [], [
  AC_MSG_ERROR([unable to find the sqrt function. Is C math library installed?])
])

AC_SEARCH_LIBS([gzopen], [z], [], [
  AC_MSG_ERROR([unable to find the gzopen function. Is zlib installed?])
])

AC_ARG_ENABLE([xdr],
              AS_HELP_STRING([--enable-xdr], [enable XDR file support]))

              AS_IF([test "x$enable_xdr" = "xyes"], [
              AC_SEARCH_LIBS([read_xtc], [xdrfile], [CFLAGS="$CFLAGS -DXDRCAP"], [
                AC_MSG_ERROR([unable to find the read_xtc function. Is xdrfile library installed?])
              ])

])

AC_ARG_ENABLE([einspline],
              AS_HELP_STRING([--enable-einspline], [enable EINSPLINE file support]))

              AS_IF([test "x$enable_einspline" = "xyes"], [
              AC_SEARCH_LIBS([create_UBspline_3d_s], [einspline], [CFLAGS="$CFLAGS -DHAVE_EINSPLINE"], [
                AC_MSG_ERROR([unable to find the create_UBspline_3d_s function. Is EINSPLINE library installed?])
              ])

])

AC_ARG_ENABLE([float],
              AS_HELP_STRING([--enable-float], [enable single precision]))

              AS_IF([test "x$enable_float" = "xyes"], [
              CFLAGS="$CFLAGS -DSINGLE"
])

AC_ARG_ENABLE([openmp],
              AS_HELP_STRING([--enable-openmp], [enable OpenMP multithreading (currently only for ILI)]))

              AS_IF([test "x$enable_openmp" = "xyes"], [
              CFLAGS="$CFLAGS -DOPENMP -fopenmp"
])

AC_ARG_ENABLE([gsl],
              AS_HELP_STRING([--enable-gsl], [enable GSL support]))

              AS_IF([test "x$enable_gsl" = "xyes"], [
#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
              AC_SEARCH_LIBS([cblas_dgemm], [gslcblas], [CFLAGS="$CFLAGS -DHAVE_GSL"], [
              AC_MSG_ERROR([unable to find the GSL CBLAS.])
              ])
              #check here, this is a bit strange, why does it work this way, but not the other way around?
              AC_SEARCH_LIBS([gsl_check_range], [gsl], [], [
              AC_MSG_ERROR([unable to find the GSL.])
              ])
])

AC_OUTPUT
