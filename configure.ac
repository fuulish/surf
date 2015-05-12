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

AC_ARG_ENABLE([depends],
              AS_HELP_STRING([--enable-depends], [enable single precision]))

              AS_IF([test "x$enable_depends" = "xyes"], [
                     LIBS="$PWD/depends/lib/libeinspline.a $PWD/depends/lib/libxdrfile.a $LIBS"
                C_INCLUDE_PATH="$C_INCLUDE_PATH:$PWD/depends/include"
                CFLAGS="$CFLAGS -I$PWD/depends/include"
                CXXFLAGS="$CXXFLAGS -I$PWD/depends/include"
                LDFLAGS="$LDFLAGS -L$PWD/depends/lib"
                LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PWD/depends/lib"
                #LDADD="$PWD/depends/lib/libeinspline.a $PWD/depends/lib/libxdrfile.a"
])

### # adding the lib to the files to link
### #LIBS="-leinspline -lxdrfile $LIBS"
###
### # introduce the optional configure parameter for a non-standard install prefix of XXX
### AC_ARG_WITH([depends],
###     [AS_HELP_STRING([--with-depends=prefix],
###         [try this for a non-standard install prefix of the libraries in <SURF>/depends])],
###     [DEPENDSPATHSET=1],
###     [DEPENDSPATHSET=0])
###
### # if optional parameter used, extend path flags for compliler and linker
### if test $DEPENDSPATHSET = 1 ; then
###     # extend the compiler and linker flags according to the path set
###     AM_C_INCLUDE_PATH="$AM_C_INCLUDE_PATH:$with_depends/include"
###     AM_CFLAGS="$AM_CFLAGS -I$with_depends/include"
###     AM_CXXFLAGS="$AM_CXXFLAGS -I$with_depends/include"
###     AM_LDFLAGS="$AM_LDFLAGS -L$with_depends/lib"
###
###     C_INCLUDE_PATH="$C_INCLUDE_PATH:$with_depends/include"
###     CFLAGS="$CFLAGS -I$with_depends/include"
###     CXXFLAGS="$CXXFLAGS -I$with_depends/include"
###     LDFLAGS="$LDFLAGS -L$with_depends/lib"
###     LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$with_depends/lib"
### fi
###
### AC_SUBST([LIBS])
### AC_SUBST([AM_C_INCLUDE_PATH])
### AC_SUBST([AM_CFLAGS])
### AC_SUBST([AM_CXXFLAGS])
### AC_SUBST([AM_LDFLAGS])

AC_ARG_ENABLE([xdr],
              AS_HELP_STRING([--enable-xdr], [enable XDR file support]))

              AS_IF([test "x$enable_xdr" = "xyes"], [
              AC_SEARCH_LIBS([read_xtc], [xdrfile], [CFLAGS="$CFLAGS -DHAVE_XDRFILE"], [
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
