#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.60])
AC_INIT([steven], [0.0], [dses5211@gmail.com])
AM_INIT_AUTOMAKE
AC_SUBST(BONMIN_CPPFLAGS)
AC_SUBST(BONMIN_LIBS)
AC_SUBST(IPOPT_CPPFLAGS)
AC_SUBST(IPOPT_LIBS)
AM_MAINTAINER_MODE([disable])
AC_CONFIG_SRCDIR([schedule/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
LT_INIT
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_CPP

AC_LANG([C++])
AC_CHECK_HEADERS([boost/program_options.hpp],,exit) 
AC_CHECK_HEADERS([boost/numeric/ublas/matrix.hpp],,exit)
AC_CHECK_HEADERS([boost/numeric/ublas/io.hpp],,exit)
AC_CHECK_HEADERS([armadillo],,exit)
AC_CHECK_HEADERS([eigen3/Eigen/Cholesky],,exit) 

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h netinet/in.h sys/ioctl.h sys/socket.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL


AC_ARG_ENABLE(debug, 
	      AC_HELP_STRING([--enable-debug], [Enable debugging]), 
	      CFLAGS="$CFLAGS -g")

AC_ARG_ENABLE(debug,
AS_HELP_STRING(
	       [--enable-debug],
	       [enable debugging, default: no]),
	       [case "${enableval}" in
		yes) debug=true ;;
		no)  debug=false ;;
		*)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
	esac],
	[debug=false])

AM_CONDITIONAL(DEBUG, test x"$debug" = x"true")

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([floor mkdir pow socket sqrt])

AC_CONFIG_FILES([Makefile
                 cluster/Makefile
		 schedule/Makefile
                 image_cluster/Makefile
                 image_schedule/Makefile
		 maxEntropy_schedule/Makefile
		 lib/Makefile])
AC_OUTPUT
