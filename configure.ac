AUTOMAKE_OPTIONS = subdir-objects

AC_PREREQ([2.69])
AC_INIT([ANTL], [1.0], [test@gmail.com])
AC_CONFIG_AUX_DIR([build])
AC_CONFIG_SRCDIR([include/ANTL/time.hpp])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S

# this macro figures out how to build with pthreads
#AX_PTHREAD(, )

# this macro figures out how to build with mpi
#AX_MPI(, )

# this macro figures out how to build with coverage
# AX_CODE_COVERAGE

# Checks for libraries. They all have a main function
AC_CHECK_LIB([gmp], [main])
AC_CHECK_LIB([m], [main])
AC_CHECK_LIB([mpfr], [main])
AC_CHECK_LIB([mpfi], [main])
AC_CHECK_LIB([ntl], [main])
AC_CHECK_LIB([gcov], [main])
AC_CHECK_LIB([lcov], [main])

# see https://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/External-Software.html
# run ./configure --help in the top level directory to list all options for configure
AC_ARG_WITH([mpfr],
  [AS_HELP_STRING([--with-mpfr],
    [support using mpfr as a multi precision backend @<:@default=check@:>@])],
  [],
  [with_mpfr=check])

LIBMPFR=
AS_IF([test "x$with_mpfr" != xno],
  [AC_CHECK_LIB([mpfr], [main],
    [AC_SUBST([LIBMPFR], ["-lmpfr"])
     AC_DEFINE([HAVE_LIBMPFR], [1],
               [Define if you have libmpfr])
    ],
    [if test "x$with_mpfr" != xcheck; then
       AC_MSG_FAILURE(
         [--with-mpfr was given, but test for mpfr failed])
     fi
    ], -lncurses)])
# see https://autotools.io/automake/conditional.html
AM_CONDITIONAL([EXTRADEPS], [test "x$with_mpfr" != xno])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([dup2 floor memset pow sqrt])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
