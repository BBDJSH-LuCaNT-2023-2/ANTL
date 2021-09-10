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
AX_PTHREAD(, )

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgmp':
AC_CHECK_LIB([gmp], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lmpfi':
AC_CHECK_LIB([mpfi], [main])
# FIXME: Replace `main' with a function in `-lmpfr':
AC_CHECK_LIB([mpfr], [main])
# FIXME: Replace `main' with a function in `-lntl':
AC_CHECK_LIB([ntl], [main])

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
