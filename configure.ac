#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)

define([svnversion], esyscmd([sh -c "svnversion -cn"]))dnl
AC_INIT([MIA],[1.0],[green@eva.mpg.de])
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/mia.c])
AC_CONFIG_HEADER([src/config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset strstr])
AC_CHECK_LIB([m],[pow])
AC_CHECK_LIB([m],[log10])

AC_CONFIG_FILES([Makefile src/Makefile man/Makefile matrices/Makefile])
AC_OUTPUT
