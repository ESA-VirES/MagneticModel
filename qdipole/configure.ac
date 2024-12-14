AC_INIT([qdipole],[0.6.0])

# get the C compiler
AC_PROG_CC

# get the Fortran compiler
AC_PROG_FC

# check if the compiler supports -fPIC
echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} $CFLAGS -fPIC -c conftest.c 2>&1`"; then
  CFLAGS="$CFLAGS -fPIC"
fi

AC_CONFIG_FILES([Makefile])
AC_CONFIG_HEADERS([qdipole_conf.h])
AC_OUTPUT
