AC_INIT(qdipole, version-0.1.0)

# get the C compiler
AC_PROG_CC

# get the Fortran compiler
AC_PROG_FC

# check if the compiler supports -fPIC
echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} $CFLAGS -fPIC -c conftest.c 2>&1`"; then
  CFLAGS="$CFLAGS -fPIC"
fi

AC_OUTPUT(Makefile)
