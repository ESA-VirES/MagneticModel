AC_INIT(wmm, version-2015)

# get the C compiler
AC_PROG_CC

# check if the compiler supports -fPIC
echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} $CFLAGS -fPIC -c conftest.c 2>&1`"; then
  CFLAGS="$CFLAGS -fPIC"
fi

AC_OUTPUT(Makefile)
AC_OUTPUT(src/lib/Makefile)
AC_OUTPUT(src/cli/Makefile)
