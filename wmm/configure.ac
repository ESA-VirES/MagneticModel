AC_INIT(wmm, version-2010)
AC_MSG_NOTICE([Hello, world.])

# get info about the build host
#AC_CANONICAL_BUILD

# get the C compiler
AC_PROG_CC

# check if the compiler support -fPIC
echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} $CFLAGS -fPIC -c conftest.c 2>&1`"; then
  CFLAGS="$CFLAGS -fPIC"
fi

AC_OUTPUT(Makefile)
AC_OUTPUT(src/lib/Makefile)
AC_OUTPUT(src/cli/Makefile)
