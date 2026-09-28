#
# CDF library conda build script
#

# TODO: add support for additional platforms

error() {
    echo "ERROR: $*" 1>&2
    exit 1
}

if [ "$OSTYPE" = 'linux-gnu' ]
then
    OS="linux"
    if [ "$ARCH" == '32' ]
    then
        ENV="gnu32"
    elif [ "$ARCH" == '64' ]
    then
        ENV="gnu"
    else
        error "Unsupported architecture $ARCH!"
    fi
else
    error "Unsupported platform $OSTYPE!"
fi

MAKE_OPTIONS="-f Makefile-classic"
BUILD_OPTIONS="SHARED=yes FORTRAN=no CURSES=no"
make OS=$OS ENV=$ENV AR=$AR RANLIBcmd=$RANLIB LD_${OS}_${ENV}=$CC CC_${OS}_${ENV}=$CC FC_${OS}=$FC $BUILD_OPTIONS all $MAKE_OPTIONS
make test $MAKE_OPTIONS
make OS=$OS ENV=$ENV CC_${OS}_${ENV}=$CC INSTALLDIR="$PREFIX" install $MAKE_OPTIONS
