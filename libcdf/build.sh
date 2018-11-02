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

BUILD_OPTIONS="SHARED=yes FORTRAN=no CURSES=no"
make OS=$OS ENV=$ENV "CC_$OS_$ENV"=$CC "FC_$OS"=$FC $BUILD_OPTIONS all
make test
make INSTALLDIR="$PREFIX" install
