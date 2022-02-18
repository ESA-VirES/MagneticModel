#
# qdipole conda build script
#
cd "$SRC_DIR/qdipole"
./configure --prefix="$PREFIX" "LDFLAGS=$LDFLAGS_USED"
make build
make install
