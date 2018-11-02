#
# qdipole conda build script
#
cd "$SRC_DIR/qdipole"
./configure --prefix="$PREFIX"
make build
make install
