#!/bin/bash
SRC_DIR="../../lib/"
g++ -O3 -I$SRC_DIR $SRC_DIR/mt19937.C $SRC_DIR/multilevel.C  $SRC_DIR/lattice.C $SRC_DIR/utils.C $SRC_DIR/interpolator.C $1
