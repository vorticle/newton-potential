#! /bin/bash
set -o xtrace;

CXX=g++-10
CXXFLAGS="-Wall -Wextra -DNDEBUG -DARMA_NO_DEBUG -O3 -march=native -mtune=native -flto -fopenmp"
#CXXFLAGS="-Wall -Wextra -O0 -ggdb"

AR=ar
ARFLAGS=cr
RANLIB=ranlib

cd geometry;
$CXX $CXXFLAGS -I.. -c interval_quadrules.cpp -o interval_quadrules.o;
$CXX $CXXFLAGS -I.. -c  cubical_quadrules.cpp -o  cubical_quadrules.o;
$CXX $CXXFLAGS -I.. -c             sphere.cpp -o             sphere.o;
cd ..;

cd math;
$CXX $CXXFLAGS -I.. -c interaction_matrices_order_1.cpp -o interaction_matrices_order_1.o;
$CXX $CXXFLAGS -I.. -c interaction_matrices_order_2.cpp -o interaction_matrices_order_2.o;
$CXX $CXXFLAGS -I.. -c interaction_matrices_order_3.cpp -o interaction_matrices_order_3.o;
$CXX $CXXFLAGS -I.. -c interaction_matrices_order_4.cpp -o interaction_matrices_order_4.o;
cd ..;

$AR $ARFLAGS libgeneral.a geometry/interval_quadrules.o geometry/cubical_quadrules.o geometry/sphere.o \
                          math/interaction_matrices_order_1.o math/interaction_matrices_order_2.o \
                          math/interaction_matrices_order_3.o math/interaction_matrices_order_4.o;
$RANLIB libgeneral.a;

$CXX $CXXFLAGS -I. general.cpp libgeneral.a -l openblas -o general;
strip -s general;

