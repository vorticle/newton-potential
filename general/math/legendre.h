/*
 * Copyright (C) 2020 Matthias Kirchhart
 *
 * This is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3, or (at your option) any later
 * version.
 *
 * This software is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this software; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */
#ifndef MATH_LEGENDRE_H
#define MATH_LEGENDRE_H

#include <armadillo>

#include <types.h>
#include <geometry/point.h>

namespace math
{

namespace legendre1d
{

// Computes the orthoNORMAL Legendre polynomials on [0,1].
template <size_t order>
arma::vec::fixed<order> N( real x );

}

namespace legendre3d
{

template <size_t order>
constexpr size_t num() { return (order*(2+3*order+order*order))/6; }
 
// Computes the orthoNORMAL Legendre polynomials on [0,1]³.
template <size_t order>
arma::vec N( geometry::point x );

// Legendre polynomials, scaled such that they are orthonormal
// on a grid of mesh-size h. x is stell expected to be in reference
// coordinates, x in [0,1]³.
template <size_t order>
arma::vec N( geometry::point x, real h );
}

}

#include <math/legendre.tpp>
#endif

