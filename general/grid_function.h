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
#ifndef GRID_FUNCTION_H
#define GRID_FUNCTION_H

#include <armadillo>
#include <unordered_map>

#include <types.h>
#include <math/legendre.h>
#include <geometry/point.h>
#include <geometry/cellid.h>

template <size_t order>
struct grid_function
{
    real operator()( geometry::point x ) const;

    real h { 1 };
    std::unordered_map<geometry::cellid,arma::vec> cells_and_coeffs;
};

template <size_t order>
real grid_function<order>::operator()( geometry::point x ) const
{
    using std::floor;
    using geometry::point;
    using geometry::cellid;
    using math::legendre3d::N;

    cellid id   { static_cast<int>( floor(x.x/h) ),
                  static_cast<int>( floor(x.y/h) ),
                  static_cast<int>( floor(x.z/h) )  }; 
    point x_ref { (x.x - id.i*h)/h,
                  (x.y - id.j*h)/h,
                  (x.z - id.k*h)/h };

    auto it  = cells_and_coeffs.find( id );
    auto end = cells_and_coeffs.end();
    real fac = 1./std::sqrt( h*h*h );
    if ( it != end ) return fac*arma::dot( it->second, N<order>( x_ref ) );
    else             return 0;
}

#endif

