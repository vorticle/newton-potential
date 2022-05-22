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
 * this software; see the file COPYING.GPL.  If not see http://www.gnu.org/licenses.
 */
#ifndef LEGENDRE_POLY_H
#define LEGENDRE_POLY_H

#include <vector>
#include <boxpot/polynomial.h>

namespace legendre1d_poly
{

template <typename number_type>
using poly6 = boxpot::polynomial<number_type,6>;

// Computes the orthogonal Legendre polynomials on the unit interval [0,1].
template <size_t order, typename number_type>
std::array< poly6<number_type>, order >
N( poly6<number_type> x )
{
    x = 2*x - 1; // Transform [0,1] to [-1,1].
    
    std::array< poly6<number_type>, order > result;
    result[0] = 1;
    if ( order == 1 )
        return result;

    result[1] = x; 
    for ( size_t n = 1; n < order - 1; ++n )
    {
        result[n+1] = (number_type(2*n+1)/number_type(n+1))*x*result[n  ] -
                      (number_type(n)    /number_type(n+1))*  result[n-1];
    }

    return result;
}

template <size_t order, typename float_type>
std::array< float_type, order >
norm_N()
{
    std::array< float_type, order > result;

    // result[ n ] = 1 / sqrt(2n+1);
    for ( size_t n = 0; n < order; ++n )
    {
        result[ n ] = float_type(1) /
                      sqrt( float_type(2)*float_type(n) + float_type(1) );
    }
    return result;
}

}

namespace legendre3d_poly
{

template <typename number_type>
using poly6 = boxpot::polynomial<number_type,6>;

template <size_t order, typename number_type>
std::vector< poly6<number_type> >
Nx( number_type x, number_type y, number_type z )
{
    poly6<number_type> x1 { 1, { 1, 0, 0, 0, 0, 0 } };
    poly6<number_type> x2 { 1, { 0, 1, 0, 0, 0, 0 } };
    poly6<number_type> x3 { 1, { 0, 0, 1, 0, 0, 0 } };
    std::array< poly6<number_type>, order > Nx { legendre1d_poly::N<order,number_type>(x1-x) };
    std::array< poly6<number_type>, order > Ny { legendre1d_poly::N<order,number_type>(x2-y) };
    std::array< poly6<number_type>, order > Nz { legendre1d_poly::N<order,number_type>(x3-z) };
  
    size_t count = 0 ;
    std::vector<poly6<number_type>> result( (order*(2 + 3*order + order*order))/6 );
    for ( size_t k = 0; k < order;         ++k )
    for ( size_t j = 0; j < order - k;     ++j )
    for ( size_t i = 0; i < order - k - j; ++i )
        result[ count++ ] = Nz[k]*Ny[j]*Nx[i];

    return result;
}

template <size_t order, typename number_type>
std::vector< poly6<number_type> >
Ny( number_type x, number_type y, number_type z )
{
    poly6<number_type> y1 { 1, { 0, 0, 0, 1, 0, 0 } };
    poly6<number_type> y2 { 1, { 0, 0, 0, 0, 1, 0 } };
    poly6<number_type> y3 { 1, { 0, 0, 0, 0, 0, 1 } };
    std::array< poly6<number_type>, order > Nx { legendre1d_poly::N<order,number_type>(y1-x) };
    std::array< poly6<number_type>, order > Ny { legendre1d_poly::N<order,number_type>(y2-y) };
    std::array< poly6<number_type>, order > Nz { legendre1d_poly::N<order,number_type>(y3-z) };
 
    size_t count = 0; 
    std::vector<poly6<number_type>> result( (order*(2 + 3*order + order*order))/6 );
    for ( size_t k = 0; k < order;         ++k )
    for ( size_t j = 0; j < order - k;     ++j )
    for ( size_t i = 0; i < order - k - j; ++i )
        result[ count++ ] = Nz[k]*Ny[j]*Nx[i];

    return result;
}

template <size_t order, typename float_type>
std::vector< float_type >
norm_N()
{
    std::array< float_type, order > norms = legendre1d_poly::norm_N<order,float_type>();

    size_t count = 0; 
    std::vector< float_type > result( (order*(2+3*order + order*order))/6 );
    for ( size_t k = 0; k < order;         ++k )
    for ( size_t j = 0; j < order - k;     ++j )
    for ( size_t i = 0; i < order - k - j; ++i )
        result[ count++ ] = norms[k]*norms[j]*norms[i];

    return result;
}

}

#endif

