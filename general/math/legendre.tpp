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

namespace math
{

namespace legendre1d
{

template <>
arma::vec::fixed<1> N<1>( real )
{
    return arma::vec::fixed<1> {{ 1.0 }};
}

template <size_t order>
arma::vec::fixed<order> N( real x )
{
    x = 2*x - 1; // Transform [0,1] to [-1,1].
    
    arma::vec::fixed<order> result;
    result(0) = 1; result(1) = x; 
    for ( size_t n = 1; n < order - 1; ++n )
    {
        result(n+1) = (real(2*n+1)/real(n+1))*x*result(n  ) -
                      (real(n)    /real(n+1))*  result(n-1);
    }

    for ( size_t n = 0; n < order; ++n )
        result( n ) *= std::sqrt( real(2*n + 1) );

    return result;
}

}

namespace legendre3d
{

template <size_t order>
arma::vec N( geometry::point x )
{
    arma::vec::fixed<order> Nx { legendre1d::N<order>(x.x) };
    arma::vec::fixed<order> Ny { legendre1d::N<order>(x.y) };
    arma::vec::fixed<order> Nz { legendre1d::N<order>(x.z) };
  
    size_t count = 0; 
    arma::vec result( num<order>() ); 
    for ( size_t k = 0; k < order;         ++k )
    for ( size_t j = 0; j < order - k;     ++j )
    for ( size_t i = 0; i < order - k - j; ++i )
        result( count++ ) = Nz(k)*Ny(j)*Nx(i);

    return result;
}

template <size_t order>
arma::vec N( geometry::point x, real h )
{
    return N<order>(x)/std::sqrt(h*h*h);
}

}

}

