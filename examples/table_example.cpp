/*
 * Copyright (C) 2020 Matthias Kirchhart and Donat Weniger
 *
 * This file is part of boxpot.
 * boxpot is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3, or (at your option) any later
 * version.
 *
 * boxpot is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with boxpot; see the files COPYING.LGPL and COPYING.GPL.  If not
 * see http://www.gnu.org/licenses.
 */
#include <cmath>
#include <iomanip>
#include <iostream>
#include <boost/multiprecision/gmp.hpp>
#include <boxpot/polynomial.h>
#include <boxpot/antiderivative.h>


using namespace boxpot;

// Cartesian id, used to identify the cell on the Cartesian grid.
struct carid { int i, j, k; };

// Implementation of the orthonormal Legendre polynomials.
namespace legendre
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

template <typename number_type, typename float_type>
float_type evaluate_antiderivative( carid id, const antiderivative<number_type> &f ) 
{
    auto y3int = [&f,&id](number_type x1, number_type x2, number_type x3,
                          number_type y1, number_type y2 ) -> float_type
    {
        std::array<number_type,6> upper { x1, x2, x3, y1, y2, id.k + 1 };
        std::array<number_type,6> lower { x1, x2, x3, y1, y2, id.k     };
    	return (f.template eval<float_type>(upper)) -
               (f.template eval<float_type>(lower));
    };

    auto y2int = [&y3int,&id]( number_type x1, number_type x2, number_type x3,
                               number_type y1 ) -> float_type
    {
      	return y3int( x1, x2, x3, y1, id.j + 1 ) - y3int( x1, x2, x3, y1, id.j );
    };

    auto y1int = [&y2int,&id]( number_type x1, number_type x2, number_type x3 ) -> float_type
    {
       	return y2int( x1, x2, x3, id.i + 1 ) - y2int( x1, x2, x3, id.i );
    };

    auto x3int = [&y1int]( number_type x1, number_type x2 ) -> float_type
    {
       	return y1int( x1, x2, 1 ) - y1int( x1, x2, 0 );
    };

    auto x2int = [&x3int]( number_type x1 ) -> float_type
    {
       	return x3int( x1, 1 ) - x3int( x1, 0 );
    };

    auto x1int = [&x2int]() -> float_type
    {
       	return x2int(1) - x2int(0);
    };

    return x1int();
}

int main()
{
    using rational = boost::multiprecision::mpq_rational;
    using real     = boost::multiprecision::mpf_float_500;
    using idx6     = multi_index<6>;
    using poly6    = polynomial<rational,6>;

    // Short-hands for convenience.
    const poly6 x1 = poly6( rational(1), idx6 { 1, 0, 0, 0, 0, 0 } );
    const poly6 x2 = poly6( rational(1), idx6 { 0, 1, 0, 0, 0, 0 } );
    const poly6 x3 = poly6( rational(1), idx6 { 0, 0, 1, 0, 0, 0 } );
    const poly6 y1 = poly6( rational(1), idx6 { 0, 0, 0, 1, 0, 0 } );
    const poly6 y2 = poly6( rational(1), idx6 { 0, 0, 0, 0, 1, 0 } );
    const poly6 y3 = poly6( rational(1), idx6 { 0, 0, 0, 0, 0, 1 } );

    constexpr size_t max_order = 7;	
    constexpr carid           j   {  -1, -1, -1 };
    constexpr multi_index<6>  idx {  1, 1, 1, 1, 1, 1 };

    poly6 Px1 = legendre::N<max_order>( x1       )[idx[0]];
    poly6 Px2 = legendre::N<max_order>( x2       )[idx[1]];
    poly6 Px3 = legendre::N<max_order>( x3       )[idx[2]];
    poly6 Py1 = legendre::N<max_order>( y1 - j.i )[idx[3]];
    poly6 Py2 = legendre::N<max_order>( y2 - j.j )[idx[4]];
    poly6 Py3 = legendre::N<max_order>( y3 - j.k )[idx[5]];

    poly6 P = Px1*Px2*Px3*Py1*Py2*Py3;
    
    real result = 0;
    size_t count = 0;
    for ( auto term: P.terms )
    {
	std::cout << "Integrating monomial: " << ++count << " "
                  << "out of: " << P.terms.size() << ". " << std::endl;
        antiderivative<rational> F;
        F.PRinv = poly6 { 1, term.first };

        F = y1_integrate(F);
        F = y2_integrate(F);
        F = y3_integrate(F);
        F = x1_integrate(F);
        F = x2_integrate(F);
        F = x3_integrate(F);

        result += static_cast<real>( term.second ) *
                  evaluate_antiderivative<rational,real>( j, F );
    }

    auto norms = legendre::norm_N<max_order,real>();
    result /= norms[ idx[0] ];
    result /= norms[ idx[1] ];
    result /= norms[ idx[2] ];
    result /= norms[ idx[3] ];
    result /= norms[ idx[4] ];
    result /= norms[ idx[5] ];

    std::cout << std::setprecision(17) << result << std::endl;
}

