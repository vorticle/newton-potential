/*
 * Copyright (C) 2020 Matthias Kirchhart and Donat Weniger
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

#include <cmath>
#include <iomanip>
#include <iostream>

#include <armadillo>

#include <boost/multiprecision/gmp.hpp>
#include <boost/math/constants/constants.hpp>

#include <boxpot/polynomial.h>
#include <boxpot/antiderivative.h>

#include "legendre_poly.h"


using namespace boxpot;

struct carid { int i, j, k; };

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

    float_type four_pi = static_cast<float_type>(4) * boost::math::constants::pi<float_type>();
    return x1int() / four_pi;
}

template <size_t n, typename number_type>
void compute_antiderivatives( antiderivative<number_type> F[n+2][n+2][n+2][n][n][n] )
{
    #pragma omp parallel
    {
        #pragma omp single
        for ( size_t i3 = 0; i3 < n+2          ; ++i3 )
        for ( size_t i2 = 0; i2 < n+2 - i3     ; ++i2 )
        for ( size_t i1 = 0; i1 < n+2 - i3 - i2; ++i1 )
        for ( size_t j3 = 0; j3 < n            ; ++j3 )
        for ( size_t j2 = 0; j2 < n   - j3     ; ++j2 )
        for ( size_t j1 = 0; j1 < n   - j3 - j2; ++j1 )
        {
            #pragma omp task
            {
                antiderivative<number_type> &f = F[i1][i2][i3][j1][j2][j3];
                f.PRinv = polynomial<number_type,6> { 1, { i1, i2, i3, j1, j2, j3 } };
                f = x1_integrate(f); f = x2_integrate(f); f = x3_integrate(f);
                f = y1_integrate(f); f = y2_integrate(f); f = y3_integrate(f);
            }
        }
    }
}

template <size_t n, typename number_type, typename float_type>
arma::mat compute_matrix( carid id, const antiderivative<number_type> F[n+2][n+2][n+2][n][n][n] )
{
    float_type ints[n+2][n+2][n+2][n][n][n];

    #pragma omp parallel
    {
        #pragma omp single
        for ( size_t i3 = 0; i3 < n+2          ; ++i3 )
        for ( size_t i2 = 0; i2 < n+2 - i3     ; ++i2 )
        for ( size_t i1 = 0; i1 < n+2 - i3 - i2; ++i1 )
        for ( size_t j3 = 0; j3 < n            ; ++j3 )
        for ( size_t j2 = 0; j2 < n   - j3     ; ++j2 )
        for ( size_t j1 = 0; j1 < n   - j3 - j2; ++j1 )
        {
            #pragma omp task
            {
                const antiderivative<number_type> &f = F[i1][i2][i3][j1][j2][j3];
                ints[i1][i2][i3][j1][j2][j3] = evaluate_antiderivative<number_type,float_type>( id, f );
            }
        }
    }

    auto x = legendre3d_poly::Nx<n+2,number_type>( 0,    0,    0    );
    auto y = legendre3d_poly::Ny<n  ,number_type>( id.i, id.j, id.k );
    auto norms_x = legendre3d_poly::norm_N<n+2,float_type>();
    auto norms_y = legendre3d_poly::norm_N<n  ,float_type>();

    arma::mat A( x.size(), y.size(), arma::fill::zeros  );

    #pragma omp parallel
    {
        #pragma omp single
        for ( size_t i = 0; i < x.size(); ++i )
        for ( size_t j = 0; j < y.size(); ++j )
        {
            #pragma omp task
            {
                auto xy = x[i]*y[j];
                float_type result = 0;
                for ( auto term: xy.terms )
                {
                    multi_index<6> idx = term.first;
                    float_type coeff = static_cast<float_type>( term.second );
                    float_type v     = ints[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]];
                    result += coeff*v;
                }
                result /= norms_x[i]*norms_y[j];

                A(i,j) = static_cast<double>(result);
                if ( std::abs(A(i,j)) < 1e-32 )
                    A(i,j) = 0;
            }
        }
    }

    return A;
}

void print_matrix( std::ostream &str, const arma::mat &A )
{
    str << std::setprecision(17) << std::scientific;
    str << "arma::mat\n{{\n";
    for ( size_t i = 0; i < A.n_rows; ++i )
    {
        str << "  { ";
        for ( size_t j = 0; j < A.n_cols; ++j )
        {
            if ( j != A.n_cols - 1 )
            {
                str << std::setw(24) << A(i,j) << ", ";
            }
            else
            {
                str << std::setw(24) << A(i,j);
            }
        }

        if ( i != A.n_rows - 1 )
        {
            str << " },\n";
        }
        else
        {
            str << " }\n";
        }
    }
    str << "}},\n\n";
    str.flush();
}

int main()
{
    constexpr size_t n = 4; // order.
    constexpr int thickness = 1;
    constexpr int N = 2*thickness + 1;
    using rational = boost::multiprecision::mpq_rational;
    using real     = boost::multiprecision::mpf_float_100;

    antiderivative<rational> F[n+2][n+2][n+2][n][n][n];
    compute_antiderivatives<n,rational>(F);

    std::cout << "#include <armadillo>\n\n\n";
    std::cout << "namespace\n{\n\n";
    std::cout << "const arma::mat A[ " << N*N*N << " ]\n{\n\n";

    std::cout.flush();

    arma::mat interaction_matrices[N][N][N];
    carid id;
    for ( int k = 0; k < N; ++k )
    for ( int j = 0; j < N; ++j )
    for ( int i = 0; i < N; ++i )
    {
        id.i = i-thickness;
        id.j = j-thickness;
        id.k = k-thickness;
        interaction_matrices[i][j][k] = compute_matrix<n,rational,real>(id,F);
        print_matrix( std::cout, interaction_matrices[i][j][k] );
    }

    std::cout << "};\n\n}\n\n";
}

