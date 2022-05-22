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

#include <unordered_map>

#include <types.h>
#include <misc/stopwatch.h>
#include <math/legendre.h>
#include <geometry/point.h>
#include <geometry/sphere.h>
#include <geometry/cellid.h>
#include <geometry/quadrules.h>

#include <grid_function.h>
#include <fmm_strategy.h>

#include <fmm/fmm.h>

          real   h;       // Mesh-size. Set by main().
constexpr size_t n =  4;  // Order of approximation for f. May be 1, 2, 3, or 4.
constexpr size_t P = 20;  // Order of multipole expansions.

real f_true( geometry::point x )
{
    if ( x.r2() < 1 )
    {
        real v = 1 - x.r2();
        real fac = 2*(-2 + 6*v - v*v)/(v*v*v*v);
        return fac * std::exp(-1/v);
    }
    else return 0;
}

real u_true( geometry::point x )
{
    return ( x.r2() < 1 ) ? std::exp(-1./(1.-x.r2())) : 0;
}


grid_function<n  > compute_f();
grid_function<n+2> compute_u_grid( grid_function<n> &f );


real compute_f_error( const grid_function<n  > &f );
real compute_u_error( const grid_function<n+2> &u );

grid_function<n> compute_f()
{
    using geometry::point;
    using geometry::sphere;
    using geometry::cellid;

    grid_function<n> result; result.h = h;

    sphere support; support.radius = 1;
    auto intersects_support = [support]( geometry::cellid id ) -> bool
    {
        point corners[8] =
        {
            point { (id.i + 0)*h, (id.j + 0)*h, (id.k + 0)*h },
            point { (id.i + 1)*h, (id.j + 0)*h, (id.k + 0)*h },
            point { (id.i + 0)*h, (id.j + 1)*h, (id.k + 0)*h },
            point { (id.i + 1)*h, (id.j + 1)*h, (id.k + 0)*h },
            point { (id.i + 0)*h, (id.j + 0)*h, (id.k + 1)*h },
            point { (id.i + 1)*h, (id.j + 0)*h, (id.k + 1)*h },
            point { (id.i + 0)*h, (id.j + 1)*h, (id.k + 1)*h },
            point { (id.i + 1)*h, (id.j + 1)*h, (id.k + 1)*h }
        };

        for ( size_t i = 0; i < 8; ++i )
        {
            if ( support.contains( corners[i] ) )
                return true;
        }
        return false;
    };

    int min_idx = std::floor(real(-1)/h) - 1;
    int max_idx = std::floor(real( 1)/h) + 1;
    for ( int i = min_idx; i <= max_idx; ++i )
    for ( int j = min_idx; j <= max_idx; ++j )
    for ( int k = min_idx; k <= max_idx; ++k )
    {
        cellid id { i, j, k };
        if ( intersects_support(id) )
            result.cells_and_coeffs[ id ] =
                arma::vec( math::legendre3d::num<n>(), arma::fill::zeros );
    }

    // Now that we initialised the grid, we can compute the coefficients.
    real J = h*h*h;
    geometry::cubical_quadrule rule = geometry::get_cubical_quadrule( 2*(n-1) );

    #pragma omp parallel
    #pragma omp single
    for ( auto it  = result.cells_and_coeffs.begin();
               it != result.cells_and_coeffs.end(); ++it )
    {
        #pragma omp task
        {
            cellid id  = it->first;
            arma::vec  &coeffs = it->second;

            geometry::point xc { id.i*h, id.j*h, id.k*h };
            for ( auto node: rule )
            {
                point x_ref = node.x;
                point x     = xc + h*x_ref;
                real  w     = J*node.w;
                real  f     = f_true(x);
                auto  N     = math::legendre3d::N<n>(x_ref,h);

                coeffs += w*f*N;
            }
        }
    }

    return result;
}

// Computes the grid, but not the coefficients for u for a given f.
grid_function<n+2>  compute_u_grid( grid_function<n> &f )
{
    grid_function<n+2> result;
    
    result.h = f.h;
    for ( auto it  = f.cells_and_coeffs.begin();
               it != f.cells_and_coeffs.end(); ++it )
    {
        geometry::cellid id = it->first;
        for ( int k = -1; k <= 1; ++k )
        for ( int j = -1; j <= 1; ++j )
        for ( int i = -1; i <= 1; ++i )
        {
            geometry::cellid myid { id.i + i, id.j + j, id.k + k };
            result.cells_and_coeffs[ myid ] =
                arma::vec( math::legendre3d::num<n+2>(), arma::fill::zeros );
        }
    }

    return result;
}

real compute_f_error( const grid_function<n  > &f )
{
    using geometry::point;
    using geometry::cellid;
    real J    = h*h*h;
    auto rule = geometry::get_cubical_quadrule( 2*(n-1) );

    real err = 0;
    #pragma omp parallel
    #pragma omp single
    for ( auto it  = f.cells_and_coeffs.begin();
               it != f.cells_and_coeffs.end(); ++it )
    {
        #pragma omp task
        {
            real my_err = 0;
            for ( auto node: rule )
            {
                cellid id   = it->first;
                real  w     = node.w*J;
                point x_ref = node.x;
                point x     = point { (id.i+x_ref.x)*h, (id.j+x_ref.y)*h, (id.k+x_ref.z)*h };
                real fh = arma::dot( math::legendre3d::N<n>( x_ref, h ), it->second );

                real f  = f_true(x);
                my_err += w*(f-fh)*(f-fh);
            }

            #pragma omp critical
            err += my_err;
        }
    }

    return std::sqrt(err);
}

real compute_u_error( const grid_function<n+2> &u )
{
    using geometry::point;
    using geometry::cellid;
    real J    = h*h*h;
    auto rule = geometry::get_cubical_quadrule( 2*(n+1) );

    real err = 0;
    #pragma omp parallel
    #pragma omp single
    for ( auto it  = u.cells_and_coeffs.begin();
               it != u.cells_and_coeffs.end(); ++it )
    {
        #pragma omp task
        {
            arma::vec c( math::legendre3d::num<n+2>(), arma::fill::zeros );
            real my_err = 0;
            for ( auto node: rule )
            {
                cellid id   = it->first;
                real  w     = node.w*J;
                point x_ref = node.x;
                point x     = point { (id.i+x_ref.x)*h, (id.j+x_ref.y)*h, (id.k+x_ref.z)*h };
                real  u     = u_true(x);
                real uh = arma::dot( math::legendre3d::N<n+2>( x_ref, h ), it->second );

                my_err += w*(u-uh)*(u-uh);
            }

            #pragma omp critical
            err += my_err;
        }
    }

    return std::sqrt(err);
}


int main( int argc, char* argv[] )
{
    std::stringstream strstr( argv[1] );
    strstr >> h;

    misc::stopwatch clock;
    grid_function<n> fh = compute_f();
    real elapsed = clock.elapsed();
    std::cout << "Time for computing fh: " << elapsed << ". " << std::endl;
    std::cout << "NDOF fh: " << fh.cells_and_coeffs.size() * math::legendre3d::num<n>() << ". " << std::endl;
           

    std::cout << u8"L²-Error of fh: " << compute_f_error(fh) << ". " << std::endl;

    clock.reset();
    fmm_strategy<P,n> S(h);
    elapsed = clock.elapsed();
    std::cout << "Time for precomputing P2M and L2P matrices: " << elapsed << ". " << std::endl;

    clock.reset();
    grid_function<n+2> uh = compute_u_grid(fh);
    elapsed = clock.elapsed();
    std::cout << "Time for setting up grid for uh: " << elapsed << ". " << std::endl;
    std::cout << "NDOF uh: " << uh.cells_and_coeffs.size() * math::legendre3d::num<n+2>() << ". " << std::endl;

    clock.reset();
    fmm::fmm( S, fh.cells_and_coeffs.begin(), fh.cells_and_coeffs.end(),
                 uh.cells_and_coeffs.begin(), uh.cells_and_coeffs.end() );
    elapsed = clock.elapsed();
    std::cout << "Time for computing uh, using FMM: " << elapsed << ". " << std::endl;
    std::cout << u8"L²-Error of uh: " << compute_u_error(uh) << ". " << std::endl;
}

