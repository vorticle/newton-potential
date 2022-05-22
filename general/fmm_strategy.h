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
 * this software; see the file COPYING. If not see http://www.gnu.org/licenses.
 */
#ifndef FMM_STRATEGY_H
#define FMM_STRATEGY_H

#include <armadillo>

#include <types.h>
#include <geometry/point.h>
#include <geometry/cellid.h>
#include <geometry/sphere.h>
#include <geometry/quadrules.h>
#include <geometry/cmplx_point.h>

#include <math/legendre.h>
#include <math/solid_harmonics.h>
#include <math/interaction_matrices.h>

#include <grid_function.h>

template <size_t expansion_order, size_t space_order>
class fmm_strategy
{
private:
    real         h;
    arma::cx_mat p2m_mat;
    arma::cx_mat l2p_mat, l2p_mat_conj;
    
    void precompute_p2m_mat();
    void precompute_l2p_mat();

public:
    using Mcoeff_t = math::solid_harmonics::cmplx_coeff<expansion_order>;
    using Lcoeff_t = math::solid_harmonics::cmplx_coeff<expansion_order>;
    using point    = geometry::point;
    using sphere   = geometry::sphere;

    // Iterator over the elements of a grid_function.
    using iterator = typename std::unordered_map<geometry::cellid,arma::vec>::iterator;

    fmm_strategy() = delete;
    fmm_strategy( real hh );

    sphere bounding_sphere( iterator it ) const;

    void p2p( iterator target, iterator source ) const;
    void p2m(       Mcoeff_t &M, point c, iterator source ) const;
    void l2p( const Lcoeff_t &L, point c, iterator target ) const;

    void m2m( Mcoeff_t &target, const Mcoeff_t &source, point c ) const;
    void m2l( Lcoeff_t &target, const Mcoeff_t &source, point c ) const;
    void l2l( Lcoeff_t &target, const Lcoeff_t &source, point c ) const;

    bool mac( sphere target, sphere source ) const;
};

#include <fmm_strategy.tpp>
#endif

