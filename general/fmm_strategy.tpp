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

template <size_t expansion_order, size_t space_order>
void fmm_strategy<expansion_order,space_order>::
p2p( iterator target, iterator source ) const
{
    using cellid = geometry::cellid;
    cellid source_id = source->first;
    cellid target_id = target->first;
    cellid delta { source_id.i - target_id.i,
                   source_id.j - target_id.j,
                   source_id.k - target_id.k };

    const arma::mat& M = math::interaction_matrix<space_order>( delta );
    target->second += (h*h)*M*(source->second);
}

template <size_t expansion_order, size_t space_order>
void fmm_strategy<expansion_order,space_order>::
p2m( Mcoeff_t &M, point, iterator source ) const
{
    constexpr size_t n_rows = math::solid_harmonics::cmplx_coeff<expansion_order>::num_coeff();
    constexpr size_t n_cols = math::legendre3d::num<space_order>();
    arma::cx_vec Nc  ( n_cols, arma::fill::zeros ); Nc.set_real( source->second );
    arma::cx_vec Mvec( M.data.data(), n_rows, false, true );
    Mvec += p2m_mat*Nc;
}

template <size_t expansion_order, size_t space_order>
void fmm_strategy<expansion_order,space_order>::
l2p( const Lcoeff_t &L, point, iterator target ) const
{
    constexpr size_t n_cols = math::solid_harmonics::cmplx_coeff<expansion_order>::num_coeff();
    const arma::cx_vec Lvec( const_cast<cmplx*>(L.data.data()), n_cols, false, true );
    target->second += arma::real(l2p_mat*Lvec + l2p_mat_conj*conj(Lvec));
}

template <size_t expansion_order, size_t space_order>
void fmm_strategy<expansion_order,space_order>::
m2m( Mcoeff_t &target, const Mcoeff_t &source, point r ) const
{
    using math::solid_harmonics::R;
    using math::solid_harmonics::idx;

    auto transfer = R<expansion_order>( -r );

    constexpr int order = expansion_order;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            for ( int nd = 0; nd <= n; ++nd )
            {
                int min_md = std::max(-nd, m - (n-nd));
                int max_md = std::min( nd, m + (n-nd));
                for ( int md = min_md; md <= max_md; ++md )
                {
                    target.data[ idx(n,m) ] += source( nd, md )*conj(transfer(n-nd,m-md));
                }
            }
        }
    }
}

template <size_t expansion_order, size_t space_order>
void fmm_strategy<expansion_order,space_order>::
m2l( Lcoeff_t &target, const Mcoeff_t &source, point r ) const
{
    using math::solid_harmonics::S;
    using math::solid_harmonics::idx;

    const auto transfer = S<expansion_order>( r );
    constexpr int order = expansion_order;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            cmplx tmp;
            for ( int nd = 0; nd < order - n; ++nd )
            {
                for ( int md = -nd; md <= nd; ++md )
                {
                    tmp += source(nd,md)*transfer(n+nd,m+md);
                }
            }
            if ( n % 2 ) target.data[ idx(n,m) ] -= tmp;
            else         target.data[ idx(n,m) ] += tmp;
        }
    }
}

template <size_t expansion_order, size_t space_order>
void fmm_strategy<expansion_order,space_order>::
l2l( Lcoeff_t &target, const Lcoeff_t &source, point r ) const
{
    using math::solid_harmonics::R;
    using math::solid_harmonics::idx;
    const auto transfer = R<expansion_order>( r );

    constexpr int order = expansion_order;
    for ( int n = 0; n < order; ++n )
    {
        for ( int m = 0; m <= n; ++m )
        {
            for ( int nd = n; nd < order; ++nd )
            {
                int min_md = std::max(-nd, m - (nd-n));
                int max_md = std::min( nd, m + (nd-n));
                for ( int md = min_md; md <= max_md; ++md )
                {
                    target.data[ idx(n,m) ] += source( nd,md )*conj(transfer(nd-n,md-m));
                }
            }
        }
    }
}

template <size_t expansion_order, size_t space_order>
bool fmm_strategy<expansion_order,space_order>::
mac( sphere target, sphere source ) const
{
    // Multipole acceptance criterion, MAC.
    // If you chose theta_crit < sqrt(3)/2 you will need
    // additional interaction matrices, as the near-field
    // is enlarged to not only directly neighbouring cells.
    constexpr real theta_crit = 0.87;

    real theta = (source.radius + target.radius)/
                 (source.centre - target.centre).r();
    return theta < theta_crit;
}

template <size_t expansion_order, size_t space_order>
void fmm_strategy<expansion_order,space_order>::
precompute_p2m_mat()
{
    constexpr real pifac = 0.079577471545947667884441881; // 1/(4pi)
    const point cref { 0.5, 0.5, 0.5 };
    const real  J    { h*h*h };
    const auto rule = geometry::get_cubical_quadrule( expansion_order + space_order - 2 );

    arma::vec    N ( math::legendre3d::num<space_order>() );
    arma::cx_vec Nc( math::legendre3d::num<space_order>(), arma::fill::zeros );
    Mcoeff_t  R;

    for ( auto node: rule )
    {
        real  w = J*node.w;
        N = math::legendre3d::N<space_order>( node.x, h );
        R = math::solid_harmonics::R<expansion_order>( h*(node.x - cref) );
        Nc.set_real(N);

        for ( int i = 0; i < R.num_coeff(); ++i )
        {
            p2m_mat.row(i) += (w*pifac*conj(R.data[i]))*Nc.t();
        }
    }
}

template <size_t expansion_order, size_t space_order>
void fmm_strategy<expansion_order,space_order>::
precompute_l2p_mat()
{
    using math::solid_harmonics::idx;
    const point cref { 0.5, 0.5, 0.5 };
    const real  J    { h*h*h };
    const auto rule = geometry::get_cubical_quadrule( expansion_order + space_order );

    arma::vec    N ( math::legendre3d::num<space_order+2>(), arma::fill::zeros );
    arma::cx_vec Nc( math::legendre3d::num<space_order+2>(), arma::fill::zeros );
    Mcoeff_t R;

    l2p_mat.zeros(); l2p_mat_conj.zeros();
    for ( auto node: rule )
    {
        real w = J*node.w;
        N = math::legendre3d::N<space_order+2>( node.x, h );
        R = math::solid_harmonics::R<expansion_order>( h*(node.x - cref) );
        Nc.zeros();
        Nc.set_real(N);

        constexpr int order { expansion_order };
        for ( int n =  0; n <  order; ++n )
        for ( int m = -n; m <= n;     ++m )
        {
            if ( m < 0 )
            {
                if ( m & 1 ) l2p_mat_conj.col( idx(n,-m) ) -= (w*conj(R(n,m)))*Nc;
                else         l2p_mat_conj.col( idx(n,-m) ) += (w*conj(R(n,m)))*Nc;
            }
            else
            {
                l2p_mat.col( idx(n,m) ) += (w*conj(R(n,m)))*Nc;
            }
        }
    }
}

template <size_t expansion_order, size_t space_order>
fmm_strategy<expansion_order,space_order>::
fmm_strategy( real hh ):
h { hh },
p2m_mat      { math::solid_harmonics::cmplx_coeff<expansion_order>::num_coeff(),
               math::legendre3d::num<space_order>(),
               arma::fill::zeros },
l2p_mat      { math::legendre3d::num<space_order+2>(),
               math::solid_harmonics::cmplx_coeff<expansion_order>::num_coeff(),
               arma::fill::zeros },
l2p_mat_conj { math::legendre3d::num<space_order+2>(),
               math::solid_harmonics::cmplx_coeff<expansion_order>::num_coeff(),
               arma::fill::zeros }
{
    precompute_p2m_mat();
    precompute_l2p_mat();
}

template <size_t expansion_order, size_t space_order>
geometry::sphere fmm_strategy<expansion_order,space_order>::
bounding_sphere( iterator it ) const
{
    geometry::sphere result;
    geometry::cellid id = it->first;
    result.centre.x = (static_cast<real>(id.i) + 0.5)*h;
    result.centre.y = (static_cast<real>(id.j) + 0.5)*h;
    result.centre.z = (static_cast<real>(id.k) + 0.5)*h;
    result.radius   = h*(std::sqrt(real(3))/2);
    return result;
}

