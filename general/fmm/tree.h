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
#ifndef FMM_TREE_H
#define FMM_TREE_H

#include <geometry/point.h>
#include <geometry/sphere.h>

#include <cmath>
#include <memory>
#include <algorithm>

namespace fmm
{

template <typename coeff_t, typename iterator>
struct tree
{
    tree() = delete;
    tree( const tree&  ) = delete;
    tree(       tree&& ) = default;
    tree& operator=( const tree&  ) = delete;
    tree& operator=(       tree&& ) = default;
   ~tree() = default;

    template <typename strategy>
    tree( iterator begin, iterator end, strategy &S );

    struct box;
    std::vector<iterator> elems;
    std::unique_ptr<box>  root;
};

template <typename coeff_t, typename iterator>
struct tree<coeff_t,iterator>::box
{
    using point  = geometry::point;
    using sphere = geometry::sphere;

    box() = delete;
    box( const box&  ) = delete;
    box(       box&& ) = delete;
    box& operator=( const box&  ) = delete;
    box& operator=(       box&& ) = delete;


    template <typename strategy> box( iterator *first, iterator *last, strategy &S );

    template <typename strategy> void compute_bounds( strategy &S );
    template <typename strategy> void split( strategy &S );

    bool is_leaf() const;

    sphere bounds;

    iterator *begin, *end;
    coeff_t   coeffs {};

    std::unique_ptr<box> children[ 8 ];
};



template <typename coeff_t, typename iterator>
template <typename strategy>
tree<coeff_t,iterator>::tree( iterator begin, iterator end,
                              strategy &S ):
elems (std::distance(begin,end)), root { nullptr }
{
    for ( size_t i = 0; i < elems.size(); ++i )
        elems[ i ] = begin++;

    root.reset( new box( elems.data(), elems.data() + elems.size(), S ) );
}

template <typename coeff_t, typename iterator>
template <typename strategy>
void tree<coeff_t,iterator>::box::split( strategy &S )
{
    if ( end - begin == 1 ) return;

    auto xless = [&S,this]( iterator p ){ return S.bounding_sphere(p).centre.x < bounds.centre.x; };
    auto yless = [&S,this]( iterator p ){ return S.bounding_sphere(p).centre.y < bounds.centre.y; };
    auto zless = [&S,this]( iterator p ){ return S.bounding_sphere(p).centre.z < bounds.centre.z; };

    iterator* oct[9];
    oct[0] = begin;
    oct[8] = end;

    using std::partition;
    oct[4] = partition( oct[0], oct[8], xless );

    oct[2] = partition( oct[0], oct[4], yless );
    oct[1] = partition( oct[0], oct[2], zless );
    oct[3] = partition( oct[2], oct[4], zless );

    oct[6] = partition( oct[4], oct[8], yless );
    oct[5] = partition( oct[4], oct[6], zless );
    oct[7] = partition( oct[6], oct[8], zless );

    for ( uint i = 0; i < 8; ++i )
    {
        if ( oct[i] != oct[i+1] )
        {
            #pragma omp task shared(S)
            children[ i ].reset( new box( oct[i], oct[i+1], S ) );
        }
    }
    #pragma omp taskwait
}


template <typename coeff_t, typename iterator>
template <typename strategy>
void tree<coeff_t,iterator>::box::compute_bounds( strategy &S )
{
    using std::min;
    using std::max;
    using geometry::point;

    std::vector< geometry::sphere > v;
    for ( iterator* p = begin; p != end; ++p )
        v.push_back( S.bounding_sphere(*p) );

    bounds = bounding_sphere(v.data(),v.data() + v.size());
}

template <typename coeff_t, typename iterator>
template <typename strategy>
tree<coeff_t,iterator>::box::box( iterator *first, iterator *last, strategy &S ):
begin {first}, end{last}
{
    compute_bounds( S );
    split ( S );
}

template <typename coeff_t, typename iterator>
bool tree<coeff_t,iterator>::box::is_leaf() const
{
    for ( const auto& v: children )
    {
        if ( v ) return false;
    }
    return true;
}

}

#endif

