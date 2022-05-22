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
#include <geometry/quadrules.h>

namespace geometry
{

const cubical_quadrule get_cubical_quadrule( size_t degree )
{
    const interval_quadrule rule( get_interval_quadrule(degree) );
    const size_t N = rule.size();
    cubical_quadrule result( N*N*N );
    for ( size_t i = 0; i < N; ++i )
    for ( size_t j = 0; j < N; ++j )
    for ( size_t k = 0; k < N; ++k )
    {
        result[ i*N*N + j*N + k ]= { point { rule[i].x,   rule[j].x,   rule[k].x  },
                                     real  { rule[i].w  * rule[j].w  * rule[k].w  } };
    }
    return result;
}

cubical_quadrule get_refined_cubical_quadrule( size_t degree, size_t level )
{
    const interval_quadrule rule( get_refined_interval_quadrule(degree,level) );
    const size_t N = rule.size();
    cubical_quadrule result( N*N*N );
    for ( size_t i = 0; i < N; ++i )
    for ( size_t j = 0; j < N; ++j )
    for ( size_t k = 0; k < N; ++k )
    {
        result[ i*N*N + j*N + k ] = { point { rule[i].x,   rule[j].x,   rule[k].x  },
                                      real  { rule[i].w  * rule[j].w  * rule[k].w  } };
    }
    return result;
}

}

