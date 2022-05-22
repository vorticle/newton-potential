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
#ifndef GEOMETRY_QUADRULES_H
#define GEOMETRY_QUADRULES_H

#include <geometry/point.h>

#include <vector>
#include <utility>

namespace geometry
{

struct interval_quad_node
{
    real x, w;
};

struct cubical_quad_node
{
    point x;
    real  w;
};

using interval_quadrule = std::vector< interval_quad_node>;
using  cubical_quadrule = std::vector<  cubical_quad_node>;

const interval_quadrule  get_interval_quadrule( size_t degree );
const  cubical_quadrule   get_cubical_quadrule( size_t degree );

interval_quadrule get_refined_interval_quadrule( size_t degree, size_t refinements );
 cubical_quadrule get_refined_cubical_quadrule ( size_t degree, size_t refinements );

}

#endif

