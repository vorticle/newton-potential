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
#ifndef GEOMETRY_SPHERE_H
#define GEOMETRY_SPHERE_H

#include <geometry/point.h>

namespace geometry
{

struct sphere
{
    constexpr bool contains( point  p ) const noexcept;
    constexpr bool contains( sphere b ) const noexcept;

    point centre;
    real  radius;
};

sphere bounding_sphere( const point  *begin, const point  *end );
sphere bounding_sphere( const sphere *begin, const sphere *end );

// Useful special case: bounding sphere of triangle in space.
sphere bounding_sphere( point x0, point x1, point x2 );

real distance( sphere A, sphere B ) noexcept;
real distance( sphere A, point  p ) noexcept;
real distance( point  p, sphere A ) noexcept;

}

#include <geometry/sphere.tpp>
#endif

