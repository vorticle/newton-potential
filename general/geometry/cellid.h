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
#ifndef GEOMETRY_CELLID_H
#define GEOMETRY_CELLID_H

#include <cmath>
#include <functional>
#include <boost/functional/hash.hpp>

namespace geometry
{

/*!
 * \brief Cell ID. Struct for identifying cells on a Cartesian grid.
 */
struct cellid
{
public:
    int    i { 0 };
    int    j { 0 };
    int    k { 0 };

    constexpr bool operator< ( cellid rhs ) const noexcept;
    constexpr bool operator> ( cellid rhs ) const noexcept;
    constexpr bool operator<=( cellid rhs ) const noexcept;
    constexpr bool operator>=( cellid rhs ) const noexcept;
    constexpr bool operator==( cellid rhs ) const noexcept;
    constexpr bool operator!=( cellid rhs ) const noexcept;
};

}

#include <geometry/cellid.tpp>
#endif

