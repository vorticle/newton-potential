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

namespace geometry
{

constexpr
bool cellid::operator< ( cellid rhs ) const noexcept
{
    return   (k!=rhs.k) ? (k<rhs.k) :
           ( (j!=rhs.j) ? (j<rhs.j) : (i<rhs.i) );
}

constexpr
bool cellid::operator> ( cellid rhs ) const noexcept
{
    return   (k!=rhs.k) ? (k>rhs.k) :
           ( (j!=rhs.j) ? (j>rhs.j) : (i>rhs.i) );
}

constexpr
bool cellid::operator<=( cellid rhs ) const noexcept
{
    return ! ( (*this) > rhs );
}

constexpr
bool cellid::operator>=( cellid rhs ) const noexcept
{
    return ! ( (*this) < rhs );
}

constexpr
bool cellid::operator==( cellid rhs ) const noexcept
{
    return i == rhs.i && j == rhs.j && k == rhs.k;
}

constexpr
bool cellid::operator!=( cellid rhs ) const noexcept
{
    return i != rhs.i || j != rhs.j || k != rhs.k;
}

}

namespace std
{

template<>
struct hash<geometry::cellid>
{
    size_t operator()( geometry::cellid id ) const noexcept
    {
        using boost::hash_combine;
        std::hash<int>  coord_hash;

        size_t result { 0 };
        hash_combine( result, coord_hash( id.i ) );
        hash_combine( result, coord_hash( id.j ) );
        hash_combine( result, coord_hash( id.k ) );

        return result;
    } 
};

}

