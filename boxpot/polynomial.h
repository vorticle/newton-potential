/*
 * Copyright (C) 2020 Matthias Kirchhart and Donat Weniger
 *
 * This file is part of boxpot.
 * boxpot is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3, or (at your option) any later
 * version.
 *
 * boxpot is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with boxpot; see the files COPYING.LGPL and COPYING.GPL.  If not
 * see http://www.gnu.org/licenses.
 */
#ifndef BOXPOT_POLYNOMIAL_H
#define BOXPOT_POLYNOMIAL_H

#include <array>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <boost/functional/hash.hpp> // For hash_combine.

#include <boxpot/rump_summation.h>

namespace boxpot
{

template <size_t dim>
using multi_index = std::array<size_t,dim>;

template <size_t dim>
bool operator==( const multi_index<dim> &lhs, const multi_index<dim> &rhs ) noexcept;

template <typename number_type, size_t dim>
class polynomial
{
public:
    using multi_index = boxpot::multi_index<dim>;

    struct index_hash
    {
        std::hash<size_t> hasher;
 
        size_t operator()( const multi_index &idx ) const noexcept
        {
            using boost::hash_combine;
            size_t result { 0 };
            for ( size_t i = 0; i < dim; ++i )
                hash_combine( result, hasher( idx[i] ) );

            return result;
        }
    };

    polynomial() = default;
    polynomial( const polynomial  &rhs ) = default;
    polynomial(       polynomial &&rhs ) = default;
    polynomial& operator=( const polynomial  &rhs ) = default;
    polynomial& operator=(       polynomial &&rhs ) = default;

    explicit polynomial( number_type coefficient, const multi_index &exponents = multi_index {} );
    explicit polynomial( const multi_index &exponents, number_type coefficient = 1 );

    number_type operator()( const std::array<number_type,dim> &pos ) const;

    polynomial& operator= ( number_type rhs );

    polynomial& operator+=( number_type rhs );
    polynomial& operator-=( number_type rhs );
    polynomial& operator*=( number_type rhs );
    polynomial& operator/=( number_type rhs );

    polynomial  operator+( number_type rhs ) const;
    polynomial  operator-( number_type rhs ) const;
    polynomial  operator*( number_type rhs ) const;
    polynomial  operator/( number_type rhs ) const;

    friend polynomial operator+( number_type lhs,       polynomial  rhs ) { rhs += lhs; return rhs; }
    friend polynomial operator*( number_type lhs,       polynomial  rhs ) { rhs *= lhs; return rhs; }
    friend polynomial operator-( number_type lhs, const polynomial &rhs )
    { polynomial result = -rhs; result += lhs; return result; }

    polynomial operator-() const;

    polynomial operator+(       polynomial  rhs ) const;
    polynomial operator-( const polynomial &rhs ) const;
    polynomial operator*( const polynomial &rhs ) const;

    polynomial& operator+=( const polynomial &rhs );
    polynomial& operator-=( const polynomial &rhs );
    polynomial& operator*=( const polynomial &rhs );

    bool operator!=(       polynomial  rhs ) const;
    bool operator==( const polynomial &rhs ) const;

    multi_index max_degrees() const noexcept;

    void compress();

    // No need to keep this private, no invariant to enforce.
    std::unordered_map<multi_index,number_type,index_hash> terms;
};

}

#include <boxpot/polynomial.tpp>
#endif

