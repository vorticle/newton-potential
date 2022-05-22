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
#ifndef BOXPOT_ANTIDERIVATIVE_H
#define BOXPOT_ANTIDERIVATIVE_H

#include <array>
#include <cmath>
#include <limits>

// For some weird reason we need to include artanh, but not arctan.
#include <boost/math/special_functions/atanh.hpp>

#include <boxpot/polynomial.h>

namespace boxpot
{

template <typename number_type>
class antiderivative
{
public:
    using multi_index = boxpot::multi_index<6>;
    using polynomial  = boxpot::polynomial<number_type,6>;

    template <typename float_type> 
    float_type eval( const std::array<number_type,6> &pos ) const;

    template <typename float_type, typename output_iterator>
    void write_summands( const std::array<number_type,6> &pos,
                         output_iterator out ) const;

    antiderivative& operator*=( number_type rhs );
    antiderivative& operator/=( number_type rhs );

    antiderivative& operator*=( const polynomial &rhs );

    antiderivative& operator+=( const antiderivative &rhs );
    antiderivative& operator-=( const antiderivative &rhs );

    antiderivative operator*( number_type rhs ) const;
    antiderivative operator/( number_type rhs ) const;
    antiderivative operator*( const polynomial     &rhs ) const;
    antiderivative operator+(       antiderivative  rhs ) const;
    antiderivative operator-( const antiderivative &rhs ) const; 

    antiderivative x1_integrate() const;
    antiderivative x2_integrate() const;
    antiderivative x3_integrate() const;
    antiderivative y1_integrate() const;
    antiderivative y2_integrate() const;
    antiderivative y3_integrate() const;

    polynomial PR, PRinv,
               Partanh_XR, Partanh_YR, Partanh_ZR,
               Parctan_RZ, Parctan_RY, Parctan_RX;

private:
    struct impl;
};

template <typename number_type>
antiderivative<number_type> operator*( number_type lhs, antiderivative<number_type> rhs );

template <typename number_type>
antiderivative<number_type> operator*( const polynomial<number_type,6> &lhs, antiderivative<number_type> rhs );

template <typename number_type> antiderivative<number_type> x1_integrate( const antiderivative<number_type> &f );
template <typename number_type> antiderivative<number_type> x2_integrate( const antiderivative<number_type> &f );
template <typename number_type> antiderivative<number_type> x3_integrate( const antiderivative<number_type> &f );
template <typename number_type> antiderivative<number_type> y1_integrate( const antiderivative<number_type> &f );
template <typename number_type> antiderivative<number_type> y2_integrate( const antiderivative<number_type> &f );
template <typename number_type> antiderivative<number_type> y3_integrate( const antiderivative<number_type> &f );

}

#include <boxpot/antiderivative.tpp>
#endif

