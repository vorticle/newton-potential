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
#ifndef BOXPOT_ANTIDERIVATIVE2D_H
#define BOXPOT_ANTIDERIVATIVE2D_H

#include <array>
#include <cmath>
#include <limits>

#include <boxpot/polynomial.h>

namespace boxpot 
{

template <typename number_type>
class antiderivative2d
{
public:
    using multi_index = boxpot::multi_index<4>;
    using polynomial  = boxpot::polynomial<number_type,4>;

    template <typename float_type> 
    float_type eval( const std::array<number_type,4> &pos ) const;

    template <typename float_type, typename output_iterator>
    void write_summands( const std::array<number_type,4> &pos,
                         output_iterator out ) const;

    antiderivative2d& operator*=( number_type rhs );
    antiderivative2d& operator/=( number_type rhs );

    antiderivative2d& operator*=( const polynomial &rhs );

    antiderivative2d& operator+=( const antiderivative2d &rhs );
    antiderivative2d& operator-=( const antiderivative2d &rhs );

    antiderivative2d operator*( number_type rhs ) const;
    antiderivative2d operator/( number_type rhs ) const;
    antiderivative2d operator*( const polynomial     &rhs ) const;
    antiderivative2d operator+(       antiderivative2d  rhs ) const;
    antiderivative2d operator-( const antiderivative2d &rhs ) const; 

    antiderivative2d x1_integrate() const;
    antiderivative2d x2_integrate() const;
    antiderivative2d y1_integrate() const;
    antiderivative2d y2_integrate() const;

    polynomial P, Plog_R2, Parctan_XY, Parctan_YX;

private:
    struct impl;
};

template <typename number_type>
antiderivative2d<number_type> operator*( number_type lhs, antiderivative2d<number_type> rhs );

template <typename number_type>
antiderivative2d<number_type> operator*( const polynomial<number_type,4> &lhs, antiderivative2d<number_type> rhs );

template <typename number_type> antiderivative2d<number_type> x1_integrate( const antiderivative2d<number_type> &f );
template <typename number_type> antiderivative2d<number_type> x2_integrate( const antiderivative2d<number_type> &f );
template <typename number_type> antiderivative2d<number_type> y1_integrate( const antiderivative2d<number_type> &f );
template <typename number_type> antiderivative2d<number_type> y2_integrate( const antiderivative2d<number_type> &f );

}

#include <boxpot/antiderivative2d.tpp>
#endif

