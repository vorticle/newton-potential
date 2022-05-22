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

namespace boxpot
{

template <typename number_type>
number_type npow( number_type x, size_t n )
{
    // Exponentiation by squaring.
    if ( n == 0 ) return 1;

    number_type y { 1 };
    while ( n > 1 )
    {
        if ( n & 1 ) // if ( is_odd(n) )
        {
            y *= x;
            x *= x;
            n  ^= 1; // n = n - 1;
            n >>= 1; // n = n / 2;
        }
        else
        {
            x *= x;
            n >>= 1; // n = n / 2;
        }
    }
    return x*y;
}

template <size_t dim> inline
bool operator==( const multi_index<dim> &lhs, const multi_index<dim> &rhs ) noexcept
{
    for ( size_t i = 0; i < dim; ++i )
        if ( lhs[i] != rhs[i] ) return false;
    return true;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>::polynomial( number_type coefficient, const multi_index &exponents )
{
    terms[ exponents ] = coefficient;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>::polynomial( const multi_index &exponents, number_type coefficient )
{
    terms[ exponents ] = coefficient;
}

template <typename number_type, size_t dim, bool is_floating_point>
struct polynomial_evaluator;

template <typename number_type, size_t dim>
struct polynomial_evaluator<number_type,dim,false>
{
    using poly = polynomial<number_type,dim>;

    static number_type eval( const poly &P, const std::array<number_type,dim> &pos )
    {
        number_type result {};
        for ( const auto &term: P.terms )
        {
            number_type summand { term.second };
            for ( size_t i = 0; i < dim; ++i )
                summand *= npow( pos[i], term.first[i] );
            result += summand;
        }

        return result;
    }
};

template <typename number_type, size_t dim>
struct polynomial_evaluator<number_type,dim,true>
{
    using poly = polynomial<number_type,dim>;

    static number_type eval( const poly &P, const std::array<number_type,dim> &pos )
    {
        std::vector<number_type> summands;
        for ( const auto &term: P.terms )
        {
            number_type summand { term.second };
            for ( size_t i = 0; i < dim; ++i )
                summand *= npow( pos[i], term.first[i] );
            summands.push_back(summand);
        }

        return rump_summation( summands.begin(), summands.end() );
    }
};

template <typename number_type, size_t dim> 
number_type polynomial<number_type,dim>::operator()( const std::array<number_type,dim> &pos ) const
{
    using eval_t = polynomial_evaluator< number_type, dim, std::is_floating_point<number_type>::value >;
    return eval_t::eval( *this, pos );
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>& polynomial<number_type,dim>::operator=( number_type rhs )
{
    terms.clear();
    terms[ multi_index() ] = rhs;
    return *this;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>& polynomial<number_type,dim>::operator+=( number_type rhs )
{
    terms[ multi_index() ] += rhs;
    return *this;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>& polynomial<number_type,dim>::operator-=( number_type rhs )
{
    terms[ multi_index() ] -= rhs;
    return *this;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>& polynomial<number_type,dim>::operator*=( number_type rhs )
{
    if ( rhs != 0 )
    {
        for ( auto &term: terms )
            term.second *= rhs;
    }
    else 
    {
        terms.clear();
    }
    return *this;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>& polynomial<number_type,dim>::operator/=( number_type rhs )
{
    for ( auto &term: terms )
        term.second /= rhs;
    return *this;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim> polynomial<number_type,dim>::operator+( number_type rhs ) const
{
    return polynomial(*this) += rhs;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim> polynomial<number_type,dim>::operator-( number_type rhs ) const
{
    return polynomial(*this) -= rhs;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim> polynomial<number_type,dim>::operator*( number_type rhs ) const
{
    return polynomial(*this) *= rhs;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim> polynomial<number_type,dim>::operator/( number_type rhs ) const
{
    return polynomial(*this) /= rhs;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim> polynomial<number_type,dim>::operator-() const
{
    polynomial result { *this };
    for ( auto &term: result.terms )
        term.second = -term.second;
    return result;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>& polynomial<number_type,dim>::operator+=( const polynomial &rhs )
{
    // Safe even if &rhs == this, no new entries are created.
    for ( const auto &term: rhs.terms )
    {
        terms[ term.first ] += term.second;
    }
    compress();
    return *this;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>& polynomial<number_type,dim>::operator-=( const polynomial &rhs )
{
    // Safe even if &rhs == this, no new entries are created.
    for ( const auto &term: rhs.terms )
    {
        terms[ term.first ] -= term.second;
    }
    compress();
    return *this;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim>& polynomial<number_type,dim>::operator*=( const polynomial &rhs )
{
    return (*this) = (*this) * rhs; 
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim> polynomial<number_type,dim>::operator+( polynomial rhs ) const
{
    rhs += *this;
    return rhs;
}

template <typename number_type, size_t dim> inline
polynomial<number_type,dim> polynomial<number_type,dim>::operator-( const polynomial &rhs ) const
{
    polynomial result { *this }; result -= rhs;
    return result;
}

template <typename number_type, size_t dim>
polynomial<number_type,dim> polynomial<number_type,dim>::operator*( const polynomial &rhs ) const
{
    polynomial result;
    for ( const auto &lterm: this->terms )
    {
        for ( const auto &rterm: rhs.terms )
        {
            multi_index exponent;
            for ( size_t i = 0; i < dim; ++i )
                exponent[i] = lterm.first[i] + rterm.first[i];

            number_type coefficient { lterm.second * rterm.second };
            result.terms[ exponent ] += coefficient;
        }
    }
    result.compress();
    return result;
}

template <typename number_type, size_t dim> inline
bool polynomial<number_type,dim>::operator!=( polynomial rhs ) const
{
    rhs -= *this;
    for ( const auto &term: rhs.terms )
        if ( term.second != 0 )
            return true;
    return false;
}

template <typename number_type, size_t dim> inline
bool polynomial<number_type,dim>::operator==( const polynomial &rhs ) const
{
    return ! ( *this != rhs );
}

template <typename number_type, size_t dim> 
multi_index<dim> polynomial<number_type,dim>::max_degrees() const noexcept
{
    multi_index result {};
    for ( const auto &term: terms )
    {
        for ( size_t i = 0; i < dim; ++i )
            result[i] = std::max(result[i],term.first[i]);
    }
    return result;
}

template <typename number_type, size_t dim> inline
void polynomial<number_type,dim>::compress()
{
    for ( auto it = terms.begin(); it != terms.end(); )
        if ( it->second == 0 ) it = terms.erase(it);
        else ++it;
}


}

