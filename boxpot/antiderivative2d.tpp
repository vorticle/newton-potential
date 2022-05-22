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

// Hide all the implementation details in a private struct. Preferred
// over a separate impl-namespace, because we do not need to repeat the
// template arguments all the time.
template <typename number_type>
struct antiderivative2d<number_type>::impl
{

    static number_type rational( number_type numerator, number_type denominator ) { return numerator / denominator; }

    template <typename float_type> static float_type log_R2   ( const std::array<number_type,4> &pos );
    template <typename float_type> static float_type arctan_XY( const std::array<number_type,4> &pos );
    template <typename float_type> static float_type arctan_YX( const std::array<number_type,4> &pos );

    static antiderivative2d x1_integrate_plain( const polynomial &P );
    static antiderivative2d x2_integrate_plain( const polynomial &P );
    static antiderivative2d y1_integrate_plain( const polynomial &P );
    static antiderivative2d y2_integrate_plain( const polynomial &P );

    static antiderivative2d x1_integrate_log_R2( const polynomial &P );
    static antiderivative2d x2_integrate_log_R2( const polynomial &P );
    static antiderivative2d y1_integrate_log_R2( const polynomial &P );
    static antiderivative2d y2_integrate_log_R2( const polynomial &P );

    static antiderivative2d x1_integrate_arctan_XY( const polynomial &P );
    static antiderivative2d x2_integrate_arctan_XY( const polynomial &P );
    static antiderivative2d y1_integrate_arctan_XY( const polynomial &P );
    static antiderivative2d y2_integrate_arctan_XY( const polynomial &P );

    static antiderivative2d x1_integrate_arctan_YX( const polynomial &P );
    static antiderivative2d x2_integrate_arctan_YX( const polynomial &P );
    static antiderivative2d y1_integrate_arctan_YX( const polynomial &P );
    static antiderivative2d y2_integrate_arctan_YX( const polynomial &P );

    static const polynomial x1; 
    static const polynomial x2; 
    static const polynomial y1; 
    static const polynomial y2; 

    static const polynomial X; 
    static const polynomial Y; 
};


// Static data members of antiderivative::impl. Quite tedious, any way to
// improve this?

template <typename number_type> const polynomial<number_type,4> antiderivative2d<number_type>::impl::x1 { multi_index {1,0,0,0} };
template <typename number_type> const polynomial<number_type,4> antiderivative2d<number_type>::impl::x2 { multi_index {0,1,0,0} };
template <typename number_type> const polynomial<number_type,4> antiderivative2d<number_type>::impl::y1 { multi_index {0,0,1,0} };
template <typename number_type> const polynomial<number_type,4> antiderivative2d<number_type>::impl::y2 { multi_index {0,0,0,1} };

template <typename number_type> const polynomial<number_type,4> antiderivative2d<number_type>::impl::X
{ antiderivative2d<number_type>::impl::x1 - antiderivative2d<number_type>::impl::y1 };

template <typename number_type> const polynomial<number_type,4> antiderivative2d<number_type>::impl::Y
{ antiderivative2d<number_type>::impl::x2 - antiderivative2d<number_type>::impl::y2 };



/////////////////////////////////////////////////
//                                             //
// Implemantation of antiderivative2d<> class. //
//                                             //
/////////////////////////////////////////////////

template <typename number_type> inline
antiderivative2d<number_type>& antiderivative2d<number_type>::operator*=( number_type rhs )
{
    P          *= rhs;
    Plog_R2    *= rhs;
    Parctan_XY *= rhs;
    Parctan_YX *= rhs;
    return *this;
}

template <typename number_type> inline
antiderivative2d<number_type>& antiderivative2d<number_type>::operator*=( const polynomial &rhs )
{
    P          *= rhs;
    Plog_R2    *= rhs;
    Parctan_XY *= rhs;
    Parctan_YX *= rhs;
    return *this;
}

template <typename number_type> inline
antiderivative2d<number_type>& antiderivative2d<number_type>::operator/=( number_type rhs )
{
    P          /= rhs;
    Plog_R2    /= rhs;
    Parctan_XY /= rhs;
    Parctan_YX /= rhs;
    return *this;
}

template <typename number_type> inline
antiderivative2d<number_type>& antiderivative2d<number_type>::operator+=( const antiderivative2d &rhs )
{
    P          += rhs.P;
    Plog_R2    += rhs.Plog_R2;
    Parctan_XY += rhs.Parctan_XY;
    Parctan_YX += rhs.Parctan_YX;
    return *this;
}

template <typename number_type> inline
antiderivative2d<number_type>& antiderivative2d<number_type>::operator-=( const antiderivative2d &rhs )
{
    P          -= rhs.P;
    Plog_R2    -= rhs.Plog_R2;
    Parctan_XY -= rhs.Parctan_XY;
    Parctan_YX -= rhs.Parctan_YX;
    return *this;
}

template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::operator*( number_type rhs ) const
{
    antiderivative2d result {*this};
    result *= rhs;
    return result; 
}

template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::operator/( number_type rhs ) const
{
    antiderivative2d result {*this};
    result /= rhs;
    return result; 
}

template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::operator*( const polynomial &rhs ) const
{
    antiderivative2d result {*this};
    result *= rhs;
    return result; 
}

template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::operator+( antiderivative2d rhs ) const
{
    rhs += *this;
    return rhs;
}

template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::operator-( const antiderivative2d &rhs ) const
{
    antiderivative2d result {*this};
    result -= rhs;
    return result;
}

template <typename number_type> inline
antiderivative2d<number_type> operator*( number_type lhs, antiderivative2d<number_type> rhs )
{
    return rhs *= lhs;
}

template <typename number_type> inline
antiderivative2d<number_type> operator*( const polynomial<number_type,4> &lhs, antiderivative2d<number_type> rhs )
{
    return rhs *= lhs;
}


template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::x1_integrate() const
{
    antiderivative2d result;
    result += impl::x1_integrate_plain( P );
    result += impl::x1_integrate_log_R2( Plog_R2 );
    result += impl::x1_integrate_arctan_XY( Parctan_XY );
    result += impl::x1_integrate_arctan_YX( Parctan_YX );
    return result;
}

template <typename number_type> inline
antiderivative2d<number_type> x1_integrate( const antiderivative2d<number_type> &f )
{
    return f.x1_integrate();
}

template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::x2_integrate() const
{
    antiderivative2d result;
    result += impl::x2_integrate_plain( P );
    result += impl::x2_integrate_log_R2( Plog_R2 );
    result += impl::x2_integrate_arctan_XY( Parctan_XY );
    result += impl::x2_integrate_arctan_YX( Parctan_YX );
    return result;
}

template <typename number_type> inline
antiderivative2d<number_type> x2_integrate( const antiderivative2d<number_type> &f )
{
    return f.x2_integrate();
}


template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::y1_integrate() const
{
    antiderivative2d result;
    result += impl::y1_integrate_plain( P );
    result += impl::y1_integrate_log_R2( Plog_R2 );
    result += impl::y1_integrate_arctan_XY( Parctan_XY );
    result += impl::y1_integrate_arctan_YX( Parctan_YX );
    return result;
}

template <typename number_type> inline
antiderivative2d<number_type> y1_integrate( const antiderivative2d<number_type> &f )
{
    return f.y1_integrate();
}

template <typename number_type> inline
antiderivative2d<number_type> antiderivative2d<number_type>::y2_integrate() const
{
    antiderivative2d result;
    result += impl::y2_integrate_plain( P );
    result += impl::y2_integrate_log_R2( Plog_R2 );
    result += impl::y2_integrate_arctan_XY( Parctan_XY );
    result += impl::y2_integrate_arctan_YX( Parctan_YX );
    return result;
}

template <typename number_type> inline
antiderivative2d<number_type> y2_integrate( const antiderivative2d<number_type> &f )
{
    return f.y2_integrate();
}


template <typename number_type> template <typename float_type> inline
float_type antiderivative2d<number_type>::eval( const std::array<number_type,4> &x ) const
{
    return static_cast<float_type>(         P(x)) +
           static_cast<float_type>(   Plog_R2(x)) * impl::template    log_R2<float_type>(x) +
           static_cast<float_type>(Parctan_XY(x)) * impl::template arctan_XY<float_type>(x) +
           static_cast<float_type>(Parctan_YX(x)) * impl::template arctan_YX<float_type>(x);
}

template <typename number_type>
template <typename float_type, typename output_iterator> inline
void antiderivative2d<number_type>::write_summands( const std::array<number_type,4> &x,
                                                  output_iterator out ) const
{
    *out++ = static_cast<float_type>(         P(x));
    *out++ = static_cast<float_type>(   Plog_R2(x)) * impl::template    log_R2<float_type>(x);
    *out++ = static_cast<float_type>(Parctan_XY(x)) * impl::template arctan_XY<float_type>(x);
    *out++ = static_cast<float_type>(Parctan_YX(x)) * impl::template arctan_YX<float_type>(x);
}



///////////////////////////////////////////////
//                                           //
// Implementation of antiderivative<>::impl  //
//                                           //
///////////////////////////////////////////////

template <typename number_type> template <typename float_type> inline
float_type antiderivative2d<number_type>::impl::log_R2( const std::array<number_type,4> &pos )
{
    number_type arg {0};
    arg += (pos[0]-pos[2])*(pos[0]-pos[2]);
    arg += (pos[1]-pos[3])*(pos[1]-pos[3]);

    if ( arg == 0 ) return 0;
    return log( static_cast<float_type>(arg) );
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative2d<number_type>::impl::arctan_XY( const std::array<number_type,4> &pos )
{
    float_type x { static_cast<float_type>( pos[0] - pos[2] ) };
    float_type y { static_cast<float_type>( pos[1] - pos[3] ) };

    if ( y == 0 ) return 0;
    return atan(x/y);
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative2d<number_type>::impl::arctan_YX( const std::array<number_type,4> &pos )
{
    float_type x { static_cast<float_type>( pos[0] - pos[2] ) };
    float_type y { static_cast<float_type>( pos[1] - pos[3] ) };

    if ( x == 0 ) return 0;
    return atan(y/x);
}

//////////////////////////////////////
// Integration of plain polynomials //
//////////////////////////////////////

template <typename number_type> auto antiderivative2d<number_type>::impl::x1_integrate_plain( const polynomial &P )
-> antiderivative2d
{
    antiderivative2d result;
    for ( auto term: P.terms )
    {
        multi_index index = term.first;
        number_type coeff = term.second;

        index[0]++; coeff /= index[0];
        result.P.terms[index] = coeff;
    } 

    return result;
}

template <typename number_type> auto antiderivative2d<number_type>::impl::x2_integrate_plain( const polynomial &P )
-> antiderivative2d
{
    antiderivative2d result;
    for ( auto term: P.terms )
    {
        auto        index = term.first;
        number_type coeff = term.second;

        index[1]++; coeff /= index[1];
        result.P.terms[index] = coeff;
    } 

    return result;
}

template <typename number_type> auto antiderivative2d<number_type>::impl::y1_integrate_plain( const polynomial &P )
-> antiderivative2d
{
    antiderivative2d result;
    for ( auto term: P.terms )
    {
        multi_index index = term.first;
        number_type coeff = term.second;

        index[2]++; coeff /= index[2];
        result.P.terms[index] = coeff;
    } 

    return result;
}

template <typename number_type> auto antiderivative2d<number_type>::impl::y2_integrate_plain( const polynomial &P )
-> antiderivative2d
{
    antiderivative2d result;
    for ( auto term: P.terms )
    {
        multi_index index = term.first;
        number_type coeff = term.second;

        index[3]++; coeff /= index[3];
        result.P.terms[index] = coeff;
    } 

    return result;
}

/////////////////////////////////////////////
// Integration of the function log(X²+Y²). //
/////////////////////////////////////////////

template <typename number_type> 
auto antiderivative2d<number_type>::impl::x1_integrate_log_R2( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.P          = -2*X;
            f.Parctan_XY =  2*Y;
            f.Plog_R2    =    X;
            break;
        case 1:
            f.P          = -(X*X + Y*Y)/2;
            f.Plog_R2    =  (X*X + Y*Y)/2;
            f += y1*cache[0];
            break;
        default:
            polynomial x1n  { multi_index {n  ,0,0,0} };
            polynomial x1n1 { multi_index {n-1,0,0,0} };
            f.P          = - x1n*x1 * rational( 2, (n+1)*(n+1) );
            f.Plog_R2    =   x1n*X  * rational( 1, n+1 );
            f.Parctan_XY =   x1n*Y  * rational( 2, n+1 );
            f += rational(n,n+1)*y1 * cache[n-1];
            f -= rational(2*n,n+1)*Y * x1_integrate_arctan_XY( x1n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[0];
        exponents[0] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative2d<number_type>::impl::x2_integrate_log_R2( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.P          = -2*Y;
            f.Parctan_YX =  2*X;
            f.Plog_R2    =    Y;
            break;
        case 1:
            f.P          = -(X*X + Y*Y)/2;
            f.Plog_R2    =  (X*X + Y*Y)/2;
            f += y2*cache[0];
            break;
        default:
            polynomial x2n  { multi_index {0,n  ,0,0} };
            polynomial x2n1 { multi_index {0,n-1,0,0} };
            f.P          = - x2n*x2 * rational( 2, (n+1)*(n+1) );
            f.Plog_R2    =   x2n*Y  * rational( 1, n+1 );
            f.Parctan_YX =   x2n*X  * rational( 2, n+1 );
            f += rational(n,n+1)*y2 * cache[n-1];
            f -= rational(2*n,n+1)*X * x2_integrate_arctan_YX( x2n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[1];
        exponents[1] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}


template <typename number_type> 
auto antiderivative2d<number_type>::impl::y1_integrate_log_R2( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.P          =  2*X;
            f.Parctan_XY = -2*Y;
            f.Plog_R2    =   -X;
            break;
        case 1:
            f.P          = -(X*X + Y*Y)/2;
            f.Plog_R2    =  (X*X + Y*Y)/2;
            f += x1*cache[0];
            break;
        default:
            polynomial y1n  { multi_index {0,0,n  ,0} };
            polynomial y1n1 { multi_index {0,0,n-1,0} };
            f.P          = - y1n*y1 * rational( 2, (n+1)*(n+1) );
            f.Plog_R2    = - y1n*X  * rational( 1, n+1 );
            f.Parctan_XY = - y1n*Y  * rational( 2, n+1 );
            f += rational(n,n+1)*x1 * cache[n-1];
            f += rational(2*n,n+1)*Y * y1_integrate_arctan_XY( y1n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[2];
        exponents[2] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative2d<number_type>::impl::y2_integrate_log_R2( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.P          =  2*Y;
            f.Parctan_YX = -2*X;
            f.Plog_R2    =   -Y;
            break;
        case 1:
            f.P          = -(X*X + Y*Y)/2;
            f.Plog_R2    =  (X*X + Y*Y)/2;
            f += x2*cache[0];
            break;
        default:
            polynomial y2n  { multi_index {0,0,0,n  } };
            polynomial y2n1 { multi_index {0,0,0,n-1} };
            f.P          = - y2n*y2 * rational(2, (n+1)*(n+1) );
            f.Plog_R2    = - y2n*Y  * rational(1, n+1 );
            f.Parctan_YX = - y2n*X  * rational(2, n+1 );
            f += rational(n,n+1)*x2 * cache[n-1];
            f += rational(2*n,n+1)*X * y2_integrate_arctan_YX( y2n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[3];
        exponents[3] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////////////
// Integration of the function arctan(X/Y). //
//////////////////////////////////////////////

template <typename number_type> 
auto antiderivative2d<number_type>::impl::x1_integrate_arctan_XY( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.Parctan_XY =  X; 
            f.Plog_R2    = -Y/2;
            break;
        case 1:
            f.P          = -X*Y/2;
            f.Parctan_XY =  (X*X + Y*Y)/2;
            f += y1*cache[0];
            break;
        default:
            polynomial x1n  { multi_index {n  ,0,0,0} };
            polynomial x1n1 { multi_index {n-1,0,0,0} };
            f.Plog_R2    = - x1n*Y  * rational(1, 2*(n+1) );
            f.Parctan_XY =   x1n*X  * rational(1, n+1 );
            f += rational(n,n+1)*y1 * cache[n-1];
            f += rational(n,2*(n+1))*Y * x1_integrate_log_R2( x1n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[0];
        exponents[0] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative2d<number_type>::impl::x2_integrate_arctan_XY( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.Parctan_XY =  Y; 
            f.Plog_R2    =  X/2;
            break;
        case 1:
            f.P          =  X*Y/2;
            f.Parctan_XY =  Y*Y/2;
            f.Parctan_YX = -X*X/2;
            f += y2*cache[0];
            break;
        default:
            polynomial x2n  { multi_index {0,n  ,0,0} };
            polynomial x2n1 { multi_index {0,n-1,0,0} };
            f.Plog_R2    =  x2n*X  * rational(  1, 2*(n+1) );
            f.Parctan_XY =  x2n*Y  * rational(  1, n+1 );
            f += rational(n,n+1)*y2 * cache[n-1];
            f -= rational(n,2*(n+1))*X * x2_integrate_log_R2( x2n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[1];
        exponents[1] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative2d<number_type>::impl::y1_integrate_arctan_XY( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.Parctan_XY = -X; 
            f.Plog_R2    =  Y/2;
            break;
        case 1:
            f.P          = -X*Y/2;
            f.Parctan_XY =  (X*X + Y*Y)/2;
            f += x1*cache[0];
            break;
        default:
            polynomial y1n  { multi_index {0,0,n  ,0} };
            polynomial y1n1 { multi_index {0,0,n-1,0} };
            f.Plog_R2    =   y1n*Y  * rational(1, 2*(n+1) );
            f.Parctan_XY = - y1n*X  * rational(1, n+1 );
            f += rational(n,n+1)*x1 * cache[n-1];
            f -= rational(n,2*(n+1))*Y * y1_integrate_log_R2( y1n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[2];
        exponents[2] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative2d<number_type>::impl::y2_integrate_arctan_XY( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.Parctan_XY = -Y; 
            f.Plog_R2    = -X/2;
            break;
        case 1:
            f.P          =  X*Y/2;
            f.Parctan_XY =  Y*Y/2;
            f.Parctan_YX = -X*X/2;
            f += x2*cache[0];
            break;
        default:
            polynomial y2n  { multi_index {0,0,0,n  } };
            polynomial y2n1 { multi_index {0,0,0,n-1} };
            f.Plog_R2    = - y2n*X  * rational(1, 2*(n+1) );
            f.Parctan_XY = - y2n*Y  * rational(1, n+1 );
            f += rational(n,n+1)*x2 * cache[n-1];
            f += rational(n,2*(n+1))*X * y2_integrate_log_R2( y2n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[3];
        exponents[3] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////////////
// Integration of the function arctan(Y/X). //
//////////////////////////////////////////////

template <typename number_type> 
auto antiderivative2d<number_type>::impl::x1_integrate_arctan_YX( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.Parctan_YX =  X; 
            f.Plog_R2    =  Y/2;
            break;
        case 1:
            f.P          =  X*Y/2;
            f.Parctan_YX =  X*X/2;
            f.Parctan_XY = -Y*Y/2;
            f += y1*cache[0];
            break;
        default:
            polynomial x1n  { multi_index {n  ,0,0,0} };
            polynomial x1n1 { multi_index {n-1,0,0,0} };
            f.Plog_R2    =  x1n*Y  * rational(1, 2*(n+1) );
            f.Parctan_YX =  x1n*X  * rational(1, n+1 );
            f += rational(n,n+1)*y1 * cache[n-1];
            f -= rational(n,2*(n+1))*Y * x1_integrate_log_R2( x1n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[0];
        exponents[0] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative2d<number_type>::impl::x2_integrate_arctan_YX( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.Parctan_YX =  Y; 
            f.Plog_R2    = -X/2;
            break;
        case 1:
            f.P          = -X*Y/2;
            f.Parctan_YX = (X*X + Y*Y)/2;
            f += y2*cache[0];
            break;
        default:
            polynomial x2n  { multi_index {0,n  ,0,0} };
            polynomial x2n1 { multi_index {0,n-1,0,0} };
            f.Plog_R2    = - x2n*X  * rational(1, 2*(n+1) );
            f.Parctan_YX =   x2n*Y  * rational(1, n+1 );
            f += rational(n,n+1)*y2 * cache[n-1];
            f += rational(n,2*(n+1))*X * x2_integrate_log_R2( x2n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[1];
        exponents[1] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative2d<number_type>::impl::y1_integrate_arctan_YX( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.Parctan_YX = -X; 
            f.Plog_R2    = -Y/2;
            break;
        case 1:
            f.P          =  X*Y/2;
            f.Parctan_YX =  X*X/2;
            f.Parctan_XY = -Y*Y/2;
            f += x1*cache[0];
            break;
        default:
            polynomial y1n  { multi_index {0,0,n  ,0} };
            polynomial y1n1 { multi_index {0,0,n-1,0} };
            f.Plog_R2    = - y1n*Y  * rational(1, 2*(n+1) );
            f.Parctan_YX = - y1n*X  * rational(1, n+1 );
            f += rational(n,n+1)*x1 * cache[n-1];
            f += rational(n,2*(n+1))*Y * y1_integrate_log_R2( y1n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[2];
        exponents[2] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative2d<number_type>::impl::y2_integrate_arctan_YX( const polynomial &P )
-> antiderivative2d
{
    thread_local std::vector<antiderivative2d> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative2d f;

        switch ( n )
        {
        case 0:
            f.Parctan_YX = -Y; 
            f.Plog_R2    =  X/2;
            break;
        case 1:
            f.P          = -X*Y/2;
            f.Parctan_YX = (X*X + Y*Y)/2;
            f += x2*cache[0];
            break;
        default:
            polynomial y2n  { multi_index {0,0,0,n  } };
            polynomial y2n1 { multi_index {0,0,0,n-1} };
            f.Plog_R2    =   y2n*X  * rational(1,2*(n+1));
            f.Parctan_YX = - y2n*Y  * rational(1,n+1);
            f += rational(n,n+1)*x2 * cache[n-1];
            f -= rational(n,2*(n+1))*X * y2_integrate_log_R2( y2n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative2d result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[3];
        exponents[3] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

}

