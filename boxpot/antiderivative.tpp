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
struct antiderivative<number_type>::impl
{

    static number_type rational( number_type numerator, number_type denominator ) { return numerator / denominator; }

    template <typename float_type> static float_type R        ( const std::array<number_type,6> &pos );
    template <typename float_type> static float_type Rinv     ( const std::array<number_type,6> &pos );
    template <typename float_type> static float_type artanh_XR( const std::array<number_type,6> &pos );
    template <typename float_type> static float_type artanh_YR( const std::array<number_type,6> &pos );
    template <typename float_type> static float_type artanh_ZR( const std::array<number_type,6> &pos );
    template <typename float_type> static float_type arctan_RZ( const std::array<number_type,6> &pos );
    template <typename float_type> static float_type arctan_RY( const std::array<number_type,6> &pos );
    template <typename float_type> static float_type arctan_RX( const std::array<number_type,6> &pos );

    static antiderivative x1_integrate_R( const polynomial &P );
    static antiderivative x2_integrate_R( const polynomial &P );
    static antiderivative x3_integrate_R( const polynomial &P );
    static antiderivative y1_integrate_R( const polynomial &P );
    static antiderivative y2_integrate_R( const polynomial &P );
    static antiderivative y3_integrate_R( const polynomial &P );

    static antiderivative x1_integrate_Rinv( const polynomial &P );
    static antiderivative x2_integrate_Rinv( const polynomial &P );
    static antiderivative x3_integrate_Rinv( const polynomial &P );
    static antiderivative y1_integrate_Rinv( const polynomial &P );
    static antiderivative y2_integrate_Rinv( const polynomial &P );
    static antiderivative y3_integrate_Rinv( const polynomial &P );

    static antiderivative x1_integrate_artanh_XR( const polynomial &P );
    static antiderivative x2_integrate_artanh_XR( const polynomial &P );
    static antiderivative x3_integrate_artanh_XR( const polynomial &P );
    static antiderivative y1_integrate_artanh_XR( const polynomial &P );
    static antiderivative y2_integrate_artanh_XR( const polynomial &P );
    static antiderivative y3_integrate_artanh_XR( const polynomial &P );

    static antiderivative x1_integrate_artanh_YR( const polynomial &P );
    static antiderivative x2_integrate_artanh_YR( const polynomial &P );
    static antiderivative x3_integrate_artanh_YR( const polynomial &P );
    static antiderivative y1_integrate_artanh_YR( const polynomial &P );
    static antiderivative y2_integrate_artanh_YR( const polynomial &P );
    static antiderivative y3_integrate_artanh_YR( const polynomial &P );

    static antiderivative x1_integrate_artanh_ZR( const polynomial &P );
    static antiderivative x2_integrate_artanh_ZR( const polynomial &P );
    static antiderivative x3_integrate_artanh_ZR( const polynomial &P );
    static antiderivative y1_integrate_artanh_ZR( const polynomial &P );
    static antiderivative y2_integrate_artanh_ZR( const polynomial &P );
    static antiderivative y3_integrate_artanh_ZR( const polynomial &P );

    static antiderivative x1_integrate_arctan_RZ( const polynomial &P );
    static antiderivative x2_integrate_arctan_RZ( const polynomial &P );
    static antiderivative x3_integrate_arctan_RZ( const polynomial &P );
    static antiderivative y1_integrate_arctan_RZ( const polynomial &P );
    static antiderivative y2_integrate_arctan_RZ( const polynomial &P );
    static antiderivative y3_integrate_arctan_RZ( const polynomial &P );

    static antiderivative x1_integrate_arctan_RY( const polynomial &P );
    static antiderivative x2_integrate_arctan_RY( const polynomial &P );
    static antiderivative x3_integrate_arctan_RY( const polynomial &P );
    static antiderivative y1_integrate_arctan_RY( const polynomial &P );
    static antiderivative y2_integrate_arctan_RY( const polynomial &P );
    static antiderivative y3_integrate_arctan_RY( const polynomial &P );

    static antiderivative x1_integrate_arctan_RX( const polynomial &P );
    static antiderivative x2_integrate_arctan_RX( const polynomial &P );
    static antiderivative x3_integrate_arctan_RX( const polynomial &P );
    static antiderivative y1_integrate_arctan_RX( const polynomial &P );
    static antiderivative y2_integrate_arctan_RX( const polynomial &P );
    static antiderivative y3_integrate_arctan_RX( const polynomial &P );

    static const polynomial x1; 
    static const polynomial x2; 
    static const polynomial x3; 
    static const polynomial y1; 
    static const polynomial y2; 
    static const polynomial y3; 

    static const polynomial X; 
    static const polynomial Y; 
    static const polynomial Z; 

    static const polynomial Xi;
    static const polynomial Upsilon;
    static const polynomial Theta;
};


// Static data members of antiderivative::impl. Quite tedious, any way to
// improve this?

template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::x1 { multi_index {1,0,0,0,0,0} };
template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::x2 { multi_index {0,1,0,0,0,0} };
template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::x3 { multi_index {0,0,1,0,0,0} };
template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::y1 { multi_index {0,0,0,1,0,0} };
template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::y2 { multi_index {0,0,0,0,1,0} };
template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::y3 { multi_index {0,0,0,0,0,1} };

template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::X
{ antiderivative<number_type>::impl::x1 - antiderivative<number_type>::impl::y1 };

template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::Y
{ antiderivative<number_type>::impl::x2 - antiderivative<number_type>::impl::y2 };

template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::Z
{ antiderivative<number_type>::impl::x3 - antiderivative<number_type>::impl::y3 };

template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::Xi
{ antiderivative<number_type>::impl::Y * antiderivative<number_type>::impl::Y  +
  antiderivative<number_type>::impl::Z * antiderivative<number_type>::impl::Z };

template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::Upsilon
{ antiderivative<number_type>::impl::X * antiderivative<number_type>::impl::X  +
  antiderivative<number_type>::impl::Z * antiderivative<number_type>::impl::Z };

template <typename number_type> const polynomial<number_type,6> antiderivative<number_type>::impl::Theta
{ antiderivative<number_type>::impl::X * antiderivative<number_type>::impl::X  +
  antiderivative<number_type>::impl::Y * antiderivative<number_type>::impl::Y };




///////////////////////////////////////////////
//                                           //
// Implemantation of antiderivative<> class. //
//                                           //
///////////////////////////////////////////////

template <typename number_type> inline
antiderivative<number_type>& antiderivative<number_type>::operator*=( number_type rhs )
{
    PR         *= rhs;
    PRinv      *= rhs;
    Partanh_XR *= rhs;
    Partanh_YR *= rhs;
    Partanh_ZR *= rhs;
    Parctan_RZ *= rhs;
    Parctan_RY *= rhs;
    Parctan_RX *= rhs;
    return *this;
}

template <typename number_type> inline
antiderivative<number_type>& antiderivative<number_type>::operator*=( const polynomial &rhs )
{
    PR         *= rhs;
    PRinv      *= rhs;
    Partanh_XR *= rhs;
    Partanh_YR *= rhs;
    Partanh_ZR *= rhs;
    Parctan_RZ *= rhs;
    Parctan_RY *= rhs;
    Parctan_RX *= rhs;
    return *this;
}

template <typename number_type> inline
antiderivative<number_type>& antiderivative<number_type>::operator/=( number_type rhs )
{
    PR         /= rhs;
    PRinv      /= rhs;
    Partanh_XR /= rhs;
    Partanh_YR /= rhs;
    Partanh_ZR /= rhs;
    Parctan_RZ /= rhs;
    Parctan_RY /= rhs;
    Parctan_RX /= rhs;
    return *this;
}

template <typename number_type> inline
antiderivative<number_type>& antiderivative<number_type>::operator+=( const antiderivative &rhs )
{
    PR         += rhs.PR;
    PRinv      += rhs.PRinv;
    Partanh_XR += rhs.Partanh_XR;
    Partanh_YR += rhs.Partanh_YR;
    Partanh_ZR += rhs.Partanh_ZR;
    Parctan_RZ += rhs.Parctan_RZ;
    Parctan_RY += rhs.Parctan_RY;
    Parctan_RX += rhs.Parctan_RX;
    return *this;
}

template <typename number_type> inline
antiderivative<number_type>& antiderivative<number_type>::operator-=( const antiderivative &rhs )
{
    PR         -= rhs.PR;
    PRinv      -= rhs.PRinv;
    Partanh_XR -= rhs.Partanh_XR;
    Partanh_YR -= rhs.Partanh_YR;
    Partanh_ZR -= rhs.Partanh_ZR;
    Parctan_RZ -= rhs.Parctan_RZ;
    Parctan_RY -= rhs.Parctan_RY;
    Parctan_RX -= rhs.Parctan_RX;
    return *this;
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::operator*( number_type rhs ) const
{
    antiderivative result {*this};
    result *= rhs;
    return result; 
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::operator/( number_type rhs ) const
{
    antiderivative result {*this};
    result /= rhs;
    return result; 
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::operator*( const polynomial &rhs ) const
{
    antiderivative result {*this};
    result *= rhs;
    return result; 
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::operator+( antiderivative rhs ) const
{
    rhs += *this;
    return rhs;
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::operator-( const antiderivative &rhs ) const
{
    antiderivative result {*this};
    result -= rhs;
    return result;
}

template <typename number_type> inline
antiderivative<number_type> operator*( number_type lhs, antiderivative<number_type> rhs )
{
    return rhs *= lhs;
}

template <typename number_type> inline
antiderivative<number_type> operator*( const polynomial<number_type,6> &lhs, antiderivative<number_type> rhs )
{
    return rhs *= lhs;
}


template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::x1_integrate() const
{
    antiderivative result;
    result += impl::x1_integrate_R( PR );
    result += impl::x1_integrate_Rinv( PRinv );
    result += impl::x1_integrate_artanh_XR( Partanh_XR );
    result += impl::x1_integrate_artanh_YR( Partanh_YR );
    result += impl::x1_integrate_artanh_ZR( Partanh_ZR );
    result += impl::x1_integrate_arctan_RZ( Parctan_RZ );
    result += impl::x1_integrate_arctan_RY( Parctan_RY );
    result += impl::x1_integrate_arctan_RX( Parctan_RX );
    return result;
}

template <typename number_type> inline
antiderivative<number_type> x1_integrate( const antiderivative<number_type> &f )
{
    return f.x1_integrate();
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::x2_integrate() const
{
    antiderivative result;
    result += impl::x2_integrate_R( PR );
    result += impl::x2_integrate_Rinv( PRinv );
    result += impl::x2_integrate_artanh_XR( Partanh_XR );
    result += impl::x2_integrate_artanh_YR( Partanh_YR );
    result += impl::x2_integrate_artanh_ZR( Partanh_ZR );
    result += impl::x2_integrate_arctan_RZ( Parctan_RZ );
    result += impl::x2_integrate_arctan_RY( Parctan_RY );
    result += impl::x2_integrate_arctan_RX( Parctan_RX );
    return result;
}

template <typename number_type> inline
antiderivative<number_type> x2_integrate( const antiderivative<number_type> &f )
{
    return f.x2_integrate();
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::x3_integrate() const
{
    antiderivative result;
    result += impl::x3_integrate_R( PR );
    result += impl::x3_integrate_Rinv( PRinv );
    result += impl::x3_integrate_artanh_XR( Partanh_XR );
    result += impl::x3_integrate_artanh_YR( Partanh_YR );
    result += impl::x3_integrate_artanh_ZR( Partanh_ZR );
    result += impl::x3_integrate_arctan_RZ( Parctan_RZ );
    result += impl::x3_integrate_arctan_RY( Parctan_RY );
    result += impl::x3_integrate_arctan_RX( Parctan_RX );
    return result;
}

template <typename number_type> inline
antiderivative<number_type> x3_integrate( const antiderivative<number_type> &f )
{
    return f.x3_integrate();
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::y1_integrate() const
{
    antiderivative result;
    result += impl::y1_integrate_R( PR );
    result += impl::y1_integrate_Rinv( PRinv );
    result += impl::y1_integrate_artanh_XR( Partanh_XR );
    result += impl::y1_integrate_artanh_YR( Partanh_YR );
    result += impl::y1_integrate_artanh_ZR( Partanh_ZR );
    result += impl::y1_integrate_arctan_RZ( Parctan_RZ );
    result += impl::y1_integrate_arctan_RY( Parctan_RY );
    result += impl::y1_integrate_arctan_RX( Parctan_RX );
    return result;
}

template <typename number_type> inline
antiderivative<number_type> y1_integrate( const antiderivative<number_type> &f )
{
    return f.y1_integrate();
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::y2_integrate() const
{
    antiderivative result;
    result += impl::y2_integrate_R( PR );
    result += impl::y2_integrate_Rinv( PRinv );
    result += impl::y2_integrate_artanh_XR( Partanh_XR );
    result += impl::y2_integrate_artanh_YR( Partanh_YR );
    result += impl::y2_integrate_artanh_ZR( Partanh_ZR );
    result += impl::y2_integrate_arctan_RZ( Parctan_RZ );
    result += impl::y2_integrate_arctan_RY( Parctan_RY );
    result += impl::y2_integrate_arctan_RX( Parctan_RX );
    return result;
}

template <typename number_type> inline
antiderivative<number_type> y2_integrate( const antiderivative<number_type> &f )
{
    return f.y2_integrate();
}

template <typename number_type> inline
antiderivative<number_type> antiderivative<number_type>::y3_integrate() const
{
    antiderivative result;
    result += impl::y3_integrate_R( PR );
    result += impl::y3_integrate_Rinv( PRinv );
    result += impl::y3_integrate_artanh_XR( Partanh_XR );
    result += impl::y3_integrate_artanh_YR( Partanh_YR );
    result += impl::y3_integrate_artanh_ZR( Partanh_ZR );
    result += impl::y3_integrate_arctan_RZ( Parctan_RZ );
    result += impl::y3_integrate_arctan_RY( Parctan_RY );
    result += impl::y3_integrate_arctan_RX( Parctan_RX );
    return result;
}

template <typename number_type> inline
antiderivative<number_type> y3_integrate( const antiderivative<number_type> &f )
{
    return f.y3_integrate();
}


template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::eval( const std::array<number_type,6> &x ) const
{
    return static_cast<float_type>(        PR(x)) * impl::template      R   <float_type>(x) +
           static_cast<float_type>(     PRinv(x)) * impl::template      Rinv<float_type>(x) +
           static_cast<float_type>(Partanh_XR(x)) * impl::template artanh_XR<float_type>(x) +
           static_cast<float_type>(Partanh_YR(x)) * impl::template artanh_YR<float_type>(x) +
           static_cast<float_type>(Partanh_ZR(x)) * impl::template artanh_ZR<float_type>(x) +
           static_cast<float_type>(Parctan_RZ(x)) * impl::template arctan_RZ<float_type>(x) + 
           static_cast<float_type>(Parctan_RY(x)) * impl::template arctan_RY<float_type>(x) +
           static_cast<float_type>(Parctan_RX(x)) * impl::template arctan_RX<float_type>(x);
}

template <typename number_type>
template <typename float_type, typename output_iterator> inline
void antiderivative<number_type>::write_summands( const std::array<number_type,6> &x,
                                                  output_iterator out ) const
{
    *out++ = static_cast<float_type>(        PR(x)) * impl::template      R   <float_type>(x);
    *out++ = static_cast<float_type>(     PRinv(x)) * impl::template      Rinv<float_type>(x);
    *out++ = static_cast<float_type>(Partanh_XR(x)) * impl::template artanh_XR<float_type>(x);
    *out++ = static_cast<float_type>(Partanh_YR(x)) * impl::template artanh_YR<float_type>(x);
    *out++ = static_cast<float_type>(Partanh_ZR(x)) * impl::template artanh_ZR<float_type>(x);
    *out++ = static_cast<float_type>(Parctan_RZ(x)) * impl::template arctan_RZ<float_type>(x);
    *out++ = static_cast<float_type>(Parctan_RY(x)) * impl::template arctan_RY<float_type>(x);
    *out++ = static_cast<float_type>(Parctan_RX(x)) * impl::template arctan_RX<float_type>(x);
}



///////////////////////////////////////////////
//                                           //
// Implementation of antiderivative<>::impl  //
//                                           //
///////////////////////////////////////////////

template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::impl::R( const std::array<number_type,6> &pos )
{
    number_type arg {0};
    arg += (pos[0]-pos[3])*(pos[0]-pos[3]);
    arg += (pos[1]-pos[4])*(pos[1]-pos[4]);
    arg += (pos[2]-pos[5])*(pos[2]-pos[5]);
    return sqrt( static_cast<float_type>(arg) );
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::impl::Rinv( const std::array<number_type,6> &pos )
{
    number_type arg {0}; 
    arg += (pos[0]-pos[3])*(pos[0]-pos[3]);
    arg += (pos[1]-pos[4])*(pos[1]-pos[4]);
    arg += (pos[2]-pos[5])*(pos[2]-pos[5]);

    if ( arg == 0 ) return 0;

    return static_cast<float_type>(1) / sqrt( static_cast<float_type>(arg) ); 
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::impl::artanh_XR( const std::array<number_type,6> &pos )
{
    float_type r { R<float_type>(pos) };
    float_type x { static_cast<float_type>( pos[0] - pos[3] ) };
    float_type y { static_cast<float_type>( pos[1] - pos[4] ) };
    float_type z { static_cast<float_type>( pos[2] - pos[5] ) };

    if ( y == 0 && z == 0 ) return 0;
    return boost::math::atanh(x/r);
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::impl::artanh_YR( const std::array<number_type,6> &pos )
{
    float_type r { R<float_type>(pos) };
    float_type x { static_cast<float_type>( pos[0] - pos[3] ) };
    float_type y { static_cast<float_type>( pos[1] - pos[4] ) };
    float_type z { static_cast<float_type>( pos[2] - pos[5] ) };

    if ( x == 0 && z == 0 ) return 0;
    return boost::math::atanh(y/r);
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::impl::artanh_ZR( const std::array<number_type,6> &pos )
{
    float_type r { R<float_type>(pos) };
    float_type x { static_cast<float_type>( pos[0] - pos[3] ) };
    float_type y { static_cast<float_type>( pos[1] - pos[4] ) };
    float_type z { static_cast<float_type>( pos[2] - pos[5] ) };

    if ( x == 0 && y == 0 ) return 0;
    return boost::math::atanh(z/r);
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::impl::arctan_RZ( const std::array<number_type,6> &pos )
{
    float_type r { R<float_type>(pos) };
    float_type x { static_cast<float_type>( pos[0] - pos[3] ) };
    float_type y { static_cast<float_type>( pos[1] - pos[4] ) };
    float_type z { static_cast<float_type>( pos[2] - pos[5] ) };

    if ( r == 0 || z == 0 ) return 0;
    return atan((x*y)/(r*z));
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::impl::arctan_RY( const std::array<number_type,6> &pos )
{
    float_type r { R<float_type>(pos) };
    float_type x { static_cast<float_type>( pos[0] - pos[3] ) };
    float_type y { static_cast<float_type>( pos[1] - pos[4] ) };
    float_type z { static_cast<float_type>( pos[2] - pos[5] ) };

    if ( r == 0 || y == 0 ) return 0;
    return atan((x*z)/(r*y));
}

template <typename number_type> template <typename float_type> inline
float_type antiderivative<number_type>::impl::arctan_RX( const std::array<number_type,6> &pos )
{
    float_type r { R<float_type>(pos) };
    float_type x { static_cast<float_type>( pos[0] - pos[3] ) };
    float_type y { static_cast<float_type>( pos[1] - pos[4] ) };
    float_type z { static_cast<float_type>( pos[2] - pos[5] ) };

    if ( r == 0 || x == 0 ) return 0;
    return atan((y*z)/(r*x));
}

/////////////////////////////////////
// Integration of the function R.  //
/////////////////////////////////////

template <typename number_type> 
auto antiderivative<number_type>::impl::x1_integrate_R( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.PR = X/2;
            f.Partanh_XR = Xi/2;
            break;
        case 1:
            f.PR = (X*X + Xi)/3 + y1*X/2;
            f.Partanh_XR = y1*Xi/2;
            break;
        default:
            polynomial x1_n1 { multi_index {n-1,0,0,0,0,0} };
            f.PR = x1_n1*(X*X + Xi)/(n+2);
            f += rational(2*n+1,n+2)*y1           * cache[n-1];
            f -= rational(  n-1,n+2)*(y1*y1 + Xi) * cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y1_integrate_R( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.PR = -X/2;
            f.Partanh_XR = -Xi/2;
            break;
        case 1:
            f.PR = (X*X + Xi)/3 - x1*X/2;
            f.Partanh_XR = -x1*Xi/2;
            break;
        default:
            polynomial y1_n1 { multi_index {0,0,0,n-1,0,0} };
            f.PR = y1_n1*(X*X + Xi)/(n+2);
            f += rational(2*n+1,n+2)*x1           * cache[n-1];
            f -= rational(  n-1,n+2)*(x1*x1 + Xi) * cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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

template <typename number_type> 
auto antiderivative<number_type>::impl::x2_integrate_R( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.PR = Y/2;
            f.Partanh_YR = Upsilon/2;
            break;
        case 1:
            f.PR = (Y*Y + Upsilon)/3 + y2*Y/2;
            f.Partanh_YR = y2*Upsilon/2;
            break;
        default:
            polynomial x2_n1 { multi_index {0,n-1,0,0,0,0} };
            f.PR = x2_n1*(Y*Y + Upsilon)/(n+2);
            f += rational(2*n+1,n+2)*y2            * cache[n-1];
            f -= rational(  n-1,n+2)*(y2*y2 + Upsilon) * cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y2_integrate_R( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[4];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.PR = -Y/2;
            f.Partanh_YR = -Upsilon/2;
            break;
        case 1:
            f.PR = (Y*Y + Upsilon)/3 - x2*Y/2;
            f.Partanh_YR = -x2*Upsilon/2;
            break;
        default:
            polynomial y2_n1 { multi_index {0,0,0,0,n-1,0} };
            f.PR = y2_n1*(Y*Y + Upsilon)/(n+2);
            f += rational(2*n+1,n+2)*x2            * cache[n-1];
            f -= rational(  n-1,n+2)*(x2*x2 + Upsilon) * cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[4];
        exponents[4] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative<number_type>::impl::x3_integrate_R( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.PR = Z/2;
            f.Partanh_ZR = Theta/2;
            break;
        case 1:
            f.PR = (Z*Z + Theta)/3 + y3*Z/2;
            f.Partanh_ZR = y3*Theta/2;
            break;
        default:
            polynomial x3_n1 { multi_index {0,0,n-1,0,0,0} };
            f.PR = x3_n1*(Z*Z + Theta)/(n+2);
            f += rational(2*n+1,n+2)*y3                * cache[n-1];
            f -= rational(  n-1,n+2)*(y3*y3 + Theta  ) * cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y3_integrate_R( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[5];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.PR = -Z/2;
            f.Partanh_ZR = -Theta/2;
            break;
        case 1:
            f.PR = (Z*Z + Theta)/3 - x3*Z/2;
            f.Partanh_ZR = -x3*Theta/2;
            break;
        default:
            polynomial y3_n1 { multi_index {0,0,0,0,0,n-1} };
            f.PR = y3_n1*(Z*Z + Theta)/(n+2);
            f += rational(2*n+1,n+2)*x3                * cache[n-1];
            f -= rational(  n-1,n+2)*(x3*x3 + Theta  ) * cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[5];
        exponents[5] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////
// Integration of the function 1/R. //
//////////////////////////////////////

template <typename number_type> 
auto antiderivative<number_type>::impl::x1_integrate_Rinv( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR = 1;
            break;
        case 1:
            f.PR = 1;
            f.Partanh_XR = y1;
            break;
        default:
            polynomial x1_n1 { multi_index {n-1,0,0,0,0,0} };
            f.PR = x1_n1/n;
            f += rational(2*n-1,n)*y1        *cache[n-1];
            f -= rational(  n-1,n)*(y1*y1+Xi)*cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y1_integrate_Rinv( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR = -1;
            break;
        case 1:
            f.PR = 1;
            f.Partanh_XR = -x1;
            break;
        default:
            polynomial y1_n1 { multi_index {0,0,0,n-1,0,0} };
            f.PR = y1_n1/n;
            f += rational(2*n-1,n)*x1        *cache[n-1];
            f -= rational(  n-1,n)*(x1*x1+Xi)*cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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

template <typename number_type> 
auto antiderivative<number_type>::impl::x2_integrate_Rinv( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = 1;
            break;
        case 1:
            f.PR = 1;
            f.Partanh_YR = y2;
            break;
        default:
            polynomial x2_n1 { multi_index {0,n-1,0,0,0,0} };
            f.PR = x2_n1/n;
            f += rational(2*n-1,n)*y2             *cache[n-1];
            f -= rational(  n-1,n)*(y2*y2+Upsilon)*cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y2_integrate_Rinv( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[4];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = -1;
            break;
        case 1:
            f.PR = 1;
            f.Partanh_YR = -x2;
            break;
        default:
            polynomial y2_n1 { multi_index {0,0,0,0,n-1,0} };
            f.PR = y2_n1/n;
            f += rational(2*n-1,n)*x2             *cache[n-1];
            f -= rational(  n-1,n)*(x2*x2+Upsilon)*cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[4];
        exponents[4] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative<number_type>::impl::x3_integrate_Rinv( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR = 1;
            break;
        case 1:
            f.PR = 1;
            f.Partanh_ZR = y3;
            break;
        default:
            polynomial x3_n1 { multi_index {0,0,n-1,0,0,0} };
            f.PR = x3_n1/n;
            f += rational(2*n-1,n)*y3           *cache[n-1];
            f -= rational(  n-1,n)*(y3*y3+Theta)*cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y3_integrate_Rinv( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[5];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR = -1;
            break;
        case 1:
            f.PR = 1;
            f.Partanh_ZR = -x3;
            break;
        default:
            polynomial y3_n1 { multi_index {0,0,0,0,0,n-1} };
            f.PR = y3_n1/n;
            f += rational(2*n-1,n)*x3           *cache[n-1];
            f -= rational(  n-1,n)*(x3*x3+Theta)*cache[n-2];
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[5];
        exponents[5] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////////////
// Integration of the function artanh(X/R). //
//////////////////////////////////////////////

template <typename number_type> 
auto antiderivative<number_type>::impl::x1_integrate_artanh_XR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR =  X;
            f.PR         = -1;
            break;
        default:
            polynomial xx { rational(1,n+1), multi_index {n+1,0,0,0,0,0} };
            f.Partanh_XR = xx;
            f -= x1_integrate_Rinv( xx );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y1_integrate_artanh_XR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR = -X;
            f.PR         =  1;
            break;
        default:
            polynomial yy { rational(1,n+1), multi_index {0,0,0,n+1,0,0} };
            f.Partanh_XR = yy;
            f += y1_integrate_Rinv( yy );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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

template <typename number_type> 
auto antiderivative<number_type>::impl::x2_integrate_artanh_XR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR =  Y;
            f.Partanh_YR =  X;
            f.Parctan_RZ = -Z;
            break;
        default:
            polynomial x2_n1 { multi_index {0,n-1,0,0,0,0} };
            f  = (x2_n1*x2/(n+1))*cache[0];
            f += (y2*rational(n,n+1))*cache[n-1];
            f -=  x2_integrate_artanh_YR( rational(n,n+1)*X*x2_n1 );
            f +=  x2_integrate_arctan_RZ( rational(n,n+1)*Z*x2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y2_integrate_artanh_XR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[4];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR = -Y;
            f.Partanh_YR = -X;
            f.Parctan_RZ =  Z;
            break;
        default:
            polynomial y2_n1 { multi_index {0,0,0,0,n-1,0} };
            f  = (y2_n1*y2/(n+1))*cache[0];
            f += (x2*rational(n,n+1))*cache[n-1];
            f +=  y2_integrate_artanh_YR( rational(n,n+1)*X*y2_n1 );
            f -=  y2_integrate_arctan_RZ( rational(n,n+1)*Z*y2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[4];
        exponents[4] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative<number_type>::impl::x3_integrate_artanh_XR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR =  Z;
            f.Partanh_ZR =  X;
            f.Parctan_RY = -Y;
            break;
        default:
            polynomial x3_n1 { multi_index {0,0,n-1,0,0,0} };
            f  = (x3_n1*x3/(n+1))*cache[0];
            f += (y3*rational(n,n+1))*cache[n-1];
            f -=  x3_integrate_artanh_ZR( rational(n,n+1)*X*x3_n1 );
            f +=  x3_integrate_arctan_RY( rational(n,n+1)*Y*x3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y3_integrate_artanh_XR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[5];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR = -Z;
            f.Partanh_ZR = -X;
            f.Parctan_RY =  Y;
            break;
        default:
            polynomial y3_n1 { multi_index {0,0,0,0,0,n-1} };
            f  = (y3_n1*y3/(n+1))*cache[0];
            f += (x3*rational(n,n+1))*cache[n-1];
            f +=  y3_integrate_artanh_ZR( rational(n,n+1)*X*y3_n1 );
            f -=  y3_integrate_arctan_RY( rational(n,n+1)*Y*y3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[5];
        exponents[5] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////////////
// Integration of the function artanh(Y/R). //
//////////////////////////////////////////////

template <typename number_type> 
auto antiderivative<number_type>::impl::x1_integrate_artanh_YR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR =  X;
            f.Partanh_XR =  Y;
            f.Parctan_RZ = -Z;
            break;
        default:
            polynomial x1_n1 { multi_index {n-1,0,0,0,0,0} };
            f  = (x1_n1*x1/(n+1))*cache[0];
            f += (y1*rational(n,n+1))*cache[n-1];
            f -=  x1_integrate_artanh_XR( rational(n,n+1)*Y*x1_n1 );
            f +=  x1_integrate_arctan_RZ( rational(n,n+1)*Z*x1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y1_integrate_artanh_YR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = -X;
            f.Partanh_XR = -Y;
            f.Parctan_RZ =  Z;
            break;
        default:
            polynomial y1_n1 { multi_index {0,0,0,n-1,0,0} };
            f  = (y1_n1*y1/(n+1))*cache[0];
            f += (x1*rational(n,n+1))*cache[n-1];
            f +=  y1_integrate_artanh_XR( rational(n,n+1)*Y*y1_n1 );
            f -=  y1_integrate_arctan_RZ( rational(n,n+1)*Z*y1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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

template <typename number_type> 
auto antiderivative<number_type>::impl::x2_integrate_artanh_YR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR =  Y;
            f.PR         = -1;
            break;
        default:
            polynomial xx { rational(1,n+1), multi_index {0,n+1,0,0,0,0} };
            f.Partanh_YR = xx;
            f -= x2_integrate_Rinv( xx );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y2_integrate_artanh_YR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[4];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = -Y;
            f.PR         =  1;
            break;
        default:
            polynomial yy { rational(1,n+1), multi_index {0,0,0,0,n+1,0} };
            f.Partanh_YR = yy;
            f += y2_integrate_Rinv( yy );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[4];
        exponents[4] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative<number_type>::impl::x3_integrate_artanh_YR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR =  Z;
            f.Partanh_ZR =  Y;
            f.Parctan_RX = -X;
            break;
        default:
            polynomial x3_n1 { multi_index {0,0,n-1,0,0,0} };
            f  = (x3_n1*x3/(n+1))*cache[0];
            f += (y3*rational(n,n+1))*cache[n-1];
            f -=  x3_integrate_artanh_ZR( rational(n,n+1)*Y*x3_n1 );
            f +=  x3_integrate_arctan_RX( rational(n,n+1)*X*x3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y3_integrate_artanh_YR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[5];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = -Z;
            f.Partanh_ZR = -Y;
            f.Parctan_RX =  X;
            break;
        default:
            polynomial y3_n1 { multi_index {0,0,0,0,0,n-1} };
            f  = (y3_n1*y3/(n+1))*cache[0];
            f += (x3*rational(n,n+1))*cache[n-1];
            f +=  y3_integrate_artanh_ZR( rational(n,n+1)*Y*y3_n1 );
            f -=  y3_integrate_arctan_RX( rational(n,n+1)*X*y3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[5];
        exponents[5] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////////////
// Integration of the function artanh(Z/R). //
//////////////////////////////////////////////

template <typename number_type> 
auto antiderivative<number_type>::impl::x1_integrate_artanh_ZR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR =  X;
            f.Partanh_XR =  Z;
            f.Parctan_RY = -Y;
            break;
        default:
            polynomial x1_n1 { multi_index {n-1,0,0,0,0,0} };
            f  = (x1_n1*x1/(n+1))*cache[0];
            f += (y1*rational(n,n+1))*cache[n-1];
            f -=  x1_integrate_artanh_XR( rational(n,n+1)*Z*x1_n1 );
            f +=  x1_integrate_arctan_RY( rational(n,n+1)*Y*x1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y1_integrate_artanh_ZR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR = -X;
            f.Partanh_XR = -Z;
            f.Parctan_RY =  Y;
            break;
        default:
            polynomial y1_n1 { multi_index {0,0,0,n-1,0,0} };
            f  = (y1_n1*y1/(n+1))*cache[0];
            f += (x1*rational(n,n+1))*cache[n-1];
            f +=  y1_integrate_artanh_XR( rational(n,n+1)*Z*y1_n1 );
            f -=  y1_integrate_arctan_RY( rational(n,n+1)*Y*y1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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

template <typename number_type> 
auto antiderivative<number_type>::impl::x2_integrate_artanh_ZR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR =  Y;
            f.Partanh_YR =  Z;
            f.Parctan_RX = -X;
            break;
        default:
            polynomial x2_n1 { multi_index {0,n-1,0,0,0,0} };
            f  = (x2_n1*x2/(n+1))*cache[0];
            f += (y2*rational(n,n+1))*cache[n-1];
            f -=  x2_integrate_artanh_YR( rational(n,n+1)*Z*x2_n1 );
            f +=  x2_integrate_arctan_RX( rational(n,n+1)*X*x2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y2_integrate_artanh_ZR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[4];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR = -Y;
            f.Partanh_YR = -Z;
            f.Parctan_RX =  X;
            break;
        default:
            polynomial y2_n1 { multi_index {0,0,0,0,n-1,0} };
            f  = (y2_n1*y2/(n+1))*cache[0];
            f += (x2*rational(n,n+1))*cache[n-1];
            f +=  y2_integrate_artanh_YR( rational(n,n+1)*Z*y2_n1 );
            f -=  y2_integrate_arctan_RX( rational(n,n+1)*X*y2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[4];
        exponents[4] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative<number_type>::impl::x3_integrate_artanh_ZR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR =  Z;
            f.PR         = -1;
            break;
        default:
            polynomial xx { rational(1,n+1), multi_index {0,0,n+1,0,0,0} };
            f.Partanh_ZR = xx;
            f -= x3_integrate_Rinv( xx );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y3_integrate_artanh_ZR( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[5];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR = -Z;
            f.PR         =  1;
            break;
        default:
            polynomial yy { rational(1,n+1), multi_index {0,0,0,0,0,n+1} };
            f.Partanh_ZR = yy;
            f += y3_integrate_Rinv( yy );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[5];
        exponents[5] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////////////////////////////
// Integration of the function arctan( (X/R)*(Y/R)*(R/Z) ). //
//////////////////////////////////////////////////////////////

template <typename number_type> 
auto antiderivative<number_type>::impl::x1_integrate_arctan_RZ( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR =  Z;
            f.Parctan_RZ =  X;
            break;
        default:
            polynomial x1_n1 { multi_index {n-1,0,0,0,0,0} };
            f  = (x1_n1*x1/(n+1))*cache[0];
            f += (y1*rational(n,n+1))*cache[n-1];
            f -= x1_integrate_artanh_YR( rational(n,n+1)*Z*x1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y1_integrate_arctan_RZ( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = -Z;
            f.Parctan_RZ = -X;
            break;
        default:
            polynomial y1_n1 { multi_index {0,0,0,n-1,0,0} };
            f  = (y1_n1*y1/(n+1))*cache[0];
            f += (x1*rational(n,n+1))*cache[n-1];
            f += y1_integrate_artanh_YR( rational(n,n+1)*Z*y1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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

template <typename number_type> 
auto antiderivative<number_type>::impl::x2_integrate_arctan_RZ( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR =  Z;
            f.Parctan_RZ =  Y;
            break;
        default:
            polynomial x2_n1 { multi_index {0,n-1,0,0,0,0} };
            f  = (x2_n1*x2/(n+1))*cache[0];
            f += (y2*rational(n,n+1))*cache[n-1];
            f -= x2_integrate_artanh_XR( rational(n,n+1)*Z*x2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y2_integrate_arctan_RZ( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[4];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR = -Z;
            f.Parctan_RZ = -Y;
            break;
        default:
            polynomial y2_n1 { multi_index {0,0,0,0,n-1,0} };
            f  = (y2_n1*y2/(n+1))*cache[0];
            f += (x2*rational(n,n+1))*cache[n-1];
            f += y2_integrate_artanh_XR( rational(n,n+1)*Z*y2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[4];
        exponents[4] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative<number_type>::impl::x3_integrate_arctan_RZ( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = -X;
            f.Partanh_XR = -Y;
            f.Parctan_RZ =  Z;
            break;
        default:
            polynomial x3_n1 { multi_index {0,0,n-1,0,0,0} };
            f  = (x3_n1*x3/(n+1))*cache[0];
            f += (y3*rational(n,n+1))*cache[n-1];
            f += x3_integrate_artanh_YR( rational(n,n+1)*X*x3_n1 );
            f += x3_integrate_artanh_XR( rational(n,n+1)*Y*x3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y3_integrate_arctan_RZ( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[5];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR =  X;
            f.Partanh_XR =  Y;
            f.Parctan_RZ = -Z;
            break;
        default:
            polynomial y3_n1 { multi_index {0,0,0,0,0,n-1} };
            f  = (y3_n1*y3/(n+1))*cache[0];
            f += (x3*rational(n,n+1))*cache[n-1];
            f -= y3_integrate_artanh_YR( rational(n,n+1)*X*y3_n1 );
            f -= y3_integrate_artanh_XR( rational(n,n+1)*Y*y3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[5];
        exponents[5] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////////////////////////////
// Integration of the function arctan( (X/R)*(Z/R)*(R/Y) ). //
//////////////////////////////////////////////////////////////

template <typename number_type> 
auto antiderivative<number_type>::impl::x1_integrate_arctan_RY( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR =  Y;
            f.Parctan_RY =  X;
            break;
        default:
            polynomial x1_n1 { multi_index {n-1,0,0,0,0,0} };
            f  = (x1_n1*x1/(n+1))*cache[0];
            f += (y1*rational(n,n+1))*cache[n-1];
            f -= x1_integrate_artanh_ZR( rational(n,n+1)*Y*x1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y1_integrate_arctan_RY( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR = -Y;
            f.Parctan_RY = -X;
            break;
        default:
            polynomial y1_n1 { multi_index {0,0,0,n-1,0,0} };
            f  = (y1_n1*y1/(n+1))*cache[0];
            f += (x1*rational(n,n+1))*cache[n-1];
            f += y1_integrate_artanh_ZR( rational(n,n+1)*Y*y1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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

template <typename number_type> 
auto antiderivative<number_type>::impl::x2_integrate_arctan_RY( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR = -Z;
            f.Partanh_ZR = -X;
            f.Parctan_RY =  Y;
            break;
        default:
            polynomial x2_n1 { multi_index {0,n-1,0,0,0,0} };
            f  = (x2_n1*x2/(n+1))*cache[0];
            f += (y2*rational(n,n+1))*cache[n-1];
            f += x2_integrate_artanh_ZR( rational(n,n+1)*X*x2_n1 );
            f += x2_integrate_artanh_XR( rational(n,n+1)*Z*x2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y2_integrate_arctan_RY( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[4];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR =  Z;
            f.Partanh_ZR =  X;
            f.Parctan_RY = -Y;
            break;
        default:
            polynomial y2_n1 { multi_index {0,0,0,0,n-1,0} };
            f  = (y2_n1*y2/(n+1))*cache[0];
            f += (x2*rational(n,n+1))*cache[n-1];
            f -= y2_integrate_artanh_ZR( rational(n,n+1)*X*y2_n1 );
            f -= y2_integrate_artanh_XR( rational(n,n+1)*Z*y2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[4];
        exponents[4] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative<number_type>::impl::x3_integrate_arctan_RY( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR =  Y;
            f.Parctan_RY =  Z;
            break;
        default:
            polynomial x3_n1 { multi_index {0,0,n-1,0,0,0} };
            f  = (x3_n1*x3/(n+1))*cache[0];
            f += (y3*rational(n,n+1))*cache[n-1];
            f -= x3_integrate_artanh_XR( rational(n,n+1)*Y*x3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y3_integrate_arctan_RY( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[5];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_XR = -Y;
            f.Parctan_RY = -Z;
            break;
        default:
            polynomial y3_n1 { multi_index {0,0,0,0,0,n-1} };
            f  = (y3_n1*y3/(n+1))*cache[0];
            f += (x3*rational(n,n+1))*cache[n-1];
            f += y3_integrate_artanh_XR( rational(n,n+1)*Y*y3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[5];
        exponents[5] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

//////////////////////////////////////////////////////////////
// Integration of the function arctan( (Y/R)*(Z/R)*(R/X) ). //
//////////////////////////////////////////////////////////////

template <typename number_type> 
auto antiderivative<number_type>::impl::x1_integrate_arctan_RX( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[0];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = -Z;
            f.Partanh_ZR = -Y;
            f.Parctan_RX =  X;
            break;
        default:
            polynomial x1_n1 { multi_index {n-1,0,0,0,0,0} };
            f  = (x1_n1*x1/(n+1))*cache[0];
            f += (y1*rational(n,n+1))*cache[n-1];
            f += x1_integrate_artanh_ZR( rational(n,n+1)*Y*x1_n1 );
            f += x1_integrate_artanh_YR( rational(n,n+1)*Z*x1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y1_integrate_arctan_RX( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[3];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR =  Z;
            f.Partanh_ZR =  Y;
            f.Parctan_RX = -X;
            break;
        default:
            polynomial y1_n1 { multi_index {0,0,0,n-1,0,0} };
            f  = (y1_n1*y1/(n+1))*cache[0];
            f += (x1*rational(n,n+1))*cache[n-1];
            f -= y1_integrate_artanh_ZR( rational(n,n+1)*Y*y1_n1 );
            f -= y1_integrate_artanh_YR( rational(n,n+1)*Z*y1_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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

template <typename number_type> 
auto antiderivative<number_type>::impl::x2_integrate_arctan_RX( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[1];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR =  X;
            f.Parctan_RX =  Y;
            break;
        default:
            polynomial x2_n1 { multi_index {0,n-1,0,0,0,0} };
            f  = (x2_n1*x2/(n+1))*cache[0];
            f += (y2*rational(n,n+1))*cache[n-1];
            f -= x2_integrate_artanh_ZR( rational(n,n+1)*X*x2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y2_integrate_arctan_RX( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[4];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_ZR = -X;
            f.Parctan_RX = -Y;
            break;
        default:
            polynomial y2_n1 { multi_index {0,0,0,0,n-1,0} };
            f  = (y2_n1*y2/(n+1))*cache[0];
            f += (x2*rational(n,n+1))*cache[n-1];
            f += y2_integrate_artanh_ZR( rational(n,n+1)*X*y2_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[4];
        exponents[4] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

template <typename number_type> 
auto antiderivative<number_type>::impl::x3_integrate_arctan_RX( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[2];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR =  X;
            f.Parctan_RX =  Z;
            break;
        default:
            polynomial x3_n1 { multi_index {0,0,n-1,0,0,0} };
            f  = (x3_n1*x3/(n+1))*cache[0];
            f += (y3*rational(n,n+1))*cache[n-1];
            f -= x3_integrate_artanh_YR( rational(n,n+1)*X*x3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
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
auto antiderivative<number_type>::impl::y3_integrate_arctan_RX( const polynomial &P )
-> antiderivative
{
    thread_local std::vector<antiderivative> cache;

    auto max_degree = P.max_degrees()[5];
    while ( cache.size() < (max_degree+1) )
    {
        auto n = cache.size();
        antiderivative f;

        switch ( n )
        {
        case 0:
            f.Partanh_YR = -X;
            f.Parctan_RX = -Z;
            break;
        default:
            polynomial y3_n1 { multi_index {0,0,0,0,0,n-1} };
            f  = (y3_n1*y3/(n+1))*cache[0];
            f += (x3*rational(n,n+1))*cache[n-1];
            f += y3_integrate_artanh_YR( rational(n,n+1)*X*y3_n1 );
            break;
        }
        cache.push_back(std::move(f));
    }

    antiderivative result;
    for ( const auto &term: P.terms )
    {
        multi_index exponents = term.first;
        number_type coeff     = term.second; 
        size_t n              = exponents[5];
        exponents[5] = 0;

        result += polynomial(exponents,coeff)*cache[n];
    }
    return result;
}

}

