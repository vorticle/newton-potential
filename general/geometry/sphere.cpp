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
#include <geometry/sphere.h>

namespace geometry
{

namespace
{

constexpr real eps = 16*std::numeric_limits<real>::epsilon();

inline
real generate_givens_rotation( real a, real b )
{
    real c, s, rho;
    if ( b == 0 )
    {
        c = 1; s = 0;
    }
    else
    {
        real rinv = 1./std::hypot(a,b);
        c = a*rinv; s = -b*rinv;
    }

    if ( c == 0 )
    {
        rho = 1;
    }
    else if ( std::abs(s) < std::abs(c) )
    {
        rho = (c>0) ? s/2 : -s/2;
    }
    else
    {
        rho = (s>0) ? 2/c : -2/c;
    }
    return rho;
}

inline
void givens_rotation_from_rho( real rho, real &c, real &s )
{
    if ( rho == 1 )
    {
        c = 0; s = 1;
    }
    else if ( std::abs(rho) < 1 )
    {
        s = 2*rho; c = std::sqrt(1-s*s);
    }
    else
    {
        c = 2/rho; s = std::sqrt(1-c*c);
    }
}

template <size_t m, size_t n>
void givens_qr( real *A )
{
    for ( size_t j = 0; j < n; ++j )
    {
        for ( size_t i = m - 1; i > j; --i )
        {
           real rho = generate_givens_rotation( A[m*j + i-1], A[m*j + i] );
           real c, s; givens_rotation_from_rho(rho,c,s);

           real t1, t2;
           for ( size_t k = j; k < n; ++k )
           {
               t1 = A[ m*k + i-1 ]; t2 = A[ m*k + i ];
               A[ m*k + i-1 ] = c*t1 - s*t2;
               A[ m*k + i   ] = s*t1 + c*t2;
           }
           A[ m*j + i ] = rho;
        }
    }
}

template <size_t m, size_t n>
real estimate_condition( const real *QR )
{
    // Given the QR decomposition of a matrix A, this method
    // computes an estimate of the condition number of A.
    // Golub, van Loan, Matrix Computations, 4th edition,
    // Algorithm 3.5.1

    // Check for actual zeros.
    for ( size_t k = 0; k < n; ++k )
    {
        if ( QR[m*k+k] == 0 )
            return 2/eps;
    }

    real p[n] {};
    real y[n] {};
    for ( size_t k = n; k-- > 0; )
    {
        real y_plus  = ( 1-p[k])/QR[m*k + k];
        real y_minus = (-1-p[k])/QR[m*k + k];
        real p_plus  = 0, p_minus = 0;
        for ( size_t j = 0; j < k; ++j )
        {
            p_plus  += std::abs( p[j] + QR[m*k+j]*y_plus  );
            p_minus += std::abs( p[j] + QR[m*k+j]*y_minus );
        }

        if ( std::abs(y_plus) + p_plus > std::abs(y_minus) + p_minus)
        {
            y[k] = y_plus;
            for ( size_t j = 0; j < k; ++j )
                p[j] = p[j] + QR[m*k+j]*y_plus;
        }
        else
        {
            y[k] = y_minus;
            for ( size_t j = 0; j < k; ++j )
                p[j] = p[j] + QR[m*k+j]*y_minus;
        }
    }

    real norm_y = 0, norm_R = 0;
    for ( size_t i = 0; i < n; ++i )
    {
        norm_y = std::max(norm_y,y[i]);
        real row_sum = 0;
        for ( size_t j = i; j < n; ++j )
            row_sum += std::abs(QR[i+j*m]);
        norm_R = std::max(norm_R,row_sum);
    }
    return norm_y*norm_R;
}

template <size_t m, size_t n>
void least_squares_solve( const real *A, real *b )
{
    // Given QR-decomposition of A, dimensions (m,n) with m >= n.
    // This method computes the least squares solution of Ax = b.

    // Compute y  = Q^T b.
    for ( size_t j = 0; j < n; ++j )
    {
        for ( size_t i = m - 1; i > j; --i )
        {
           real rho = A[ m*j + i ];
           real c, s; givens_rotation_from_rho(rho,c,s);

           real t1 = b[i-1], t2 = b[i];
           b[ i-1 ] = c*t1 - s*t2;
           b[ i   ] = s*t1 + c*t2;
        }
    }

    // Solve Rx = y.
    for ( size_t i = n; i-- > 0; )
    {
        for ( size_t j = i + 1; j < n; ++j )
        {
            b[i] -= A[ m*j + i ]*b[j];
        }
        b[i] /= A[ m*i + i ];
    }
}

template <size_t m, size_t n>
void transposed_least_squares_solve( const real *A, real *a, real *b )
{
    // Given QR-decomposition of A, dimensions (m,n) with m >= n.
    // This method computes the minimum norm solution of A^T x = a, A^T x = b.
    // Golub van Loan, Algorithm 5.6.2.

    // L := R^T;
    // Solve Ly = a and Ly = b respectively.
    for ( size_t i = 0; i < n; ++i )
    {
        for ( size_t j = 0; j < i; ++j )
        {
            a[i] -= A[m*i + j]*a[j];
            b[i] -= A[m*i + j]*b[j];
        }
        a[i] /= A[m*i + i];
        b[i] /= A[m*i + i];
    }
    for ( size_t i = n; i < m; ++i ) a[i] = b[i] = 0;

    // Compute x = Qy.
    for ( size_t j = n; j-- > 0; )
    {
        for ( size_t i = j + 1; i < m; ++i )
        {
           real rho = A[ m*j + i ];
           real c, s; givens_rotation_from_rho(rho,c,s);

           real t1 = a[i-1], t2 = a[i];
           a[ i-1 ] =  c*t1 + s*t2;
           a[ i   ] = -s*t1 + c*t2;

           t1 = b[i-1]; t2 = b[i];
           b[ i-1 ] =  c*t1 + s*t2;
           b[ i   ] = -s*t1 + c*t2;
        }
    }
}

class basis
{
public:
    basis() = delete;
    basis( sphere S );
    basis( sphere S0, sphere S1 );
    basis( sphere S0, sphere S1, sphere S2 );
    basis( sphere S0, sphere S1, sphere S2, sphere S3 );

    sphere operator()( size_t i ) noexcept { return V[i]; }

    sphere ball()   { update(); return bounds;        }
    real   radius() { update(); return bounds.radius; }
    point  centre() { update(); return bounds.centre; }

    size_t size() const noexcept { return n;    }
    bool is_violated_by( sphere S ) noexcept;
    bool is_violated_by( const basis &B ) noexcept;

private:
    size_t   n;
    sphere   V[ 4 ];


    void   update();
    bool   clean { false };
    sphere bounds;
};

inline
basis::basis( sphere S )
{
    V[0] = S; n = 1;
    bounds = S;
    clean = true;
}

inline
basis::basis( sphere S0, sphere S1 )
{
    V[0] = S0;
    V[1] = S1;
    n = 2;
    clean = false;
}

inline
basis::basis( sphere S0, sphere S1, sphere S2 )
{
    V[0] = S0;
    V[1] = S1;
    V[2] = S2;
    n = 3;
    clean = false;
}

inline
basis::basis( sphere S0, sphere S1, sphere S2, sphere S3 )
{
    V[0] = S0;
    V[1] = S1;
    V[2] = S2;
    V[3] = S3;
    n = 4;
    clean = false;
}

void basis::update()
{
    if ( clean == false )
    {
        if ( n == 1 )
        {
            bounds = V[0];
            clean  = true;
            return;
        }

        const real r0 = V[0].radius;
        real max_rho = r0;
        real AQR[3*3], a[3], b[3];
        for ( size_t i = 1; i < n; ++i )
        {
            AQR[ 0 + 3*(i-1) ] = (V[i].centre.x - V[0].centre.x);
            AQR[ 1 + 3*(i-1) ] = (V[i].centre.y - V[0].centre.y);
            AQR[ 2 + 3*(i-1) ] = (V[i].centre.z - V[0].centre.z);

            real ri = V[i].radius;
            a[i-1] = ri - r0;
            b[i-1] = ((V[i].centre-V[0].centre).r2() + r0*r0 - ri*ri)/2;
            max_rho = std::max( max_rho, ri );
        }

        real condition = 0;
        switch (n)
        {
        case 2: givens_qr<3,1>(AQR);
                condition = estimate_condition<3,1>(AQR);
                transposed_least_squares_solve<3,1>(AQR,a,b);
                break;
        case 3: givens_qr<3,2>(AQR);
                condition = estimate_condition<3,2>(AQR);
                transposed_least_squares_solve<3,2>(AQR,a,b);
                break;
        case 4: givens_qr<3,3>(AQR);
                condition = estimate_condition<3,3>(AQR);
                transposed_least_squares_solve<3,3>(AQR,a,b);
                break;
        }
        if ( condition*eps > 1 )
        {
            n--; update(); return;
        }
        point e { a[0], a[1], a[2] };
        point f { b[0], b[1], b[2] };

        // At this point we know the centre of the desired circle is
        // V[0].centre + f + rho*e. It remains to find rho. Rho is the
        // solution of a quadratic equation, there are thus two possible
        // solutions rho_minus and rho_plus.
        const real p = (r0 + scal_prod(e,f))/(1-e.r2());
        const real q = (f.r2() - r0*r0)/(1-e.r2());
        const real root = std::sqrt(std::abs(p*p+q));
        const real rho_minus = p - root;
        const real rho_plus  = p + root;

        const point z_minus = rho_minus*e + f;
        const point z_plus  = rho_plus *e + f;

        auto in_convex_hull = [this,AQR]( point z ) -> bool
        {
            // Solve Al = z.
            real l[4] = { 1, z.x, z.y, z.z };
            switch ( n )
            {
            case 2: least_squares_solve<3,1>(AQR,l+1); break;
            case 3: least_squares_solve<3,2>(AQR,l+1); break;
            case 4: least_squares_solve<3,3>(AQR,l+1); break;
            }

            for ( size_t i = 1; i <= n; ++i )
                l[0] -= l[i];

            for ( size_t i = 0; i <= n; ++i )
            {
                if ( l[i] < -eps || l[i] > 1+eps )
                    return false;
            }
            return true;
        };

        // Try if one of the solutions is feasible.
        // Otherwise shrink and try again.
        if ( rho_minus >= max_rho && in_convex_hull(z_minus) )
        {
            bounds.radius = rho_minus;
            bounds.centre = V[0].centre + z_minus;
            clean = true;
        }
        else if ( rho_plus >= max_rho && in_convex_hull(z_plus) )
        {
            bounds.radius = rho_plus;
            bounds.centre = V[0].centre + z_plus;
            clean = true;
        }
        else
        {
            n--; update();
        }
    }
}

inline
bool basis::is_violated_by( sphere S ) noexcept
{
    for ( size_t i = 0; i < n; ++i )
    {
        if ( V[i].centre == S.centre && V[i].radius == S.radius )
            return false;
    }

    update();
    real excess = std::max(0., (S.centre-bounds.centre).r() + S.radius - bounds.radius );
    return ( excess*excess > eps*bounds.radius*bounds.radius );
}

inline
bool basis::is_violated_by( const basis &B ) noexcept
{
    for ( size_t i = 0; i < B.size(); ++i )
    {
        if ( is_violated_by(B.V[i]) )
            return true;
    }
    return false;
}

basis basis_computation( basis V, sphere B )
{
    std::vector<basis> candidates;
    switch ( V.size() )
    {
    case 1:
        candidates.push_back( basis { B } );
        candidates.push_back( basis { B, V(0) } );
        break;
    case 2:
        candidates.push_back( basis { B } );
        candidates.push_back( basis { B, V(0) } );
        candidates.push_back( basis { B, V(1) } );
        candidates.push_back( basis { B, V(1), V(0) } );
        break;
    case 3:
        candidates.push_back( basis { B } );
        candidates.push_back( basis { B, V(0) } );
        candidates.push_back( basis { B, V(1) } );
        candidates.push_back( basis { B, V(2) } );
        candidates.push_back( basis { B, V(0), V(1) } );
        candidates.push_back( basis { B, V(0), V(2) } );
        candidates.push_back( basis { B, V(1), V(2) } );
        candidates.push_back( basis { B, V(0), V(1), V(2) } );
        candidates.push_back( basis { B, V(1), V(0), V(2) } );
        candidates.push_back( basis { B, V(2), V(0), V(1) } );
        break;
    case 4:
        candidates.push_back( basis { B } );
        candidates.push_back( basis { B, V(0) } );
        candidates.push_back( basis { B, V(1) } );
        candidates.push_back( basis { B, V(2) } );
        candidates.push_back( basis { B, V(3) } );
        candidates.push_back( basis { B, V(0), V(1) } );
        candidates.push_back( basis { B, V(0), V(2) } );
        candidates.push_back( basis { B, V(0), V(3) } );
        candidates.push_back( basis { B, V(1), V(2) } );
        candidates.push_back( basis { B, V(1), V(3) } );
        candidates.push_back( basis { B, V(2), V(3) } );
        candidates.push_back( basis { B, V(1), V(2), V(3) } );
        candidates.push_back( basis { B, V(0), V(2), V(3) } );
        candidates.push_back( basis { B, V(0), V(1), V(3) } );
        candidates.push_back( basis { B, V(0), V(1), V(2) } );
        break;
    }

    for ( size_t i = 0; i < candidates.size(); ++i )
    {
        if ( ! candidates[i].is_violated_by(V) )
        {
            return candidates[i];
        }
    }

    throw std::runtime_error { "bounding_sphere: Could not compute a new basis." };
}

}

/*!
 * \brief Computes the minimal bounding sphere of a set of points.
 * \see Kaspar Fischer and Bernd Gärtner, "The Smallest Enclosing Ball of Balls:
 *      Combinatorial Structure and Algorithms".
 * \return The minimal bounding sphere of the given set of points.
 */
sphere bounding_sphere( const point *begin, const point *end )
{
    if ( begin == end )
        throw std::logic_error { "Cannot compute bounding sphere of an empty set." };

    basis V { sphere { *begin, 0 } };
    while ( true )
    {
        real max_dist = length(*begin - V.centre());
        const point* farthest = begin;
        for ( const point *i = begin + 1; i != end; ++i )
        {
            real dist = length(*i-V.centre());
            if ( dist > max_dist )
            {
                max_dist = dist;
                farthest = i;
            }
        }

        if ( V.is_violated_by( sphere { *farthest, 0 } ) )
        {
            V = basis_computation( V, sphere { *farthest, 0 } );
        }
        else break;
    }
    return V.ball();
}

/*!
 * \brief Computes the minimal bounding sphere of a set of spheres.
 * \see Kaspar Fischer and Bernd Gärtner, "The Smallest Enclosing Ball of Balls:
 *      Combinatorial Structure and Algorithms".
 * \return The minimal bounding sphere of the given set of spheres.
 */
sphere bounding_sphere( const sphere *begin, const sphere *end )
{
    if ( begin == end )
        throw std::logic_error { "Cannot compute bounding sphere of an empty set." };

    basis V { *begin };
    while ( true )
    {
        real max_dist = length(begin->centre - V.centre()) + begin->radius;
        const sphere* farthest = begin;
        for ( const sphere *i = begin + 1; i != end; ++i )
        {
            real dist = length(i->centre - V.centre()) + i->radius;
            if ( dist > max_dist )
            {
                max_dist = dist;
                farthest = i;
            }
        }

        if ( V.is_violated_by(*farthest) )
        {
            V = basis_computation(V,*farthest);
        }
        else break;
    }
    return V.ball();
}

sphere bounding_sphere( point x0, point x1, point x2 )
{
    real e1 = (x1-x0).r2();
    real e2 = (x2-x0).r2();
    if ( e1 < e2 )
    {
        std::swap(x1,x2);
        std::swap(e1,e2);
    }

    // x0 -- x1 is longest edge of triangle
    // Take bounding sphere of that edge.
    real  radius_squared = 0.25 * e1;
    point centre         = 0.5  * (x0+x1);

    if ( (x2 - centre).r2() <= radius_squared )
    {
        // If this sphere already contains x2, we are done.
        return sphere { centre, std::sqrt(radius_squared) };
    }
    else
    {
        // Otherwise we know that it is the minimal sphere through the
        // three points, the rotated circumcircle in that plane.
        // Explicit formula for that is here.
        point E1 = x1 - x0;
        point E2 = x2 - x0;
        point cross = cross_prod(E1,E2);
        point vec   = cross_prod(e1*E2 - e2*E1, cross) / ( 2*cross.r2() );
        return sphere { x0 + vec, vec.r() };
    }
}

}

