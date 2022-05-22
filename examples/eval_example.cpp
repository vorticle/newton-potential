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
#include <cmath>
#include <iostream>
#include <boost/multiprecision/gmp.hpp>
#include <boxpot/antiderivative.h>

using namespace boxpot;

int main()
{
    using rational = boost::multiprecision::mpq_rational;
    using idx6     = multi_index<6>;
    using poly6    = polynomial<rational,6>;

    // Short-hands for convenience.
    const poly6 x1 = poly6( rational(1), idx6 { 1, 0, 0, 0, 0, 0 } );
    const poly6 x2 = poly6( rational(1), idx6 { 0, 1, 0, 0, 0, 0 } );
    const poly6 x3 = poly6( rational(1), idx6 { 0, 0, 1, 0, 0, 0 } );
    const poly6 y1 = poly6( rational(1), idx6 { 0, 0, 0, 1, 0, 0 } );
    const poly6 y2 = poly6( rational(1), idx6 { 0, 0, 0, 0, 1, 0 } );
    const poly6 y3 = poly6( rational(1), idx6 { 0, 0, 0, 0, 0, 1 } );

    poly6 P = 16*x1*y2*y2*y2 + rational(2,3)*x1*x2*x2;

    antiderivative<rational> F; 
    F.PRinv = P; 

    F = y1_integrate(F);
    F = y2_integrate(F);
    F = y3_integrate(F);
    F = x1_integrate(F);
    F = x2_integrate(F);
    F = x3_integrate(F);

    using float100 = boost::multiprecision::mpf_float_100;
    std::cout << std::setprecision(50)
              << F.eval<float100>( { rational(1234,1000), 1, 2, 3, 4, 5 } );
    return 0;
}

