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
#include <boxpot/rump_summation.h>

using namespace boxpot;

int main()
{
    using rational = boost::multiprecision::mpq_rational;
    using real     = boost::multiprecision::mpf_float_100;
 
    // To see the effect of round off errors,
    // you can for example instead define the types as:
    // using rational = double;
    // using real     = double;
 
    antiderivative<rational> F;
    F.PRinv = 1;

    F = x1_integrate(F);
    F = x2_integrate(F);
    F = x3_integrate(F);
    F = y1_integrate(F);
    F = y2_integrate(F);
    F = y3_integrate(F);

    // The following lines evaluate F at the 64 corners of the
    // hyper-cube [0,1]³x[0,1]³, manage the correct sign, and write the
    // resulting terms of the sum into the vector summands.

    std::vector<real> summands;
    auto y3int = [&F,&summands]( rational sign,
                                 rational x1, rational x2, rational x3,
                                 rational y1, rational y2 )
    {
        std::array<rational,6> upper { x1, x2, x3, y1, y2, 1 };
        std::array<rational,6> lower { x1, x2, x3, y1, y2, 0 };

        size_t begin = summands.size();
        F.write_summands<real>( upper, std::back_inserter(summands) );
        if ( sign < 0 )
        {
            for ( size_t i = begin; i < summands.size(); ++i )
                summands[i] = -summands[i];
        }

        begin = summands.size();
        F.write_summands<real>( lower, std::back_inserter(summands) );
        if ( sign > 0 )
        {
            for ( size_t i = begin; i < summands.size(); ++i )
                summands[i] = -summands[i];
        }
    };

    auto y2int = [&y3int]( rational sign,
                           rational x1, rational x2, rational x3,
                           rational y1 )
    {
        y3int( sign, x1, x2, x3, y1, 100 );
        y3int(-sign, x1, x2, x3, y1, 0 );
    };

    auto y1int = [&y2int]( rational sign,
                           rational x1, rational x2, rational x3 )
    {
        y2int( sign, x1, x2, x3, 1 );
        y2int(-sign, x1, x2, x3, 0 );
    };

    auto x3int = [&y1int]( rational sign,
                           rational x1, rational x2 )
    {
        y1int( sign, x1, x2, 1 );
        y1int(-sign, x1, x2, 0 );
    };

    auto x2int = [&x3int]( rational sign, rational x1 )
    {
        x3int( sign, x1, 1 );
        x3int(-sign, x1, 0 );
    };

    auto x1int = [&x2int]()
    {
        // F(100) - F(0), thus once sign 1 and once sign - 1.
        x2int( 1, 100 );
        x2int(-1,   0 );
    };

    x1int();

    real  sum = 0, l1norm = 0;
    for ( size_t i = 0; i < summands.size(); ++i )
    {
        l1norm += abs(summands[i]);
        sum    +=     summands[i] ;
    }

    // If you are using real = float, double, or long double
    // uncomment the following line to use Rump's summation algorithm.
    // sum = rump_summation( summands.begin(), summands.end() );

    real cond = l1norm/abs(sum);
    std::cout << "Condition: "  << std::setprecision(6)  << cond << ".\n";
    std::cout << "Result:    "  << std::setprecision(50) << sum  << ".\n";
    return 0;
}

