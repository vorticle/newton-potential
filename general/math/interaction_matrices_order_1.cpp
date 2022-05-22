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

#include <stdexcept>
#include <armadillo>
#include <geometry/cellid.h>

namespace
{

const arma::mat A[ 27 ]
{

arma::mat
{{
  {  4.60592019399116009e-02 },
  { -4.57797273965582230e-03 },
  {  5.52226285673207320e-05 },
  { -4.57797273965582230e-03 },
  {  1.36672834168741946e-03 },
  {  5.52226285673207320e-05 },
  { -4.57797273965582230e-03 },
  {  1.36672834168741946e-03 },
  {  1.36672834168741946e-03 },
  {  5.52226285673207320e-05 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 },
  { -8.28205569892252746e-03 },
  {  0.00000000000000000e+00 },
  {  6.55032436422000805e-04 },
  { -8.28205569892252746e-03 },
  {  0.00000000000000000e+00 },
  {  3.48346882347923400e-03 },
  {  6.55032436422000805e-04 }
}},

arma::mat
{{
  {  4.60592019399116009e-02 },
  {  4.57797273965582230e-03 },
  {  5.52226285673207320e-05 },
  { -4.57797273965582230e-03 },
  { -1.36672834168741946e-03 },
  {  5.52226285673207320e-05 },
  { -4.57797273965582230e-03 },
  { -1.36672834168741946e-03 },
  {  1.36672834168741946e-03 },
  {  5.52226285673207320e-05 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  { -8.28205569892252746e-03 },
  {  6.55032436422000805e-04 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 },
  { -8.28205569892252746e-03 },
  {  3.48346882347923400e-03 },
  {  0.00000000000000000e+00 },
  {  6.55032436422000805e-04 }
}},

arma::mat
{{
  {  7.80563627878484989e-02 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 },
  { -1.99442695377312906e-02 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  4.21241087025435793e-03 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  {  8.28205569892252746e-03 },
  {  6.55032436422000805e-04 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 },
  { -8.28205569892252746e-03 },
  { -3.48346882347923400e-03 },
  {  0.00000000000000000e+00 },
  {  6.55032436422000805e-04 }
}},

arma::mat
{{
  {  4.60592019399116009e-02 },
  { -4.57797273965582230e-03 },
  {  5.52226285673207320e-05 },
  {  4.57797273965582230e-03 },
  { -1.36672834168741946e-03 },
  {  5.52226285673207320e-05 },
  { -4.57797273965582230e-03 },
  {  1.36672834168741946e-03 },
  { -1.36672834168741946e-03 },
  {  5.52226285673207320e-05 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 },
  {  8.28205569892252746e-03 },
  {  0.00000000000000000e+00 },
  {  6.55032436422000805e-04 },
  { -8.28205569892252746e-03 },
  {  0.00000000000000000e+00 },
  { -3.48346882347923400e-03 },
  {  6.55032436422000805e-04 }
}},

arma::mat
{{
  {  4.60592019399116009e-02 },
  {  4.57797273965582230e-03 },
  {  5.52226285673207320e-05 },
  {  4.57797273965582230e-03 },
  {  1.36672834168741946e-03 },
  {  5.52226285673207320e-05 },
  { -4.57797273965582230e-03 },
  { -1.36672834168741946e-03 },
  { -1.36672834168741946e-03 },
  {  5.52226285673207320e-05 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  { -8.28205569892252746e-03 },
  {  6.55032436422000805e-04 },
  { -8.28205569892252746e-03 },
  {  3.48346882347923400e-03 },
  {  6.55032436422000805e-04 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 }
}},

arma::mat
{{
  {  7.80563627878484989e-02 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 },
  { -1.99442695377312906e-02 },
  {  0.00000000000000000e+00 },
  {  4.21241087025435793e-03 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  {  8.28205569892252746e-03 },
  {  6.55032436422000805e-04 },
  { -8.28205569892252746e-03 },
  { -3.48346882347923400e-03 },
  {  6.55032436422000805e-04 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 }
}},

arma::mat
{{
  {  7.80563627878484989e-02 },
  { -1.99442695377312906e-02 },
  {  4.21241087025435793e-03 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 }
}},

arma::mat
{{
  {  1.49789680899495681e-01 },
  {  0.00000000000000000e+00 },
  { -1.06825659728629584e-02 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.06825659728629584e-02 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.06825659728629584e-02 }
}},

arma::mat
{{
  {  7.80563627878484989e-02 },
  {  1.99442695377312906e-02 },
  {  4.21241087025435793e-03 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  { -8.28205569892252746e-03 },
  {  6.55032436422000805e-04 },
  {  8.28205569892252746e-03 },
  { -3.48346882347923400e-03 },
  {  6.55032436422000805e-04 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 }
}},

arma::mat
{{
  {  7.80563627878484989e-02 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 },
  {  1.99442695377312906e-02 },
  {  0.00000000000000000e+00 },
  {  4.21241087025435793e-03 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  {  8.28205569892252746e-03 },
  {  6.55032436422000805e-04 },
  {  8.28205569892252746e-03 },
  {  3.48346882347923400e-03 },
  {  6.55032436422000805e-04 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 }
}},

arma::mat
{{
  {  4.60592019399116009e-02 },
  { -4.57797273965582230e-03 },
  {  5.52226285673207320e-05 },
  { -4.57797273965582230e-03 },
  {  1.36672834168741946e-03 },
  {  5.52226285673207320e-05 },
  {  4.57797273965582230e-03 },
  { -1.36672834168741946e-03 },
  { -1.36672834168741946e-03 },
  {  5.52226285673207320e-05 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 },
  { -8.28205569892252746e-03 },
  {  0.00000000000000000e+00 },
  {  6.55032436422000805e-04 },
  {  8.28205569892252746e-03 },
  {  0.00000000000000000e+00 },
  { -3.48346882347923400e-03 },
  {  6.55032436422000805e-04 }
}},

arma::mat
{{
  {  4.60592019399116009e-02 },
  {  4.57797273965582230e-03 },
  {  5.52226285673207320e-05 },
  { -4.57797273965582230e-03 },
  { -1.36672834168741946e-03 },
  {  5.52226285673207320e-05 },
  {  4.57797273965582230e-03 },
  {  1.36672834168741946e-03 },
  { -1.36672834168741946e-03 },
  {  5.52226285673207320e-05 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  { -8.28205569892252746e-03 },
  {  6.55032436422000805e-04 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 },
  {  8.28205569892252746e-03 },
  { -3.48346882347923400e-03 },
  {  0.00000000000000000e+00 },
  {  6.55032436422000805e-04 }
}},

arma::mat
{{
  {  7.80563627878484989e-02 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -2.77581632989939639e-03 },
  {  1.99442695377312906e-02 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  {  4.21241087025435793e-03 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  {  8.28205569892252746e-03 },
  {  6.55032436422000805e-04 },
  {  0.00000000000000000e+00 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 },
  {  8.28205569892252746e-03 },
  {  3.48346882347923400e-03 },
  {  0.00000000000000000e+00 },
  {  6.55032436422000805e-04 }
}},

arma::mat
{{
  {  4.60592019399116009e-02 },
  { -4.57797273965582230e-03 },
  {  5.52226285673207320e-05 },
  {  4.57797273965582230e-03 },
  { -1.36672834168741946e-03 },
  {  5.52226285673207320e-05 },
  {  4.57797273965582230e-03 },
  { -1.36672834168741946e-03 },
  {  1.36672834168741946e-03 },
  {  5.52226285673207320e-05 }
}},

arma::mat
{{
  {  5.63802507983433218e-02 },
  {  0.00000000000000000e+00 },
  { -1.13779176424918193e-03 },
  {  8.28205569892252746e-03 },
  {  0.00000000000000000e+00 },
  {  6.55032436422000805e-04 },
  {  8.28205569892252746e-03 },
  {  0.00000000000000000e+00 },
  {  3.48346882347923400e-03 },
  {  6.55032436422000805e-04 }
}},

arma::mat
{{
  {  4.60592019399116009e-02 },
  {  4.57797273965582230e-03 },
  {  5.52226285673207320e-05 },
  {  4.57797273965582230e-03 },
  {  1.36672834168741946e-03 },
  {  5.52226285673207320e-05 },
  {  4.57797273965582230e-03 },
  {  1.36672834168741946e-03 },
  {  1.36672834168741946e-03 },
  {  5.52226285673207320e-05 }
}}

};

}

namespace math
{

template <size_t order>
const arma::mat& interaction_matrix( geometry::cellid delta );

template <>
const arma::mat& interaction_matrix<1>( geometry::cellid delta )
{
    if ( -1 > delta.i || delta.i > 1 ||
         -1 > delta.j || delta.j > 1 ||
         -1 > delta.k || delta.k > 1 )
    {
        throw std::out_of_range { "interaction_matrix<1>: id out of bounds." };
    }

    delta.i += 1; delta.j += 1; delta.k += 1;
    size_t idx = 3*3*delta.k + 3*delta.j + delta.i;
    return A[idx];
}

}
