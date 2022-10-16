// SF-Power.mjs
// ----------------------------------------------------------------------------
// Copyright (C) 2007 Brian Gough
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2004 Gerard Jungman
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
// ----------------------------------------------------------------------------

// Author:  G. Jungman
// Adaptation to JavaScript: Igor Izvarin

import { GSL_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_IS_ODD }       from "./SF-Math.mjs";

import { EVAL_RESULT_DI }   from "./SF-Evaluate.mjs";


//*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error handling *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_pow_int_e(x0, n0)
{
    var n     = n0;
    var x     = x0;
    var u     = 0.0;
    var value = 1.0;
    var count = 0;
    var r     = { val: 0.0, err: 0.0 }; // Result;
   
    if (n < 0)
    {
        n = -n;
       
        if (x == 0.0)
        {
            u = 1.0 / x;
            if (n % 2 != 0) // correct sign of infinity
            {
                r.val = u;
            }
            else
            {
                r.val = u * u;
            }
            //result.err = 0.0; --GSL_POSINF;
            throw "SF.OverflowException";
        }
       
        x = 1.0 / x;
    }
   
    // repeated squaring method 
    // returns 0.0^0 = 1.0, so continuous in x
    //
    while (true)
    {
        if (GSL_IS_ODD(n))
        {
            value = value * x;
        }
        n = Math.trunc(n / 2);
        x = x * x;
        count = count + 1;
        if (n == 0)
        {
            break;
        }
    }
   
    r.val = value;
    r.err = 2.0 * GSL_DBL_EPSILON * (count + 1) * Math.abs(value); 
   
    return r;

} // gsl_sf_pow_int_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-

export function gsl_sf_pow_int( x, n )
{ // gsl_sf_pow_int
    return EVAL_RESULT_DI( gsl_sf_pow_int_e, { x: x, i: n }, "gsl_sf_pow_int" );
} // gsl_sf_pow_int

// ----------------------------------------------------------------------------
// EOF SF-Power.mjs
