// SF-BesselJn.mjs
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

import { gsl_sf_bessel_J0_e }              from "./SF-BesselJ0.mjs";
import { gsl_sf_bessel_J1_e }              from "./SF-BesselJ1.mjs";
import { GSL_IS_ODD }                      from "./SF-Math.mjs";
import { GSL_SQRT_DBL_MIN }                from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }                 from "./SF-Machine.mjs";
import { GSL_ROOT4_DBL_EPSILON }           from "./SF-Machine.mjs";
import { GSL_ROOT5_DBL_EPSILON }           from "./SF-Machine.mjs";
import { gsl_sf_bessel_J_CF1 }             from "./SF-Bessel.mjs";
import { gsl_sf_bessel_IJ_taylor_e }       from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Jnu_asympx_e }      from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Jnu_asymp_Olver_e } from "./SF-BesselOlver.mjs";

import { EVAL_RESULT_DD }                  from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Jn_e(n0, x0)
{
    var n    = n0;
    var sign = 1;

    var x = x0;

    var b  = { val: 0.0, err: 0.0 }; // Result;
    var b0 = { val: 0.0, err: 0.0 }; // Result;
    var b1 = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;
   
    if (n < 0)
    {
        // reduce to case n >= 0
        n = -n;
        if (GSL_IS_ODD(n))
        {
            sign = -sign;
        }
    }
   
    if (x < 0.0)
    {
        // reduce to case x >= 0.
        x = -x;
        if (GSL_IS_ODD(n))
        {
            sign = -sign;
        }
    }
   
    if (n == 0)
    {//console.log("*** 1");
        b0 = gsl_sf_bessel_J0_e(x);
        r.val = (sign) * b0.val;
        r.err = b0.err;
    }
    else if (n == 1)
    {//console.log("*** 2");
        b1 = gsl_sf_bessel_J1_e(x);
        r.val = (sign) * b1.val;
        r.err = b1.err;
    }
    else
    {//console.log("*** 3");
        if (x == 0.0)
        {//console.log("*** 3.1");
            r.val = 0.0;
            r.err = 0.0;
        }
        else if (x * x < 10.0 * ((n) + 1.0) * GSL_ROOT5_DBL_EPSILON)
        {//console.log("*** 3.2");
            b = gsl_sf_bessel_IJ_taylor_e((n), x, -1, 50, GSL_DBL_EPSILON);
            r.val = (sign) * b.val;
            r.err = b.err;
            r.err = r.err + GSL_DBL_EPSILON * Math.abs(r.val);
        }
        else if (GSL_ROOT4_DBL_EPSILON * x > (n * n + 1))
        {//console.log("*** 3.3");
            r = gsl_sf_bessel_Jnu_asympx_e((n), x);
            r.val = r.val * (sign);
        }
        else if (n > 50)
        {//console.log("*** 3.4");
            r = gsl_sf_bessel_Jnu_asymp_Olver_e((n), x);
            r.val = r.val * (sign);
        }
        else if (x > 1000.0)
        {//console.log("*** 3.5");
            // We need this to avoid feeding large x to CF1; note that
            // due to the above check, we know that n <= 50.
            //
            r = gsl_sf_bessel_Jnu_asympx_e((n), x);
            r.val = r.val * (sign);
        }
        else
        {//console.log("*** 3.6 | n=" + n);
            var ans   = 0.0;
            var err   = 0.0;
            var ratio = 0.0;
            var sgn   = 0.0;
            var Jkp1  = 0.0;
            var Jk    = 0.0;
            var Jkm1  = 0.0;
            var rs = {};

            rs = gsl_sf_bessel_J_CF1((n), x); //, ratio, sgn);
            ratio = rs.ratio;
            sgn = rs.sgn;
            
            // backward recurrence
            Jkp1 = GSL_SQRT_DBL_MIN * ratio;
            Jk   = GSL_SQRT_DBL_MIN;
            
            for (let k = n; k >= 1; k--)
            {
                Jkm1 = 2.0 * (k) / x * Jk - Jkp1;
                Jkp1 = Jk;
                Jk   = Jkm1;
            }
            
            if (Math.abs(Jkp1) > Math.abs(Jk))
            {
                b1 = gsl_sf_bessel_J1_e(x);
                ans = b1.val / Jkp1 * GSL_SQRT_DBL_MIN;
                err = b1.err / Jkp1 * GSL_SQRT_DBL_MIN;
            }
            else
            {
                b0 = gsl_sf_bessel_J0_e(x);
                ans = b0.val / Jk * GSL_SQRT_DBL_MIN;
                err = b0.err / Jk * GSL_SQRT_DBL_MIN;
            }
            
            r.val = (sign) * ans;
            r.err = Math.abs(err);
        }
    }

    return r;

} // gsl_sf_bessel_Jn_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Jn_array( nmin, nmax, x, result_array )
{
  
    if ( nmin < 0 || nmax < nmin )
    {
        //FOR n IN REVERSE nmin .. nmax LOOP
        //    result_array(n-nmin) = 0.0;
        //END LOOP;
        throw "SF.DomainException"; // WITH "domain error";
    }
    else if ( x == 0.0 )
    {
        for ( let n = nmax; n >= nmin; n-- )
        {
            result_array[n-nmin] = 0.0;
        }
        if ( nmin == 0 )
        {
            result_array[0] = 1.0;
        }
    }
    else
    {
        var r_Jnp1 = { val: 0.0, err: 0.0 }; // Result;
        var r_Jn   = { val: 0.0, err: 0.0 }; // Result;
        var Jnp1   = 0.0;
        var Jn     = 0.0;
        var Jnm1   = 0.0;

        r_Jnp1 = gsl_sf_bessel_Jn_e( nmax + 1, x );
        r_Jn = gsl_sf_bessel_Jn_e( nmax, x );
        Jnp1 = r_Jnp1.val;
        Jn   = r_Jn.val;
        for ( let n = nmax; n >= nmin; n-- )
        {
            result_array[n-nmin] = Jn;
            Jnm1 = -Jnp1 + 2.0 * (n) / x * Jn;
            Jnp1 = Jn;
            Jn   = Jnm1;
        }
    }

} // gsl_sf_bessel_Jn_array

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Jn( n, x )
{ // gsl_sf_bessel_Jn
    return EVAL_RESULT_DD( gsl_sf_bessel_Jn_e, { x: n, y: x }, "gsl_sf_bessel_Jn" );
} // gsl_sf_bessel_Jn

// ----------------------------------------------------------------------------
// EOF SF-BesselJn.mjs
