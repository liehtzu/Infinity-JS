// SF-BesselIn.mjs
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

import { GSL_IS_ODD }                      from "./SF-Math.mjs";
import { M_E }                             from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }                 from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }                 from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }                from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON }           from "./SF-Machine.mjs";
import { GSL_ROOT6_DBL_EPSILON }           from "./SF-Machine.mjs";
import { gsl_sf_bessel_IJ_taylor_e }       from "./SF-Bessel.mjs";
import { gsl_sf_bessel_I_CF1_ser }         from "./SF-Bessel.mjs";
import { gsl_sf_bessel_I0_scaled_e }       from "./SF-BesselI0.mjs";
import { gsl_sf_bessel_I1_scaled_e }       from "./SF-BesselI1.mjs";
import { gsl_sf_bessel_Inu_scaled_asymp_unif_e } from "./SF-Bessel.mjs";

import { EVAL_RESULT_DD }                  from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_In_scaled_e(n0, x)
{
    var ax = Math.abs(x);

    var n = Math.abs(n0);  // I(-n, z) = I(n, z)

    var r = { val: 0.0, err: 0.0 }; // Result;
   
    if (n == 0)
    {
        r = gsl_sf_bessel_I0_scaled_e(x);
    }
    else if (n == 1)
    {
        r = gsl_sf_bessel_I1_scaled_e(x);
    }
    else if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (x * x < 10.0 * (n + 1) / M_E)
    {
        let t  = { val: 0.0, err: 0.0 }; // Result;
        let ex = Math.exp(-ax);

        t = gsl_sf_bessel_IJ_taylor_e((n), ax, 1, 50, GSL_DBL_EPSILON);
        r.val = t.val * ex;
        r.err = t.err * ex;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        if (x < 0.0 && GSL_IS_ODD(n))
        {
            r.val = -r.val;
        }
    }
    else if (n < 150 && ax < 1.0e7)
    {
        let I0_scaled = { val: 0.0, err: 0.0 }; // Result;
        let rat  = 0.0;
        let Ikp1 = 0.0;
        let Ik   = 0.0;
        let Ikm1 = 0.0;

        I0_scaled = gsl_sf_bessel_I0_scaled_e(ax);
        rat = gsl_sf_bessel_I_CF1_ser((n), ax);
        Ikp1 = rat * GSL_SQRT_DBL_MIN;
        Ik   = GSL_SQRT_DBL_MIN;
        for (let k = n; k >= 1; k--)
        {
            Ikm1 = Ikp1 + 2.0 * (k) / ax * Ik;
            Ikp1 = Ik;
            Ik   = Ikm1;
        }
        r.val = I0_scaled.val * (GSL_SQRT_DBL_MIN / Ik);
        r.err = I0_scaled.err * (GSL_SQRT_DBL_MIN / Ik);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        if (x < 0.0 && GSL_IS_ODD(n))
        {
            r.val = -r.val;
        }
    }
    else if ( Math.min(0.29 / ((n) * (n)), 0.5 / ((n) * (n) + x * x)) < 0.5 * GSL_ROOT3_DBL_EPSILON)
    {
        r = gsl_sf_bessel_Inu_scaled_asymp_unif_e((n), ax);
        if (x < 0.0 && GSL_IS_ODD(n))
        {
            r.val = -r.val;
        }
    }
    else
    {
        let nhi = 2 + Math.trunc(1.2 / GSL_ROOT6_DBL_EPSILON);
        let r_Ikp1 = { val: 0.0, err: 0.0 }; // Result;
        let r_Ik   = { val: 0.0, err: 0.0 }; // Result;
        let Ikp1 = 0.0;
        let Ik   = 0.0;
        let Ikm1 = 0.0;

        r_Ikp1 = gsl_sf_bessel_Inu_scaled_asymp_unif_e((nhi + 1), ax);
        r_Ik = gsl_sf_bessel_Inu_scaled_asymp_unif_e((nhi), ax);
        Ikp1 = r_Ikp1.val;
        Ik   = r_Ik.val;
        for (let k = nhi; k >= n - 1; k--)
        {
            Ikm1 = Ikp1 + 2.0 * (k) / ax * Ik;
            Ikp1 = Ik;
            Ik   = Ikm1;
        }
        r.val = Ik;
        r.err = Ik * (r_Ikp1.err / r_Ikp1.val + r_Ik.err / r_Ik.val);
        if (x < 0.0 && GSL_IS_ODD(n))
        {
            r.val = -r.val;
        }
    }

    return r;

} // gsl_sf_bessel_In_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_In_scaled_array( nmin, nmax, x, result_array )
{
  
    if ( nmax < nmin || nmin < 0 )
    {
        //FOR j IN 0 .. nmax - nmin LOOP
        //    result_array(j) = 0.0;
        //END LOOP;
        throw "SF.DomainException";
    }
    else if ( x == 0.0 )
    {
        for (let j = 0; j <= nmax - nmin; j++ )
        {
            result_array[j] = 0.0;
        }
        if ( nmin == 0 )
        {
            result_array[0] = 1.0;
        }
    }
    else if ( nmax == 0 )
    {
        var I0_scaled = { val: 0.0, err: 0.0 }; // Result;

        I0_scaled = gsl_sf_bessel_I0_scaled_e( x );
        result_array[0] = I0_scaled.val;
    }
    else
    {
        var ax = Math.abs( x );
        var two_over_x = 2.0 / ax;
        var r_Inp1 = { val: 0.0, err: 0.0 }; // Result;
        var r_In   = { val: 0.0, err: 0.0 }; // Result;
        var Inp1 = 0.0;
        var In1  = 0.0;
        var Inm1 = 0.0;

        // starting values
        r_Inp1 = gsl_sf_bessel_In_scaled_e( nmax + 1, ax );
        r_In = gsl_sf_bessel_In_scaled_e( nmax,     ax );
        Inp1 = r_Inp1.val;
        In1  = r_In.val;
        for ( let n = nmax; n >= nmin; n-- )
        {
            result_array[n-nmin] = In1;
            Inm1 = Inp1 + (n) * two_over_x * In1;
            Inp1 = In1;
            In1  = Inm1;
        }
        
        // deal with signs
        if ( x < 0.0 )
        {
            for ( let n = nmin; n <= nmax; n++ )
            {
                if ( GSL_IS_ODD( n ) )
                {
                    result_array[n-nmin] = -result_array[n-nmin];
                }
            }
        }
    }

} // gsl_sf_bessel_In_scaled_array

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_In_e(n_in, x)
{
    var ax = Math.abs(x);
    var ex = 0.0;
    var n  = Math.abs(n_in);  // I(-n, z) = I(n, z)
    var In_scaled = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    In_scaled = gsl_sf_bessel_In_scaled_e(n, ax);
  
    // In_scaled is always less than 1,
    // so this overflow check is conservative.
    //
    if (ax > GSL_LOG_DBL_MAX - 1.0)
    {
        throw "SF.OverflowException";
    }
    else
    {
        ex = Math.exp(ax);
        r.val = ex * In_scaled.val;
        r.err = ex * In_scaled.err;
        r.err = r.err + ax * GSL_DBL_EPSILON * Math.abs(r.val);
        if (x < 0.0 && GSL_IS_ODD(n))
        {
            r.val = -r.val;
        }
    }

    return r;

} // gsl_sf_bessel_In_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_In_array( nmin, nmax, x, result_array )
{
    var ax  = Math.abs( x );
    var eax = Math.exp( ax );

    if ( ax > GSL_LOG_DBL_MAX - 1.0 )
    {
        //FOR j IN 0 .. nmax - nmin LOOP
        //    result_array(j) = 0.0; // FIXME: should be Inf
        //END LOOP;
        throw "SF.OverflowException";
    }
    else
    {
        gsl_sf_bessel_In_scaled_array( nmin, nmax, x, result_array );
        for ( let j = 0; j <= nmax - nmin; j++ )
        {
            result_array[j] = result_array[j] * eax;
        }
    }

} // gsl_sf_bessel_In_array

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_In_scaled( n, x )
{ // gsl_sf_bessel_In_scaled
    return EVAL_RESULT_DD( gsl_sf_bessel_In_scaled_e, { x: n, y: x }, "gsl_sf_bessel_In_scaled" );
} // gsl_sf_bessel_In_scaled;

export function gsl_sf_bessel_In( n, x )
{ // gsl_sf_bessel_In
    return EVAL_RESULT_DD( gsl_sf_bessel_In_e, { x: n, y: x }, "gsl_sf_bessel_In" );
} // gsl_sf_bessel_In;

// ----------------------------------------------------------------------------
// EOF SF-BesselIn.mjs
