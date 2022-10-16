// SF-BesselKn.mjs
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

import { M_EULER }                         from "./SF-Math.mjs";
import { GSL_IS_ODD }                      from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }                 from "./SF-Machine.mjs";
import { GSL_DBL_MAX }                     from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }                 from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON }           from "./SF-Machine.mjs";
import { gsl_sf_bessel_K0_scaled_e }       from "./SF-BesselK0.mjs";
import { gsl_sf_bessel_K1_scaled_e }       from "./SF-BesselK1.mjs";
import { gsl_sf_fact_e }                   from "./SF-Gamma.mjs";
import { gsl_sf_lnfact_e }                 from "./SF-Gamma.mjs";
import { gsl_sf_psi_int_e }                from "./SF-Psi.mjs";
import { gsl_sf_bessel_Knu_scaled_asympx_e } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Knu_scaled_asymp_unif_e } from "./SF-Bessel.mjs";

import { EVAL_RESULT_ID }                  from "./SF-Evaluate.mjs";

//*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// [Abramowitz+Stegun, 9.6.11]
// assumes n >= 1
//
function bessel_Kn_scaled_small_x(n, x)
{
    var y       = 0.25 * x * x;
    var ln_x_2  = Math.log(0.5 * x);
    var ex      = Math.exp(x);
    var k_term  = 0.0;
    var term1   = 0.0;
    var sum1    = 0.0;
    var ln_pre1 = 0.0;
    var term2   = 0.0;
    var sum2    = 0.0;
    var pre2    = 0.0;
    var k       = 0;
    var ln_nm1_fact = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    ln_nm1_fact = gsl_sf_lnfact_e(n - 1);
   
    ln_pre1 = -(n) * ln_x_2 + ln_nm1_fact.val;
    if (ln_pre1 > GSL_LOG_DBL_MAX - 3.0)
    {
        throw "SF.OverflowException";
    }
   
    sum1 = 1.0;
    k_term = 1.0;
    for (k = 1; k <= n - 1; k++)
    {
        k_term = k_term * (-y / (k * (n - k)));
        sum1 = sum1 + k_term;
    }
    term1 = 0.5 * Math.exp(ln_pre1) * sum1;
   
    pre2 = 0.5 * Math.exp((n) * ln_x_2);
    if (pre2 > 0.0)
    {
        const KMAX = 20;
        var psi_n = { val: 0.0, err: 0.0 }; // Result;
        var npk_fact = { val: 0.0, err: 0.0 }; // Result;
        var yk        = 1.0;
        var k_fact    = 1.0;
        var psi_kp1   = -M_EULER;
        var psi_npkp1 = 0.0;

        psi_n = gsl_sf_psi_int_e(n);
        npk_fact = gsl_sf_fact_e(n);
        psi_npkp1 = psi_n.val + 1.0 / (n);
        sum2 = (psi_kp1 + psi_npkp1 - 2.0 * ln_x_2) / npk_fact.val;
        for (k = 1; k <= KMAX - 1; k++)
        {
            psi_kp1   = psi_kp1 + 1.0 / (k);
            psi_npkp1 = psi_npkp1 + 1.0 / (n + k);
            k_fact    = k_fact * (k);
            npk_fact.val = npk_fact.val * (n + k);
            yk = yk * y;
            k_term = yk * (psi_kp1 + psi_npkp1 - 2.0 * ln_x_2) / (k_fact * npk_fact.val);
            sum2 = sum2 + k_term;
        }
        term2 = pre2 * sum2;
        if (GSL_IS_ODD(n))
        {
            term2 = -term2;
        }
    }
    else
    {
        term2 = 0.0;
    }
   
    r.val = ex * (term1 + term2);
    r.err = ex * GSL_DBL_EPSILON * (Math.abs(ln_pre1) * Math.abs(term1) + Math.abs(term2));
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
   
    return r;

} // bessel_Kn_scaled_small_x

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Kn_scaled_e(n0, x)
{
    var n = Math.abs(n0); // K(-n, z) = K(n, z)
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (n == 0)
    {
        r = gsl_sf_bessel_K0_scaled_e(x);
    }
    else if (n == 1)
    {
        r = gsl_sf_bessel_K1_scaled_e(x);
    }
    else if (x <= 5.0)
    {
        r = bessel_Kn_scaled_small_x(n, x);
    }
    else if (GSL_ROOT3_DBL_EPSILON * x > 0.25 * (n * n + 1))
    {
        r = gsl_sf_bessel_Knu_scaled_asympx_e((n), x);
    }
    else if (Math.min(0.29 / (n * n), 0.5 / ((n * n) + x * x)) < GSL_ROOT3_DBL_EPSILON)
    {
        r = gsl_sf_bessel_Knu_scaled_asymp_unif_e((n), x);
    }
    else
    {
        // Upward recurrence. [Gradshteyn + Ryzhik, 8.471.1]
        var two_over_x = 2.0 / x;
        var r_b_jm1 = { val: 0.0, err: 0.0 }; // Result;
        var r_b_j   = { val: 0.0, err: 0.0 }; // Result;
        var b_jm1   = 0.0;
        var b_j     = 0.0;
        var b_jp1   = 0.0;
        var j       = 0;

        r_b_jm1 = gsl_sf_bessel_K0_scaled_e(x);
        r_b_j = gsl_sf_bessel_K1_scaled_e(x);
        b_jm1 = r_b_jm1.val;
        b_j   = r_b_j.val;
        for (j = 1; j <= n - 1; j++)
        {
            b_jp1 = b_jm1 + (j) * two_over_x * b_j;
            b_jm1 = b_j;
            b_j   = b_jp1; 
        }
        
        r.val = b_j;
        r.err = (n) * (Math.abs(b_j) * (Math.abs(r_b_jm1.err / r_b_jm1.val) + Math.abs(r_b_j.err / r_b_j.val)));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_Kn_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Kn_e( n, x )
{
    var ex = 0.0;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_bessel_Kn_scaled_e( n, x );
    ex = Math.exp( -x );
    r.val = r.val * ex;
    r.err = r.err * ex;
    r.err = r.err + x * GSL_DBL_EPSILON * Math.abs( r.val );
    return r;

} // gsl_sf_bessel_Kn_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Kn_scaled_array( nmin, nmax, x, result_array )
{
  
    if ( nmin < 0 || nmax < nmin || x <= 0.0 )
    {
        //FOR j IN 0 .. nmax - nmin LOOP
        //    result_array(j) = 0.0;
        //END LOOP;
        throw "SF.DomainException";
    }
    else if ( nmax == 0 )
    {
        var b = { val: 0.0, err: 0.0 }; // Result;

        b = gsl_sf_bessel_K0_scaled_e( x );
        result_array[0] = b.val;
    }
    else
    {
        var two_over_x = 0.0;
        var r_Knm1 = { val: 0.0, err: 0.0 }; // Result;
        var r_Kn   = { val: 0.0, err: 0.0 }; // Result;
        var Knp1   = 0.0;
        var Kn     = 0.0;
        var Knm1   = 0.0;

        two_over_x = 2.0 / x;
        r_Knm1 = gsl_sf_bessel_Kn_scaled_e( nmin,     x );
        r_Kn = gsl_sf_bessel_Kn_scaled_e( nmin + 1, x );
        Kn   = r_Kn.val;
        Knm1 = r_Knm1.val;
        for ( let n = nmin + 1; n <= nmax + 1; n++ )
        {
            if ( Knm1 < GSL_DBL_MAX )
            {
                result_array[n-1-nmin] = Knm1;
                Knp1 = Knm1 + (n) * two_over_x * Kn;
                Knm1 = Kn;
                Kn   = Knp1;
            }
            else
            {
                // Overflow. Set the rest of the elements to
                // zero and bug out.
                // FIXME: Note: this relies on the convention
                // that the test x < DBL_MIN fails for x not
                // a number. This may be only an IEEE convention,
                // so the portability is unclear.
                //
                for ( let j = n; j <= nmax + 1; j++ )
                {
                    result_array[j-1-nmin] = 0.0;
                }
                throw "SF.OverflowException";
            }
        }
    }

} // gsl_sf_bessel_Kn_scaled_array

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Kn_array( nmin, nmax, x, result_array )
{
    var ex = 0.0;

    gsl_sf_bessel_Kn_scaled_array( nmin, nmax, x, result_array );
    ex = Math.exp( -x );
    for ( let i = 0; i <= nmax - nmin; i++ )
    {
        result_array[i] = result_array[i] * ex;
    }

} // gsl_sf_bessel_Kn_array

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Kn_scaled( n, x )
{ // gsl_sf_bessel_Kn_scaled
    return EVAL_RESULT_ID( gsl_sf_bessel_Kn_scaled_e, { i: n, x: x }, "gsl_sf_bessel_Kn_scaled" );
} // gsl_sf_bessel_Kn_scaled;

export function gsl_sf_bessel_Kn( n, x )
{ // gsl_sf_bessel_Kn
    return EVAL_RESULT_ID( gsl_sf_bessel_Kn_e, { i: n, x: x }, "gsl_sf_bessel_Kn" );
} // gsl_sf_bessel_Kn;

// ----------------------------------------------------------------------------
// EOF SF-BesselKn.mjs
