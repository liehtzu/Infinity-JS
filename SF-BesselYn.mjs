// SF-BesselYn.mjs
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

import { gsl_sf_bessel_Y0_e }              from "./SF-BesselY0.mjs";
import { gsl_sf_bessel_Y1_e }              from "./SF-BesselY1.mjs";
import { GSL_IS_ODD }                      from "./SF-Math.mjs";
import { M_EULER }                         from "./SF-Math.mjs";
import { M_PI }                            from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }                 from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }                 from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON }           from "./SF-Machine.mjs";
import { gsl_sf_bessel_Ynu_asympx_e }      from "./SF-Bessel.mjs";
import { gsl_sf_lnfact_e }                 from "./SF-Gamma.mjs";
import { gsl_sf_fact_e }                   from "./SF-Gamma.mjs";
import { gsl_sf_psi_int_e }                from "./SF-Psi.mjs";
import { gsl_sf_bessel_Ynu_asymp_Olver_e } from "./SF-BesselOlver.mjs";

import { EVAL_RESULT_DD }                  from "./SF-Evaluate.mjs";

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// assumes n >= 1
function bessel_Yn_small_x(n, x)
{
    var y       = 0.25 * x * x;
    var ln_x_2  = Math.log(0.5 * x);
    var k_term  = 0.0;
    var term1   = 0.0;
    var sum1    = 0.0;
    var ln_pre1 = 0.0;
    var term2   = 0.0;
    var sum2    = 0.0;
    var pre2    = 0.0;
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
    for (let k = 1; k <= n - 1; k++)
    {
        k_term = k_term * (y / (k * (n - k)));
        sum1 = sum1 + k_term;
    }
    term1 = -Math.exp(ln_pre1) * sum1 / M_PI;
    
    pre2 = -Math.exp((n) * ln_x_2) / M_PI;
    if (Math.abs(pre2) > 0.0)
    {
        const KMAX      = 20;
        var psi_n     = { val: 0.0, err: 0.0 }; // Result;
        var npk_fact  = { val: 0.0, err: 0.0 }; // Result;
        var yk        = 1.0;
        var k_fact    = 1.0;
        var psi_kp1   = -M_EULER;
        var psi_npkp1 = 0.0;

        psi_n = gsl_sf_psi_int_e(n);
        npk_fact = gsl_sf_fact_e(n);
        psi_npkp1 = psi_n.val + 1.0 / (n);
        sum2 = (psi_kp1 + psi_npkp1 - 2.0 * ln_x_2) / npk_fact.val;
        for (let k = 1; k <= KMAX - 1; k++)
        {
            psi_kp1   = psi_kp1 + 1.0 / (k);
            psi_npkp1 = psi_npkp1 + 1.0/ (n + k); 
            k_fact    = k_fact * (k);
            npk_fact.val = npk_fact.val * (n + k);
            yk = -yk * y;
            k_term = yk * (psi_kp1 + psi_npkp1 - 2.0 * ln_x_2) / (k_fact * npk_fact.val);
            sum2 = sum2 + k_term;
        }
        term2 = pre2 * sum2;
    }
    else
    {
        term2 = 0.0;
    }

    r.val = term1 + term2;
    r.err = GSL_DBL_EPSILON * (Math.abs(ln_pre1) * Math.abs(term1) + Math.abs(term2));
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
   
    return r;

} // bessel_Yn_small_x

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Yn_e(n0, x)
{
    var sign = 1;
    var n    = n0;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    if (n < 0)
    {
        // reduce to case n >= 0
        n = -n;
        if (GSL_IS_ODD(n))
        {
            sign = -1;
        }
    }
   
    if (n == 0)
    {
        r = gsl_sf_bessel_Y0_e(x);
        r.val = r.val * (sign);
    }
    else if (n == 1)
    {
        r = gsl_sf_bessel_Y1_e(x);
        r.val = r.val * (sign);
    }
    else
    {
        if (x <= 0.0)
        {
            throw "SF.DomainException";
        }
        if (x < 5.0)
        {
            r = bessel_Yn_small_x(n, x);
            r.val = r.val * (sign);
        }
        else if (GSL_ROOT3_DBL_EPSILON * x > (n * n + 1))
        {
            r = gsl_sf_bessel_Ynu_asympx_e((n), x);
            r.val = r.val * (sign);
        }
        else if (n > 50)
        {
            r = gsl_sf_bessel_Ynu_asymp_Olver_e((n), x);
            r.val = r.val * (sign);
        }
        else
        {
            var two_over_x = 2.0 /x;
            var r_by  = { val: 0.0, err: 0.0 }; // Result; 
            var r_bym = { val: 0.0, err: 0.0 }; // Result;
            var bym = 0.0;
            var by  = 0.0;
            var byp = 0.0;
            var j = 0;

            r_by  = gsl_sf_bessel_Y1_e(x);
            r_bym = gsl_sf_bessel_Y0_e(x);
            bym = r_bym.val;
            by  = r_by.val;
            for (j = 1; j <= n - 1; j++)
            {
                byp = (j) * two_over_x * by - bym;
                bym = by;
                by  = byp;
            }
            r.val = (sign) * by;
            r.err = Math.abs(r.val) * (Math.abs(r_by.err / r_by.val) + Math.abs(r_bym.err / r_bym.val));
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
    }

    return r;

} // gsl_sf_bessel_Yn_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Yn_array( nmin, nmax, x, result_array )
{
  
    if ( nmin < 0 || nmax < nmin || x <= 0.0 )
    {
        //FOR j IN 0 .. nmax - nmin LOOP
        //    result_array(j) = 0.0;
        //END LOOP;
        throw "SF.DomainException";
    }
    else
    {
        var r_Ynm1   = { val: 0.0, err: 0.0 }; // Result;
        var r_Yn     = { val: 0.0, err: 0.0 }; // Result;
        var Ynp1     = 0.0;
        var Yn       = 0.0;
        var Ynm1     = 0.0;

        r_Ynm1 = gsl_sf_bessel_Yn_e( nmin,   x );
        r_Yn = gsl_sf_bessel_Yn_e( nmin + 1, x );
        Yn   = r_Yn.val;
        Ynm1 = r_Ynm1.val;
        for ( let n = nmin + 1; n <= nmax + 1; n++ )
        {
            result_array[n-nmin-1] = Ynm1;
            Ynp1 = -Ynm1 + 2.0 * (n) / x * Yn;
            Ynm1 = Yn;
            Yn   = Ynp1;
        }
    }

} // gsl_sf_bessel_Yn_array

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Yn( n, x )
{ // gsl_sf_bessel_Yn
    return EVAL_RESULT_DD( gsl_sf_bessel_Yn_e, { x: n, y: x }, "gsl_sf_bessel_Yn" );
 } gsl_sf_bessel_Yn

// ----------------------------------------------------------------------------
// EOF SF-BesselYn.mjs
