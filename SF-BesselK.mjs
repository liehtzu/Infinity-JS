// SF-BesselK.mjs
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

import { M_PI }                  from "./SF-Math.mjs";
import { M_SQRTPI }              from "./SF-Math.mjs";
import { M_SQRT2 }               from "./SF-Math.mjs";
import { GSL_IS_ODD }            from "./SF-Math.mjs";
import { GSL_SQRT_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_MAX }     from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON } from "./SF-Machine.mjs";
import { gsl_sf_doublefact_e }   from "./SF-Gamma.mjs";
import { gsl_sf_bessel_il_scaled_e }   from "./SF-BesselI.mjs";
import { gsl_sf_bessel_Knu_scaled_asympx_e } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Knu_scaled_asymp_unif_e } from "./SF-Bessel.mjs";

import { EVAL_RESULT_D }  from "./SF-Evaluate.mjs";
import { EVAL_RESULT_ID } from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

//*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// [Abramowitz+Stegun, 10.2.4 + 10.2.6]
// with lmax=15, precision ~ 15D for x < 3
//
// assumes l >= 1
//
function bessel_kl_scaled_small_x(l, x)
{
    var den      = 0.0;
    var num_fact = { val: 0.0, err: 0.0 }; // Result;
    var r        = { val: 0.0, err: 0.0 }; // Result;

    den = Math.pow(x, l + 1); // gsl_sf_pow_int(x, l+1);
    num_fact = gsl_sf_doublefact_e(2 * l - 1);
   
    if (den == 0.0)
    {
        throw "SF.OverflowException";
    }
    else
    {
        const lmax = 50;
        var ipos_term = { val: 0.0, err: 0.0 }; // Result;
        var ineg_term = 0.0;
        var sgn       = 0.0;
        var ex        = 0.0;
        var t         = 0.0;
        var sum       = 0.0;
        var t_coeff   = 0.0;
        var t_power   = 0.0;
        var delta1    = 0.0;

        if (GSL_IS_ODD(l))
        {
            sgn = -1.0;
        }
        else
        {
            sgn = 1.0;
        }
        ex  = Math.exp(x);
        t   = 0.5 * x * x;
        sum = 1.0;
        t_coeff = 1.0;
        t_power = 1.0;
        for (let i = 1; i <= lmax - 1; i++)
        {
            t_coeff = t_coeff / (i * (2 * (i - l) - 1));
            t_power = t_power * t;
            delta1 = t_power * t_coeff;
            sum = sum + delta1;
            if (Math.abs(delta1 / sum) < GSL_DBL_EPSILON) break;
        }
        
        ipos_term = gsl_sf_bessel_il_scaled_e(l, x);
        ineg_term  = sgn * num_fact.val / den * sum;
        r.val = -sgn * 0.5 * M_PI * (ex * ipos_term.val - ineg_term);
        r.val = r.val * ex;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // bessel_kl_scaled_small_x

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_k0s_scaled_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else
    {
        r.val = M_PI / (2.0 * x);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        //CHECK_UNDERFLOW(result);
    }

    return r;

} // gsl_sf_bessel_k0s_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_k1s_scaled_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < (M_SQRTPI + 1.0) / (M_SQRT2 * GSL_SQRT_DBL_MAX))
    {
        throw "SF.OverflowException";
    }
    else
    {
        r.val = M_PI / (2.0 * x) * (1.0 + 1.0 / x);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        //CHECK_UNDERFLOW(result);
    }

    return r;

 } // gsl_sf_bessel_k1s_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_k2s_scaled_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 2.0 / GSL_ROOT3_DBL_MAX)
    {
        throw "SF.OverflowException";
    }
    else
    {
        r.val = M_PI / (2.0 * x) * (1.0 + 3.0 / x * (1.0 + 1.0 / x));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        //CHECK_UNDERFLOW(result);
    }

    return r;

} // gsl_sf_bessel_k2s_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_kl_scaled_e(l, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    var pre = 0.0;

    if (l < 0 || x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (l == 0)
    {
        r = gsl_sf_bessel_k0s_scaled_e(x);
    }
    else if (l == 1)
    {
        r = gsl_sf_bessel_k1s_scaled_e(x);
    }
    else if (l == 2)
    {
        r = gsl_sf_bessel_k2s_scaled_e(x);
    }
    else if (x < 3.0)
    {
        r = bessel_kl_scaled_small_x(l, x);
    }
    else if (GSL_ROOT3_DBL_EPSILON * x > (l * l + l + 1))
    {
        r = gsl_sf_bessel_Knu_scaled_asympx_e((l) + 0.5, x);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = r.val * pre;
        r.err = r.err * pre;
    }
    else if (Math.min(0.29 / (l * l + 1), 0.5 / ((l * l + 1) + x * x)) < GSL_ROOT3_DBL_EPSILON)
    {
        r = gsl_sf_bessel_Knu_scaled_asymp_unif_e((l) + 0.5, x);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = r.val * pre;
        r.err = r.err * pre;
    }
    else
    {
        // recurse upward
        var r_bk  = { val: 0.0, err: 0.0 }; // Result;
        var r_bkm = { val: 0.0, err: 0.0 }; // Result;
        var bkp   = 0.0;
        var bk    = 0.0;
        var bkm   = 0.0;

        r_bk  = gsl_sf_bessel_k1s_scaled_e(x);
        r_bkm = gsl_sf_bessel_k0s_scaled_e(x);
        bk  = r_bk.val;
        bkm = r_bkm.val;
        for (let j = 1; j <= l - 1; j++)
        {
            bkp = (2 * j + 1) / x * bk + bkm;
            bkm = bk;
            bk  = bkp;
        }
        r.val = bk;
        r.err = Math.abs(bk) * (Math.abs(r_bk.err / r_bk.val) + Math.abs(r_bkm.err / r_bkm.val));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_kl_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_kl_scaled_array(lmax, x, result_array)
{

    if (lmax < 0 || x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (lmax == 0)
    {
        var r = { val: 0.0, err: 0.0 }; // Result;

        r = gsl_sf_bessel_k0s_scaled_e(x);
        result_array[0] = r.val;
    }
    else
    {
        var kellp1   = 0.0;
        var kell     = 0.0;
        var kellm1   = 0.0;
        var r_kell   = { val: 0.0, err: 0.0 }; // Result;
        var r_kellm1 = { val: 0.0, err: 0.0 }; // Result;

        r_kell   = gsl_sf_bessel_k1s_scaled_e(x);
        r_kellm1 = gsl_sf_bessel_k0s_scaled_e(x);
        kell   = r_kell.val;
        kellm1 = r_kellm1.val;
        result_array[0] = kellm1;
        result_array[1] = kell;
        for (let ell = 1; ell <= lmax - 1; ell++)
        {
            kellp1 = (2 * ell + 1) / x * kell + kellm1;
            result_array[ell+1] = kellp1;
            kellm1 = kell;
            kell   = kellp1;
        }
    }

} // gsl_sf_bessel_kl_scaled_array

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_k0s_scaled(x)
{ // gsl_sf_bessel_k0s_scaled
    return EVAL_RESULT_D(gsl_sf_bessel_k0s_scaled_e, x, "gsl_sf_bessel_k0s_scaled");
} // gsl_sf_bessel_k0s_scaled

export function gsl_sf_bessel_k1s_scaled(x)
{ // gsl_sf_bessel_k1s_scaled
    return EVAL_RESULT_D(gsl_sf_bessel_k1s_scaled_e, x, "gsl_sf_bessel_k1s_scaled");
} // gsl_sf_bessel_k1s_scaled

export function gsl_sf_bessel_k2s_scaled(x)
{ // gsl_sf_bessel_k2s_scaled
    return EVAL_RESULT_D(gsl_sf_bessel_k2s_scaled_e, x, "gsl_sf_bessel_k2s_scaled");
} // gsl_sf_bessel_k2s_scaled

export function gsl_sf_bessel_kl_scaled(l, x)
{ // gsl_sf_bessel_kl_scaled
    return EVAL_RESULT_ID(gsl_sf_bessel_kl_scaled_e, { i: l, x: x }, "gsl_sf_bessel_kl_scaled");
} // gsl_sf_bessel_kl_scaled

// ----------------------------------------------------------------------------
// EOF SF-BesselK.mjs
