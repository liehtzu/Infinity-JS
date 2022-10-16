// SF-BesselI.mjs
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

import { M_E }                   from "./SF-Math.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { GSL_IS_ODD }            from "./SF-Math.mjs";
import { GSL_DBL_MIN }           from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }      from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON }   from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_ROOT6_DBL_EPSILON } from "./SF-Machine.mjs";
import { gsl_sf_bessel_IJ_taylor_e } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Inu_scaled_asymp_unif_e } from "./SF-Bessel.mjs";

import { EVAL_RESULT_D }  from "./SF-Evaluate.mjs";
import { EVAL_RESULT_ID } from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// i_{l+1}/i_l
//
function bessel_il_CF1(l, x, threshold)
{
    const kmax = 2000;
    var k    = 0;
    var ak   = 0.0;
    var tk   = 0.0;
    var sum  = 0.0;
    var rhok = 0.0;
    var r    = 0.0;

    tk   = 1.0;
    sum  = 1.0;
    rhok = 0.0;
    k    = 1;
    while (k <= kmax)
    {
        ak = (x / (2 * l + 1 + 2 * k)) * (x / (2 * l + 3 + 2 * k));
        rhok = -ak * (1.0 + rhok) / (1.0 + ak * (1.0 + rhok));
        tk  = tk * rhok;
        sum = sum + tk;
        if (Math.abs(tk / sum) < threshold) break;
        k = k + 1;
    }

    r = x / (2.0 * (l) + 3.0) * sum;

    if (k >= kmax)
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // bessel_il_CF1

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_i0s_scaled_e(x)
{
    var ax = Math.abs(x);
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0)
    {
        r.val = 1.0;
        r.err = 0.0;
    }
    else if (ax < 0.2)
    {
        const eax = Math.exp(-ax);
        const y   = ax * ax;
        const c1  = 1.0 / 6.0;
        const c2  = 1.0 / 120.0;
        const c3  = 1.0 / 5040.0;
        const c4  = 1.0 / 362880.0;
        const c5  = 1.0 / 39916800.0;
        const sum = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * c5))));

        r.val = eax * sum;
        r.err = 2.0 * GSL_DBL_EPSILON * r.val;
    }
    else if (ax < -0.5 * GSL_LOG_DBL_EPSILON)
    {
        r.val = (1.0 - Math.exp(-2.0 * ax)) / (2.0 * ax);
        r.err = 2.0 * GSL_DBL_EPSILON * r.val;
    }
    else
    {
        r.val = 1.0 / (2.0 * ax);
        r.err = 2.0 * GSL_DBL_EPSILON * r.val;
    }

    return r;

} // gsl_sf_bessel_i0s_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_i1s_scaled_e(x)
{
    var ax = Math.abs(x);
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (ax < 3.0 * GSL_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (ax < 0.25)
    {
        const eax = Math.exp(-ax);
        const y   = x * x;
        const c1  = 1.0 / 10.0;
        const c2  = 1.0 / 280.0;
        const c3  = 1.0 / 15120.0;
        const c4  = 1.0 / 1330560.0;
        const c5  = 1.0 / 172972800.0;
        const sum = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * c5))));

        r.val = eax * x / 3.0 * sum;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var ex = Math.exp(-2.0 * ax);

        r.val = 0.5 * (ax * (1.0 + ex) - (1.0 - ex)) / (ax * ax);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        if (x < 0.0)
        {
            r.val = -r.val;
        }
    }

    return r;

} // gsl_sf_bessel_i1s_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_i2s_scaled_e(x)
{
    var ax = Math.abs(x);
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (ax < 4.0 * GSL_SQRT_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (ax < 0.25)
    {
        const y   = x * x;
        const c1  = 1.0 / 14.0;
        const c2  = 1.0 / 504.0;
        const c3  = 1.0 / 33264.0;
        const c4  = 1.0 / 3459456.0;
        const c5  = 1.0 / 518918400.0;
        const sum = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * c5))));
        const pre = Math.exp(-ax) * x * x / 15.0;

        r.val = pre * sum;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var ex = Math.exp(-2.0 * ax);
        var x2 = x * x;

        r.val = 0.5 * ((3.0 + x2) * (1.0 - ex) - 3.0 * ax * (1.0 + ex)) / (ax * ax * ax);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_i2s_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_il_scaled_e(l, x0)
{
    var sgn = 1.0;
    var x   = x0;
    var ax  = Math.abs(x);
    var pre = x0;
    var il  = { val: 0.0, err: 0.0 }; // Result;
    var b   = { val: 0.0, err: 0.0 }; // Result;
    var r   = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0)
    {
        // i_l(-x) = (-1)^l i_l(x)
        if (GSL_IS_ODD(l))
        {
            sgn = -1.0;
        }
        else
        {
            sgn =  1.0;
        }
        x = -x;
    }

    if (l < 0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        if (l == 0)
        {
            r.val = 1.0;
        }
        else
        {
            r.val = 0.0;
        }
        r.err = 0.0;
    }
    else if (l == 0)
    {
        il = gsl_sf_bessel_i0s_scaled_e(x);
        r.val = sgn * il.val;
        r.err = il.err;
    }
    else if (l == 1)
    {
        il = gsl_sf_bessel_i1s_scaled_e(x);
        r.val = sgn * il.val;
        r.err = il.err;
    }
    else if (l == 2)
    {
        il = gsl_sf_bessel_i2s_scaled_e(x);
        r.val = sgn * il.val;
        r.err = il.err;
    }
    else if (x * x < 10.0 * ((l) + 1.5) / M_E)
    {
        b = gsl_sf_bessel_IJ_taylor_e((l) + 0.5, x, 1, 50, GSL_DBL_EPSILON);
        pre = Math.exp(-ax) * Math.sqrt((0.5 * M_PI) / x);
        r.val = sgn * pre * b.val;
        r.err = pre * b.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (l < 150)
    {
        let i0_scaled = { val: 0.0, err: 0.0 }; // Result;
        let rat       = 0.0;
        let iellp1    = 0.0;
        let iell      = 0.0;
        let iellm1    = 0.0;

        i0_scaled = gsl_sf_bessel_i0s_scaled_e(ax);
        rat = bessel_il_CF1(l, ax, GSL_DBL_EPSILON);
        iellp1 = rat * GSL_SQRT_DBL_MIN;
        iell   = GSL_SQRT_DBL_MIN;
        for (let ell = l; ell >= 1; ell--)
        {
            iellm1 = iellp1 + (2 * ell + 1) / x * iell;
            iellp1 = iell;
            iell   = iellm1;
        }
        r.val = sgn * i0_scaled.val * (GSL_SQRT_DBL_MIN / iell);
        r.err = i0_scaled.err * (GSL_SQRT_DBL_MIN / iell);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (Math.min(0.29 / (l * l + 1), 0.5 / ((l * l + 1) + x * x)) < 0.5 * GSL_ROOT3_DBL_EPSILON)
    {
        r = gsl_sf_bessel_Inu_scaled_asymp_unif_e((l) + 0.5, x);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = r.val * sgn * pre;
        r.err = r.err * pre;
    }
    else
    {
        // recurse down from safe values
        let rt_term  = Math.sqrt((0.5 * M_PI) / x);
        const LMAX   = 2 + Math.trunc(1.2 / GSL_ROOT6_DBL_EPSILON);
        let r_iellp1 = { val: 0.0, err: 0.0 }; // Result;
        let r_iell   = { val: 0.0, err: 0.0 }; // Result;
        let iellp1   = 0.0;
        let iell     = 0.0;
        let iellm1   = 0.0;

        r_iellp1 = gsl_sf_bessel_Inu_scaled_asymp_unif_e((LMAX + 1) + 0.5, x);
        r_iell = gsl_sf_bessel_Inu_scaled_asymp_unif_e((LMAX    ) + 0.5, x);
        iellp1 = r_iellp1.val;
        iell   = r_iell.val;
        iellm1 = 0.0;
        iellp1 = iellp1 *rt_term;
        iell   = iell * rt_term;
        for (let ell = LMAX; ell >= l + 1; ell--)
        {
            iellm1 = iellp1 + (2 * ell + 1) / x * iell;
            iellp1 = iell;
            iell   = iellm1;
        }
        r.val = sgn * iellm1;
        r.err = Math.abs(r.val) * (GSL_DBL_EPSILON + Math.abs(r_iellp1.err / r_iellp1.val) + Math.abs(r_iell.err / r_iell.val));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_il_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_il_scaled_array(lmax, x, result_array)
{

    if (x == 0.0)
    {
        result_array[0] = 1.0;
        for (let ell = lmax; ell >= 1; ell--)
        {
            result_array[ell] = 0.0;
        }
    }
    else
    {
        var r_iellp1 = { val: 0.0, err: 0.0 }; // Result;
        var r_iell   = { val: 0.0, err: 0.0 }; // Result;
        var iellp1   = 0.0;
        var iell     = 0.0;
        var iellm1   = 0.0;

        r_iellp1 = gsl_sf_bessel_il_scaled_e(lmax + 1, x);
        r_iell   = gsl_sf_bessel_il_scaled_e(lmax,     x);
        iellp1 = r_iellp1.val;
        iell   = r_iell.val;
        result_array[lmax] = iell;
        for (let ell = lmax; ell >= 1; ell--)
        {
            iellm1 = iellp1 + (2 * ell + 1) / x * iell;
            iellp1 = iell;
            iell   = iellm1;
            result_array[ell-1] = iellm1;
        }
    }

} // gsl_sf_bessel_il_scaled_array

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_i0s_scaled(x)
{ // gsl_sf_bessel_i0s_scaled
    return EVAL_RESULT_D(gsl_sf_bessel_i0s_scaled_e, x, "gsl_sf_bessel_i0s_scaled");
} // gsl_sf_bessel_i0s_scaled

export function gsl_sf_bessel_i1s_scaled(x)
{ // gsl_sf_bessel_i1s_scaled
    return EVAL_RESULT_D(gsl_sf_bessel_i1s_scaled_e, x, "gsl_sf_bessel_i1s_scaled");
} // gsl_sf_bessel_i1s_scaled

export function gsl_sf_bessel_i2s_scaled(x)
{ // gsl_sf_bessel_i2s_scaled
    return EVAL_RESULT_D(gsl_sf_bessel_i2s_scaled_e, x, "gsl_sf_bessel_i2s_scaled");
} // gsl_sf_bessel_i2s_scaled

export function gsl_sf_bessel_il_scaled(l, x)
{ // gsl_sf_bessel_il_scaled
    return EVAL_RESULT_ID(gsl_sf_bessel_il_scaled_e, { i: l, x: x }, "gsl_sf_bessel_il_scaled");
} // gsl_sf_bessel_il_scaled

// ----------------------------------------------------------------------------
// EOF SF-BesselI.mjs
