// SF-BesselY.mjs
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
import { GSL_DBL_MAX }           from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_MAX }     from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON } from "./SF-Machine.mjs";
import { gsl_sf_sin_e }          from "./SF-Trigonometric.mjs";
import { gsl_sf_cos_e }          from "./SF-Trigonometric.mjs";
import { gsl_sf_doublefact_e }   from "./SF-Gamma.mjs";
import { gsl_sf_bessel_Ynu_asympx_e } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Ynu_asymp_Olver_e } from "./SF-BesselOlver.mjs";

import { EVAL_RESULT_D }  from "./SF-Evaluate.mjs";
import { EVAL_RESULT_ID } from "./SF-Evaluate.mjs";

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// [Abramowitz+Stegun, 10.1.3]
// with lmax=15, precision ~ 15D for x < 3
//
// checked OK [GJ] Wed May 13 15:41:25 MDT 1998 
//
function bessel_yl_small_x(l, x)
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
        const lmax  = 200;
        var t       = -0.5 * x * x;
        var sum     = 1.0;
        var t_coeff = 1.0;
        var t_power = 1.0;
        var delta1  = 0.0;

        for (let i = 1; i <= lmax; i++)
        {
            t_coeff = t_coeff / (i * (2 * (i - l) - 1));
            t_power = t_power * t;
            delta1  = t_power * t_coeff;
            sum = sum + delta1;
            if (Math.abs(delta1 / sum) < 0.5 * GSL_DBL_EPSILON) break;
        }
        r.val = -num_fact.val / den * sum;
        r.err = GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // bessel_yl_small_x

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_y0s_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (1.0 / GSL_DBL_MAX > 0.0 && x < 1.0 / GSL_DBL_MAX) // ???
    {
        throw "SF.OverflowException";
    }
    else
    {
        var cos_result = { val: 0.0, err: 0.0 }; // Result;

        cos_result = gsl_sf_cos_e(x);
        r.val = -cos_result.val / x;
        r.err = Math.abs(cos_result.err / x);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_y0s_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_y1s_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 1.0 / GSL_SQRT_DBL_MAX)
    {
        throw "SF.OverflowException";
    }
    else if (x < 0.25)
    {
        const y  = x * x;
        const c1 =  1.0 / 2.0;
        const c2 = -1.0 / 8.0;
        const c3 =  1.0 / 144.0;
        const c4 = -1.0 / 5760.0;
        const c5 =  1.0 / 403200.0;
        const c6 = -1.0 / 43545600.0;
        const sum = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * c6)))));

        r.val = -sum / y;
        r.err = GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var cos_result = { val: 0.0, err: 0.0 }; // Result;
        var sin_result = { val: 0.0, err: 0.0 }; // Result;
        var cx = 0.0;
        var sx = 0.0;

        cos_result = gsl_sf_cos_e(x);
        sin_result = gsl_sf_sin_e(x);
        cx = cos_result.val;
        sx = sin_result.val;
        r.val = -(cx / x + sx) / x;
        r.err = (Math.abs(cos_result.err / x) + sin_result.err) / Math.abs(x);
        r.err = r.err + GSL_DBL_EPSILON * (Math.abs(sx / x) + Math.abs(cx / (x * x)));
    }

    return r;

} // gsl_sf_bessel_y1s_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_y2s_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 1.0 / GSL_ROOT3_DBL_MAX)
    {
        throw "SF.OverflowException";
    }
    else if (x < 0.5)
    {
        var y  = x * x;
        var c1 =  1.0 / 6.0;
        var c2 =  1.0 / 24.0;
        var c3 = -1.0 / 144.0;
        var c4 =  1.0 / 3456.0;
        var c5 = -1.0 / 172800.0;
        var c6 =  1.0 / 14515200.0;
        var c7 = -1.0 / 1828915200.0;
        const sum = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * c7))))));

        r.val = -3.0 / (x * x * x) * sum;
        r.err = GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var cos_result = { val: 0.0, err: 0.0 }; // Result;
        var sin_result = { val: 0.0, err: 0.0 }; // Result;
        var sx = 0.0;
        var cx = 0.0;
        var a  = 0.0;

        cos_result = gsl_sf_cos_e(x);
        sin_result = gsl_sf_sin_e(x);
        sx = sin_result.val;
        cx = cos_result.val;
        a  = 3.0 / (x * x);
        r.val = (1.0 - a) / x * cx - a * sx;
        r.err = cos_result.err * Math.abs((1.0 - a) / x) + sin_result.err * Math.abs(a);
        r.err = r.err + GSL_DBL_EPSILON * (Math.abs(cx / x) + Math.abs(sx / (x * x)));
    }

    return r;

} // gsl_sf_bessel_y2s_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_yl_e(l, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    var pre = 0.0;

    if (l < 0 || x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (l == 0)
    {
        r = gsl_sf_bessel_y0s_e(x);
    }
    else if (l == 1)
    {
        r = gsl_sf_bessel_y1s_e(x);
    }
    else if (l == 2)
    {
        r = gsl_sf_bessel_y2s_e(x);
    }
    else if (x < 3.0)
    {
        r = bessel_yl_small_x(l, x);
    }
    else if (GSL_ROOT3_DBL_EPSILON * x > (l * l + l + 1))
    {
        r = gsl_sf_bessel_Ynu_asympx_e((l) + 0.5, x);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = r.val * pre;
        r.err = r.err * pre;
    }
    else if (l > 40)
    {
        r = gsl_sf_bessel_Ynu_asymp_Olver_e((l) + 0.5, x);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = r.val * pre;
        r.err = r.err * pre;
    }
    else
    {
        // recurse upward
        var r_by  = { val: 0.0, err: 0.0 }; // Result;
        var r_bym = { val: 0.0, err: 0.0 }; // Result;
        var bym = 0.0;
        var by  = 0.0;
        var byp = 0.0;

        r_by  = gsl_sf_bessel_y1s_e(x);
        r_bym = gsl_sf_bessel_y0s_e(x);
        bym = r_bym.val;
        by  = r_by.val;
        for (let j = 1; j <= l - 1; j++)
        {
            byp = (2 * j + 1) / x * by - bym;
            bym = by;
            by  = byp;
        }
        r.val = by;
        r.err = Math.abs(r.val) * (GSL_DBL_EPSILON + Math.abs(r_by.err / r_by.val) + Math.abs(r_bym.err / r_bym.val));
    }

    return r;

} // gsl_sf_bessel_yl_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_yl_array(lmax, x, result_array)
{
  
    if (lmax < 0 || x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (lmax == 0)
    {
        var r = { val: 0.0, err: 0.0 }; // Result;

        r = gsl_sf_bessel_y0s_e(x);
        result_array[0] = r.val;
    }
    else
    {
        var r_yell   = { val: 0.0, err: 0.0 }; // Result;
        var r_yellm1 = { val: 0.0, err: 0.0 }; // Result;
        var yellp1   = 0.0;
        var yell     = 0.0;
        var yellm1   = 0.0;

        r_yell   = gsl_sf_bessel_y1s_e(x);
        r_yellm1 = gsl_sf_bessel_y0s_e(x);
        yell   = r_yell.val;
        yellm1 = r_yellm1.val;
        result_array[0] = yellm1;
        result_array[1] = yell;
        
        for (let ell = 1; ell <= lmax - 1; ell++)
        {
            yellp1 = (2 * ell + 1) / x * yell - yellm1;
            result_array[ell+1] = yellp1;
            yellm1 = yell;
            yell   = yellp1;
        }
    }

} // gsl_sf_bessel_yl_array

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_y0s(x)
{
    return EVAL_RESULT_D(gsl_sf_bessel_y0s_e, x, "gsl_sf_bessel_y0s");
} // gsl_sf_bessel_y0s

export function gsl_sf_bessel_y1s(x)
{
    return EVAL_RESULT_D(gsl_sf_bessel_y1s_e, x, "gsl_sf_bessel_y1s");
} // gsl_sf_bessel_y1s

export function gsl_sf_bessel_y2s(x)
{
    return EVAL_RESULT_D(gsl_sf_bessel_y2s_e, x, "gsl_sf_bessel_y2s");
} // gsl_sf_bessel_y2s

export function gsl_sf_bessel_yl(l, x)
{
    return EVAL_RESULT_ID(gsl_sf_bessel_yl_e, { i: l, x: x }, "gsl_sf_bessel_yl");
} // gsl_sf_bessel_yl

// ----------------------------------------------------------------------------
// EOF SF-BesselY.mjs
