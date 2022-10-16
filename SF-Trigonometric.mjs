// SF-Trigonometric.mjs
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

import { GSL_DBL_MAX }           from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }       from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON }   from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_ROOT4_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_IS_ODD }            from "./SF-Math.mjs";
import { GSL_SIGN }              from "./SF-Math.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { M_LN2 }                 from "./SF-Math.mjs";

import { cheb_eval_e }           from "./SF-Chebyshev.mjs";
import { gsl_sf_log_1plusx_e }   from "./SF-Logarithmic.mjs";
import { gsl_sf_complex_log_e }  from "./SF-Logarithmic.mjs";

import { EVAL_RESULT_D }         from "./SF-Evaluate.mjs";
import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

function ldexp(x, e)
{
    return x * Math.pow(2.0, e);
} // ldexp

// ----------------------------------------------------------------------------

// sinh(x) series
// double-precision for |x| < 1.0
function sinh_series(x)
{
    const y  = x * x;
    const c0 = 1.0 / 6.0;
    const c1 = 1.0 / 120.0;
    const c2 = 1.0 / 5040.0;
    const c3 = 1.0 / 362880.0;
    const c4 = 1.0 / 39916800.0;
    const c5 = 1.0 / 6227020800.0;
    const c6 = 1.0 / 1307674368000.0;
    const c7 = 1.0 / 355687428096000.0;

    return x * (1.0 + y * (c0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * c7))))))));
} // sinh_series

// ----------------------------------------------------------------------------

// cosh(x)-1 series
// double-precision for |x| < 1.0
function cosh_m1_series(x)
{
    const y  = x * x;
    const c0 = 0.5;
    const c1 = 1.0 / 24.0;
    const c2 = 1.0 / 720.0;
    const c3 = 1.0 / 40320.0;
    const c4 = 1.0 / 3628800.0;
    const c5 = 1.0 / 479001600.0;
    const c6 = 1.0 / 87178291200.0;
    const c7 = 1.0 / 20922789888000.0;
    const c8 = 1.0 / 6402373705728000.0;

    return y * (c0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * (c7 + y * c8))))))));
} // cosh_m1_series

// ----------------------------------------------------------------------------

// Chebyshev expansion for f(t) = sinc((t+1)/2), -1 < t < 1
const sinc_data =//: CONSTANT Series(0..16) =
    [
    1.133648177811747875422,
    -0.532677564732557348781,
    -0.068293048346633177859,
    0.033403684226353715020,
    0.001485679893925747818,
    -0.000734421305768455295,
    -0.000016837282388837229,
    0.000008359950146618018,
    0.000000117382095601192,
    -0.000000058413665922724,
    -0.000000000554763755743,
    0.000000000276434190426,
    0.000000000001895374892,
    -0.000000000000945237101,
    -0.000000000000004900690,
    0.000000000000002445383,
    0.000000000000000009925
    ];
const sinc_cs = { length: 16, c: sinc_data, order: 16, a: -1.0, b: 1.0, order_sp: 10 };


// Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
// g(x) = (sin(x)/x - 1)/(x*x)
const sin_data =//: CONSTANT Series(0..11) =
    [
    -0.3295190160663511504173,
    0.0025374284671667991990,
    0.0006261928782647355874,
    -4.6495547521854042157541e-06,
    -5.6917531549379706526677e-07,
    3.7283335140973803627866e-09,
    3.0267376484747473727186e-10,
    -1.7400875016436622322022e-12,
    -1.0554678305790849834462e-13,
    5.3701981409132410797062e-16,
    2.5984137983099020336115e-17,
    -1.1821555255364833468288e-19
    ];
const sin_cs = { length: 11, c: sin_data, order: 11, a: -1.0, b: 1.0, order_sp: 11 };

// Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
// g(x) = (2(cos(x) - 1)/(x^2) + 1) / x^2
const cos_data =//: CONSTANT Series(0..10) =
    [
    0.165391825637921473505668118136,
    -0.00084852883845000173671196530195,
    -0.000210086507222940730213625768083,
    1.16582269619760204299639757584e-6,
    1.43319375856259870334412701165e-7,
    -7.4770883429007141617951330184e-10,
    -6.0969994944584252706997438007e-11,
    2.90748249201909353949854872638e-13,
    1.77126739876261435667156490461e-14,
    -7.6896421502815579078577263149e-17,
    -3.7363121133079412079201377318e-18
    ];
const  cos_cs = { length: 10, c: cos_data, order: 10, a: -1.0, b: 1.0, order_sp: 10 };


// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

// I would have prefered just using the library sin() function.
// But after some experimentation I decided that there was
// no good way to understand the error; library sin() is just a black box.
// So we have to roll our own.
export function gsl_sf_sin_e(x)
{
    const P1 = 7.85398125648498535156e-1;
    const P2 = 3.77489470793079817668e-8;
    const P3 = 2.69515142907905952645e-15;

    const sgn_x = GSL_SIGN(x);
    const abs_x = Math.abs(x);

    var r = { val: 0.0, err: 0.0 };   // Result;

    if (abs_x < GSL_ROOT4_DBL_EPSILON)
    {
        const x2 = x * x;
        r.val = x * (1.0 - x2 / 6.0);
        r.err = Math.abs(x * x2 * x2 / 100.0);
    }
    else
    {
        var sgn_result = sgn_x;
        var y          = Math.floor(abs_x / (0.25 * M_PI));
        var octant     = Math.trunc(y - ldexp(Math.floor(ldexp(y,-3)),3));
        var z          = 0.0;

        if (GSL_IS_ODD(octant))
        {
            octant = octant + 1;
            octant = octant & 0x07;
            y = y + 1.0;
        }

        if (octant > 3)
        {
            octant = octant - 4;
            sgn_result = -sgn_result;
        }

        z = ((abs_x - y * P1) - y * P2) - y * P3;

        if (octant == 0)
        {
            var sin_cs_result = { val: 0.0, err: 0.0 };   // Result;
            const t = 8.0 * Math.abs(z) / M_PI - 1.0;
            sin_cs_result = cheb_eval_e(sin_cs, t);
            r.val = z * (1.0 + z * z * sin_cs_result.val);
        }
        else // octant == 2
        {
            var cos_cs_result = { val: 0.0, err: 0.0 };   // Result;
            const t = 8.0 * Math.abs(z) / M_PI - 1.0;
            cos_cs_result = cheb_eval_e(cos_cs, t);
            r.val = 1.0 - 0.5 * z * z * (1.0 - z * z * cos_cs_result.val);
        }

        r.val = r.val * sgn_result;

        if (abs_x > 1.0 / Number.EPSILON)
        {
            r.err = Math.abs(r.val);
        }
        else if (abs_x > 100.0 / GSL_SQRT_DBL_EPSILON)
        {
            r.err = 2.0 * abs_x * Number.EPSILON * Math.abs(r.val);
        }
        else if (abs_x > 0.1 / GSL_SQRT_DBL_EPSILON)
        {
            r.err = 2.0 * GSL_SQRT_DBL_EPSILON * Math.abs(r.val);
        }
        else
        {
            r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
        }

    }

    return r;

}// gsl_sf_sin_e

// ----------------------------------------------------------------------------

export function gsl_sf_cos_e(x)
{
    const P1 = 7.85398125648498535156e-1;
    const P2 = 3.77489470793079817668e-8;
    const P3 = 2.69515142907905952645e-15;

    var abs_x = Math.abs(x);

    var r = { val: 0.0, err: 0.0 }; // Result;
    var x2 = 0.0;
    var t  = 0.0;

    if (abs_x < GSL_ROOT4_DBL_EPSILON)
    {
        x2 = x * x;
        r.val = 1.0 - 0.5 * x2;
        r.err = Math.abs(x2 * x2 / 12.0);
    }
    else
    {
        var sgn_result = 1.0;
        var y          = Math.floor(abs_x / (0.25 * M_PI));
        var octant     = Math.trunc(y - ldexp(Math.floor(ldexp(y,-3)),3)); // Math.trunc == INTEGER
        var z          = 0.0;

        if (GSL_IS_ODD(octant))
        {
            octant = octant + 1;
            octant = octant & 0x07; // 07;
            y = y + 1.0;
        }

        if (octant > 3)
        {
            octant = octant - 4;
            sgn_result = -sgn_result;
        }

        if (octant > 1)
        {
            sgn_result = -sgn_result;
        }

        z = ((abs_x - y * P1) - y * P2) - y * P3;

        if (octant == 0)
        {
            var cos_cs_result = { val: 0.0, err: 0.0 }; // Result;
            t = 8.0 * Math.abs(z) / M_PI - 1.0;
            cos_cs_result = cheb_eval_e(cos_cs, t);
            r.val = 1.0 - 0.5 * z * z * (1.0 - z * z * cos_cs_result.val);
        }
        else // octant == 2
        {
            var sin_cs_result = { val: 0.0, err: 0.0 }; // Result;
            t = 8.0 * Math.abs(z) / M_PI - 1.0;
            sin_cs_result = cheb_eval_e(sin_cs, t);
            r.val = z * (1.0 + z * z * sin_cs_result.val);
        }

        r.val = r.val * sgn_result;

        if (abs_x > 1.0 / Number.EPSILON)
        {
            r.err = Math.abs(r.val);
        }
        else if (abs_x > 100.0 / Math.sqrt(Number.EPSILON))
        {
            r.err = 2.0 * abs_x * Number.EPSILON * Math.abs(r.val);
        }
        else if (abs_x > 0.1 / Math.sqrt(Number.EPSILON))
        {
            r.err = 2.0 * Math.sqrt(Number.EPSILON) * Math.abs(r.val);
        }
        else
        {
            r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
        }
    }

   return r;

} // gsl_sf_cos_e

// ----------------------------------------------------------------------------

export function gsl_sf_hypot_e(x, y)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0 && y == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else
    {
        const a = Math.abs(x);
        const b = Math.abs(y);
        const min = Math.min(a, b);
        const max = Math.max(a, b);
        const rat = min / max;
        const root_term = Math.sqrt(1.0 + rat * rat);

        if (max < GSL_DBL_MAX / root_term)
        {
            r.val = max * root_term;
            r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
        else
        {
            throw "SF.OverflowException";
        }
    }

    return r;

} // gsl_sf_hypot_e

// ----------------------------------------------------------------------------

export function gsl_sf_complex_sin_e(zr, zi, szr, szi)
{
    var ch_m1 = 0.0;
    var sh = 0.0;
    var ch = 0.0;
    var ex = 0.0;

    if (Math.abs(zi) < 1.0)
    {
        sh = sinh_series(zi);
        ch_m1 = cosh_m1_series(zi);
        szr.val = Math.sin(zr) * (ch_m1 + 1.0);
        szi.val = Math.cos(zr) * sh;
        szr.err = 2.0 * GSL_DBL_EPSILON * Math.abs(szr.val);
        szi.err = 2.0 * GSL_DBL_EPSILON * Math.abs(szi.val);
    }
    else if (Math.abs(zi) < GSL_LOG_DBL_MAX)
    {
        ex = Math.exp(zi);
        ch = 0.5 * (ex + 1.0 / ex);
        sh = 0.5 * (ex - 1.0 / ex);
        szr.val = Math.sin(zr) * ch;
        szi.val = Math.cos(zr) * sh;
        szr.err = 2.0 * GSL_DBL_EPSILON * Math.abs(szr.val);
        szi.err = 2.0 * GSL_DBL_EPSILON * Math.abs(szi.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

} // gsl_sf_complex_sin_e

// ----------------------------------------------------------------------------

export function gsl_sf_complex_cos_e(zr, zi, czr, czi)
{
    var ch_m1 = 0.0;
    var sh = 0.0;
    var ch = 0.0;
    var ex = 0.0;

    if (Math.abs(zi) < 1.0)
    {
        sh = sinh_series(zi);
        ch_m1 = cosh_m1_series(zi);
        czr.val =  Math.cos(zr) * (ch_m1 + 1.0);
        czi.val = -Math.sin(zr) * sh;
        czr.err = 2.0 * GSL_DBL_EPSILON * Math.abs(czr.val);
        czi.err = 2.0 * GSL_DBL_EPSILON * Math.abs(czi.val);
    }
    else if (Math.abs(zi) < GSL_LOG_DBL_MAX)
    {
        ex = Math.exp(zi);
        ch = 0.5 * (ex + 1.0 / ex);
        sh = 0.5 * (ex - 1.0 / ex);
        czr.val =  Math.cos(zr)*ch;
        czi.val = -Math.sin(zr)*sh;
        czr.err = 2.0 * GSL_DBL_EPSILON * Math.abs(czr.val);
        czi.err = 2.0 * GSL_DBL_EPSILON * Math.abs(czi.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

} // gsl_sf_complex_cos_e

// ----------------------------------------------------------------------------

export function gsl_sf_complex_logsin_e( zr, zi, lszr, lszi )
{
    var sin_r = { val: 0.0, err: 0.0 }; // Result;
    var sin_i = { val: 0.0, err: 0.0 }; // Result;
    var log_r = { val: 0.0, err: 0.0 }; // Result;
    var log_i = { val: 0.0, err: 0.0 }; // Result;

    if (zi > 60.0)
    {
        lszr.val = -M_LN2 + zi;
        lszi.val =  0.5 * M_PI - zr;
        lszr.err = 2.0 * GSL_DBL_EPSILON * Math.abs( lszr.val );
        lszi.err = 2.0 * GSL_DBL_EPSILON * Math.abs( lszi.val );
    }
    else if (zi < -60.0)
    {
        lszr.val = -M_LN2 - zi;
        lszi.val = -0.5 * M_PI + zr;
        lszr.err = 2.0 * GSL_DBL_EPSILON * Math.abs( lszr.val );
        lszi.err = 2.0 * GSL_DBL_EPSILON * Math.abs( lszi.val );
    }
    else
    {
        gsl_sf_complex_sin_e( zr, zi, sin_r, sin_i ); // ok by construction
        gsl_sf_complex_log_e( sin_r.val, sin_i.val, log_r, log_i ); //lszr, lszi);
        lszr.val = log_r.val;
        lszr.err = log_r.err;
        lszi.val = log_i.val;
        lszi.err = log_i.err;
    }
    lszi.val = gsl_sf_angle_restrict_symm_e( lszi.val );

} // gsl_sf_complex_logsin_e

// ----------------------------------------------------------------------------

export function gsl_sf_lnsinh_e(x)
{
    var eps = 0.0;
    var r   = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (Math.abs(x) < 1.0)
    {
        eps = sinh_series(x);
        r.val = Math.log(eps);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < -0.5 * GSL_LOG_DBL_EPSILON)
    {
        r.val = x + Math.log(0.5 * (1.0 - Math.exp(-2.0 * x)));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        r.val = -M_LN2 + x;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_lnsinh_e

// ----------------------------------------------------------------------------

export function gsl_sf_lncosh_e(x)
{
    var eps = 0.0;
    var r   = { val: 0.0, err: 0.0 }; // Result;

    if (Math.abs(x) < 1.0)
    {
        eps = cosh_m1_series(x);
        r = gsl_sf_log_1plusx_e(eps);
    }
    else if (x < -0.5 * GSL_LOG_DBL_EPSILON)
    {
        r.val = x + Math.log(0.5 * (1.0 + Math.exp(-2.0 * x)));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        r.val = -M_LN2 + x;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_lncosh_e

// ----------------------------------------------------------------------------

export function gsl_sf_polar_to_rect(r, theta, x, y)
{
    var t = theta;
    var c = 0.0;
    var s = 0.0;

    t = gsl_sf_angle_restrict_symm_e(t);
    c = Math.cos(t);
    s = Math.sin(t);
    x.val = r * Math.cos(t);
    y.val = r * Math.sin(t);
    x.err = r * Math.abs(s * GSL_DBL_EPSILON * t);
    x.err = x.err + 2.0 * GSL_DBL_EPSILON * Math.abs(x.val);
    y.err = r * Math.abs(c * GSL_DBL_EPSILON * t);
    y.err = y.err + 2.0 * GSL_DBL_EPSILON * Math.abs(y.val);
} // gsl_sf_polar_to_rect

// ----------------------------------------------------------------------------

// PROCEDURE gsl_sf_rect_to_polar(x: LONG_FLOAT; y: LONG_FLOAT; r: IN OUT Result; theta: IN OUT Result) IS
// BEGIN -- gsl_sf_rect_to_polar
//     r = gsl_sf_hypot_e(x, y);
//     IF (r.val > 0.0)
//         theta.val = Arctan(y, x);
//         theta.err = 2.0 * Number.EPSILON * Math.abs(theta.val);
//     else
//         RAISE SF.DomainException;
//     END IF;
// END gsl_sf_rect_to_polar;

// ----------------------------------------------------------------------------

export function gsl_sf_angle_restrict_symm_err_e(theta)
{
    // synthetic extended precision constants
    const P1 = 4.0 * 7.8539812564849853515625e-01;
    const P2 = 4.0 * 3.7748947079307981766760e-08;
    const P3 = 4.0 * 2.6951514290790594840552e-15;
    const TwoPi = 2.0 * (P1 + P2 + P3);
  
    var y = GSL_SIGN(theta) * 2.0 * Math.floor(Math.abs(theta) / TwoPi);
    var r = ((theta - y * P1) - y * P2) - y * P3;
    var d = 0.0;
    var s = { val: 0.0, err: 0.0 }; // Result;

    if (r >  M_PI)
    {
        r = (((r - 2.0 * P1) - 2.0 * P2) - 2.0 * P3); // r-TwoPi
    }
    else if (r < -M_PI)
    {
        r = (((r + 2.0 * P1) + 2.0 * P2) + 2.0 * P3); // r+TwoPi
    }

    s.val = r;

    if (Math.abs(theta) > 0.0625 / GSL_DBL_EPSILON)
    {
        //result.val = 0.0; //GSL_NAN;
        //result.err = 0.0; //GSL_NAN;
        //GSL_ERROR ("error", GSL_ELOSS);
        throw "SF.AccuracyLossException";
    }
    else if (Math.abs(theta) > 0.0625 / GSL_SQRT_DBL_EPSILON)
    {
        s.err = 2.0 * GSL_DBL_EPSILON * Math.abs(s.val - theta);
    }
    else
    {
        d = Math.abs(s.val - theta);
        if (d < M_PI)
        {
            s.err = 2.0 * GSL_DBL_EPSILON * d;
        }
        else
        {
            s.err = 2.0 * GSL_DBL_EPSILON * M_PI;
        }
    }

    return s;

} // gsl_sf_angle_restrict_symm_err_e

// ----------------------------------------------------------------------------

export function gsl_sf_angle_restrict_pos_err_e(theta)
{
    // synthetic extended precision constants
    const P1 = 4.0 * 7.85398125648498535156e-01;
    const P2 = 4.0 * 3.77489470793079817668e-08;
    const P3 = 4.0 * 2.69515142907905952645e-15;
    const TwoPi = 2.0 * (P1 + P2 + P3);
  
    const y = 2.0 * Math.floor(theta / TwoPi);
  
    var r = ((theta - y * P1) - y * P2) - y * P3;
    var d = 0.0;
    var s = { val: 0.0, err: 0.0 }; // Result;

    if (r > TwoPi)
    {
        r = (((r - 2.0 * P1) - 2.0 * P2) - 2.0 * P3); // r-TwoPi
    }
    else if (r < 0.0)
    {
        // may happen due to FP rounding
        r = (((r + 2.0 * P1) + 2.0 * P2) + 2.0 * P3); // r+TwoPi
    }

    s.val = r;

    if (Math.abs(theta) > 0.0625 / Number.EPSILON)
    {
        //s.val = 0.0; --GSL_NAN;
        //s.err = Math.abs(result.val);
        throw "SF.AccuracyLossException";
    }
    else if (Math.abs(theta) > 0.0625 / GSL_SQRT_DBL_EPSILON)
    {
        s.err = Number.EPSILON * Math.abs(s.val - theta);
    }
    else
    {
        d = Math.abs(s.val - theta);
        if (d < M_PI)
        {
            s.err = 2.0 * Number.EPSILON * d;
        }
        else
        {
            s.err = 2.0 * Number.EPSILON * M_PI;
        }
    }

    return s;

} // gsl_sf_angle_restrict_pos_err_e

// ----------------------------------------------------------------------------

export function gsl_sf_angle_restrict_symm_e(theta)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_angle_restrict_symm_err_e(theta);
    return r.val;

} // gsl_sf_angle_restrict_symm_e

// ----------------------------------------------------------------------------

export function gsl_sf_angle_restrict_pos_e(theta)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_angle_restrict_pos_err_e(theta);
    return r.val;

} // gsl_sf_angle_restrict_pos_e

// ----------------------------------------------------------------------------

export function gsl_sf_sin_err_e(x, dx)
{
    var r = { val: 0.0, err: 0.0 };   // Result;

    r = gsl_sf_sin_e(x);
    r.err = r.err + Math.abs(Math.cos(x) * dx);
    r.err = r.err + Number.EPSILON * Math.abs(r.val);
    return r;
} // gsl_sf_sin_err_e

// ----------------------------------------------------------------------------

export function gsl_sf_cos_err_e(x, dx)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    r = gsl_sf_cos_e(x);
    r.err = r.err + Math.abs(Math.sin(x) * dx);
    r.err = r.err + Number.EPSILON * Math.abs(r.val);
    return r;
} // gsl_sf_cos_err_e

// ----------------------------------------------------------------------------

export function gsl_sf_sinc_e(x)
{
    const ax = Math.abs(x);
    const q  = M_PI * ax;
    var s = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (ax < 0.8)
    {
        // Do not go to the limit of the fit since
        // there is a zero there and the Chebyshev
        // accuracy will go to zero.
        return cheb_eval_e(sinc_cs, 2.0 * ax - 1.0);
    }
    else if (ax < 100.0)
    {
        // Small arguments are no problem.
        // We trust the library sin() to
        // roughly machine precision.
        r.val = Math.sin(M_PI * ax) / (M_PI * ax);
        r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
        return r;
    }
    else
    {
        // Large arguments must be handled separately.
        s = gsl_sf_sin_e(q);
        r.val = s.val / q;
        r.err = s.err / q + 2.0 * Number.EPSILON * Math.abs(r.val);
        return r;
    }

} // gsl_sf_sinc_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_sin(x)
{ // gsl_sf_sin
    return EVAL_RESULT_D(gsl_sf_sin_e, x, "gsl_sf_sin");
} // gsl_sf_sin

export function gsl_sf_cos(x)
{ // gsl_sf_cos
    return EVAL_RESULT_D(gsl_sf_cos_e, x, "gsl_sf_cos");
} // gsl_sf_cos

export function gsl_sf_hypot(x, y)
{ // gsl_sf_hypot
    return EVAL_RESULT_DD(gsl_sf_hypot_e, { x: x, y: y }, "gsl_sf_hypot");
} // gsl_sf_hypot

export function gsl_sf_lnsinh(x)
{ // gsl_sf_lnsinh
    return EVAL_RESULT_D(gsl_sf_lnsinh_e, x, "gsl_sf_lnsinh");
} // gsl_sf_lnsinh

export function gsl_sf_lncosh(x)
{ // gsl_sf_lncosh
    return EVAL_RESULT_D(gsl_sf_lncosh_e, x, "gsl_sf_lncosh");
} // gsl_sf_lncosh

export function gsl_sf_angle_restrict_symm( theta )
{ // gsl_sf_angle_restrict_symm
    return gsl_sf_angle_restrict_symm_e( theta );
} // gsl_sf_angle_restrict_symm

export function gsl_sf_angle_restrict_pos( theta )
{
    return gsl_sf_angle_restrict_pos_e( theta );
}

export function gsl_sf_sinc(x)
{ // gsl_sf_sinc
    return EVAL_RESULT_D(gsl_sf_sinc_e, x, "gsl_sf_sinc");
} // gsl_sf_sinc

// ----------------------------------------------------------------------------
// EOF Trigonometric.mjs
