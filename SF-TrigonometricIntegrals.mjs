// SF-TrigonometricIntegrals.mjs
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

import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_DBL_MIN }          from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }     from "./SF-Machine.mjs";
import { cheb_eval_e }          from "./SF-Chebyshev.mjs";
import { M_PI }                 from "./SF-Math.mjs";
import { GSL_SIGN }             from "./SF-Math.mjs";
import { gsl_sf_sin_e }         from "./SF-Trigonometric.mjs";
import { gsl_sf_cos_e }         from "./SF-Trigonometric.mjs";


// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// based on SLATEC r9sifg.f, W. Fullerton
//
// series for f1   on the interval  2.00000e-02 to  6.25000e-02
//                                        with weighted error   2.82e-17
//                                         log weighted error  16.55
//                               significant figures required  15.36
//                                    decimal places required  17.20

const f1_data =//: CONSTANT Series(0..19) := --[20]
    [
   -0.1191081969051363610,
   -0.0247823144996236248,
    0.0011910281453357821,
   -0.0000927027714388562,
    0.0000093373141568271,
   -0.0000011058287820557,
    0.0000001464772071460,
   -0.0000000210694496288,
    0.0000000032293492367,
   -0.0000000005206529618,
    0.0000000000874878885,
   -0.0000000000152176187,
    0.0000000000027257192,
   -0.0000000000005007053,
    0.0000000000000940241,
   -0.0000000000000180014,
    0.0000000000000035063,
   -0.0000000000000006935,
    0.0000000000000001391,
   -0.0000000000000000282
    ];
const f1_cs = { length: 19, c: f1_data, order: 19, a: -1.0, b: 1.0, order_sp: 10 };

// ----------------------------------------------------------------------------

// series for f2   on the interval  0.00000e+00 to  2.00000e-02
//                                        with weighted error   4.32e-17
//                                         log weighted error  16.36
//                               significant figures required  14.75
//                                    decimal places required  17.10
//
const f2_data =//: CONSTANT Series(0..28) := --[29]
    [
   -0.0348409253897013234,
   -0.0166842205677959686,
    0.0006752901241237738,
   -0.0000535066622544701,
    0.0000062693421779007,
   -0.0000009526638801991,
    0.0000001745629224251,
   -0.0000000368795403065,
    0.0000000087202677705,
   -0.0000000022601970392,
    0.0000000006324624977,
   -0.0000000001888911889,
    0.0000000000596774674,
   -0.0000000000198044313,
    0.0000000000068641396,
   -0.0000000000024731020,
    0.0000000000009226360,
   -0.0000000000003552364,
    0.0000000000001407606,
   -0.0000000000000572623,
    0.0000000000000238654,
   -0.0000000000000101714,
    0.0000000000000044259,
   -0.0000000000000019634,
    0.0000000000000008868,
   -0.0000000000000004074,
    0.0000000000000001901,
   -0.0000000000000000900,
    0.0000000000000000432
    ];
const f2_cs = { length: 28, c: f2_data, order: 28, a: -1.0, b: 1.0, order_sp: 14 };

// ----------------------------------------------------------------------------

// series for g1   on the interval  2.00000e-02 to  6.25000e-02
//                                        with weighted error   5.48e-17
//                                         log weighted error  16.26
//                               significant figures required  15.47
//                                    decimal places required  16.92
//
const g1_data =//: CONSTANT Series(0..20) := --[21]
    [
   -0.3040578798253495954,
   -0.0566890984597120588,
    0.0039046158173275644,
   -0.0003746075959202261,
    0.0000435431556559844,
   -0.0000057417294453025,
    0.0000008282552104503,
   -0.0000001278245892595,
    0.0000000207978352949,
   -0.0000000035313205922,
    0.0000000006210824236,
   -0.0000000001125215474,
    0.0000000000209088918,
   -0.0000000000039715832,
    0.0000000000007690431,
   -0.0000000000001514697,
    0.0000000000000302892,
   -0.0000000000000061400,
    0.0000000000000012601,
   -0.0000000000000002615,
    0.0000000000000000548
    ];
const g1_cs = { length: 20, c: g1_data, order: 20, a: -1.0, b: 1.0, order_sp: 13 };

// ----------------------------------------------------------------------------

// series for g2   on the interval  0.00000e+00 to  2.00000e-02
//                                        with weighted error   5.01e-17
//                                         log weighted error  16.30
//                               significant figures required  15.12
//                                    decimal places required  17.07
//
const g2_data =//: Series(0..33) := --[34]
    [
   -0.0967329367532432218,
   -0.0452077907957459871,
    0.0028190005352706523,
   -0.0002899167740759160,
    0.0000407444664601121,
   -0.0000071056382192354,
    0.0000014534723163019,
   -0.0000003364116512503,
    0.0000000859774367886,
   -0.0000000238437656302,
    0.0000000070831906340,
   -0.0000000022318068154,
    0.0000000007401087359,
   -0.0000000002567171162,
    0.0000000000926707021,
   -0.0000000000346693311,
    0.0000000000133950573,
   -0.0000000000053290754,
    0.0000000000021775312,
   -0.0000000000009118621,
    0.0000000000003905864,
   -0.0000000000001708459,
    0.0000000000000762015,
   -0.0000000000000346151,
    0.0000000000000159996,
   -0.0000000000000075213,
    0.0000000000000035970,
   -0.0000000000000017530,
    0.0000000000000008738,
   -0.0000000000000004487,
    0.0000000000000002397,
   -0.0000000000000001347,
    0.0000000000000000801,
   -0.0000000000000000501
    ];
const g2_cs = { length: 33, c: g2_data, order: 33, a: -1.0, b: 1.0, order_sp: 20 };

// ----------------------------------------------------------------------------

// x >= 4.0
function fg_asymp(x, f, g)
{
    const xbig  = 1.0 / GSL_SQRT_DBL_EPSILON;
    const xmaxf = 1.0 / GSL_DBL_MIN;
    const xmaxg = 1.0 / GSL_SQRT_DBL_MIN;
    const xbnd  = 7.07106781187;
    const x2    = x * x;
    var c1 = { val: 0.0, err: 0.0 }; // Result;
    var c2 = { val: 0.0, err: 0.0 }; // Result;

    if (x <= xbnd)
    {
        c1 = cheb_eval_e(f1_cs, (1.0 / x2 - 0.04125) / 0.02125);
        c2 = cheb_eval_e(g1_cs, (1.0 / x2 - 0.04125) / 0.02125);
        f.val = (1.0 + c1.val) / x;
        g.val = (1.0 + c2.val) / x2;
        f.err = c1.err / x  + 2.0 * Number.EPSILON * Math.abs(f.val);
        g.err = c2.err / x2 + 2.0 * Number.EPSILON * Math.abs(g.val);
    }
    else if (x <= xbig)
    {
        c1 = cheb_eval_e(f2_cs, 100.0 / x2 - 1.0);
        c2 = cheb_eval_e(g2_cs, 100.0 / x2 - 1.0);
        f.val = (1.0 + c1.val) / x;
        g.val = (1.0 + c2.val) / x2;
        f.err = c1.err / x  + 2.0 * Number.EPSILON * Math.abs(f.val);
        g.err = c2.err / x2 + 2.0 * Number.EPSILON * Math.abs(g.val);
    }
    else
    {
        if (x < xmaxf)
        {
            f.val = 1.0 / x;
        }
        else
        {
            f.val = 0.0;
        }
        if (x < xmaxg)
        {
            g.val = 1.0 / x2;
        }
        else
        {
            g.val = 0.0;
        }
        f.err = 2.0 * Number.EPSILON * Math.abs(f.val);
        g.err = 2.0 * Number.EPSILON * Math.abs(g.val);
    }

} // fg_asymp

// ----------------------------------------------------------------------------

// based on SLATEC si.f, W. Fullerton
//
// series for si   on the interval  0.00000e+00 to  1.60000e+01
//                                        with weighted error   1.22e-17
//                                         log weighted error  16.91
//                               significant figures required  16.37
//                                    decimal places required  17.45

const si_data =//: CONSTANT Series(0..11) := --[12]
    [
   -0.1315646598184841929,
   -0.2776578526973601892,
    0.0354414054866659180,
   -0.0025631631447933978,
    0.0001162365390497009,
   -0.0000035904327241606,
    0.0000000802342123706,
   -0.0000000013562997693,
    0.0000000000179440722,
   -0.0000000000001908387,
    0.0000000000000016670,
   -0.0000000000000000122
    ];
const si_cs = { length: 11, c: si_data, order: 11, a: -1.0, b: 1.0, order_sp: 9 };

// ----------------------------------------------------------------------------
//
// series for ci   on the interval  0.00000e+00 to  1.60000e+01
//                                        with weighted error   1.94e-18
//                                         log weighted error  17.71
//                               significant figures required  17.74
//                                    decimal places required  18.27

const ci_data =//: CONSTANT Series(0..12) := --[13]
    [
   -0.34004281856055363156,
   -1.03302166401177456807,
    0.19388222659917082877,
   -0.01918260436019865894,
    0.00110789252584784967,
   -0.00004157234558247209,
    0.00000109278524300229,
   -0.00000002123285954183,
    0.00000000031733482164,
   -0.00000000000376141548,
    0.00000000000003622653,
   -0.00000000000000028912,
    0.00000000000000000194
    ];
const ci_cs = { length: 12, c: ci_data, order: 12, a: -1.0, b: 1.0, order_sp: 9 };


// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_Si_e(x)
{
    const ax = Math.abs(x);
    var c = { val: 0.0, err: 0.0 }; // Result;
    var f = { val: 0.0, err: 0.0 }; // Result;
    var g = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (ax < GSL_SQRT_DBL_EPSILON)
    {
        r.val = x;
        r.err = 0.0;
    }
    else if (ax <= 4.0)
    {
        c = cheb_eval_e(si_cs, (x * x - 8.0) * 0.125);
        r.val =  x * (0.75 + c.val);
        r.err = ax * c.err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else
    {
        // Note there is no loss of precision
        // here bcause of the leading constant.
        fg_asymp(ax, f, g);
        r.val = 0.5 * M_PI - f.val * Math.cos(ax) - g.val * Math.sin(ax);
        r.err = f.err + g.err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
        if (x < 0.0)
        {
            r.val = -r.val;
        }
    }

    return r;

} // gsl_sf_Si_e

// ----------------------------------------------------------------------------

export function gsl_sf_Ci_e(x)
{
    var lx = 0.0;
    var y  = 0.0;
    var sin_result = { val: 0.0, err: 0.0 }; // Result;
    var cos_result = { val: 0.0, err: 0.0 }; // Result;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var f = { val: 0.0, err: 0.0 }; // Result;
    var g = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x <= 4.0)
    {
        lx = Math.log(x);
        y  = (x * x - 8.0) * 0.125;
        c = cheb_eval_e(ci_cs, y);
        r.val = lx - 0.5 + c.val;
        r.err = 2.0 * Number.EPSILON * (Math.abs(lx) + 0.5) + c.err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else
    {
        sin_result = gsl_sf_sin_e(x);
        cos_result = gsl_sf_cos_e(x);
        fg_asymp(x, f, g);
        r.val = f.val * sin_result.val - g.val * cos_result.val;
        r.err = Math.abs(f.err * sin_result.val);
        r.err = r.err + Math.abs(g.err * cos_result.val);
        r.err = r.err + Math.abs(f.val * sin_result.err);
        r.err = r.err + Math.abs(g.val * cos_result.err);
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_Ci_e

// ----------------------------------------------------------------------------

const atanint_data =//: CONSTANT Series(0..20) := --[21]
    [
    1.91040361296235937512,
   -0.4176351437656746940e-01,
    0.275392550786367434e-02,
   -0.25051809526248881e-03,
    0.2666981285121171e-04,
   -0.311890514107001e-05,
    0.38833853132249e-06,
   -0.5057274584964e-07,
    0.681225282949e-08,
   -0.94212561654e-09,
    0.13307878816e-09,
   -0.1912678075e-10,
    0.278912620e-11,
   -0.41174820e-12,
    0.6142987e-13,
   -0.924929e-14,
    0.140387e-14,
   -0.21460e-15,
    0.3301e-16,
   -0.511e-17,
    0.79e-18
    ];
const atanint_cs = { length: 20, c: atanint_data, order: 20, a: -1.0, b: 1.0, order_sp: 10 };

// ----------------------------------------------------------------------------

export function gsl_sf_atanint_e(x)
{
    const ax  = Math.abs(x);
    const sgn = GSL_SIGN(x);
    var t   = 0.0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (ax == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (ax < 0.5 * GSL_SQRT_DBL_EPSILON)
    {
        r.val = x;
        r.err = 0.0;
    }
    else if (ax <= 1.0)
    {
        t = 2.0 * (x * x - 0.5);
        c = cheb_eval_e(atanint_cs, t);
        r.val = x * c.val;
        r.err = x * c.err;
        r.err = r.err + Number.EPSILON * Math.abs(r.val);
    }
    else if (ax < 1.0 / GSL_SQRT_DBL_EPSILON)
    {
        t = 2.0 * (1.0 / (x * x) - 0.5);
        c = cheb_eval_e(atanint_cs, t);
        r.val = sgn * (0.5 * M_PI * Math.log(ax) + c.val / ax);
        r.err = c.err / ax + Math.abs(r.val * Number.EPSILON);
        r.err = r.err + Number.EPSILON * Math.abs(r.val);
    }
    else
    {
        r.val = sgn * 0.5 * M_PI * Math.log(ax);
        r.err = 2.0 * Math.abs(r.val * Number.EPSILON);
    }

    return r;

} // gsl_sf_atanint_e

// // -*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

// FUNCTION gsl_sf_Si(x: LONG_FLOAT) RETURN LONG_FLOAT IS
// BEGIN -- gsl_sf_Si
//     RETURN EVAL_RESULT(gsl_sf_Si_e'Access, x, "gsl_sf_Si_e");
// END gsl_sf_Si;

// FUNCTION gsl_sf_Ci(x: LONG_FLOAT) RETURN LONG_FLOAT IS
// BEGIN -- gsl_sf_Ci
//     RETURN EVAL_RESULT(gsl_sf_Ci_e'Access, x, "gsl_sf_Ci_e");
// END gsl_sf_Ci;

// FUNCTION gsl_sf_atanint(x: LONG_FLOAT) RETURN LONG_FLOAT IS
// BEGIN -- gsl_sf_atanint
//     RETURN EVAL_RESULT(gsl_sf_atanint_e'Access, x, "gsl_sf_atanint_e");
// END gsl_sf_atanint;

// END SF.TrigonometricIntegrals;

// // ----------------------------------------------------------------------------
// // SF-TrigonometricIntegrals.mjs
