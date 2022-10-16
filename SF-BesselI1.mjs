// SF-BesselI1.mjs
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

import { M_SQRT2 }              from "./SF-Math.mjs";
import { GSL_DBL_MIN }          from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }      from "./SF-Machine.mjs";
import { cheb_eval_e }          from "./SF-Chebyshev.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const ROOT_EIGHT = 2.0 * M_SQRT2;

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// based on SLATEC besi1(), besi1e()
//
// chebyshev expansions
//
// series for bi1        on the interval  0.          to  9.00000d+00
//                                        with weighted error   2.40e-17
//                                         log weighted error  16.62
//                               significant figures required  16.23
//                                    decimal places required  17.14
//
// series for ai1        on the interval  1.25000d-01 to  3.33333d-01
//                                        with weighted error   6.98e-17
//                                         log weighted error  16.16
//                               significant figures required  14.53
//                                    decimal places required  16.82
//
// series for ai12       on the interval  0.          to  1.25000d-01
//                                       with weighted error   3.55e-17
//                                        log weighted error  16.45
//                              significant figures required  14.69
//                                   decimal places required  17.12
//

const bi1_data =//: CONSTANT Series(0..10) = --[11]
    [
   -0.001971713261099859,
    0.407348876675464810,
    0.034838994299959456,
    0.001545394556300123,
    0.000041888521098377,
    0.000000764902676483,
    0.000000010042493924,
    0.000000000099322077,
    0.000000000000766380,
    0.000000000000004741,
    0.000000000000000024
    ];
const bi1_cs = { length: 10, c: bi1_data, order: 10, a: -1.0, b: 1.0, order_sp: 10 };

const ai1_data =//: CONSTANT Series(0..20) = --[21]
    [
   -0.02846744181881479,
   -0.01922953231443221,
   -0.00061151858579437,
   -0.00002069971253350,
    0.00000858561914581,
    0.00000104949824671,
   -0.00000029183389184,
   -0.00000001559378146,
    0.00000001318012367,
   -0.00000000144842341,
   -0.00000000029085122,
    0.00000000012663889,
   -0.00000000001664947,
   -0.00000000000166665,
    0.00000000000124260,
   -0.00000000000027315,
    0.00000000000002023,
    0.00000000000000730,
   -0.00000000000000333,
    0.00000000000000071,
   -0.00000000000000006
    ];
const ai1_cs = { length: 20, c: ai1_data, order: 20, a: -1.0, b: 1.0, order_sp: 11 };

const ai12_data =//: CONSTANT Series(0..21) = --[22]
    [
    0.02857623501828014,
   -0.00976109749136147,
   -0.00011058893876263,
   -0.00000388256480887,
   -0.00000025122362377,
   -0.00000002631468847,
   -0.00000000383538039,
   -0.00000000055897433,
   -0.00000000001897495,
    0.00000000003252602,
    0.00000000001412580,
    0.00000000000203564,
   -0.00000000000071985,
   -0.00000000000040836,
   -0.00000000000002101,
    0.00000000000004273,
    0.00000000000001041,
   -0.00000000000000382,
   -0.00000000000000186,
    0.00000000000000033,
    0.00000000000000028,
   -0.00000000000000003
    ];
const ai12_cs = { length: 21, c: ai12_data, order: 21, a: -1.0, b: 1.0, order_sp: 9 };

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_I1_scaled_e(x)
{
    var xmin    = 2.0 * GSL_DBL_MIN;
    var x_small = ROOT_EIGHT * GSL_SQRT_DBL_EPSILON;
    var y  = Math.abs(x);
    var ey = 0.0;
    var sy = 0.0;
    var b  = 0.0;
    var s  = 0.0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (y == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (y < xmin)
    {
        throw "SF.UnderflowException";
    }
    else if (y < x_small)
    {
        r.val = 0.5 * x;
        r.err = 0.0;
    }
    else if (y <= 3.0)
    {
        ey = Math.exp(-y);
        c = cheb_eval_e(bi1_cs, y * y / 4.5 - 1.0);
        r.val = x * ey * (0.875 + c.val);
        r.err = ey * c.err + y * GSL_DBL_EPSILON * Math.abs(r.val);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (y <= 8.0)
    {
        sy = Math.sqrt(y);
        c = cheb_eval_e(ai1_cs, (48.0 / y - 11.0) / 5.0);
        b = (0.375 + c.val) / sy;
        if (x > 0.0)
        {
            s =  1.0;
        }
        else
        {
            s = -1.0;
        }
        r.val = s * b;
        r.err = c.err / sy;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        sy = Math.sqrt(y);
        c = cheb_eval_e(ai12_cs, 16.0 / y - 1.0);
        b = (0.375 + c.val) / sy;
        if (x > 0.0)
        {
            s =  1.0;
        }
        else
        {
            s = -1.0;
        }
        r.val = s * b;
        r.err = c.err / sy;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_I1_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_I1_e(x)
{
    var xmin    = 2.0 * GSL_DBL_MIN;
    var x_small = ROOT_EIGHT * GSL_SQRT_DBL_EPSILON;
    var y  = Math.abs(x);
    var ey = 0.0;
    var c  = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (y == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (y < xmin)
    {
        throw "SF.UnderflowException";
    }
    else if (y < x_small)
    {
        r.val = 0.5 * x;
        r.err = 0.0;
    }
    else if (y <= 3.0)
    {
        c = cheb_eval_e(bi1_cs, y * y / 4.5 - 1.0);
        r.val = x * (0.875 + c.val);
        r.err = y * c.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (y < GSL_LOG_DBL_MAX)
    {
        ey = Math.exp(y);
        c = gsl_sf_bessel_I1_scaled_e(x);
        r.val = ey * c.val;
        r.err = ey * c.err + y * GSL_DBL_EPSILON * Math.abs(r.val);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

    return r;

} // gsl_sf_bessel_I1_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_I1_scaled( x )
{ // gsl_sf_bessel_I1_scaled
    return EVAL_RESULT_D( gsl_sf_bessel_I1_scaled_e, x, "gsl_sf_bessel_I1_scaled" );
} // gsl_sf_bessel_I1_scaled

export function gsl_sf_bessel_I1( x )
{ // gsl_sf_bessel_I1
    return EVAL_RESULT_D( gsl_sf_bessel_I1_e, x, "gsl_sf_bessel_I1" );
} // gsl_sf_bessel_I1

// ----------------------------------------------------------------------------
// EOF SF-BesselI1.mjs
