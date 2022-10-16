// SF-BesselI0.mjs
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

import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }      from "./SF-Machine.mjs";
import { cheb_eval_e }          from "./SF-Chebyshev.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*


// based on SLATEC besi0

// chebyshev expansions
//
// series for bi0        on the interval  0.          to  9.00000d+00
//                                        with weighted error   2.46e-18
//                                         log weighted error  17.61
//                               significant figures required  17.90
//                                    decimal places required  18.15
//
// series for ai0        on the interval  1.25000d-01 to  3.33333d-01
//                                        with weighted error   7.87e-17
//                                         log weighted error  16.10
//                               significant figures required  14.69
//                                    decimal places required  16.76
//
//
// series for ai02       on the interval  0.          to  1.25000d-01
//                                        with weighted error   3.79e-17
//                                         log weighted error  16.42
//                               significant figures required  14.86
//                                    decimal places required  17.09
//

const bi0_data =
    [
   -0.07660547252839144951,
   1.92733795399380827000,
    0.22826445869203013390, 
    0.01304891466707290428,
    0.00043442709008164874,
    0.00000942265768600193,
    0.00000014340062895106,
    0.00000000161384906966,
    0.00000000001396650044,
    0.00000000000009579451,
    0.00000000000000053339,
    0.00000000000000000245
    ];
const bi0_cs = { length: 11, c: bi0_data, order: 11, a: -1.0, b: 1.0, order_sp: 11 };

const ai0_data =
    [
    0.07575994494023796, 
    0.00759138081082334,
    0.00041531313389237,
    0.00001070076463439,
   -0.00000790117997921,
   -0.00000078261435014,
    0.00000027838499429,
    0.00000000825247260,
   -0.00000001204463945,
    0.00000000155964859,
    0.00000000022925563,
   -0.00000000011916228,
    0.00000000001757854,
    0.00000000000112822,
   -0.00000000000114684,
    0.00000000000027155,
   -0.00000000000002415,
   -0.00000000000000608,
    0.00000000000000314,
   -0.00000000000000071,
    0.00000000000000007
    ];
const ai0_cs = { length: 20, c: ai0_data, order: 20, a: -1.0, b: 1.0, order_sp: 13 };

const ai02_data =
    [
    0.05449041101410882,
    0.00336911647825569,
    0.00006889758346918,
    0.00000289137052082,
    0.00000020489185893,
    0.00000002266668991,
    0.00000000339623203,
    0.00000000049406022,
    0.00000000001188914,
   -0.00000000003149915,
   -0.00000000001321580,
   -0.00000000000179419,
    0.00000000000071801,
    0.00000000000038529,
    0.00000000000001539,
   -0.00000000000004151,
   -0.00000000000000954,
    0.00000000000000382,
    0.00000000000000176,
   -0.00000000000000034,
   -0.00000000000000027,
    0.00000000000000003
    ];
const ai02_cs = { length: 21, c: ai02_data, order: 21, a: -1.0, b: 1.0, order_sp: 11 };

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_I0_scaled_e(x)
{
    var y  = Math.abs(x);
    var ey = 0.0;
    var sy = 0.0;
    var c  = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (y < 2.0 * GSL_SQRT_DBL_EPSILON)
    {
        r.val = 1.0 - y;
        r.err = 0.5 * y * y;
    }
    else if (y <= 3.0)
    {
        ey = Math.exp(-y);
        c = cheb_eval_e(bi0_cs, y * y / 4.5 - 1.0);
        r.val = ey * (2.75 + c.val);
        r.err = GSL_DBL_EPSILON * Math.abs(r.val) + ey * c.err;
    }
    else if (y <= 8.0)
    {
        sy = Math.sqrt(y);
        c = cheb_eval_e(ai0_cs, (48.0 / y - 11.0) / 5.0);
        r.val = (0.375 + c.val) / sy;
        r.err = 2.0 * GSL_DBL_EPSILON * (0.375 + Math.asb(c.val)) / sy;
        r.err = r.err + c.err / sy;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        sy = Math.sqrt(y);
        c = cheb_eval_e(ai02_cs, 16.0 / y - 1.0);
        r.val = (0.375 + c.val) / sy;
        r.err = 2.0 * GSL_DBL_EPSILON * (0.375 + Math.abs(c.val)) / sy;
        r.err = r.err + c.err / sy;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_I0_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_I0_e(x)
{
    var y  = Math.abs(x);
    var ey = 0.0;
    var c  = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (y < 2.0 * GSL_SQRT_DBL_EPSILON)
    {
        r.val = 1.0;
        r.err = 0.5 * y * y;
    }
    else if (y <= 3.0)
    {
        c = cheb_eval_e(bi0_cs, y * y / 4.5 - 1.0);
        r.val = 2.75 + c.val;
        r.err = GSL_DBL_EPSILON * (2.75 + Math.abs(c.val));
        r.err = r.err + c.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (y < GSL_LOG_DBL_MAX - 1.0)
    {
        ey = Math.exp(y);
        c = gsl_sf_bessel_I0_scaled_e(x);
        r.val = ey * c.val;
        r.err = ey * c.err + y * GSL_DBL_EPSILON * Math.abs(r.val);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

    return r;

} // gsl_sf_bessel_I0_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_I0_scaled( x )
{ // gsl_sf_bessel_I0_scaled
    return EVAL_RESULT_D( gsl_sf_bessel_I0_scaled_e, x, "gsl_sf_bessel_I0_scaled" );
} // gsl_sf_bessel_I0_scaled;

export function gsl_sf_bessel_I0( x )
{ // gsl_sf_bessel_I0
    return EVAL_RESULT_D( gsl_sf_bessel_I0_e, x, "gsl_sf_bessel_I0" );
} // gsl_sf_bessel_I0;

// ----------------------------------------------------------------------------
// EOF SF-BesselI0.mjs
