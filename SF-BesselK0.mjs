// SF-BesselK0.mjs
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

import { M_LN2 }                 from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { cheb_eval_e }           from "./SF-Chebyshev.mjs";
import { gsl_sf_bessel_I0_e }    from "./SF-BesselI0.mjs";
import { gsl_sf_exp_mult_err_e } from "./SF-Exponential.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// based on SLATEC bk0(), bk0e()

// chebyshev expansions 
//
// series for bk0        on the interval  0.          to  4.00000d+00
//                                        with weighted error   3.57e-19
//                                         log weighted error  18.45
//                               significant figures required  17.99
//                                    decimal places required  18.97
//
// series for ak0        on the interval  1.25000d-01 to  5.00000d-01
//                                        with weighted error   5.34e-17
//                                         log weighted error  16.27
//                               significant figures required  14.92
//                                    decimal places required  16.89
//
// series for ak02       on the interval  0.          to  1.25000d-01
//                                        with weighted error   2.34e-17
//                                         log weighted error  16.63
//                               significant figures required  14.67
//                                    decimal places required  17.20
//

const bk0_data =
    [
   -0.03532739323390276872,
    0.3442898999246284869, 
    0.03597993651536150163,
    0.00126461541144692592,
    0.00002286212103119451,
    0.00000025347910790261,
    0.00000000190451637722,
    0.00000000001034969525,
    0.00000000000004259816,
    0.00000000000000013744,
    0.00000000000000000035
    ];
const bk0_cs = { length: 10, c: bk0_data, order: 10, a: -1.0, b: 1.0, order_sp: 10 };

const ak0_data =
    [
   -0.07643947903327941,
   -0.02235652605699819,
    0.00077341811546938,
   -0.00004281006688886,
    0.00000308170017386,
   -0.00000026393672220,
    0.00000002563713036,
   -0.00000000274270554,
    0.00000000031694296,
   -0.00000000003902353,
    0.00000000000506804,
   -0.00000000000068895,
    0.00000000000009744,
   -0.00000000000001427,
    0.00000000000000215,
   -0.00000000000000033,
    0.00000000000000005
    ];
const ak0_cs = { length: 16, c: ak0_data, order: 16, a: -1.0, b: 1.0, order_sp: 10 };

const ak02_data =
    [
   -0.01201869826307592,
   -0.00917485269102569,
    0.00014445509317750,
   -0.00000401361417543,
    0.00000015678318108,
   -0.00000000777011043,
    0.00000000046111825,
   -0.00000000003158592,
    0.00000000000243501,
   -0.00000000000020743,
    0.00000000000001925,
   -0.00000000000000192,
    0.00000000000000020,
   -0.00000000000000002
    ];
const ak02_cs = { length: 13, c: ak02_data, order: 13, a: -1.0, b: 1.0, order_sp: 8 };

// -*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_K0_scaled_e(x)
{
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x <= 2.0)
    {
        let lx = Math.log(x);
        let ex = Math.exp(x);
        let I0 = { val: 0.0, err: 0.0 }; // Result;

        c = cheb_eval_e(bk0_cs, 0.5 * x * x - 1.0);
        I0 = gsl_sf_bessel_I0_e(x);
        r.val = ex * ((-lx + M_LN2) * I0.val - 0.25 + c.val);
        r.err = ex * ((M_LN2 + Math.abs(lx)) * I0.err + c.err);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x <= 8.0)
    {
        let sx = Math.sqrt(x);

        c = cheb_eval_e(ak0_cs, (16.0 / x - 5.0) / 3.0);
        r.val = (1.25 + c.val) / sx;
        r.err = c.err / sx;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        let sx = Math.sqrt(x);

        c = cheb_eval_e(ak02_cs, 16.0 / x - 1.0);
        r.val = (1.25 + c.val) / sx;
        r.err = (c.err + GSL_DBL_EPSILON) / sx;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_K0_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_K0_e(x)
{
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x <= 2.0)
    {
        let lx = Math.log(x);
        let I0 = { val: 0.0, err: 0.0 }; // Result;

        c = cheb_eval_e(bk0_cs, 0.5 * x * x - 1.0);
        I0 = gsl_sf_bessel_I0_e(x);
        r.val = (-lx + M_LN2) * I0.val - 0.25 + c.val;
        r.err = (Math.abs(lx) + M_LN2) * I0.err + c.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        let K0_scaled = { val: 0.0, err: 0.0 }; // Result;

        K0_scaled = gsl_sf_bessel_K0_scaled_e(x);
        r = gsl_sf_exp_mult_err_e(-x, GSL_DBL_EPSILON * Math.abs(x), K0_scaled.val, K0_scaled.err);
    }

    return r;

} // gsl_sf_bessel_K0_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_K0_scaled( x )
{ // gsl_sf_bessel_K0_scaled
    return EVAL_RESULT_D( gsl_sf_bessel_K0_scaled_e, x, "gsl_sf_bessel_K0_scaled" );
} // gsl_sf_bessel_K0_scaled

export function gsl_sf_bessel_K0( x )
{ // gsl_sf_bessel_K0
    return EVAL_RESULT_D( gsl_sf_bessel_K0_e, x, "gsl_sf_bessel_K0" );
} // gsl_sf_bessel_K0

// ----------------------------------------------------------------------------
// EOF SF-BesselK0.mjs
