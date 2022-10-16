// SF-BesselK1.mjs
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
import { GSL_DBL_MIN }           from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { cheb_eval_e }           from "./SF-Chebyshev.mjs";
import { gsl_sf_bessel_I1_e }    from "./SF-BesselI1.mjs";
import { gsl_sf_exp_mult_err_e } from "./SF-Exponential.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// based on SLATEC besk1(), besk1e()

// chebyshev expansions 
//
// series for bk1        on the interval  0.          to  4.00000d+00
//                                        with weighted error   7.02e-18
//                                         log weighted error  17.15
//                               significant figures required  16.73
//                                    decimal places required  17.67
//
// series for ak1        on the interval  1.25000d-01 to  5.00000d-01
//                                        with weighted error   6.06e-17
//                                         log weighted error  16.22
//                               significant figures required  15.41
//                                    decimal places required  16.83
//
// series for ak12       on the interval  0.          to  1.25000d-01
//                                        with weighted error   2.58e-17
//                                         log weighted error  16.59
//                               significant figures required  15.22
//                                    decimal places required  17.16
//

const bk1_data =
    [
    0.0253002273389477705,
   -0.3531559607765448760, 
   -0.1226111808226571480, 
   -0.0069757238596398643,
   -0.0001730288957513052,
   -0.0000024334061415659,
   -0.0000000221338763073,
   -0.0000000001411488392,
   -0.0000000000006666901,
   -0.0000000000000024274,
   -0.0000000000000000070
    ];
const bk1_cs = { length: 10, c: bk1_data, order: 10, a: -1.0, b: 1.0, order_sp: 8 };

const ak1_data =
    [
    0.27443134069738830, 
    0.07571989953199368,
   -0.00144105155647540,
    0.00006650116955125,
   -0.00000436998470952,
    0.00000035402774997,
   -0.00000003311163779,
    0.00000000344597758,
   -0.00000000038989323,
    0.00000000004720819,
   -0.00000000000604783,
    0.00000000000081284,
   -0.00000000000011386,
    0.00000000000001654,
   -0.00000000000000248,
    0.00000000000000038,
   -0.00000000000000006
    ];
const ak1_cs = { length: 16, c: ak1_data, order: 16, a: -1.0, b: 1.0, order_sp: 9 };

const ak12_data =
    [
    0.06379308343739001,
    0.02832887813049721,
   -0.00024753706739052,
    0.00000577197245160,
   -0.00000020689392195,
    0.00000000973998344,
   -0.00000000055853361,
    0.00000000003732996,
   -0.00000000000282505,
    0.00000000000023720,
   -0.00000000000002176,
    0.00000000000000215,
   -0.00000000000000022,
    0.00000000000000002
    ];
const ak12_cs = { length: 13, c: ak12_data, order: 13, a: -1.0, b: 1.0, order_sp: 7 };

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_K1_scaled_e(x)
{
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 2.0 * GSL_DBL_MIN)
    {
        throw "SF.OverflowException";
    }
    else if (x <= 2.0)
    {
        let lx = Math.log(x);
        let ex = Math.exp(x);
        let I1 = { val: 0.0, err: 0.0 }; // Result;

        c = cheb_eval_e(bk1_cs, 0.5 * x * x - 1.0);
        I1 = gsl_sf_bessel_I1_e(x);
        r.val = ex * ((lx - M_LN2) * I1.val + (0.75 + c.val) / x);
        r.err = ex * (c.err / x + Math.abs(lx) * I1.err);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x <= 8.0)
    {
        let sx = Math.sqrt(x);

        c = cheb_eval_e(ak1_cs, (16.0 / x - 5.0) / 3.0);
        r.val = (1.25 + c.val) / sx;
        r.err = c.err / sx;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        let sx = Math.sqrt(x);

        c = cheb_eval_e(ak12_cs, 16.0 / x - 1.0);
        r.val = (1.25 + c.val) / sx;
        r.err = c.err / sx;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_K1_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_K1_e(x)
{
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 2.0*GSL_DBL_MIN)
    {
        throw "SF.OverflowException";
    }
    else if (x <= 2.0)
    {
        let lx = Math.log(x);
        let I1 = { val: 0.0, err: 0.0 }; // Result;

        c = cheb_eval_e(bk1_cs, 0.5 * x * x - 1.0);
        I1 = gsl_sf_bessel_I1_e(x);
        r.val = (lx - M_LN2) * I1.val + (0.75 + c.val) / x;
        r.err = c.err / x + Math.abs(lx) * I1.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        let K1_scaled = { val: 0.0, err: 0.0 }; // Result;

        K1_scaled = gsl_sf_bessel_K1_scaled_e(x);
        r = gsl_sf_exp_mult_err_e(-x, 0.0, K1_scaled.val, K1_scaled.err);
        r.err = Math.abs(r.val) * (GSL_DBL_EPSILON * Math.abs(x) + K1_scaled.err / K1_scaled.val);
    }

    return r;

} // gsl_sf_bessel_K1_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_K1_scaled( x )
{ // gsl_sf_bessel_K1_scaled
    return EVAL_RESULT_D( gsl_sf_bessel_K1_scaled_e, x, "gsl_sf_bessel_K1_scaled" );
} // gsl_sf_bessel_K1_scaled

export function gsl_sf_bessel_K1( x )
{ // gsl_sf_bessel_K1
    return EVAL_RESULT_D( gsl_sf_bessel_K1_e, x, "gsl_sf_bessel_K1" );
} // gsl_sf_bessel_K1

// ----------------------------------------------------------------------------
// EOF SF-BesselK1.mjs
