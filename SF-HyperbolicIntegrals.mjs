// SF-HyperbolicIntegrals.mjs
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

import { gsl_sf_expint_E1_e } from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_Ei_e } from "./SF-ExponentialIntegral.mjs";
import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { cheb_eval_e } from "./SF-Chebyshev.mjs";

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// based on SLATEC shi.f, W. Fullerton
//
// series for shi  on the interval  0.00000e+00 to  1.40625e-01
//                                        with weighted error   4.67e-20
//                                         log weighted error  19.33
//                               significant figures required  17.07
//                                    decimal places required  19.75

const shi_data =//: Series(0..6) := --[7]
    [
    0.0078372685688900950695,
    0.0039227664934234563973,
    0.0000041346787887617267,
    0.0000000024707480372883,
    0.0000000000009379295591,
    0.0000000000000002451817,
    0.0000000000000000000467
    ];
const shi_cs = { length: 6, c: shi_data, order: 6, a: -1.0, b: 1.0, order_sp: 6 };

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_Shi_e(x)
{
    const xsml = GSL_SQRT_DBL_EPSILON; // sqrt (r1mach(3))
    const ax   = Math.abs(x);
    var c    = { val: 0.0, err: 0.0 }; // Result;
    var r    = { val: 0.0, err: 0.0 }; // Result;
    var Ei   = { val: 0.0, err: 0.0 }; // Result;
    var E1   = { val: 0.0, err: 0.0 }; // Result;

    if (ax < xsml)
    {
        r.val = x;
        r.err = 0.0;
    }
    else if (ax <= 0.375)
    {
        c = cheb_eval_e(shi_cs, 128.0 * x * x / 9.0 - 1.0);
        r.val = x * (1.0 + c.val);
        r.err = x * c.err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else
    {
        Ei = gsl_sf_expint_Ei_e(x);
        E1 = gsl_sf_expint_E1_e(x);
        r.val = 0.5 * (Ei.val + E1.val);
        r.err = 0.5 * (Ei.err + E1.err);
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_Shi_e

// ----------------------------------------------------------------------------

export function gsl_sf_Chi_e(x)
{
    var r  = { val: 0.0, err: 0.0 }; // Result;
    var Ei = { val: 0.0, err: 0.0 }; // Result;
    var E1 = { val: 0.0, err: 0.0 }; // Result;

    Ei = gsl_sf_expint_Ei_e(x);
    E1 = gsl_sf_expint_E1_e(x);

    //IF ((status_Ei = GSL_EDOM) OR (status_E1 = GSL_EDOM)) THEN
    //    DOMAIN_ERROR(result);
    //ELSIF ((status_Ei = GSL_EUNDRFLW) AND (status_E1 = GSL_EUNDRFLW)) THEN
    //    UNDERFLOW_ERROR(result);
    //ELSIF ((status_Ei = GSL_EOVRFLW) OR (status_E1 = GSL_EOVRFLW)) THEN
    //    OVERFLOW_ERROR(result);
    //ELSE
    r.val = 0.5 * (Ei.val - E1.val);
    r.err = 0.5 * (Ei.err + E1.err);
    r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    //END IF;

    return r;

} // gsl_sf_Chi_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

// FUNCTION gsl_sf_Shi(x: LONG_FLOAT) RETURN LONG_FLOAT IS
// BEGIN -- gsl_sf_Shi
//     RETURN EVAL_RESULT(gsl_sf_Shi_e'Access, x, "gsl_sf_Shi");
// END gsl_sf_Shi;

// FUNCTION gsl_sf_Chi(x: LONG_FLOAT) RETURN LONG_FLOAT IS
// BEGIN -- gsl_sf_Chi
//     RETURN EVAL_RESULT(gsl_sf_Chi_e'Access, x, "gsl_sf_Chi");
// END gsl_sf_Chi;

// END SF.HyperbolicIntegrals;

// ----------------------------------------------------------------------------
// SF-HyperbolicIntegrals.mjs
