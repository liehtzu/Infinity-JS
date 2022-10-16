// SF-BesselY0.mjs
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

import { GSL_DBL_EPSILON }    from "./SF-Machine.mjs";
import { M_PI }               from "./SF-Math.mjs";
import { M_LN2 }              from "./SF-Math.mjs";
import { cheb_eval_e }        from "./SF-Chebyshev.mjs";
import { gsl_sf_bessel_J0_e } from "./SF-BesselJ0.mjs";
import { gsl_sf_bessel_sin_pi4_e } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_amp_phase_bm0_cs } from "./SF-BesselAF.mjs";
import { gsl_sf_bessel_amp_phase_bth0_cs } from "./SF-BesselAF.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// based on SLATEC besy0, 1980 version, w. fullerton

// chebyshev expansions
//
// series for by0        on the interval  0.          to  1.60000d+01
//                                        with weighted error   1.20e-17
//                                         log weighted error  16.92
//                               significant figures required  16.15
//                                    decimal places required  17.48
//

const by0_data =
    [
   -0.011277839392865573,
   -0.128345237560420350,
   -0.104378847997942490,
    0.023662749183969695,
   -0.002090391647700486,
    0.000103975453939057,
   -0.000003369747162423,
    0.000000077293842676,
   -0.000000001324976772,
    0.000000000017648232,
   -0.000000000000188105,
    0.000000000000001641,
   -0.000000000000000011
    ];
const by0_cs = { length: 12, c: by0_data, order: 12, a: -1.0, b: 1.0, order_sp: 8 };


// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Y0_e(x)
{
    const two_over_pi = 2.0 / M_PI;
    const xmax        = 1.0 / GSL_DBL_EPSILON;

    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 4.0)
    {
        var J0 = { val: 0.0, err: 0.0 }; // Result;
        var c  = { val: 0.0, err: 0.0 }; // Result;

        J0 = gsl_sf_bessel_J0_e(x);
        c = cheb_eval_e(by0_cs, 0.125 * x * x - 1.0);
        r.val = two_over_pi * (-M_LN2 + Math.log(x)) * J0.val + 0.375 + c.val;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val) + c.err;
    }
    else if (x < xmax)
    {
        // Leading behaviour of phase is x, which is exact,
        // so the error is bounded.
        //
        var z = 32.0 / (x * x) - 1.0;
        var sqrtx = 0.0;
        var ampl  = 0.0;
        var c1 = { val: 0.0, err: 0.0 }; // Result;
        var c2 = { val: 0.0, err: 0.0 }; // Result;
        var sp = { val: 0.0, err: 0.0 }; // Result;

        c1 = cheb_eval_e(gsl_sf_bessel_amp_phase_bm0_cs,  z);
        c2 = cheb_eval_e(gsl_sf_bessel_amp_phase_bth0_cs, z);
        sp = gsl_sf_bessel_sin_pi4_e(x, c2.val / x);
        sqrtx = Math.sqrt(x);
        ampl  = (0.75 + c1.val) / sqrtx;
        r.val = ampl * sp.val;
        r.err = Math.abs(sp.val) * c1.err / sqrtx + Math.abs(ampl) * sp.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.UnderflowException";
    }

    return r;

} // gsl_sf_bessel_Y0_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Y0( x )
{ // gsl_sf_bessel_Y0
    return EVAL_RESULT_D( gsl_sf_bessel_Y0_e, x, "gsl_sf_bessel_Y0" );
} // gsl_sf_bessel_Y0

// ----------------------------------------------------------------------------
// EOF SF-BesselY0.mjs
