// SF-BesselY1.mjs
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

import { GSL_DBL_MIN }          from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { M_PI }                 from "./SF-Math.mjs";
import { cheb_eval_e }          from "./SF-Chebyshev.mjs";
import { gsl_sf_bessel_J1_e }   from "./SF-BesselJ1.mjs";
import { gsl_sf_bessel_cos_pi4_e } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_amp_phase_bm1_cs } from "./SF-BesselAF.mjs";
import { gsl_sf_bessel_amp_phase_bth1_cs } from "./SF-BesselAF.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

// based on SLATEC besy1, 1977 version, w. fullerton

// chebyshev expansions
//
// series for by1        on the interval  0.          to  1.60000d+01
//                                        with weighted error   1.87e-18
//                                         log weighted error  17.73
//                               significant figures required  17.83
//                                    decimal places required  18.30
//

const by1_data =
    [
    0.03208047100611908629,
    1.262707897433500450,
    0.00649996189992317500,
   -0.08936164528860504117,
    0.01325088122175709545,
   -0.00089790591196483523,
    0.00003647361487958306,
   -0.00000100137438166600,
    0.00000001994539657390,
   -0.00000000030230656018,
    0.00000000000360987815,
   -0.00000000000003487488,
    0.00000000000000027838,
   -0.00000000000000000186
    ];
const by1_cs = { length: 13, c: by1_data, order: 13, a: -1.0, b: 1.0, order_sp: 10 };

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Y1_e(x)
{
    const two_over_pi = 2.0 / M_PI;
    const xmin        = 1.571 * GSL_DBL_MIN; //exp ( amax1(alog(r1mach(1)), -alog(r1mach(2)))+.01)
    const x_small     = 2.0 * GSL_SQRT_DBL_EPSILON;
    const xmax        = 1.0 / GSL_DBL_EPSILON;

    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < xmin)
    {
        throw "SF.OverflowException";
    }
    else if (x < x_small)
    {
        const lnterm = Math.log(0.5 * x);
        let J1 = { val: 0.0, err: 0.0 }; // Result;
        let c  = { val: 0.0, err: 0.0 }; // Result;

        J1 = gsl_sf_bessel_J1_e(x);
        c = cheb_eval_e(by1_cs, -1.0);
        r.val = two_over_pi * lnterm * J1.val + (0.5 + c.val) / x;
        r.err = Math.abs(lnterm) * (Math.abs(GSL_DBL_EPSILON * J1.val) + J1.err) + c.err / x;
    }
    else if (x < 4.0)
    {
        const lnterm = Math.log(0.5 * x);
        let J1 = { val: 0.0, err: 0.0 }; // Result;
        let c  = { val: 0.0, err: 0.0 }; // Result;

        c = cheb_eval_e(by1_cs, 0.125 * x * x - 1.0);
        J1 = gsl_sf_bessel_J1_e(x);
        r.val = two_over_pi * lnterm * J1.val + (0.5 + c.val) / x;
        r.err = Math.abs(lnterm) * (Math.abs(GSL_DBL_EPSILON * J1.val) + J1.err) + c.err / x;
    }
    else if (x < xmax)
    {
        const z = 32.0 / (x * x) - 1.0;
        let sqrtx = 0.0;
        let ampl  = 0.0;
        let ca = { val: 0.0, err: 0.0 }; // Result;
        let ct = { val: 0.0, err: 0.0 }; // Result;
        let cp = { val: 0.0, err: 0.0 }; // Result;

        ca = cheb_eval_e(gsl_sf_bessel_amp_phase_bm1_cs,  z);
        ct = cheb_eval_e(gsl_sf_bessel_amp_phase_bth1_cs, z);
        cp = gsl_sf_bessel_cos_pi4_e(x, ct.val / x);
        sqrtx = Math.sqrt(x);
        ampl  = (0.75 + ca.val) / sqrtx;
        r.val = -ampl * cp.val;
        r.err = Math.abs(cp.val) * ca.err / sqrtx + Math.abs(ampl) * cp.err;
        r.err = r.err + GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.UnderflowException";
    }

    return r;

} // gsl_sf_bessel_Y1_e


// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Y1( x )
{ // gsl_sf_bessel_Y1
    return EVAL_RESULT_D( gsl_sf_bessel_Y1_e, x, "gsl_sf_bessel_Y1" );
} // gsl_sf_bessel_Y1

// ----------------------------------------------------------------------------
// EOF SF-BesselY1.mjs
