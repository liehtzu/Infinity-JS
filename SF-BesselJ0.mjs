// SF-BesselJ0.mjs
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

import { GSL_DBL_EPSILON }                 from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }            from "./SF-Machine.mjs";
import { gsl_sf_bessel_cos_pi4_e }         from "./SF-Bessel.mjs";
import { cheb_eval_e }                     from "./SF-Chebyshev.mjs";
import { gsl_sf_bessel_amp_phase_bm0_cs }  from "./SF-BesselAF.mjs";
import { gsl_sf_bessel_amp_phase_bth0_cs } from "./SF-BesselAF.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*


// based on SLATEC besj0, 1977 version, w. fullerton

// chebyshev expansions for Bessel functions
//
// series for bj0        on the interval  0.          to  1.60000d+01
//                                        with weighted error   7.47e-18
//                                         log weighted error  17.13
//                               significant figures required  16.98
//                                    decimal places required  17.68
//

const bj0_data =//: CONSTANT Series(0..12) := --[13]
    [
    0.100254161968939137, 
   -0.665223007764405132, 
    0.248983703498281314, 
   -0.0332527231700357697,
    0.0023114179304694015,
   -0.0000991127741995080,
    0.0000028916708643998,
   -0.0000000612108586630,
    0.0000000009838650793,
   -0.0000000000124235515,
    0.0000000000001265433,
   -0.0000000000000010619,
    0.0000000000000000074
    ];
const bj0_cs = { length: 12, c: bj0_data, order: 12, a: -1.0, b: 1.0, order_sp: 9 };


// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_J0_e(x)
// ****************************************************************************
{
    var y = Math.abs(x);

    var r = { val: 0.0, err: 0.0 }; // Result;

    if (y < 2.0 * GSL_SQRT_DBL_EPSILON)
    {
        r.val = 1.0;
        r.err = y * y;
    }
    else if (y <= 4.0)
    {
        r = cheb_eval_e(bj0_cs, 0.125 * y * y - 1.0);
    }
    else
    {
        var z     = 0.0;
        var sqrty = 0.0;
        var ampl  = 0.0;
        var ca = { val: 0.0, err: 0.0 }; // Result;
        var ct = { val: 0.0, err: 0.0 }; // Result;
        var cp = { val: 0.0, err: 0.0 }; // Result;

        z = 32.0 / (y * y) - 1.0;
        ca = cheb_eval_e(gsl_sf_bessel_amp_phase_bm0_cs,  z);
        ct = cheb_eval_e(gsl_sf_bessel_amp_phase_bth0_cs, z);
        cp = gsl_sf_bessel_cos_pi4_e(y, ct.val / y);
        sqrty = Math.sqrt(y);
        ampl  = (0.75 + ca.val) / sqrty;
        r.val = ampl * cp.val;
        r.err = Math.abs(cp.val) * ca.err / sqrty + Math.abs(ampl) * cp.err;
        r.err = r.err + GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_J0_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_J0( x )
{ // gsl_sf_bessel_J0
    return EVAL_RESULT_D( gsl_sf_bessel_J0_e, x, "gsl_sf_bessel_J0" );
} // gsl_sf_bessel_J0

// ----------------------------------------------------------------------------
// EOF SF-BesselJ0.mjs
