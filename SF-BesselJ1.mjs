// SF-BesselJ1.mjs
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

import { M_SQRT2 }                         from "./SF-Math.mjs";
import { GSL_DBL_MIN }                     from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }                 from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }            from "./SF-Machine.mjs";
import { gsl_sf_bessel_sin_pi4_e }         from "./SF-Bessel.mjs";
import { cheb_eval_e }                     from "./SF-Chebyshev.mjs";
import { gsl_sf_bessel_amp_phase_bm1_cs }  from "./SF-BesselAF.mjs";
import { gsl_sf_bessel_amp_phase_bth1_cs } from "./SF-BesselAF.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const ROOT_EIGHT = 2.0 * M_SQRT2;

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*


// based on SLATEC besj1, 1983 version, w. fullerton

// chebyshev expansions
//
// series for bj1        on the interval  0.          to  1.60000d+01
//                                        with weighted error   4.48e-17
//                                         log weighted error  16.35
//                               significant figures required  15.77
//                                    decimal places required  16.89
//

const bj1_data =//: CONSTANT Series(0..11) = --[12]
    [
   -0.11726141513332787,
   -0.25361521830790640,
    0.050127080984469569,
   -0.004631514809625081,
    0.000247996229415914,
   -0.000008678948686278,
    0.000000214293917143,
   -0.000000003936093079,
    0.000000000055911823,
   -0.000000000000632761,
    0.000000000000005840,
   -0.000000000000000044
    ];
const bj1_cs = { length: 11, c: bj1_data, order: 11, a: -1.0, b: 1.0, order_sp: 8 };

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_J1_e(x)
{
    var y = Math.abs(x);
    var z = Math.abs(x);
    var sqrty = Math.abs(x);
    var ampl  = Math.abs(x);

    var c  = { val: 0.0, err: 0.0 }; // Result;
    var ca = { val: 0.0, err: 0.0 }; // Result;
    var ct = { val: 0.0, err: 0.0 }; // Result;
    var sp = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (y == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (y < 2.0 * GSL_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (y < ROOT_EIGHT * GSL_SQRT_DBL_EPSILON)
    {
        r.val = 0.5 * x;
        r.err = 0.0;
    }
    else if (y < 4.0)
    {
        c = cheb_eval_e(bj1_cs, 0.125 * y * y - 1.0);
        r.val = x * (0.25 + c.val);
        r.err = Math.abs(x * c.err);
    }
    else
    {
        // Because the leading term in the phase is y,
        // which we assume is exactly known, the error
        // in the cos() evaluation is bounded.
        //
        z = 32.0 / (y * y) - 1.0;
        ca = cheb_eval_e(gsl_sf_bessel_amp_phase_bm1_cs,  z);
        ct = cheb_eval_e(gsl_sf_bessel_amp_phase_bth1_cs, z);
        sp = gsl_sf_bessel_sin_pi4_e(y, ct.val / y);
        sqrty = Math.sqrt(y);
        ampl  = (0.75 + ca.val) / sqrty;
        if (x < 0.0)
        {
            r.val = -ampl * sp.val;
        }
        else
        {
            r.val =  ampl * sp.val;
        }
        r.err = Math.abs(sp.val) * ca.err / sqrty + Math.abs(ampl) * sp.err;
        r.err = r.err + GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

 } // gsl_sf_bessel_J1_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_J1( x )
{ // gsl_sf_bessel_J1
    return EVAL_RESULT_D( gsl_sf_bessel_J1_e, x, "gsl_sf_bessel_J1" );
} // gsl_sf_bessel_J1

// ----------------------------------------------------------------------------
// EOF SF-BesselJ1.mjs
