// SF-Clausen.mjs
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
import { M_PI }                 from "./SF-Math.mjs";
import { cheb_eval_e }          from "./SF-Chebyshev.mjs";
import { gsl_sf_angle_restrict_pos_e } from "./SF-Trigonometric.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const aclaus_data =//: CONSTANT Series(0..14) := --[15]
    [
    2.142694363766688447e+00,
    0.723324281221257925e-01,
    0.101642475021151164e-02,
    0.3245250328531645e-04,
    0.133315187571472e-05,
    0.6213240591653e-07,
    0.313004135337e-08,
    0.16635723056e-09,
    0.919659293e-11,
    0.52400462e-12,
    0.3058040e-13,
    0.18197e-14,
    0.1100e-15,
    0.68e-17,
    0.4e-18
    ];
const aclaus_cs = { length: 14, c: aclaus_data, order: 14, a: -1.0, b: 1.0, order_sp: 8 };
// FIXME:  this is a guess, correct value needed here BJG


// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_clausen_e(x0)
{
    const p0  = 6.28125;
    const p1  = 0.19353071795864769253e-02;
    var x_cut = M_PI * GSL_SQRT_DBL_EPSILON;
    var sgn   = 1.0;
    var t     = 0.0;
    var x     = x0;
    var c     = { val: 0.0, err: 0.0 }; // Result;
    var r     = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0)
    {
        x   = -x;
        sgn = -1.0;
    }
   
    // Argument reduction to [0, 2pi)
    x = gsl_sf_angle_restrict_pos_e(x);
   
    // Further reduction to [0,pi)
    if (x > M_PI)
    {
        // simulated extra precision: 2PI := p0 + p1
        x = (p0 - x) + p1;
        sgn = -sgn;
    }
   
    if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (x < x_cut)
    {
        r.val = x * (1.0 - Math.log(x));
        r.err = x * GSL_DBL_EPSILON;
    }
    else
    {
        t = 2.0 * (x * x / (M_PI * M_PI) - 0.5);
        c = cheb_eval_e(aclaus_cs, t);
        r.val = x * (c.val - Math.log(x));
        r.err = x * (c.err + GSL_DBL_EPSILON);
    }
   
    r.val = r.val * sgn;
   
    return r;

} // gsl_sf_clausen_e


// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_clausen( x )
{ // gsl_sf_clausen
    return EVAL_RESULT_D( gsl_sf_clausen_e, x, "gsl_sf_clausen" );
} // gsl_sf_clausen

// ----------------------------------------------------------------------------
// EOF SF-Clausen.mjs
