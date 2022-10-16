// SF-Elementary.mjs
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

import { GSL_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MAX } from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }  from "./SF-Machine.mjs";

import { EVAL_RESULT_DD }   from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

export function gsl_sf_multiply_e(x, y)
{
    var ax = Math.abs(x);
    var ay = Math.abs(y);

    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0 || y == 0.0)
    {
        // It is necessary to eliminate this immediately.
        r.val = 0.0;
        r.err = 0.0;
    }
    else if ((ax <= 1.0 && ay >= 1.0) || (ay <= 1.0 && ax >= 1.0))
    {
        // Straddling 1.0 is always safe.
        r.val = x * y;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var f   = 1.0 - 2.0 * GSL_DBL_EPSILON;
        var min = Math.min(Math.abs(x), Math.abs(y));
        var max = Math.max(Math.abs(x), Math.abs(y));

        if (max < 0.9 * GSL_SQRT_DBL_MAX || min < (f * GSL_DBL_MAX) / max)
        {
            r.val = x * y; //GSL_COERCE_DBL(x * y);
            r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            //CHECK_UNDERFLOW(result);
        }
        else
        {
            throw "SF.OverflowException";
        }
    }

    return r;

} // gsl_sf_multiply_e

// ----------------------------------------------------------------------------

export function gsl_sf_multiply_err_e(x, dx, y, dy)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_multiply_e(x, y);
    r.err = r.err + Math.abs(dx * y) + Math.abs(dy * x);
    return r;

} // gsl_sf_multiply_err_e


// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_multiply( x, y )
{ // gsl_sf_multiply
    return EVAL_RESULT_DD( gsl_sf_multiply_e, { x: x, y: y }, "gsl_sf_multiply" );
} // gsl_sf_multiply

// ----------------------------------------------------------------------------
// EOF SF-Elementary.mjs
