// SF-Chebyshev.mjs
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
// Translation to JavaScript: Igor Izvarin

import { GSL_PREC_DOUBLE } from "./SF-Mode.mjs";

export function cheb_eval_e(cs, x)
{
    var d      = 0.0;
    var dd     = 0.0;
    var y      = 0.0;
    var y2     = 0.0;
    var e      = 0.0;
    var temp   = 0.0;
    var j      =0;
    var r = { val: 0.0, err: 0.0 };   // Result;

    y  = (2.0 * x - cs.a - cs.b) / (cs.b - cs.a);
    y2 = 2.0 * y;

    e = 0.0;

    for (j = cs.order; j >= 1; j--)
    {
        temp = d;
        d = y2 * d - dd + cs.c[j];
        e = e + Math.abs(y2 * temp) + Math.abs(dd) + Math.abs(cs.c[j]);
        dd = temp;
    }

    temp = d;
    d = y * d - dd + 0.5 * cs.c[0];
    e = e + Math.abs(y * temp) + Math.abs(dd) + 0.5 * Math.abs(cs.c[0]);

    r.val = d;
    r.err = Number.EPSILON * e + Math.abs(cs.c[cs.order]);

    return r;

} // cheb_eval_e

// ----------------------------------------------------------------------------

export function cheb_eval_mode_e(cs, x, mode)
{
    var d          = 0.0;
    var dd         = 0.0;
    var y          = 0.0;
    var y2         = 0.0;
    var temp       = 0.0;
    var eval_order = 0;
    var j      =0;
    var r          = { val: 0.0, err: 0.0 };   // Result;

    y  = (2.0 * x - cs.a - cs.b) / (cs.b - cs.a);
    y2 = 2.0 * y;

    if (mode == GSL_PREC_DOUBLE)
    {
        eval_order = cs.order;
    }
    else
    {
        eval_order = cs.order_sp;
    }

    for (j = eval_order; j >= 1; j--)
    {
        temp = d;
        d = y2 * d - dd + cs.c[j];
        dd = temp;
    }

    r.val = y * d - dd + 0.5 * cs.c[0];
    r.err = Number.EPSILON * Math.abs(r.val) + Math.abs(cs.c[eval_order]);

    return r;

} // cheb_eval_mode_e

// ----------------------------------------------------------------------------
// EOF SF-Chebyshev.mjs
