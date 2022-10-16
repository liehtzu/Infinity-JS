// SF-Gegenbauer.mjs
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

import { GSL_DBL_EPSILON } from "./SF-Machine.mjs";

import { EVAL_RESULT_DD }  from "./SF-Evaluate.mjs";
import { EVAL_RESULT_IDD } from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// See: [Thompson, Atlas for Computing Mathematical Functions]

export function gsl_sf_gegenpoly_1_e(lambda, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
  
    if (lambda == 0.0)
    {
        r.val = 2.0 * x;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        r.val = 2.0 * lambda * x;
        r.err = 4.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_gegenpoly_1_e

// ----------------------------------------------------------------------------

export function gsl_sf_gegenpoly_2_e(lambda, x)
{
    var txx = 0.0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (lambda == 0.0)
    {
        txx = 2.0 * x * x;
        r.val = -1.0 + txx;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(txx);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        r.val = lambda * (-1.0 + 2.0 * (1.0 + lambda) * x * x);
        r.err = GSL_DBL_EPSILON * (2.0 * Math.abs(r.val) + Math.abs(lambda));
    }

    return r;

} // gsl_sf_gegenpoly_2_e

// ----------------------------------------------------------------------------

export function gsl_sf_gegenpoly_3_e(lambda, x)
{
    var c = 0.0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (lambda == 0.0)
    {
        r.val = x * (-2.0 + 4.0 / 3.0 * x * x);
        r.err = GSL_DBL_EPSILON * (2.0 * Math.abs(r.val) + Math.abs(x));
    }
    else
    {
        c = 4.0 + lambda * (6.0 + 2.0 * lambda);
        r.val = 2.0 * lambda * x * ( -1.0 - lambda + c * x * x / 3.0);
        r.err = GSL_DBL_EPSILON * (2.0 * Math.abs(r.val) + Math.abs(lambda * x));
    }

    return r;

} // gsl_sf_gegenpoly_3_e

// ----------------------------------------------------------------------------

export function gsl_sf_gegenpoly_n_e(n, lambda, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (lambda <= -0.5 || n < 0)
    {
        throw "SF.DomainException";
    }
    else if (n == 0)
    {
        r.val = 1.0;
        r.err = 0.0;
    }
    else if (n == 1)
    {
        r = gsl_sf_gegenpoly_1_e(lambda, x);
    }
    else if (n == 2)
    {
        r = gsl_sf_gegenpoly_2_e(lambda, x);
    }
    else if (n == 3)
    {
        r = gsl_sf_gegenpoly_3_e(lambda, x);
    }
    else
    {
        if (lambda == 0.0 && (x >= -1.0 || x <= 1.0))
        {
            // 2 T_n(x)/n
            var z = (n) * Math.acos(x);

            r.val = 2.0 * Math.cos(z) / (n);
            r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(z * r.val);
        }
        else
        {
            var g2 = { val: 0.0, err: 0.0 }; // Result;
            var g3 = { val: 0.0, err: 0.0 }; // Result;
            var gkm2 = 0.0;
            var gkm1 = 0.0;
            var gk   = 0.0;
            var k    = 0;

            g2 = gsl_sf_gegenpoly_2_e(lambda, x);
            g3 = gsl_sf_gegenpoly_3_e(lambda, x);
            gkm2 = g2.val;
            gkm1 = g3.val;
            for (k = 4; k <= n; k++)
            {
                gk = (2.0 * ((k) + lambda - 1.0) * x * gkm1 - ((k) + 2.0 * lambda - 2.0) * gkm2) / (k);
                gkm2 = gkm1;
                gkm1 = gk;
            }
            r.val = gk;
            r.err = 2.0 * GSL_DBL_EPSILON * 0.5 * (n) * Math.abs(gk);
        }
    }

    return r;

} // gsl_sf_gegenpoly_n_e

// ----------------------------------------------------------------------------

export function gsl_sf_gegenpoly_array(nmax, lambda, x, result_array)
{
    var term1 = 0.0;
    var term2 = 0.0;

    if (lambda <= -0.5 || nmax < 0)
    {
        throw "SF.DomainException";
    }
   
    // n = 0
    result_array[0] = 1.0;
    if (nmax == 0)
    {
        return;
    }
   
    // n = 1
    if (lambda == 0.0)
    {
        result_array[1] = 2.0 * x;
    }
    else
    {
        result_array[1] = 2.0 * lambda * x;
    }
   
    // n <= nmax
    for (let k = 2; k <= nmax; k++)
    {
        term1 = 2.0 * ((k) + lambda - 1.0) * x * result_array[k-1];
        term2 = ((k) + 2.0 * lambda - 2.0)     * result_array[k-2];
        result_array[k] = (term1 - term2) / (k);
    }

} // gsl_sf_gegenpoly_array

// //*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_gegenpoly_1( lambda, x )
{ //sl_sf_gegenpoly_1
    return EVAL_RESULT_DD( gsl_sf_gegenpoly_1_e, { x: lambda, y: x }, "gsl_sf_gegenpoly_1" );
} // gsl_sf_gegenpoly_1

export function gsl_sf_gegenpoly_2( lambda, x )
{ //sl_sf_gegenpoly_2
    return EVAL_RESULT_DD( gsl_sf_gegenpoly_2_e, { x: lambda, y: x }, "gsl_sf_gegenpoly_2" );
} // gsl_sf_gegenpoly_2

export function gsl_sf_gegenpoly_3( lambda, x )
{ //sl_sf_gegenpoly_3
    return EVAL_RESULT_DD( gsl_sf_gegenpoly_3_e, { x: lambda, y: x }, "gsl_sf_gegenpoly_3" );
} // gsl_sf_gegenpoly_3

export function gsl_sf_gegenpoly_n( n, lambda, x )
{ //sl_sf_gegenpoly_n
    return EVAL_RESULT_IDD( gsl_sf_gegenpoly_n_e, { i: n, x: lambda, y: x }, "gsl_sf_gegenpoly_n" );
} // gsl_sf_gegenpoly_n

// ----------------------------------------------------------------------------
// EOF SF-Gegenbauer.mjs
