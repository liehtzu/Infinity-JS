// SF-Lambert.mjs
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

import { GSL_DBL_EPSILON }  from "./SF-Machine.mjs";
import { M_E }              from "./SF-Math.mjs";

import { EVAL_RESULT_D }     from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// Started with code donated by K. Briggs; added
// error estimates, GSL foo, and minor tweaks.
// Some Lambert-ology from
//  [Corless, Gonnet, Hare, and Jeffrey, "On Lambert's W Function".]
//

// Halley iteration (eqn. 5.12, Corless et al)
function halley_iteration(x, w_initial, max_iters)
// *************************************************************************************************
{
    var w   = 0.0;
    var e   = 0.0;
    var p   = 0.0;
    var t   = 0.0;
    var tol = 0.0;
    var i = 0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    w = w_initial;
   
    for (i = 0; i <= max_iters - 1; i++)
    {
        e = Math.exp(w);
        p = w + 1.0;
        t = w * e - x;
        // printf("FOO: %20.16g  %20.16g\n", w, t);
       
        if (w > 0.0)
        {
            t = (t / p) / e;  // Newton iteration
        }
        else
        {
            t = t / (e * p - 0.5 * (p + 1.0) * t / p);  // Halley iteration
        }
       
        w = w - t;
       
        tol = 10.0 * GSL_DBL_EPSILON * Math.max(Math.abs(w), 1.0 / (Math.abs(p) * e));
       
        if (Math.abs(t) < tol)
        {
            r.val = w;
            r.err = 2.0 * tol;
            return r;
        }
    }
   
    // should never get here
    r.val = w;
    r.err = Math.abs(w);

    return r;

} // halley_iteration

// ----------------------------------------------------------------------------

// series which appears for q near zero;
// only the argument is different for the different branches
//
function series_eval(r)
// *****************************************************
{

    const c =//: CONSTANT ARRAY(0..11) OF LONG_FLOAT =
        [
       -1.0,
        2.331643981597124203363536062168,
       -1.812187885639363490240191647568,
        1.936631114492359755363277457668,
       -2.353551201881614516821543561516,
        3.066858901050631912893148922704,
       -4.175335600258177138854984177460,
        5.858023729874774148815053846119,
       -8.401032217523977370984161688514,
        12.250753501314460424,
       -18.100697012472442755,
        27.029044799010561650
        ];
    var t_8 = 0.0;
    var t_5 = 0.0;
    var t_1 = 0.0;

    t_8 = c[8] + r*(c[9] + r*(c[10] + r*c[11]));
    t_5 = c[5] + r*(c[6] + r*(c[7]  + r*t_8));
    t_1 = c[1] + r*(c[2] + r*(c[3]  + r*(c[4] + r*t_5)));
    return c[0] + r * t_1;

} // series_eval

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_lambert_W0_e(x)
// ****************************************************************************
{
    const MAX_ITERS = 10;
    const one_over_E = 1.0 / M_E;
    var q  = x + one_over_E;
    var s  = 0.0;
    var w  = 0.0;
    var p  = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (q < 0.0)
    {
        // Strictly speaking this is an error. But because of the
        // arithmetic operation connecting x and q, I am a little
        // lenient in case of some epsilon overshoot. The following
        // answer is quite accurate in that case. Anyway, we have
        // to code = GSL_EDOM.
        //
        r.val = -1.0;
        r.err = Math.sqrt(-q);
    }
    else if (q == 0.0)
    {
        r.val = -1.0;
        r.err = GSL_DBL_EPSILON; // cannot error is zero, maybe q = 0 by "accident"
    }
    else if (q < 1.0e-03)
    {
        // series near -1/E in sqrt(q)
        s = Math.sqrt(q);
        r.val = series_eval(s);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        if (x < 1.0)
        {
            // obtain initial approximation from series near x=0;
            // no need for extra care, since the Halley iteration
            // converges nicely on this branch
            //
            p = Math.sqrt(2.0 * M_E * q);
            w = -1.0 + p * (1.0 + p * (-1.0 / 3.0 + p * 11.0 / 72.0)); 
        }
        else
        {
            // obtain initial approximation from rough asymptotic */
            w = Math.log(x);
            if (x > 3.0)
            {
                w = w - Math.log(w);
            }
        }
       
        r = halley_iteration(x, w, MAX_ITERS);
    }

    return r;

} // gsl_sf_lambert_W0_e

// ----------------------------------------------------------------------------

export function gsl_sf_lambert_Wm1_e(x)
// ****************************************************************************
{
    const MAX_ITERS = 32;
    var one_over_E = 0.0;
    var q   = 0.0;
    var w   = 0.0;
    var s   = 0.0;
    var L_1 = 0.0;
    var L_2 = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x > 0.0)
    {
        r = gsl_sf_lambert_W0_e(x);
    }
    else if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else
    {
        one_over_E = 1.0 / M_E;
        q = x + one_over_E;
       
        if (q < 0.0)
        {
            // As in the W0 branch above, code = some reasonable answer anyway.
            r.val = -1.0; 
            r.err =  Math.sqrt(-q);
            return r;
        }
       
        if (x < -1.0e-6)
        {
            // Obtain initial approximation from series about q = 0,
            // as long as we're not very close to x = 0.
            // Use full series and try to bail out if q is too small,
            // since the Halley iteration has bad convergence properties
            // in finite arithmetic for q very small, because the
            // increment alternates and p is near zero.
            //
            s = -Math.sqrt(q);
            w = series_eval(s);
            if (q < 3.0e-3)
            {
                // this approximation is good enough
                r.val = w;
                r.err = 5.0 * GSL_DBL_EPSILON * Math.abs(w);
                return r;
            }
        }
        else
        {
            // Obtain initial approximation from asymptotic near zero.
            L_1 = Math.log(-x);
            L_2 = Math.log(-L_1);
            w = L_1 - L_2 + L_2 / L_1;
        }
       
        r = halley_iteration(x, w, MAX_ITERS);
    }

    return r;

} // gsl_sf_lambert_Wm1_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_lambert_W0( x )
{ // gsl_sf_lambert_W0
    return EVAL_RESULT_D( gsl_sf_lambert_W0_e, x, "gsl_sf_lambert_W0" );
} // gsl_sf_lambert_W0;

export function gsl_sf_lambert_Wm1( x )
{ // gsl_sf_lambert_Wm1
    return EVAL_RESULT_D( gsl_sf_lambert_Wm1_e, x, "gsl_sf_lambert_Wm1" );
} // gsl_sf_lambert_Wm1

// ----------------------------------------------------------------------------
// SF-Lambert.mjs
