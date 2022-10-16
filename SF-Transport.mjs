// SF-Transport.mjs
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

import { GSL_DBL_EPSILON }        from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }   from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON }    from "./SF-Machine.mjs";
import { cheb_eval_e }            from "./SF-Chebyshev.mjs";

import { EVAL_RESULT_D }          from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const transport2_data =
    [
    1.671760446434538503,
   -0.147735359946794490,
    0.148213819946936338e-01,
   -0.14195330326305613e-02,
    0.1306541324415708e-03,
   -0.117155795867579e-04,
    0.10333498445756e-05,
   -0.901911304223e-07,
    0.78177169833e-08,
   -0.6744565684e-09,
    0.579946394e-10,
   -0.49747619e-11,
    0.425961e-12,
   -0.36422e-13,
    0.3111e-14,
   -0.265e-15,
    0.23e-16,
   -0.19e-17
    ];
const transport2_cs = { length: 17, c: transport2_data, order: 17, a: -1.0, b: 1.0, order_sp: 9 };

const transport3_data =
    [
    0.762012543243872007,
   -0.105674387705058533,
    0.119778084819657810e-01,
   -0.12144015203698307e-02,
    0.1155099769392855e-03,
   -0.105815992124423e-04,
    0.9474663385302e-06,
   -0.836221212858e-07,
    0.73109099278e-08,
   -0.6350594779e-09,
    0.549118282e-10,
   -0.47321395e-11,
    0.4067695e-12,
   -0.348971e-13,
    0.29892e-14,
   -0.256e-15,
    0.219e-16,
   -0.19e-17
    ];
const transport3_cs = { length: 17, c: transport3_data, order: 17, a: -1.0, b: 1.0, order_sp: 9 };

const transport4_data =
    [
    0.4807570994615110579,
   -0.8175378810321083956e-01,
    0.1002700665975162973e-01,
   -0.10599339359820151e-02,
    0.1034506245030405e-03,
   -0.96442705485899e-05,
    0.8745544408515e-06,
   -0.779321207981e-07,
    0.68649886141e-08,
   -0.5999571076e-09,
    0.521366241e-10,
   -0.45118382e-11,
    0.3892159e-12,
   -0.334936e-13,
    0.28767e-14,
   -0.2467e-15,
    0.211e-16,
   -0.18e-17
    ];
const transport4_cs = { length: 17, c: transport4_data, order: 17, a: -1.0, b: 1.0, order_sp: 9 };

const transport5_data =
    [
    0.347777777133910789,
   -0.66456988976050428e-01,
    0.8611072656883309e-02,
   -0.9396682223755538e-03,
    0.936324806081513e-04,
   -0.88571319340833e-05,
    0.811914989145e-06,
   -0.72957654233e-07,
    0.646971455e-08,
   -0.568490283e-09,
    0.49625598e-10,
   -0.4310940e-11,
    0.373100e-12,
   -0.32198e-13,
    0.2772e-14,
   -0.238e-15,
    0.21e-16,
   -0.18e-17
    ];
const transport5_cs = { length: 17, c: transport5_data, order: 17, a: -1.0, b: 1.0, order_sp: 9 };

// ----------------------------------------------------------------------------

function transport_sumexp(numexp, order, t, x)
{
    var rk     = 0.0;
    var sumexp = 0.0;
    var sum2   = 0.0;
    var xk     = 0.0;
    var xk1    = 0.0;
    var k = 0;
    var j = 0;

    rk = (numexp);
    sumexp = 0.0;
    for (k = 1; k <= numexp; k++)
    {
        sum2 = 1.0;
        xk   = 1.0 / (rk * x);
        xk1  = 1.0;
        for (j = 1; j <= order; j++)
        {
            sum2 = sum2 * xk1 * xk + 1.0;
            xk1  = xk1 + 1.0;
        }
        sumexp = sumexp * t;
        sumexp = sumexp + sum2;
        rk = rk - 1.0;
    }
    return sumexp;

} // transport_sumexp

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_transport_2_e(x)
{
    const val_infinity = 3.289868133696452873;
    var t      = 0.0;
    var et     = 0.0;
    var sumexp = 0.0;
    var numexp = 0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 3.0*GSL_SQRT_DBL_EPSILON)
    {
        r.val = x;
        r.err = GSL_DBL_EPSILON * Math.abs(x) + x * x / 2.0;
        return r;
    }
    else if (x <= 4.0)
    {
        t = (x * x / 8.0 - 0.5) - 0.5;
        c = cheb_eval_e(transport2_cs, t);
        r.val = x * c.val;
        r.err = x * c.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }
    else if (x < -GSL_LOG_DBL_EPSILON)
    {
        numexp = Math.trunc((-GSL_LOG_DBL_EPSILON) / x) + 1;
        sumexp = transport_sumexp(numexp, 2, Math.exp(-x), x);
        t = 2.0 * Math.log(x) - x + Math.log(sumexp);
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + Math.abs(t) * et);
        }
        return r;
    }
    else if (x < 2.0 / GSL_DBL_EPSILON)
    {
        numexp = 1;
        sumexp = transport_sumexp(numexp, 2, 1.0, x);
        t = 2.0 * Math.log(x) - x + Math.log(sumexp);
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t)+1.0) * et);
        }
        return r;
    }
    else
    {
        t = 2.0 * Math.log(x) - x;
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }

} // gsl_sf_transport_2_e

// ----------------------------------------------------------------------------

export function gsl_sf_transport_3_e(x)
{
    const val_infinity = 7.212341418957565712;
    var t      = 0.0;
    var et     = 0.0;
    var x2     = 0.0;
    var sumexp = 0.0;
    var numexp = 0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else if (x < 3.0 * GSL_SQRT_DBL_EPSILON)
    {
        r.val = 0.5 * x * x;
        r.err = 2.0 * GSL_DBL_EPSILON * r.val;
        //CHECK_UNDERFLOW(result);
        return r;
    }
    else if (x <= 4.0)
    {
        x2 = x * x;
        t = (x2 / 8.0 - 0.5) - 0.5;
        c = cheb_eval_e(transport3_cs, t);
        r.val = x2 * c.val;
        r.err = x2 * c.err;
        r.err = r.err + GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }
    else if (x < -GSL_LOG_DBL_EPSILON)
    {
        numexp = Math.trunc((-GSL_LOG_DBL_EPSILON)/x) + 1;
        sumexp = transport_sumexp(numexp, 3, Math.exp(-x), x);
        t = 3.0 * Math.log(x) - x + Math.log(sumexp);
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + Math.abs(t) * et);
        }
        return r;
    }
    else if (x < 3.0/GSL_DBL_EPSILON)
    {
        numexp = 1;
        sumexp = transport_sumexp(numexp, 3, 1.0, x);
        t = 3.0 * Math.log(x) - x + Math.log(sumexp);
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }
    else
    {
        t = 3.0 * Math.log(x) - x;
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }

} // gsl_sf_transport_3_e

// ----------------------------------------------------------------------------

export function gsl_sf_transport_4_e(x)
{
    const val_infinity = 25.97575760906731660;
    var t      = 0.0;
    var et     = 0.0;
    var x2     = 0.0;
    var sumexp = 0.0;
    var numexp = 0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else if (x < 3.0 * GSL_SQRT_DBL_EPSILON)
    {
        r.val = x * x * x / 3.0;
        r.err = 3.0 * GSL_DBL_EPSILON * r.val;
        //CHECK_UNDERFLOW(result);
        return r;
    }
    else if (x <= 4.0)
    {
        x2 = x * x;
        t = (x2 / 8.0 - 0.5) - 0.5;
        c = cheb_eval_e(transport4_cs, t);
        r.val = x2 * x * c.val;
        r.err = x2 * x * c.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }
    else if (x < -GSL_LOG_DBL_EPSILON)
    {
        numexp = Math.trunc((-GSL_LOG_DBL_EPSILON) / x) + 1;
        sumexp = transport_sumexp(numexp, 4, Math.exp(-x), x);
        t = 4.0 * Math.log(x) - x + Math.log(sumexp);
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }
    else if (x < 3.0 / GSL_DBL_EPSILON)
    {
        numexp = 1;
        sumexp = transport_sumexp(numexp, 4, 1.0, x);
        t = 4.0 * Math.log(x) - x + Math.log(sumexp);
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }
    else
    {
        t = 4.0 * Math.log(x) - x;
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }

} // gsl_sf_transport_4_e

// ----------------------------------------------------------------------------

export function gsl_sf_transport_5_e(x)
{
    const val_infinity = 124.4313306172043912;
    var t      = 0.0;
    var et     = 0.0;
    var x2     = 0.0;
    var sumexp = 0.0;
    var numexp = 0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else if (x < 3.0 * GSL_SQRT_DBL_EPSILON)
    {
        r.val = x * x * x * x / 4.0;
        r.err = 4.0 * GSL_DBL_EPSILON * r.val;
        //CHECK_UNDERFLOW(result);
        return r;
    }
    else if (x <= 4.0)
    {
        x2 = x * x;
        t = (x2 / 8.0 - 0.5) - 0.5;
        c = cheb_eval_e(transport5_cs, t);
        r.val = x2 * x2 * c.val;
        r.err = x2 * x2 * c.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }
    else if (x < -GSL_LOG_DBL_EPSILON)
    {
        numexp = Math.trunc((-GSL_LOG_DBL_EPSILON) / x) + 1;
        sumexp = transport_sumexp(numexp, 5, Math.exp(-x), x);
        t = 5.0 * Math.log(x) - x + Math.log(sumexp);
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }
    else if (x < 3.0 / GSL_DBL_EPSILON)
    {
        numexp = 1;
        sumexp = transport_sumexp(numexp, 5, 1.0, x);
        t = 5.0 * Math.log(x) - x + Math.log(sumexp);
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }
    else
    {
        t = 5.0 * Math.log(x) - x;
        if (t < GSL_LOG_DBL_EPSILON)
        {
            r.val = val_infinity;
            r.err = 2.0 * GSL_DBL_EPSILON * val_infinity;
        }
        else
        {
            et = Math.exp(t);
            r.val = val_infinity - et;
            r.err = 2.0 * GSL_DBL_EPSILON * (val_infinity + (Math.abs(t) + 1.0) * et);
        }
        return r;
    }

} // gsl_sf_transport_5_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_transport_2( x )
{ // gsl_sf_transport_2
     return EVAL_RESULT_D( gsl_sf_transport_2_e, x, "gsl_sf_transport_2" );
} // gsl_sf_transport_2

export function gsl_sf_transport_3( x )
{ // gsl_sf_transport_3
     return EVAL_RESULT_D( gsl_sf_transport_3_e, x, "gsl_sf_transport_3" );
} // gsl_sf_transport_3

export function gsl_sf_transport_4( x )
{ // gsl_sf_transport_4
     return EVAL_RESULT_D( gsl_sf_transport_4_e, x, "gsl_sf_transport_4" );
} // gsl_sf_transport_4

export function gsl_sf_transport_5( x )
{ // gsl_sf_transport_5
     return EVAL_RESULT_D( gsl_sf_transport_5_e, x, "gsl_sf_transport_5" );
} // gsl_sf_transport_5

// ----------------------------------------------------------------------------
// EOF SF-Transport.mjs
