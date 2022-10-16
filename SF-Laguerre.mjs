// SF-Laguerre.mjs
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

import { GSL_DBL_MAX }           from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_IS_ODD }            from "./SF-Math.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { gsl_sf_lngamma_e }      from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_sgn_e }  from "./SF-Gamma.mjs";
import { gsl_sf_lnfact_e }       from "./SF-Gamma.mjs";
import { gsl_sf_taylorcoeff_e }  from "./SF-Gamma.mjs";
import { gsl_sf_exp_mult_err_e } from "./SF-Exponential.mjs";
import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";
import { EVAL_RESULT_IDD }       from "./SF-Evaluate.mjs";

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*


// based on the large 2b-4a asymptotic for 1F1
// [Abramowitz+Stegun, 13.5.21]
// L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
//
// The second term (ser_term2) is from Slater,"The Confluent
// Hypergeometric Function" p.73.  I think there may be an error in
// the first term of the expression given there, comparing with AS
// 13.5.21 (cf sin(a\pi+\Theta) vs sin(a\pi) + sin(\Theta)) - but the
// second term appears correct.
//
function laguerre_large_n( n, alpha, x )
{
    var a      = 0.0;
    var b      = 0.0;
    var eta    = 0.0;
    var cos2th = 0.0;
    var sin2th = 0.0;
    var eps    = 0.0;
    var pre_h  = 0.0;
    var pre_term1 = 0.0;
    var pre_term2 = 0.0;
    var lnpre_val = 0.0;
    var lnpre_err = 0.0;
    var ser_term1 = 0.0;
    var ser_term2 = 0.0;
    var phi1      = 0.0;
    var A1        = 0.0;
    var ser_val   = 0.0;
    var ser_err   = 0.0;

    var lnfact = { val: 0.0, err: 0.0 }; // Result;
    var lg_b   = { val: 0.0, err: 0.0 }; // Result;
    var r      = { val: 0.0, err: 0.0 }; // Result;

    a         = -(n);
    b         = alpha + 1.0;
    eta       = 2.0 * b - 4.0 * a;
    cos2th    = x / eta;
    sin2th    = 1.0 - cos2th;
    eps       = Math.asin(Math.sqrt(cos2th)); // theta = pi/2 - eps
    pre_h     = 0.25 * M_PI * M_PI * eta * eta * cos2th * sin2th;
    lg_b      = gsl_sf_lngamma_e(b + (n));
    lnfact    = gsl_sf_lnfact_e(n);
    pre_term1 = 0.5 * (1.0 - b) * Math.log(0.25 * x * eta);
    pre_term2 = 0.25 * Math.log(pre_h);
    lnpre_val = lg_b.val - lnfact.val + 0.5 * x + pre_term1 - pre_term2;
    lnpre_err = lg_b.err + lnfact.err + GSL_DBL_EPSILON * (Math.abs(pre_term1) + Math.abs(pre_term2));

    phi1      = 0.25 * eta * (2.0 * eps + Math.sin(2.0 * eps));
    ser_term1 = -Math.sin(phi1);

    A1        = (1.0 / 12.0) * (5.0 / (4.0 * sin2th) + (3.0 * b * b - 6.0 * b + 2.0) * sin2th - 1.0);
    ser_term2 = -A1 * Math.cos(phi1) / (0.25 * eta * Math.sin(2.0 * eps));

    ser_val   = ser_term1 + ser_term2;
    ser_err   = ser_term2 * ser_term2 + GSL_DBL_EPSILON * (Math.abs(ser_term1) + Math.abs(ser_term2));

    r = gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err, ser_val, ser_err);
    r.err = r.err + 2.0 * GSL_SQRT_DBL_EPSILON * Math.abs(r.val);

    return r;

} // laguerre_large_n

// ----------------------------------------------------------------------------

// Evaluate polynomial based on confluent hypergeometric representation.
//
// L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
//
// assumes n > 0 and a != negative integer greater than -n
//
function laguerre_n_cp(n, a, x)
{
    var poly_1F1_val = 0.0;
    var poly_1F1_err = 0.0;
    var lnpre_val    = 0.0;
    var lnpre_err    = 0.0;
    var t            = 0.0;
    var u            = 0.0;
    var k = 0;

    var lnfact = { val: 0.0, err: 0.0 }; // Result;
    var lg1    = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var lg2    = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var r      = { val: 0.0, err: 0.0 }; // Result;

    lnfact = gsl_sf_lnfact_e(n);
    lg1 = gsl_sf_lngamma_sgn_e(a + 1.0 + (n)); //, lg1, s1);
    lg2 = gsl_sf_lngamma_sgn_e(a + 1.0); //, lg2, s2);
    poly_1F1_val = 1.0;
    poly_1F1_err = 0.0;

    lnpre_val = (lg1.val - lg2.val) - lnfact.val;
    lnpre_err = lg1.err + lg2.err + lnfact.err + 2.0 * GSL_DBL_EPSILON * Math.abs(lnpre_val);

    for (k = n - 1; k >= 0; k--)
    {
        t = (-n + k) / (a + 1.0 + (k)) * (x / (k + 1));
        u = t + 1.0 / poly_1F1_val;
        if (u > 0.9 * GSL_DBL_MAX / poly_1F1_val)
        {
            // internal error only, don't call the error handler
            // INTERNAL_OVERFLOW_ERROR(result);
            throw "SF.OverflowException";
        }
        else
        {
            // Collect the Horner terms.
            poly_1F1_val = 1.0 + t * poly_1F1_val;
            poly_1F1_err = poly_1F1_err + GSL_DBL_EPSILON + Math.abs(t) * poly_1F1_err;
        }
    }

    r = gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err, poly_1F1_val, poly_1F1_err);

    return r;

} // laguerre_n_cp

// ----------------------------------------------------------------------------

// Evaluate the polynomial based on the confluent hypergeometric
// function in a safe way, with no restriction on the arguments.
//
// assumes x != 0
//
function laguerre_n_poly_safe(n, a, x)
{
    var b  = 0.0;
    var mx = 0.0;
    var tc_sgn  = 0.0;
    var term    = 0.0;
    var sum_val = 0.0;
    var sum_err = 0.0;
    var k = 0;

    var tc = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        b  = a + 1.0;
        mx = -x;
        if (x < 0.0)
        {
            tc_sgn = 1.0;
        }
        else
        {
            if (GSL_IS_ODD(n))
            {
                tc_sgn = -1.0;
            }
            else
            {
                tc_sgn = 1.0;
            }
        }
        tc = gsl_sf_taylorcoeff_e(n, Math.abs(x));

        term = tc.val * tc_sgn;
        sum_val = term;
        sum_err = tc.err;
        for (k = n - 1; k >= 0; k--)
        {
            term    = term * ((b + (k)) / (n - k)) * (k + 1) / mx;
            sum_val = sum_val + term;
            sum_err = sum_err + 4.0 * GSL_DBL_EPSILON * Math.abs(term);
        }
        r.val = sum_val;
        r.err = sum_err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }
    catch (e)
    {
        r.val = 0.0; // FIXME: should be Inf
        r.err = 0.0;
        return r;
    }

} // laguerre_n_poly_safe

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_laguerre_1_e(a, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r.val = 1.0 + a - x;
    r.err = 2.0 * GSL_DBL_EPSILON * (1.0 + Math.abs(a) + Math.abs(x));

    return r;

} // gsl_sf_laguerre_1_e

// ----------------------------------------------------------------------------

export function gsl_sf_laguerre_2_e(a, x)
{
    var c0 = 0.0;
    var c1 = 0.0;
    var c2 = 0.0;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (a == -2.0)
    {
        r.val = 0.5 * x * x;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        c0 = 0.5 * (2.0 + a) * (1.0 + a);
        c1 = -(2.0 + a);
        c2 = -0.5 / (2.0 + a);
        r.val = c0 + c1 * x * (1.0 + c2 * x);
        r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs(c0) + 2.0 * Math.abs(c1 * x) * (1.0 + 2.0 * Math.abs(c2 * x)));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_laguerre_2_e

// ----------------------------------------------------------------------------

export function gsl_sf_laguerre_3_e(a, x)
{
    var c0 = 0.0;
    var c1 = 0.0;
    var c2 = 0.0;
    var c3 = 0.0;
    var x2_6 = 0.0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (a == -2.0)
    {
        x2_6 = x * x / 6.0;
        r.val = x2_6 * (3.0 - x);
        r.err = x2_6 * (3.0 + Math.abs(x)) * 2.0 * GSL_DBL_EPSILON;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (a == -3.0)
    {
        r.val = -x * x / 6.0;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        c0 = (3.0 + a) * (2.0 + a) * (1.0 + a) / 6.0;
        c1 = -c0 * 3.0 / (1.0 + a);
        c2 = -1.0 / (2.0 + a);
        c3 = -1.0 / (3.0 * (3.0 + a));
        r.val = c0 + c1 * x * (1.0 + c2 * x * (1.0 + c3 * x));
        r.err = 1.0 + 2.0 * Math.abs(c3 * x);
        r.err = 1.0 + 2.0 * Math.abs(c2 * x) * r.err;
        r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs(c0) + 2.0 * Math.abs(c1 * x) * r.err);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_laguerre_3_e

// ----------------------------------------------------------------------------

export function gsl_sf_laguerre_n_e(n, a, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (n < 0)
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
        r.val = 1.0 + a - x;
        r.err = 2.0 * GSL_DBL_EPSILON * (1.0 + Math.abs(a) + Math.abs(x));
    }
    else if (x == 0.0)
    {
        var product = a + 1.0;
        var k = 0;

        for (k = 2; k <= n; k++)
        {
            product = product * (a + (k)) / (k);
        }
        r.val = product;
        r.err = 2.0 * (n + 1) * GSL_DBL_EPSILON * Math.abs(product) + GSL_DBL_EPSILON;
    }
    else if (x < 0.0 && a > -1.0)
    {
        // In this case all the terms in the polynomial
        // are of the same sign. Note that this also
        // catches overflows correctly.
        //
        r = laguerre_n_cp(n, a, x);
    }
    else if (n < 5 || (x > 0.0 && a < (-n-1)))
    {
        // Either the polynomial will not lose too much accuracy
        // or all the terms are negative. In any case,
        // the error estimate here is good. We try both
        // explicit summation methods, as they have different
        // characteristics. One may underflow/overflow while the
        // other does not.
        //
        try
        {
            r = laguerre_n_cp(n, a, x);
        }
        catch (e)
        {
            r = laguerre_n_poly_safe(n, a, x);
        }
    }
    else if (n > 10000000 && x > 0.0 && a > -1.0 && x < 2.0 * (a + 1.0) + 4.0 * (n))
    {
        r = laguerre_large_n(n, a, x);
    }
    else if (a >= 0.0 || (x > 0.0 && a < (-n-1)))
    {
        var lg2  = { val: 0.0, err: 0.0 }; // Result;
        var Lkm1 = 0.0;
        var Lk   = 0.0;
        var Lkp1 = 0.0;

        lg2  = gsl_sf_laguerre_2_e(a, x);
        Lkm1 = 1.0 + a - x;
        Lk   = lg2.val;

        for (let k = 2; k <= n - 1; k++)
        {
            Lkp1 = (-((k) + a) * Lkm1 + (2.0 * (k) + a + 1.0 - x) * Lk) / (k + 1);
            Lkm1 = Lk;
            Lk   = Lkp1;
        }
        r.val = Lk;
        r.err = (Math.abs(lg2.err / lg2.val) + GSL_DBL_EPSILON) * (n) * Math.abs(Lk);
    }
    else
    {
        // Despair... or magic?
        r = laguerre_n_poly_safe(n, a, x);
    }

    return r;

} // gsl_sf_laguerre_n_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_laguerre_1( a, x )
{ // gsl_sf_laguerre_1
    return EVAL_RESULT_DD( gsl_sf_laguerre_1_e, { x: a, y: x }, "gsl_sf_laguerre_1" );
} // gsl_sf_laguerre_1

export function gsl_sf_laguerre_2( a, x )
{ // gsl_sf_laguerre_2
    return EVAL_RESULT_DD( gsl_sf_laguerre_2_e, { x: a, y: x }, "gsl_sf_laguerre_2" );
} // gsl_sf_laguerre_2

export function gsl_sf_laguerre_3( a, x )
{ // gsl_sf_laguerre_3
    return EVAL_RESULT_DD( gsl_sf_laguerre_3_e, { x: a, y: x }, "gsl_sf_laguerre_3" );
} // gsl_sf_laguerre_3

export function gsl_sf_laguerre_n( n, a, x )
{ // gsl_sf_laguerre_n
    return EVAL_RESULT_IDD( gsl_sf_laguerre_n_e, { i: n, x: a, y: x }, "gsl_sf_laguerre_n" );
} // gsl_sf_laguerre_n

// ----------------------------------------------------------------------------
// SF-Laguerre.mjs
