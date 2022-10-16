// SF-Exponential.mjs
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

import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON }   from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MIN }       from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }       from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }      from "./SF-Machine.mjs";
import { GSL_SIGN }              from "./SF-Math.mjs";
import { GSL_IS_ODD }            from "./SF-Math.mjs";
import { M_LN10 }                from "./SF-Math.mjs";
import { gsl_sf_lnfact_e }       from "./SF-Gamma.mjs";

import { EVAL_RESULT_D }  from "./SF-Evaluate.mjs";
import { EVAL_RESULT_DD } from "./SF-Evaluate.mjs";
import { EVAL_RESULT_ID } from "./SF-Evaluate.mjs";

const INT32_MAX = 2147483647;
const INT32_MIN = -2147483647 - 1;

// Evaluate the continued fraction for exprel.
// [Abramowitz+Stegun, 4.2.41]
function exprel_n_CF(n, x)
{
    const RECUR_BIG = GSL_SQRT_DBL_MAX;
    const maxiter   = 5000;

    var nn   = 0;
    var Anm2 = 0.0;
    var Bnm2 = 0.0;
    var Anm1 = 0.0;
    var Bnm1 = 0.0;
    var a1   = 0.0;
    var b1   = 0.0;
    var a2   = 0.0;
    var b2   = 0.0;
    var ann  = 0.0;
    var bnn  = 0.0;
    var fn   = 0.0;
    var An   = 0.0;
    var Bn   = 0.0;
    var old_fn = 0.0;
    var del   = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    nn = 1;
    Anm2 = 1.0;
    Bnm2 = 0.0;
    Anm1 = 0.0;
    Bnm1 = 1.0;
    a1 = 1.0;
    b1 = 1.0;
    a2 = -x;
    b2 = n + 1.0;
  
    An = b1 * Anm1 + a1 * Anm2;   // A1
    Bn = b1 * Bnm1 + a1 * Bnm2;   // B1
    
    // One explicit step, before we get to the main pattern.
    nn = nn + 1;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    An = b2 * Anm1 + a2 * Anm2;   // A2
    Bn = b2 * Bnm1 + a2 * Bnm2;   // B2
  
    fn = An / Bn;
  
    while (nn < maxiter)
    {
        nn = nn + 1;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        if (GSL_IS_ODD(nn))
        {
            ann = Math.trunc((nn - 1) / 2) * x;
        }
        else
        {
            ann = -(n + Math.trunc(nn / 2) - 1.0) * x;
        }
        bnn = n + (nn - 1);
        An = bnn * Anm1 + ann * Anm2;
        Bn = bnn * Bnm1 + ann * Bnm2;
  
        if ((Math.abs(An) > RECUR_BIG) || (Math.abs(Bn) > RECUR_BIG))
        {
            An = An / RECUR_BIG;
            Bn = Bn / RECUR_BIG;
            Anm1 = Anm1 / RECUR_BIG;
            Bnm1 = Bnm1 / RECUR_BIG;
            Anm2 = Anm2 / RECUR_BIG;
            Bnm2 = Bnm2 / RECUR_BIG;
        }
  
        old_fn = fn;
        fn = An / Bn;
        del = old_fn / fn;
          
        if (Math.abs(del - 1.0) < 2.0 * GSL_DBL_EPSILON)
        {
            break;
        }
    }
  
    r.val = fn;
    r.err = 4.0 * ((nn) + 1.0) * GSL_DBL_EPSILON * Math.abs(fn);

    if (nn >= maxiter)
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // exprel_n_CF

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-**)

export function gsl_sf_exp_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x > GSL_LOG_DBL_MAX)
    {
        throw "SF.OverflowException";
    }
    else if (x < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else
    {
        r.val = Math.exp(x);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }

} // gsl_sf_exp_e

// ----------------------------------------------------------------------------

export function gsl_sf_exp_e10_e(x)
{
    var N = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x > (INT32_MAX - 1))
    {
        throw "SF.OverflowException";
    }
    else if (x < (INT32_MIN + 1))
    {
        throw "SF.UnderflowException";
    }
    else
    {
        if ((x > GSL_LOG_DBL_MAX) || (x < GSL_LOG_DBL_MIN))
        {
            N = Math.trunc(Math.floor(x / M_LN10));
        }
        else
        {
            N = 0;
        }
        r.val = Math.exp(x - (N) * M_LN10);
        r.err = 2.0 * (Math.abs(x) + 1.0) * GSL_DBL_EPSILON * Math.abs(r.val);
        r.e10 = N;
    }

    return r;

} // gsl_sf_exp_e10_e

// ----------------------------------------------------------------------------

export function gsl_sf_exp_mult_e(x, y)
{
    var ay   = 0.0;
    var ex   = 0.0;
    var ly   = 0.0;
    var lnr  = 0.0;
    var sy   = 0.0;
    var M    = 0.0;
    var N    = 0.0;
    var a    = 0.0;
    var b    = 0.0;
    var berr = 0.0;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    ay = Math.abs(y);
  
    if (y == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if ((x < 0.5 * GSL_LOG_DBL_MAX && x > 0.5 * GSL_LOG_DBL_MIN)
        && (ay < 0.8 * GSL_SQRT_DBL_MAX && ay > 1.2 * GSL_SQRT_DBL_MIN))
    {
        ex = Math.exp(x);
        r.val = y * ex;
        r.err = (2.0 + Math.abs(x)) * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        ly  = Math.log(ay);
        lnr = x + ly;
  
        if (lnr > GSL_LOG_DBL_MAX - 0.01)
        {
            throw "SF.OverflowException";
        }
        else if (lnr < GSL_LOG_DBL_MIN + 0.01)
        {
            throw "SF.UnderflowException";
        }
        else
        {
            sy   = GSL_SIGN(y);
            M    = Math.floor(x);
            N    = Math.floor(ly);
            a    = x  - M;
            b    = ly - N;
            berr = 2.0 * GSL_DBL_EPSILON * (Math.abs(ly) + Math.abs(N));
            r.val = sy * Math.exp(M + N) * Math.exp(a + b);
            r.err = berr * Math.abs(r.val);
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * (M + N + 1.0) * Math.abs(r.val);
        }
    }

    return r;

} // gsl_sf_exp_mult_e

// ----------------------------------------------------------------------------

export function gsl_sf_exp_mult_e10_e(x, y)
{
    var ay = Math.abs(y);
    var ex = 0.0;
    var ly = 0.0;
    var sy = 0.0;
    var l10_val = 0.0;
    var arg_val = 0.0;
    var arg_err = 0.0;
    var n = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;
  
    if (y == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
        r.e10 = 0;
    }
    else if (((x < 0.5 * GSL_LOG_DBL_MAX) && (x > 0.5 * GSL_LOG_DBL_MIN))
        && ((ay < 0.8 * GSL_SQRT_DBL_MAX) && (ay > 1.2 * GSL_SQRT_DBL_MIN)))
    {
        ex = Math.exp(x);
        r.val = y * ex;
        r.err = (2.0 + Math.abs(x)) * GSL_DBL_EPSILON * Math.abs(r.val);
        r.e10 = 0;
    }
    else
    {
        ly = Math.log(ay);
        l10_val = (x + ly) / M_LN10;
  
        if (l10_val > (INT32_MAX - 1))
        {
            throw "SF.OverflowException";
        }
        else if (l10_val < (INT32_MIN + 1))
        {
            throw "SF.UnderflowException";
        }
        else
        {
            sy = GSL_SIGN(y);
            n = Math.trunc(Math.floor(l10_val));
            arg_val = (l10_val - (n)) * M_LN10;
            arg_err = 2.0 * GSL_DBL_EPSILON * (Math.abs(x) + Math.abs(ly) + M_LN10 * Math.abs((n)));
  
            r.val = sy * Math.exp(arg_val);
            r.err = arg_err * Math.abs(r.val);
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            r.e10 = n;
        }
    }

    return r;

} // gsl_sf_exp_mult_e10_e

// ----------------------------------------------------------------------------

export function gsl_sf_exp_mult_err_e(x, dx, y, dy)
{
    var ay  = 0.0;
    var ex  = 0.0;
    var ly  = 0.0;
    var lnr = 0.0;
    var sy  = 0.0;
    var M   = 0.0;
    var N   = 0.0;
    var a   = 0.0;
    var b   = 0.0;
    var eMN = 0.0;
    var eab = 0.0;
    var r   = { val: 0.0, err: 0.0 }; // Result;


    ay = Math.abs(y);

    if (y == 0.0)
    {
        r.val = 0.0;
        r.err = Math.abs(dy * Math.exp(x));
    }
    else if (((x < 0.5 * GSL_LOG_DBL_MAX) && (x > 0.5 * GSL_LOG_DBL_MIN))
        && ((ay < 0.8 * GSL_SQRT_DBL_MAX) && (ay > 1.2 * GSL_SQRT_DBL_MIN)))
    {
       
        ex = Math.exp(x);
        r.val = y * ex;
        r.err = ex * (Math.abs(dy) + Math.abs(y * dx));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        ly  = Math.log(ay);
        lnr = x + ly;
  
        if (lnr > GSL_LOG_DBL_MAX - 0.01)
        {
            throw "SF.OverflowException";
        }
        else if (lnr < GSL_LOG_DBL_MIN + 0.01)
        {
            throw "SF.UnderflowException";
        }
        else
        {
            sy  = GSL_SIGN(y);
            M   = Math.floor(x);
            N   = Math.floor(ly);
            a   = x  - M;
            b   = ly - N;
            eMN = Math.exp(M + N);
            eab = Math.exp(a + b);
            r.val = sy * eMN * eab;
            r.err = eMN * eab * 2.0*GSL_DBL_EPSILON;
            r.err = r.err + eMN * eab * Math.abs(dy / y);
            r.err = r.err + eMN * eab * Math.abs(dx);
        }
    }

    return r;

} // gsl_sf_exp_mult_err_e

// ----------------------------------------------------------------------------

export function gsl_sf_exp_mult_err_e10_e(x, dx, y, dy)
{
    var ay = Math.abs(y);
    var ex = 0.0;
    var ly = 0.0;
    var sy = 0.0;
    var l10_val = 0.0;
    var arg_val = 0.0;
    var arg_err = 0.0;
    var n = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;
  
    if (y == 0.0)
    {
        r.val = 0.0;
        r.err = Math.abs(dy * Math.exp(x));
        r.e10 = 0;
    }
    else if (((x < 0.5 * GSL_LOG_DBL_MAX) && (x > 0.5 * GSL_LOG_DBL_MIN))
        && ((ay < 0.8 * GSL_SQRT_DBL_MAX) && (ay > 1.2 * GSL_SQRT_DBL_MIN)))
    {   
        ex = Math.exp(x);
        r.val = y * ex;
        r.err = ex * (Math.abs(dy) + Math.abs(y * dx));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        r.e10 = 0;
    }
    else
    {
        ly = Math.log(ay);
        l10_val = (x + ly) / M_LN10;
  
        if (l10_val > (INT32_MAX - 1))
        {
            throw "SF.OverflowException";
        }
        else if (l10_val < (INT32_MIN + 1))
        {
            throw "SF.UnderflowException";
        }
        else
        {
            sy  = GSL_SIGN(y);
            n   = Math.trunc(Math.floor(l10_val));
            arg_val = (l10_val - (n)) * M_LN10;
            arg_err = dy / Math.abs(y) + dx + 2.0 * GSL_DBL_EPSILON * Math.abs(arg_val);
  
            r.val = sy * Math.exp(arg_val);
            r.err = arg_err * Math.abs(r.val);
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            r.e10 = n;
        }
    }

    return r;

} // gsl_sf_exp_mult_err_e10_e

// ----------------------------------------------------------------------------

export function gsl_sf_expm1_e(x)
{
    const cut = 0.02;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        r.val = -1.0;
        r.err = GSL_DBL_EPSILON;
    }
    else if (x < -cut)
    {
        r.val = Math.exp(x) - 1.0;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < cut)
    {
        r.val = x * (1.0 + 0.5 * x * (1.0 + x / 3.0 * (1.0 + 0.25 * x * (1.0 + 0.2 * x))));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < GSL_LOG_DBL_MAX)
    {
        r.val = Math.exp(x) - 1.0;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

    return r;

} // gsl_sf_expm1_e

// ----------------------------------------------------------------------------

export function gsl_sf_exprel_e(x)
{
    const cut = 0.002;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        r.val = -1.0 / x;
        r.err = GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < -cut)
    {
        r.val = (Math.exp(x) - 1.0) / x;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < cut)
    {
        r.val = (1.0 + 0.5 * x * (1.0 + x / 3.0 * (1.0 + 0.25 * x * (1.0 + 0.2 * x))));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < GSL_LOG_DBL_MAX)
    {
        r.val = (Math.exp(x) - 1.0) / x;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

    return r;

} // gsl_sf_exprel_e

// ----------------------------------------------------------------------------

export function gsl_sf_exprel_2_e(x)
{
    const cut = 0.002;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        r.val = -2.0 / x * (1.0 + 1.0 / x);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < -cut)
    {
        r.val = 2.0 * (Math.exp(x) - 1.0 - x) / (x * x);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < cut)
    {
        r.val = (1.0 + 1.0 / 3.0 * x * (1.0 + 0.25 * x * (1.0 + 0.2 * x * (1.0 + 1.0 / 6.0 * x))));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < GSL_LOG_DBL_MAX)
    {
        r.val = 2.0 * (Math.exp(x) - 1.0 - x) / (x * x);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

    return r;

} // gsl_sf_exprel_2_e

// ----------------------------------------------------------------------------

export function gsl_sf_exprel_n_CF_e(n, x)
{
    return exprel_n_CF(n, x);
} // gsl_sf_exprel_n_CF_e

// ----------------------------------------------------------------------------

export function gsl_sf_exprel_n_e(n, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (n < 0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        r.val = 1.0;
        r.err = 0.0;
    }
    else if (Math.abs(x) < GSL_ROOT3_DBL_EPSILON * (n))
    {
        r.val = 1.0 + x / (n + 1) * (1.0 + x / (n + 2));
        r.err = 2.0 * GSL_DBL_EPSILON;
    }
    else if (n == 0)
    {
        r = gsl_sf_exp_e(x);
    }
    else if (n == 1)
    {
        r = gsl_sf_exprel_e(x);
    }
    else if (n == 2)
    {
        r = gsl_sf_exprel_2_e(x);
    }
    else
    {
        if ((x > (n)) && ((-x + (n) * (1.0 + Math.log(x / (n))) < GSL_LOG_DBL_EPSILON)))
        {//console.log(1);
            // x is much larger than n.
            // Ignore polynomial part, so
            // exprel_N(x) ~= e^x N!/x^N
            var lnf_N   = { val: 0.0, err: 0.0 }; // Result;
            var lnr_val = 0.0;
            var lnr_err = 0.0;
            var lnterm  = 0.0;

            lnf_N   = gsl_sf_lnfact_e(n);
            lnterm  = (n) * Math.log(x);
            lnr_val = x + lnf_N.val - lnterm;
            lnr_err = GSL_DBL_EPSILON * (Math.abs(x) + Math.abs(lnf_N.val) + Math.abs(lnterm));
            lnr_err = lnr_err + lnf_N.err;
            r = gsl_sf_exp_err_e(lnr_val, lnr_err);
        }
        else if (x > (n))
        {//console.log(2);
            // Write the identity
            //   exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
            // then use the asymptotic expansion
            // Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)
            let ln_x      = Math.log(x);
            let lg_N      = 0.0;
            let lnpre_val = 0.0;
            let lnpre_err = 0.0;
            let lnf_N     = { val: 0.0, err: 0.0 }; // Result;

            lnf_N     = gsl_sf_lnfact_e(n); // log(N!)
            lg_N      = lnf_N.val - Math.log((n)); // log(Gamma(N))
            lnpre_val = x + lnf_N.val - (n) * ln_x;
            lnpre_err = GSL_DBL_EPSILON * (Math.abs(x) + Math.abs(lnf_N.val) + Math.abs((n) * ln_x));
            lnpre_err = lnpre_err + lnf_N.err;
            if (lnpre_val < GSL_LOG_DBL_MAX - 5.0)
            {
                let bigG_ratio = { val: 0.0, err: 0.0 }; // Result;
                let pre        = { val: 0.0, err: 0.0 }; // Result;
                let ln_bigG_ratio_pre = 0.0;
                let bigGsum    = 0.0;
                let term       = 0.0;
                let k = 0;

                pre = gsl_sf_exp_err_e(lnpre_val, lnpre_err);
                ln_bigG_ratio_pre = -x + (n - 1) * ln_x - lg_N;
                bigGsum = 1.0;
                term = 1.0;
                for (k = 1; k <= n - 1; k++)
                {
                    term = term * (n - k) / x;
                    bigGsum = bigGsum + term;
                }
                bigG_ratio = gsl_sf_exp_mult_e(ln_bigG_ratio_pre, bigGsum);
                r.val = pre.val * (1.0 - bigG_ratio.val);
                r.err = pre.val * (2.0 * GSL_DBL_EPSILON + bigG_ratio.err);
                r.err = r.err + pre.err * Math.abs(1.0 - bigG_ratio.val);
                r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            }
            else
            {
                throw "SF.OverflowException";
            }
        }
        else if (x > -10.0 * (n))
        {//console.log(3);
            r = exprel_n_CF((n), x);
        }
        else
        {//console.log(4);
            // x -> -Inf asymptotic:
            // exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
            //             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
            var sum  = 1.0;
            var term = 1.0;
            var k    = 0;

            for (k = 1; k <= n - 1; k++)
            {
                term = term * (n - k) / x;
                sum  = sum  + term;
            }
            r.val = -(n) / x * sum;
            r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
    }

    return r;

} // gsl_sf_exprel_n_e


// ----------------------------------------------------------------------------

export function gsl_sf_exp_err_e(x, dx)
{
    var adx = 0.0;
    var edx = 0.0;
    var ex  = 0.0;
    var r   = { val: 0.0, err: 0.0 }; // Result;

    adx = Math.abs(dx);
    if (x + adx > GSL_LOG_DBL_MAX)
    {
        throw "SF.OverflowException";
    }
    else if (x - adx < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else
    {
        ex  = Math.exp(x);
        edx = Math.exp(adx);
        r.val = ex;
        r.err = ex * Math.max(GSL_DBL_EPSILON, edx - 1.0 / edx);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_exp_err_e

// ----------------------------------------------------------------------------

export function gsl_sf_exp_err_e10_e(x, dx)
{
    var adx = Math.abs(dx);
    var ex  = 0.0;
    var n   = 0;
    var r   = { val: 0.0, err: 0.0 }; // Result;
  
    if (x + adx > (INT32_MAX - 1))
    {
        throw "SF.OverflowException";
    }
    else if (x - adx < (INT32_MIN + 1))
    {
        throw "SF.UnderflowException";
    }
    else
    {
        n = Math.trunc(Math.floor(x / M_LN10));
        ex = Math.exp(x - (n) * M_LN10);
        r.val = ex;
        r.err = ex * (2.0 * GSL_DBL_EPSILON * (Math.abs(x) + 1.0) + adx);
        r.e10 = n;
    }

    return r;

} // gsl_sf_exp_err_e10_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_exp( x )
{ // gsl_sf_exp
    return EVAL_RESULT_D( gsl_sf_exp_e, x, "gsl_sf_exp_e" );
} // gsl_sf_exp;

export function gsl_sf_exp_mult( x, y )
{ // gsl_sf_expm1
    return EVAL_RESULT_DD( gsl_sf_exp_mult_e, { x: x, y: y }, "gsl_sf_exp_mult_e" );
} // gsl_sf_exp_mult;

export function gsl_sf_expm1( x )
{ // gsl_sf_expm1
    return EVAL_RESULT_D( gsl_sf_expm1_e, x, "gsl_sf_expm1_e" );
} // gsl_sf_expm1;

export function gsl_sf_exprel( x )
{ // gsl_sf_expm1
    return EVAL_RESULT_D(gsl_sf_exprel_e, x, "gsl_sf_exprel_e" );
} // gsl_sf_exprel;

export function gsl_sf_exprel_2( x )
{ // gsl_sf_expm1
    return EVAL_RESULT_D( gsl_sf_exprel_2_e, x, "gsl_sf_exprel_2_e" );
} // gsl_sf_exprel_2;

export function gsl_sf_exprel_n( n, x )
{ // gsl_sf_expm1
    return EVAL_RESULT_ID( gsl_sf_exprel_n_e, { i: n, x: x }, "gsl_sf_exprel_n_e" );
} // gsl_sf_exprel_n;

// ----------------------------------------------------------------------------
// EOF SF-Exponential.mjs
