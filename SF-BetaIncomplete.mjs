// SF-BetaIncomplete.mjs
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
import { GSL_DBL_MIN }          from "./SF-Machine.mjs";
import { gsl_sf_log_1plusx_e }  from "./SF-Logarithmic.mjs";
import { gsl_sf_exp_err_e }     from "./SF-Exponential.mjs";
import { gsl_sf_beta_e }        from "./SF-Beta.mjs";
import { gsl_sf_lnbeta_e }      from "./SF-Beta.mjs";
import { gsl_sf_log_e }         from "./SF-Logarithmic.mjs";
import { gsl_sf_hyperg_2F1_e }  from "./SF-Hypergeometric2F1.mjs";

import { EVAL_RESULT_3D }       from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

function isnegint(x)
{
    return (x < 0.0) && (x == Math.floor(x));
} // isnegint

function beta_cont_frac(a, b, x)
{
    const max_iter   = 512; // control iterations
    var iter_count = 0;
    var k          = 0;
    const cutoff     = 2.0 * GSL_DBL_MIN; // control the zero cutoff
    var cf         = 0.0;
    var num_term   = 0.0;
    var den_term   = 0.0;
    var delta_frac = 0.0;
    var coeff      = 0.0;
    var r          = { val: 0.0, err: 0.0 }; // Result;

    // standard initialization for continued fraction
    num_term = 1.0;
    den_term = 1.0 - (a + b) * x / (a + 1.0);
    if (Math.abs(den_term) < cutoff)
    {
        den_term = cutoff;
    }
    den_term = 1.0 / den_term;
    cf = den_term;

    iter_count = 0;
    while (iter_count < max_iter)
    {
        k = iter_count + 1;
        coeff = (k) * (b - (k)) * x / (((a - 1.0) + (2 * k)) * (a + (2 * k)));

        // first step
        den_term = 1.0 + coeff * den_term;
        num_term = 1.0 + coeff / num_term;
        if (Math.abs(den_term) < cutoff)
        {
            den_term = cutoff;
        }
        if (Math.abs(num_term) < cutoff)
        {
            num_term = cutoff;
        }
        den_term = 1.0 / den_term;

        delta_frac = den_term * num_term;
        cf = cf * delta_frac;

        coeff = -(a + (k)) * (a + b + (k)) * x / ((a + (2 * k)) * (a + (2 * k + 1)));

        // second step
        den_term = 1.0 + coeff * den_term;
        num_term = 1.0 + coeff / num_term;
        if (Math.abs(den_term) < cutoff)
        {
            den_term = cutoff;
        }
        if (Math.abs(num_term) < cutoff)
        {
            num_term = cutoff;
        }
        den_term = 1.0 / den_term;

        delta_frac = den_term * num_term;
        cf = cf * delta_frac;

        if (Math.abs(delta_frac - 1.0) < 2.0 * GSL_DBL_EPSILON)
        {
            break;
        }

        iter_count = iter_count + 1;
    }

    r.val = cf;
    r.err = (iter_count) * 4.0 * GSL_DBL_EPSILON * Math.abs(cf);

    if (iter_count >= max_iter)
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // beta_cont_frac

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_beta_inc_e(a, b, x)
{
    var prefactor = 0.0;
    var f         = { val: 0.0, err: 0.0 }; // Result;
    var beta      = { val: 0.0, err: 0.0 }; // Result;
    var r         = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0 || x > 1.0)
    {
        throw "SF.DomainException";
    }
    else if (isnegint(a) || isnegint(b))
    {
        throw "SF.DomainException";
    }
    else if (isnegint(a + b))
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else if (x == 1.0)
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if (a <= 0.0 || b <= 0.0)
    {
        f = gsl_sf_hyperg_2F1_e(a, 1.0 - b, a + 1.0, x);
        beta = gsl_sf_beta_e(a, b);
        prefactor = (x ** a) / a;
        r.val = prefactor * f.val / beta.val;
        r.err = Math.abs(prefactor) * f.err / Math.abs(beta.val) + Math.abs(r.val / beta.val) * beta.err;
        //CHECK_UNDERFLOW(result);
        return r;
    }
    else
    {
        var ln_beta    = { val: 0.0, err: 0.0 }; // Result;
        var ln_x       = { val: 0.0, err: 0.0 }; // Result;
        var ln_1mx     = { val: 0.0, err: 0.0 }; // Result;
        //var prefactor  = { val: 0.0, err: 0.0 }; // Result;
        var cf         = { val: 0.0, err: 0.0 }; // Result;
        var ln_pre_val = 0.0;
        var ln_pre_err = 0.0;
        var term       = 0.0;

        try
        {
            ln_beta = gsl_sf_lnbeta_e(a, b);
            ln_1mx  = gsl_sf_log_1plusx_e(-x);
            ln_x    = gsl_sf_log_e(x);
        }
        catch (e)
        {
            throw "SF.SanityException";
        }

        ln_pre_val = -ln_beta.val + a * ln_x.val + b * ln_1mx.val;
        ln_pre_err =  ln_beta.err + Math.abs(a * ln_x.err) + Math.abs(b * ln_1mx.err);
        prefactor  = gsl_sf_exp_err_e(ln_pre_val, ln_pre_err);

        if (x < (a + 1.0) / (a + b + 2.0))
        {
            // Apply continued fraction directly.
            cf = beta_cont_frac(a, b, x);
            r.val = prefactor.val * cf.val / a;
            r.err = (Math.abs(prefactor.err * cf.val) + Math.abs(prefactor.val * cf.err)) / a;
            //CHECK_UNDERFLOW(result);
            return r;
        }
        else
        {
            // Apply continued fraction after hypergeometric transformation.
            cf = beta_cont_frac(b, a, 1.0 - x);
            term = prefactor.val * cf.val / b;
            r.val = 1.0 - term;
            r.err = Math.abs(prefactor.err * cf.val) / b;
            r.err = r.err + Math.abs(prefactor.val * cf.err) / b;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * (1.0 + Math.abs(term));
            //CHECK_UNDERFLOW(result);
            return r;
        }
    }

} // gsl_sf_beta_inc_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_beta_inc( a, b, x )
{ // gsl_sf_beta_inc
    return EVAL_RESULT_3D( gsl_sf_beta_inc_e, { x: a, y: b, z: x }, "gsl_sf_beta_inc" );
} // gsl_sf_beta_inc

// ----------------------------------------------------------------------------
// EOF SF-BetaIncomplete.mjs
