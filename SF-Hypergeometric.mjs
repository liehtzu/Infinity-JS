// SF-Hypergeometric.mjs
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

import { M_PI }                   from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }        from "./SF-Machine.mjs";
import { GSL_DBL_MAX }            from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }        from "./SF-Machine.mjs";
import { gsl_sf_lngamma_sgn_e }   from "./SF-Gamma.mjs";
import { gsl_sf_exp_mult_err_e }  from "./SF-Exponential.mjs";
import { gsl_sf_bessel_Inu_scaled_e } from "./SF-BesselInu.mjs";
import { gsl_sf_bessel_Knu_scaled_e } from "./SF-BesselKnu.mjs";
import { gsl_sf_bessel_Jnu_e }   from "./SF-BesselJnu.mjs";
import { gsl_sf_bessel_Ynu_e }   from "./SF-BesselYnu.mjs";

import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const SUM_LARGE = 1.0e-5 * GSL_DBL_MAX;


export function gsl_sf_hyperg_1F1_series_e( a, b, x )
{
    var an          = a;
    var bn          = b;
    var n           = 1.0;
    var del         = 1.0;
    var abs_del     = 1.0;
    var max_abs_del = 1.0;
    var sum_val     = 1.0;
    var sum_err     = 0.0;
    var u           = 0.0;
    var abs_u       = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;
  
    while ( abs_del / Math.abs( sum_val ) > 0.25 * GSL_DBL_EPSILON )
    {
        if ( bn == 0.0 )
        {
            throw "SF.DomainException";
        }
       
        if ( an == 0.0 )
        {
            r.val = sum_val;
            r.err = sum_err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * n * Math.abs( sum_val );
            return r;
        }
       
        if ( n > 10000.0 )
        {
            r.val = sum_val;
            r.err = sum_err;
            throw "SF.GenericException"; // WITH "hypergeometric series failed to converge"; -- GSL_EFAILED
        }
       
        u = x * (an / (bn * n));
        abs_u = Math.abs( u );
        if ( abs_u > 1.0 && max_abs_del > GSL_DBL_MAX / abs_u )
        {
            r.val = sum_val;
            r.err = Math.abs( sum_val );
            throw "SF.OverflowException";
        }
        del = del * u;
        sum_val = sum_val + del;
        if ( Math.abs( sum_val ) > SUM_LARGE )
        {
            r.val = sum_val;
            r.err = Math.abs( sum_val );
            throw "SF.OverflowException";
        }
       
        abs_del = Math.abs( del );
        max_abs_del = Math.max( abs_del, max_abs_del );
        sum_err = sum_err + 2.0 * GSL_DBL_EPSILON * abs_del;
       
        an = an + 1.0;
        bn = bn + 1.0;
        n  = n  + 1.0;
    }
  
    r.val = sum_val;
    r.err = sum_err;
    r.err = r.err + abs_del;
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * n * Math.abs( sum_val );
  
    return r;

} // gsl_sf_hyperg_1F1_series_e

// ----------------------------------------------------------------------------

function gsl_sf_hyperg_1F1_large_b_e( a, b, x )
{

    if ( Math.abs( x / b ) < 1.0 )
    {
        var u   = x / b;
        var v   = 1.0 / (1.0 - u);
        var pre = (v ** a);
        var uv  = u * v;
        var uv2 = uv * uv;
        var t1  = a * (a + 1.0) / (2.0 * b) * uv2;
        var t2a = a * (a + 1.0) / (24.0 * b * b) * uv2;
        var t2b = 12.0 + 16.0 * (a + 2.0) * uv + 3.0 * (a + 2.0) * (a + 3.0) * uv2;
        var t2  = t2a * t2b;
        var r   = { val: 0.0, err: 0.0 }; // Result;

        r.val = pre * (1.0 - t1 + t2);
        r.err = pre * GSL_DBL_EPSILON * (1.0 + Math.abs( t1 ) + Math.abs( t2 ));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else
    {
        throw "SF.DomainException";
    }

} // gsl_sf_hyperg_1F1_large_b_e

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_U_large_b_e( a, b, x, r, ln_multiplier) //: IN OUT LONG_FLOAT )
{
    var N   = Math.floor( b ); // b = N + eps
    var eps = b - N;

    if ( Math.abs(eps) < GSL_SQRT_DBL_EPSILON )
    {
        var lnpre_val = 0.0;
        var lnpre_err = 0.0;
        var M = { val: 0.0, err: 0.0 }; // Result;

        if ( b > 1.0 )
        {
            var tmp    = (1.0 - b) * Math.log( x );
            var lg_bm1 = { val: 0.0, err: 0.0 }; // Result;
            var lg_a   = { val: 0.0, err: 0.0 }; // Result;

            lg_bm1 = gsl_sf_lngamma_e( b - 1.0 );
            lg_a   = gsl_sf_lngamma_e( a );
            lnpre_val = tmp + x + lg_bm1.val - lg_a.val;
            lnpre_err = lg_bm1.err + lg_a.err + GSL_DBL_EPSILON * (Math.abs( x ) + Math.abs( tmp ));
            M = gsl_sf_hyperg_1F1_large_b_e( 1.0 - a, 2.0 - b, -x );
        }
        else
        {
            var lg_1mb   = { val: 0.0, err: 0.0 }; // Result;
            var lg_1pamb = { val: 0.0, err: 0.0 }; // Result;

            lg_1mb    = gsl_sf_lngamma_e( 1.0 - b );
            lg_1pamb  = gsl_sf_lngamma_e( 1.0 + a - b );
            lnpre_val = lg_1mb.val - lg_1pamb.val;
            lnpre_err = lg_1mb.err + lg_1pamb.err;
            M = gsl_sf_hyperg_1F1_large_b_e( a, b, x );
        }
        
        if ( lnpre_val > GSL_LOG_DBL_MAX - 10.0 )
        {
            result.val  = M.val;
            result.err  = M.err;
            ln_multiplier.Double = lnpre_val;
            throw "SF.OverflowException";
        }
        else
        {
            var epre = { val: 0.0, err: 0.0 }; // Result;

            epre = gsl_sf_exp_err_e( lnpre_val, lnpre_err );
            r.val = epre.val * M.val;
            r.err = epre.val * M.err + epre.err * Math.abs( M.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            ln_multiplier.Double = 0.0;
            return;

        }
    }
    else
    {
        var omb_lnx    = (1.0 - b) * Math.log( x );
        var sgn_1mb    = 0.0;
        var sgn_1pamb  = 0.0;
        var sgn_bm1    = 0.0;
        var sgn_a      = 0.0;
        var lnpre1_val = 0.0;
        var lnpre2_val = 0.0;
        var lnpre1_err = 0.0;
        var lnpre2_err = 0.0;
        var sgpre1     = 0.0;
        var sgpre2     = 0.0;
        var lg_1mb     = { val: 0.0, err: 0.0 }; // Result;
        var lg_1pamb   = { val: 0.0, err: 0.0 }; // Result;
        var lg_bm1     = { val: 0.0, err: 0.0 }; // Result;
        var lg_a       = { val: 0.0, err: 0.0 }; // Result;
        var M1         = { val: 0.0, err: 0.0 }; // Result;
        var M2         = { val: 0.0, err: 0.0 }; // Result;

        M1 = gsl_sf_hyperg_1F1_large_b_e( a, b, x );
        M2 = gsl_sf_hyperg_1F1_large_b_e( 1.0 - a, 2.0 - b, x );
        
        gsl_sf_lngamma_sgn_e( 1.0 - b,     lg_1mb,   sgn_1mb );
        gsl_sf_lngamma_sgn_e( 1.0 + a - b, lg_1pamb, sgn_1pamb );
        
        gsl_sf_lngamma_sgn_e( b - 1.0, lg_bm1, sgn_bm1 );
        gsl_sf_lngamma_sgn_e( a,       lg_a,   sgn_a );

        lnpre1_val = lg_1mb.val - lg_1pamb.val;
        lnpre1_err = lg_1mb.err + lg_1pamb.err;
        lnpre2_val = lg_bm1.val - lg_a.val - omb_lnx - x;
        lnpre2_err = lg_bm1.err + lg_a.err + GSL_DBL_EPSILON * (Math.abs( omb_lnx ) + Math.abs( x ));
        sgpre1 = sgn_1mb * sgn_1pamb;
        sgpre2 = sgn_bm1 * sgn_a;
        
        if ( lnpre1_val > GSL_LOG_DBL_MAX - 10.0 || lnpre2_val > GSL_LOG_DBL_MAX - 10.0 )
        {
            //max_lnpre_val : LONG_FLOAT = LONG_FLOAT'Max(lnpre1_val, lnpre2_val);
            //max_lnpre_err : LONG_FLOAT = LONG_FLOAT'Max(lnpre1_err, lnpre2_err);
            //lp1 : LONG_FLOAT = lnpre1_val - max_lnpre_val;
            //lp2 : LONG_FLOAT = lnpre2_val - max_lnpre_val;
            //t1  : LONG_FLOAT = sgpre1 * Exp(lp1);
            //t2  : LONG_FLOAT = sgpre2 * Exp(lp2);

            //r.val = t1 * M1.val + t2 * M2.val;
            //r.err = Math.abs(t1) * M1.err + Math.abs(t2) * M2.err;
            //r.err = r.err + GSL_DBL_EPSILON * Exp(max_lnpre_err) * (Math.abs(t1 * M1.val) + Math.abs(t2 * M2.val));
            //r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            //ln_multiplier = max_lnpre_val;
            throw "SF.OverflowException";
        }
        else
        {
            var t1 = sgpre1 * Math.exp( lnpre1_val );
            var t2 = sgpre2 * Math.exp( lnpre2_val );

            r.val = t1 * M1.val + t2 * M2.val;
            r.err = Math.abs( t1 ) * M1.err + Math.abs( t2 ) * M2.err;
            r.err = r.err + GSL_DBL_EPSILON * (Math.exp( lnpre1_err ) * Math.abs( t1 * M1.val ) + Math.exp( lnpre2_err ) * Math.abs( t2 * M2.val ));
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            ln_multiplier.Double = 0.0;
            return;
        }
    }

} // gsl_sf_hyperg_U_large_b_e

// ----------------------------------------------------------------------------

// [Carlson, p.109] says the error in truncating this asymptotic series
// is less than the absolute value of the first neglected term.
//
// A termination argument is provided, so that the series will
// be summed at most up to n=n_trunc. If n_trunc is set negative,
// then the series is summed until it appears to start diverging.
//
export function gsl_sf_hyperg_2F0_series_e( a, b, x, n_trunc )
{
    const maxiter = 2000;
    var an           = a;
    var bn           = b;  
    var n            = 1.0;
    var sum          = 1.0;
    var del          = 1.0;
    var abs_del      = 1.0;
    var max_abs_del  = 1.0;
    var last_abs_del = 1.0;
    var u            = 0.0;
    var abs_u        = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;
    
    while ( abs_del / Math.abs( sum ) > GSL_DBL_EPSILON && Math.trunc( n ) < maxiter )
    {
  
        u = an * (bn/n * x);
        abs_u = Math.abs( u );
       
        if ( abs_u > 1.0 && (max_abs_del > GSL_DBL_MAX / abs_u) )
        {
            r.val = sum;
            r.err = Math.abs( sum );
            throw "SF.OverflowException";
        }
       
        del = del * u;
        sum = sum + del;
       
        abs_del = Math.abs( del );
       
        if ( abs_del > last_abs_del ) break; // series is probably starting to grow
       
        last_abs_del = abs_del;
        max_abs_del  = Math.max( abs_del, max_abs_del );
       
        an = an + 1.0;
        bn = bn + 1.0;
        n  = n  + 1.0;
        
        if ( an == 0.0 || bn == 0.0 ) break; // series terminated
        
        if ( n_trunc >= 0 && Math.trunc( n ) >= n_trunc ) break; // reached requested timeout
    }
  
    r.val = sum;
    r.err = GSL_DBL_EPSILON * n + abs_del;
    if ( Math.trunc( n ) >= maxiter )
    {
        throw "SF.MaxIterationsException";
    }
    else
    {
        return r;
    }

} // gsl_sf_hyperg_2F0_series_e

// ----------------------------------------------------------------------------
// EOF SF-Hypergeometric.mjs
