// SF-CoulombBound.mjs
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
// Translation to Ada: Igor Izvarin

import { GSL_DBL_EPSILON }     from "./SF-Machine.mjs";
import { gsl_sf_lnfact_e }     from "./SF-Gamma.mjs";
import { gsl_sf_exp_err_e }    from "./SF-Exponential.mjs";
import { gsl_sf_laguerre_n_e } from "./SF-Laguerre.mjs";

import { EVAL_RESULT_DD }      from "./SF-Evaluate.mjs";
import { EVAL_RESULT_IIDD }    from "./SF-Evaluate.mjs";

// normalization for hydrogenic wave functions
function R_norm( n, l, Z )
{
    var A        = 0.0;
    var pre      = 0.0;
    var diff_val = 0.0;
    var diff_err = 0.0;
    var ln_a     = { val: 0.0, err: 0.0 }; // Result;
    var ln_b     = { val: 0.0, err: 0.0 }; // Result;
    var ex       = { val: 0.0, err: 0.0 }; // Result;
    var r        = { val: 0.0, err: 0.0 }; // Result;

    A   = 2.0 * Z / (n);
    pre = Math.sqrt( A * A * A / (2.0 * (n)) );
    ln_a = gsl_sf_lnfact_e( n + l );
    ln_b = gsl_sf_lnfact_e( n - l - 1 );
    diff_val = 0.5 * (ln_b.val - ln_a.val);
    diff_err = 0.5 * (ln_b.err + ln_a.err) + GSL_DBL_EPSILON * Math.abs( diff_val );
    ex = gsl_sf_exp_err_e( diff_val, diff_err );
    r.val = pre * ex.val;
    r.err = pre * ex.err;
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
    return r;

} // R_norm

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hydrogenicR_1_e( Z, s )
{
    var A    = 0.0;
    var norm = 0.0;
    var ea   = 0.0;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    if ( Z > 0.0 && s >= 0.0 )
    {
        A     = 2.0 * Z;
        norm  = A * Math.sqrt( Z );
        ea    = Math.exp( -Z * s );
        r.val = norm * ea;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val ) * Math.abs( Z * s );
        //CHECK_UNDERFLOW(result);
        return r;
    }
    else
    {
        throw "SF.DomainException";
    }

} // gsl_sf_hydrogenicR_1_e

// ----------------------------------------------------------------------------

export function gsl_sf_hydrogenicR_e( n, l, Z, s )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( n < 1 || l > n - 1 || Z <= 0.0 || s < 0.0 )
    {
        throw "SF.DomainException";
    }
    else
    {
        var A     = 0.0;
        var rho   = 0.0;
        var ea    = 0.0;
        var pp    = 0.0;
        var W_val = 0.0;
        var W_err = 0.0;
        var norm  = { val: 0.0, err: 0.0 }; // Result;
        var lag   = { val: 0.0, err: 0.0 }; // Result;

        A = 2.0 * Z / (n);
        norm = R_norm( n, l, Z );
        rho = A * s;
        ea = Math.exp( -0.5 * rho );
        pp = Math.pow( rho, l ); // gsl_sf_pow_int(rho, l);
        lag = gsl_sf_laguerre_n_e( n - l - 1, (2 * l + 1), rho );
        W_val = norm.val * ea * pp;
        W_err = norm.err * ea * pp;
        W_err = W_err + norm.val * ((0.5 * rho + 1.0) * GSL_DBL_EPSILON) * ea * pp;
        W_err = W_err + norm.val * ea * ((l + 1) * GSL_DBL_EPSILON) * pp;
        r.val = W_val * lag.val;
        r.err = W_val * lag.err + W_err * Math.abs( lag.val );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        //IF ((l = 0 OR (r > 0.0 AND l > 0)) AND lag.val /= 0.0 
        //    AND stat_lag = GSL_SUCCESS AND stat_norm = GSL_SUCCESS)
        //THEN
        //    CHECK_UNDERFLOW(result);
        //END IF;
        return r;
    }

} // gsl_sf_hydrogenicR_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hydrogenicR_1( Z, s )
{ // gsl_sf_hydrogenicR_1
    return EVAL_RESULT_DD( gsl_sf_hydrogenicR_1_e, { x: Z, y: s }, "gsl_sf_hydrogenicR_1" );
} // gsl_sf_hydrogenicR_1;

export function gsl_sf_hydrogenicR( n, l, Z, s )
{ // gsl_sf_hydrogenicR
    return EVAL_RESULT_IIDD( gsl_sf_hydrogenicR_e, { i1: n, i2: l, x: Z, y: s }, "gsl_sf_hydrogenicR" );
} // gsl_sf_hydrogenicR;

// ----------------------------------------------------------------------------
// EOF SF-CoulombBound.mjs
