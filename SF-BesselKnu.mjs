// SF-BesselKnu.mjs
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

import { GSL_LOG_DBL_MAX }       from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { M_LN2 }                 from "./SF-Math.mjs";
import { gsl_sf_exp_mult_err_e } from "./SF-Exponential.mjs";
import { gsl_sf_lngamma_e }      from "./SF-Gamma.mjs";
import { gsl_sf_bessel_K0_scaled_e } from "./SF-BesselK0.mjs";
import { gsl_sf_bessel_K_scaled_temme } from "./SF-BesselTemme.mjs";
import { gsl_sf_bessel_K_scaled_steed_temme_CF2 } from "./SF-Bessel.mjs";

import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Knu_scaled_e( nu, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x <= 0.0 || nu < 0.0 )
    {
        throw "SF.DomainException";
    }
    else
    {
        var N      = Math.trunc( nu + 0.5 );
        var mu     = nu - (N); // -1/2 <= mu <= 1/2
        var K_mu   = 0.0;
        var K_mup1 = 0.0;
        var Kp_mu  = 0.0;
        var K_nu   = 0.0;
        var K_nup1 = 0.0;
        var K_num1 = 0.0;
        var rk = { K_nu: 0.0, K_nup1: 0.0, Kp_nu: 0.0 };

        if ( x < 2.0 )
        {
            rk = gsl_sf_bessel_K_scaled_temme( mu, x ); //, K_mu, K_mup1, Kp_mu);
        }
        else
        {
            rk = gsl_sf_bessel_K_scaled_steed_temme_CF2( mu, x ); //, K_mu, K_mup1, Kp_mu);
        }
        K_mu = rk.K_nu;
        K_mup1 = rk.K_nup1;
        Kp_mu = rk.Kp_nu;
        
        // recurse forward to obtain K_num1, K_nu
        K_nu   = K_mu;
        K_nup1 = K_mup1;
        
        for ( let i = 0; i <= N - 1; i++ )
        {
            K_num1 = K_nu;
            K_nu   = K_nup1;
            K_nup1 = 2.0 * (mu + (i + 1)) / x * K_nu + K_num1;
        }
        
        r.val = K_nu;
        r.err = 2.0 * GSL_DBL_EPSILON * (N + 4) * Math.abs( r.val );
    }

    return r;

} // gsl_sf_bessel_Knu_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Knu_e( nu, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_bessel_Knu_scaled_e( nu, x );
    r = gsl_sf_exp_mult_err_e( -x, 0.0, r.val, r.err );
    return r;

} // gsl_sf_bessel_Knu_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_lnKnu_e( nu, x )
{
    var K_scaled = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x <= 0.0 || nu < 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( nu == 0.0 )
    {
        // This cannot underflow, and
        // it will not throw GSL_EDOM
        // since that is already checked.
        //
        K_scaled = gsl_sf_bessel_K0_scaled_e( x );
        r.val = -x + Math.log( Math.abs( K_scaled.val ) );
        r.err = GSL_DBL_EPSILON * Math.abs( x ) + Math.abs( K_scaled.err / K_scaled.val );
        r.err = r.err + GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( x < 2.0 && nu > 1.0 )
    {
        // Make use of the inequality
        // Knu(x) <= 1/2 (2/x)^nu Gamma(nu),
        // which follows from the integral representation
        // [Abramowitz+Stegun, 9.6.23 (2)]. With this
        // we decide whether or not there is an overflow
        // problem because x is small.
        //
        var ln_bound = 0.0;
        var xi    = 0.0;
        var sum   = 0.0;
        var lg_nu = { val: 0.0, err: 0.0 }; // Result;

        lg_nu = gsl_sf_lngamma_e( nu );
        ln_bound = -M_LN2 - nu * Math.log( 0.5 * x ) + lg_nu.val;
        if ( ln_bound > GSL_LOG_DBL_MAX - 20.0 )
        {
            // x must be very small or nu very large (or both).
            //
            xi  = 0.25 * x * x;
            sum = 1.0 - xi / (nu - 1.0);
            if ( nu > 2.0 )
            {
                sum = sum + (xi / (nu - 1.0)) * (xi / (nu - 2.0));
            }
            r.val = ln_bound + Math.log( sum );
            r.err = lg_nu.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            return r;
        }
        // can drop-through here
    }
   
   
    // We passed the above tests, so no problem.
    // Evaluate as usual. Note the possible drop-through
    // in the above code!
    //
    K_scaled = gsl_sf_bessel_Knu_scaled_e( nu, x );
    r.val = -x + Math.log( Math.abs( K_scaled.val ) );
    r.err = GSL_DBL_EPSILON * Math.abs( x ) + Math.abs( K_scaled.err / K_scaled.val );
    r.err = r.err + GSL_DBL_EPSILON * Math.abs( r.val );
    return r;

} // gsl_sf_bessel_lnKnu_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Knu_scaled( nu, x )
{ // gsl_sf_bessel_Knu_scaled
    return EVAL_RESULT_DD( gsl_sf_bessel_Knu_scaled_e, { x: nu, y: x }, "gsl_sf_bessel_Knu_scaled" );
} // gsl_sf_bessel_Knu_scaled;

export function gsl_sf_bessel_Knu( nu, x )
{ // gsl_sf_bessel_Knu
    return EVAL_RESULT_DD( gsl_sf_bessel_Knu_e, { x: nu, y: x }, "gsl_sf_bessel_Knu" );
} // gsl_sf_bessel_Knu;

export function gsl_sf_bessel_lnKnu( nu, x )
{ // gsl_sf_bessel_lnKnu
    return EVAL_RESULT_DD( gsl_sf_bessel_lnKnu_e, { x: nu, y: x }, "gsl_sf_bessel_lnKnu" );
} // gsl_sf_bessel_lnKnu;

// ----------------------------------------------------------------------------
// EOF SF-BesselKnu.mjs
