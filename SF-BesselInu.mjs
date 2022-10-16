// SF-BesselInu.mjs
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

import { GSL_DBL_EPSILON }              from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON }        from "./SF-Machine.mjs";
import { gsl_sf_bessel_IJ_taylor_e }    from "./SF-Bessel.mjs";
import { gsl_sf_bessel_I_CF1_ser }      from "./SF-Bessel.mjs";
import { gsl_sf_exp_mult_err_e }        from "./SF-Exponential.mjs";
import { gsl_sf_bessel_K_scaled_temme } from "./SF-BesselTemme.mjs";
import { gsl_sf_bessel_K_scaled_steed_temme_CF2 } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Inu_scaled_asymp_unif_e }  from "./SF-Bessel.mjs";

import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Inu_scaled_e( nu, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x < 0.0 || nu < 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( x * x < 10.0 * (nu + 1.0) )
    {
        var b  = { val: 0.0, err: 0.0 }; // Result;
        var ex = 0.0;

        ex = Math.exp( -x );
        b = gsl_sf_bessel_IJ_taylor_e( nu, x, 1, 100, GSL_DBL_EPSILON );
        r.val = b.val * ex;
        r.err = b.err * ex;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
    }
    else if ( 0.5 / (nu * nu + x * x) < GSL_ROOT3_DBL_EPSILON )
    {
        r = gsl_sf_bessel_Inu_scaled_asymp_unif_e( nu, x );
    }
    else
    {
        var N          = Math.trunc( nu + 0.5 );
        var mu         = nu - (N); // -1/2 <= mu <= 1/2 
        var K_mu       = 0.0;
        var K_mup1     = 0.0;
        var Kp_mu      = 0.0;
        var K_nu       = 0.0;
        var K_nup1     = 0.0;
        var K_num1     = 0.0;
        var I_nu_ratio = 0.0;

        var rk = { K_nu: 0.0, K_nup1: 0.0, Kp_nu: 0.0 };

        // obtain K_mu, K_mup1
        if ( x < 2.0 )
        {
            rk = gsl_sf_bessel_K_scaled_temme( mu, x );//, K_mu, K_mup1, Kp_mu);
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
        
        // calculate I_{nu+1}/I_nu
        I_nu_ratio = gsl_sf_bessel_I_CF1_ser(nu, x);
        
        // solve for I_nu
        r.val = 1.0 / (x * (K_nup1 + I_nu_ratio * K_nu));
        r.err = GSL_DBL_EPSILON * (0.5 * (N) + 2.0) * Math.abs( r.val );
    }

    return r;

} // gsl_sf_bessel_Inu_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Inu_e( nu, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_bessel_Inu_scaled_e( nu, x );
    r = gsl_sf_exp_mult_err_e( x, Math.abs( x * GSL_DBL_EPSILON ), r.val, r.err );
    return r;

} // gsl_sf_bessel_Inu_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Inu_scaled( nu, x )
{
    return EVAL_RESULT_DD( gsl_sf_bessel_Inu_scaled_e, { x: nu, y: x }, "gsl_sf_bessel_Inu_scaled" );
}// gsl_sf_bessel_Inu_scaled

export function gsl_sf_bessel_Inu( nu, x )
{
    return EVAL_RESULT_DD( gsl_sf_bessel_Inu_e, { x: nu, y: x }, "gsl_sf_bessel_Inu" );
} // gsl_sf_bessel_Inu

// ----------------------------------------------------------------------------
// EOF SF-BesselInu.mjs
