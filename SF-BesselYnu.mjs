// SF-BesselYnu.mjs
// ----------------------------------------------------------------------------
// 
// Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, 2017 Konrad Griessinger
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
import { gsl_sf_sin_pi_e }       from "./SF-SinCosPi.mjs";
import { gsl_sf_cos_pi_e }       from "./SF-SinCosPi.mjs";
import { gsl_sf_bessel_Jnupos_e }          from "./SF-BesselJnu.mjs";
import { gsl_sf_bessel_JY_mu_restricted }  from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Y_temme }           from "./SF-BesselTemme.mjs";
import { gsl_sf_bessel_Ynu_asymp_Olver_e } from "./SF-BesselOlver.mjs";

import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Ynupos_e( nu, x )
{
    var r = { val: 0.0, err: 0.0 };
    /* CHECK_POINTER(result) */

    if ( nu > 50.0 )
    {
        return gsl_sf_bessel_Ynu_asymp_Olver_e( nu, x );
    }
    else
    {
        // -1/2 <= mu <= 1/2
        var N = Math.trunc( nu + 0.5 );
        var mu = nu - N;

        var Y_mu   = { val: 0.0, err: 0.0 };
        var Y_mup1 = { val: 0.0, err: 0.0 };
        var Ynm1 = 0.0;
        var Yn   = 0.0;
        var Ynp1 = 0.0;
        var n = 0;

        if ( x < 2.0 )
        {
            // Determine Ymu, Ymup1 directly. This is really
            // an optimization since this case could as well
            // be handled by a call to gsl_sf_bessel_JY_mu_restricted(),
            // as below.
            //
            gsl_sf_bessel_Y_temme( mu, x, Y_mu, Y_mup1 );
        }
        else
        {
            // Determine Ymu, Ymup1 and Jmu, Jmup1.
            //
            var J_mu   = { val: 0.0, err: 0.0 };
            var J_mup1 = { val: 0.0, err: 0.0 };
            gsl_sf_bessel_JY_mu_restricted( mu, x, J_mu, J_mup1, Y_mu, Y_mup1 );
        }

        // Forward recursion to get Ynu, Ynup1.
        //
        Ynm1 = Y_mu.val;
        Yn   = Y_mup1.val;
        for ( n = 1; n <= N; n++ )
        {
            Ynp1 = 2.0 * (mu + n) / x * Yn - Ynm1;
            Ynm1 = Yn;
            Yn   = Ynp1;
        }

        r.val  = Ynm1; // Y_nu
        r.err  = (N + 1.0) * Math.abs( Ynm1 ) * (Math.abs( Y_mu.err / Y_mu.val ) + Math.abs( Y_mup1.err / Y_mup1.val ));
        r.err += 2.0 * GSL_DBL_EPSILON * Math.abs( Ynm1 );

        return r;
    }
}

export function gsl_sf_bessel_Ynu_e( nu, x )
{
    var r = { val: 0.0, err: 0.0 };
    /* CHECK_POINTER(result) */

    if ( x <= 0.0 )
    {
      throw "SF.DomainException";
    }
    else if ( nu < 0.0 )
    {
        r = gsl_sf_bessel_Jnupos_e( -nu, x );
        var Jval = r.val;
        var Jerr = r.err;
        r = gsl_sf_bessel_Ynupos_e( -nu, x );
        var Yval = r.val;
        var Yerr = r.err;
        // var s = sin(M_PI*nu), c = cos(M_PI*nu);
        r = gsl_sf_sin_pi_e( nu );
        var s = r.val;
        var serr = r.err;
        r = gsl_sf_cos_pi_e( nu );
        var c = r.val;
        var cerr = r.err;
        r.val = c * Yval - s * Jval;
        r.err = Math.abs( c * Yerr ) + Math.abs( s * Jerr ) + Math.abs( cerr * Yval ) + Math.abs( serr * Jval );
        return r; //GSL_ERROR_SELECT_4(Jstatus, Ystatus, sinstatus, cosstatus);
    }
    else
    {
        return gsl_sf_bessel_Ynupos_e( nu, x );
    }
}

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Ynu( nu, x )
{
    return EVAL_RESULT_DD( gsl_sf_bessel_Ynu_e, { x: nu, y: x }, "gsl_sf_bessel_Ynu_e" );
}

// ----------------------------------------------------------------------------
// EOF SF-BesselYnu.mjs
