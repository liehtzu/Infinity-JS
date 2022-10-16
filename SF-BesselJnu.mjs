// SF-BesselJnu.mjs
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

import { M_PI }                  from "./SF-Math.mjs";
import { GSL_SIGN }              from "./SF-Math.mjs";
import { GSL_SQRT_DBL_MIN }      from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { gsl_sf_sin_pi_e }       from "./SF-SinCosPi.mjs";
import { gsl_sf_cos_pi_e }       from "./SF-SinCosPi.mjs";
import { gsl_sf_bessel_IJ_taylor_e }  from "./SF-Bessel.mjs";
import { gsl_sf_bessel_J_CF1 }        from "./SF-Bessel.mjs";
import { gsl_sf_bessel_JY_steed_CF2 } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Jnu_asympx_e } from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Ynupos_e }     from "./SF-BesselYnu.mjs";
import { gsl_sf_bessel_Y_temme }      from "./SF-BesselTemme.mjs";
import { gsl_sf_bessel_Jnu_asymp_Olver_e } from "./SF-BesselOlver.mjs";

import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

export function gsl_sf_bessel_Jnupos_e( nu, x )
{

    var r = { val: 0.0, err: 0.0 };
    // CHECK_POINTER(result)

    if ( x == 0.0 )
    {
        if ( nu == 0.0 )
        {
            r.val = 1.0;
            r.err = 0.0;
        }
        else
        {
            r.val = 0.0;
            r.err = 0.0;
        }
        return r;
    }
    else if ( x * x < 10.0 * (nu + 1.0) )
    {
        return gsl_sf_bessel_IJ_taylor_e( nu, x, -1, 100, GSL_DBL_EPSILON);
    }
    else if (nu > 50.0)
    {
        return gsl_sf_bessel_Jnu_asymp_Olver_e( nu, x );
    }
    else if (x > 1000.0)
    {
        // We need this to avoid feeding large x to CF1; note that
        // due to the above check, we know that n <= 50. See similar
        // block in bessel_Jn.c.
        //
        return gsl_sf_bessel_Jnu_asympx_e( nu, x );
    }
    else
    {
        // -1/2 <= mu <= 1/2
        var N = Math.trunc( nu + 0.5 );
        var mu = nu - N;

        // Determine the J ratio at nu.
        var Jnup1_Jnu = 0.0;
        var sgn_Jnu   = 0.0;
        var r2 = { ratio: 0.0, sgn: 0.0 };
        r2 = gsl_sf_bessel_J_CF1( nu, x ); //, &Jnup1_Jnu, &sgn_Jnu );
        Jnup1_Jnu = r2.ratio;
        sgn_Jnu   = r2.sgn;

        if ( x < 2.0 )
        {
            // Determine Y_mu, Y_mup1 directly and recurse forward to nu.
            // Then use the CF1 information to solve for J_nu and J_nup1.
            //
            var Y_mu   = { val: 0.0, err: 0.0 };
            var Y_mup1 = { val: 0.0, err: 0.0 };
            gsl_sf_bessel_Y_temme( mu, x, Y_mu, Y_mup1 );

            var Ynm1 = Y_mu.val;
            var Yn   = Y_mup1.val;
            var Ynp1 = 0.0;

            for ( let n = 1; n < N; n++ )
            {
                Ynp1 = 2.0 * (mu + n) / x * Yn - Ynm1;
                Ynm1 = Yn;
                Yn   = Ynp1;
            }

            r.val = 2.0 / (M_PI * x) / (Jnup1_Jnu * Yn - Ynp1);
            r.err = GSL_DBL_EPSILON * (N + 2.0) * Math.abs( r.val );
            return r; // GSL_ERROR_SELECT_2(stat_mu, stat_CF1);
        }
        else
        {
            // Recurse backward from nu to mu, determining the J ratio
            // at mu. Use this together with a Steed method CF2 to
            // determine the actual J_mu, and thus obtain the normalization.
            //
            var Jmu = 0.0;
            var Jmup1_Jmu = 0.0;
            var sgn_Jmu = 0.0;
            var Jmuprime_Jmu = 0.0;
            var P = 0.0;
            var Q = 0.0;
            var r1 = { P: 0.0, Q: 0.0 };
            var gamma = 0.0;

            r1 = gsl_sf_bessel_JY_steed_CF2( mu, x ); //, &P, &Q);
            P = r1.P;
            Q = r1.Q;

            var Jnp1 = sgn_Jnu * GSL_SQRT_DBL_MIN * Jnup1_Jnu;
            var Jn   = sgn_Jnu * GSL_SQRT_DBL_MIN;
            var Jnm1 = 0.0;

            for ( let n = N; n > 0; n-- )
            {
                Jnm1 = 2.0 * (mu + n) / x * Jn - Jnp1;
                Jnp1 = Jn;
                Jn   = Jnm1;
            }
            Jmup1_Jmu = Jnp1 / Jn;
            sgn_Jmu   = GSL_SIGN( Jn );
            Jmuprime_Jmu = mu / x - Jmup1_Jmu;

            gamma = (P - Jmuprime_Jmu) / Q;
            Jmu   = sgn_Jmu * Math.sqrt( 2.0 / (M_PI * x) / (Q + gamma * (P - Jmuprime_Jmu)) );

            r.val = Jmu * (sgn_Jnu * GSL_SQRT_DBL_MIN) / Jn;
            r.err = 2.0 * GSL_DBL_EPSILON * (N + 2.0) * Math.abs( r.val );

            return r; // GSL_ERROR_SELECT_2(stat_CF2, stat_CF1);
        }
    }
}

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Jnu_e(nu, x)
{

    var r = { val: 0.0, err: 0.0 };
    // CHECK_POINTER(result)

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else if (nu < 0.0)
    {
        r = gsl_sf_bessel_Jnupos_e(-nu, x);
        var Jval = r.val;
        var Jerr = r.err;
        r = gsl_sf_bessel_Ynupos_e(-nu, x);
        var Yval = r.val;
        var Yerr = r.err;
        // var s = sin(M_PI*nu), c = cos(M_PI*nu);
        r = gsl_sf_sin_pi_e(nu);
        var s = r.val;
        var serr = r.err;
        r = gsl_sf_cos_pi_e(nu);
        var c = r.val;
        var cerr = r.err;
        r.val = s * Yval + c * Jval;
        r.err = Math.abs(c * Yerr) + Math.abs(s * Jerr) + Math.abs(cerr * Yval) + Math.abs(serr * Jval);
        return r; // GSL_ERROR_SELECT_4(Jstatus, Ystatus, sinstatus, cosstatus);
    }
    else
    {
        return gsl_sf_bessel_Jnupos_e(nu, x);
    }
}

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_Jnu( nu, x )
{ // gsl_sf_bessel_Jnu
    return EVAL_RESULT_DD( gsl_sf_bessel_Jnu_e, { x: nu, y: x }, "gsl_sf_bessel_Jnu" );
} // gsl_sf_bessel_Jnu

// ----------------------------------------------------------------------------
// EOF SF-BesselJnu.mjs
