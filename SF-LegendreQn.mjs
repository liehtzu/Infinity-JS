// SF-LegendreQn.mjs
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
import { GSL_ROOT4_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_ROOT6_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_DBL_MIN      }      from "./SF-Machine.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { M_SQRT3 }               from "./SF-Math.mjs";
import { gsl_sf_exp_mult_e }     from "./SF-Exponential.mjs";
import { gsl_sf_multiply_e }     from "./SF-Elementary.mjs";
import { gsl_sf_bessel_Y0_e }    from "./SF-BesselY0.mjs";
import { gsl_sf_bessel_Y1_e }    from "./SF-BesselY1.mjs";
import { gsl_sf_bessel_K0_scaled_e } from "./SF-BesselK0.mjs";
import { gsl_sf_bessel_K1_scaled_e } from "./SF-BesselK1.mjs";

import { EVAL_RESULT_D }  from "./SF-Evaluate.mjs";
import { EVAL_RESULT_ID } from "./SF-Evaluate.mjs";

// -- Evaluate f_{ell+1}/f_ell
// -- f_ell = Q^{b}_{a+ell}(x)
// -- x > 1
// --
function legendreQ_CF1_xgt1( ell, a, b, x )
{
    const RECUR_BIG = GSL_SQRT_DBL_MAX;
    const maxiter   = 5000;
    var n      = 0;
    var Anm2   = 0.0;
    var Bnm2   = 0.0;
    var Anm1   = 0.0;
    var Bnm1   = 0.0;
    var a1     = 0.0;
    var b1     = 0.0;
    var An     = 0.0;
    var Bn     = 0.0;
    var an1    = 0.0;
    var bn1    = 0.0;
    var fn     = 0.0;
    var old_fn = 0.0;
    var del    = 0.0;
    var lna    = 0.0;
    var r      = 0.0;

    Anm2 = 1.0;
    Bnm2 = 0.0;
    Anm1 = 0.0;
    Bnm1 = 1.0;
    a1 = (ell) + 1.0 + a + b;
    b1 = (2.0 * ((ell + 1) + a) + 1.0) * x;
    An = b1 * Anm1 + a1 * Anm2;
    Bn = b1 * Bnm1 + a1 * Bnm2;
    fn = An / Bn;

    n = 1;
    while (n < maxiter)
    {
        n = n + 1;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        lna  = (ell + n) + a;
        an1  = b * b - lna * lna;
        bn1  = (2.0 * lna + 1.0) * x;
        An   = bn1 * Anm1 + an1 * Anm2;
        Bn   = bn1 * Bnm1 + an1 * Bnm2;

        if ( Math.abs( An ) > RECUR_BIG || Math.abs( Bn ) > RECUR_BIG )
        {
            An = An / RECUR_BIG;
            Bn = Bn / RECUR_BIG;
            Anm1 = Anm1 / RECUR_BIG;
            Bnm1 = Bnm1 / RECUR_BIG;
            Anm2 = Anm2 / RECUR_BIG;
            Bnm2 = Bnm2 / RECUR_BIG;
        }

        old_fn = fn;
        fn     = An / Bn;
        del    = old_fn / fn;

        if ( Math.abs( del - 1.0) < 4.0 * GSL_DBL_EPSILON ) break;
    }

    r = fn;

    if ( n >= maxiter )
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // legendreQ_CF1_xgt1

// ----------------------------------------------------------------------------

// Uniform asymptotic for Q_l(x).
// Assumes x > -1.0 and x != 1.0.
// Discards second order and higher terms.
//
function legendre_Ql_asymp_unif( ell, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x < 1.0 )
    {
        let u   = 0.0;
        let th  = 0.0;
        let pre = 0.0;
        let B00 = 0.0;
        let sum = 0.0;
        let sin_th = 0.0;
        let cot_th = 0.0;
        let Y0  = { val: 0.0, err: 0.0 }; // Result;
        let Y1  = { val: 0.0, err: 0.0 }; // Result;

        u   = ell + 0.5;
        th  = Math.acos( x );

        // B00 = 1/8 (1 - th cot(th) / th^2
        // pre = sqrt(th/sin(th))
        //
        if ( th < GSL_ROOT4_DBL_EPSILON )
        {
            B00 = (1.0 + th * th / 15.0) / 24.0;
            pre = 1.0 + th * th / 12.0;
        }
        else
        {
            sin_th = Math.sqrt( 1.0 - x * x );
            cot_th = x / sin_th;
            B00 = 1.0 / 8.0 * (1.0 - th * cot_th) / (th * th);
            pre = Math.sqrt( th / sin_th );
        }

        Y0 = gsl_sf_bessel_Y0_e( u * th );
        Y1 = gsl_sf_bessel_Y1_e( u * th );
        
        sum = -0.5 * M_PI * (Y0.val + th / u * Y1.val * B00);
        
        r = gsl_sf_multiply_e( pre, sum );
        r.err = r.err + 0.5 * M_PI * Math.abs( pre ) * (Y0.err + Math.abs( th / u * B00 ) * Y1.err);
        r.err = r.err + GSL_DBL_EPSILON * Math.abs( r.val );
        
        return r;
    }
    else
    {
        let u         = 0.0;
        let xi        = 0.0;
        let pre       = 0.0;
        let B00       = 0.0;
        let sum       = 0.0;
        let sinh_xi   = 0.0;
        let coth_xi   = 0.0;
        let K0_scaled = { val: 0.0, err: 0.0 }; // Result;
        let K1_scaled = { val: 0.0, err: 0.0 }; // Result;

        u  = ell + 0.5;
        xi = Math.acosh( x );

        // B00 = -1/8 (1 - xi coth(xi) / xi^2
        // pre = sqrt(xi/sinh(xi))
        //
        if ( xi < GSL_ROOT4_DBL_EPSILON )
        {
            B00 = (1.0 - xi * xi / 15.0) / 24.0;
            pre = 1.0 - xi * xi / 12.0;
        }
        else
        {
            sinh_xi = Math.sqrt( x * x - 1.0 );
            coth_xi = x / sinh_xi;
            B00 = -1.0 / 8.0 * (1.0 - xi * coth_xi) / (xi * xi);
            pre = Math.sqrt( xi / sinh_xi );
        }

        K0_scaled = gsl_sf_bessel_K0_scaled_e( u * xi );
        K1_scaled = gsl_sf_bessel_K1_scaled_e( u * xi );
        
        sum = K0_scaled.val - xi / u * K1_scaled.val * B00;
        
        r = gsl_sf_exp_mult_e( -u * xi, pre * sum );
        r.err = GSL_DBL_EPSILON * Math.abs( r.val ) * Math.abs( u * xi );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        
        return r;
    }

} // legendre_Ql_asymp_unif

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_legendre_Q0_e( x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x <= -1.0 || x == 1.0 )
    {
        throw "SF.DomainException";
    }
    else if ( x * x < GSL_ROOT6_DBL_EPSILON ) // |x| <~ 0.05
    {
        let c3  = 1.0 / 3.0;
        let c5  = 1.0 / 5.0;
        let c7  = 1.0 / 7.0;
        let c9  = 1.0 / 9.0;
        let c11 = 1.0 / 11.0;
        let y   = x * x;
        let series = 1.0 + y * (c3 + y * (c5 + y * (c7 + y * (c9 + y * c11))));

        r.val = x * series;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( x );
        return r;
    }
    else if ( x < 1.0 )
    {
        r.val = 0.5 * Math.log( (1.0 + x) / (1.0 - x) );
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( x < 10.0 )
    {
        r.val = 0.5 * Math.log( (x + 1.0) / (x - 1.0) );
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( x * GSL_DBL_MIN < 2.0 )
    {
        let y  = 1.0 / (x * x);
        let c1 = 1.0 / 3.0;
        let c2 = 1.0 / 5.0;
        let c3 = 1.0 / 7.0;
        let c4 = 1.0 / 9.0;
        let c5 = 1.0 / 11.0;
        let c6 = 1.0 / 13.0;
        let c7 = 1.0 / 15.0;

        r.val = (1.0 / x) * (1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * c7)))))));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else
    {
        throw "SF.UnderflowException";
    }

} // gsl_sf_legendre_Q0_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_Q1_e( x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x <= -1.0 || x == 1.0 )
    {
        throw "SF.DomainException";
    }
    else if ( x * x < GSL_ROOT6_DBL_EPSILON ) // |x| <~ 0.05
    {
        let c3  = 1.0 / 3.0;
        let c5  = 1.0 / 5.0;
        let c7  = 1.0 / 7.0;
        let c9  = 1.0 / 9.0;
        let c11 = 1.0 / 11.0;
        let y   = x * x;
        let series = 1.0 + y * (c3 + y * (c5 + y * (c7 + y * (c9 + y * c11))));

        r.val = x * x * series - 1.0;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( x < 1.0 )
    {
        r.val = 0.5 * x * (Math.log( (1.0 + x) / (1.0 - x) )) - 1.0;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( x < 6.0 )
    {
        r.val = 0.5 * x * Math.log( (x + 1.0) / (x - 1.0) ) - 1.0;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( x * GSL_SQRT_DBL_MIN < 0.99 / M_SQRT3 )
    {
        let y  = 1.0 / (x * x);
        let c1 = 3.0 / 5.0;
        let c2 = 3.0 / 7.0;
        let c3 = 3.0 / 9.0;
        let c4 = 3.0 / 11.0;
        let c5 = 3.0 / 13.0;
        let c6 = 3.0 / 15.0;
        let c7 = 3.0 / 17.0;
        let c8 = 3.0 / 19.0;
        let sum = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * (c7 + y * c8)))))));

        r.val = sum / (3.0 * x * x);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else
    {
        throw "SF.UnderflowException";
    }

} // gsl_sf_legendre_Q1_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_Ql_e( l, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x <= -1.0 || x == 1.0 || l < 0 )
    {
        throw "SF.DomainException";
    }
    else if ( l == 0 )
    {
        return gsl_sf_legendre_Q0_e( x );
    }
    else if ( l == 1 )
    {
        return gsl_sf_legendre_Q1_e( x );
    }
    else if ( l > 100000 )
    {
        return legendre_Ql_asymp_unif( (l), x );
    }
    else if ( x < 1.0 )
    {
        // Forward recurrence.
        //
        let Qellm1 = 0.0;
        let Qell   = 0.0;
        let Qellp1 = 0.0;
        let Q0 = { val: 0.0, err: 0.0 }; // Result;
        let Q1 = { val: 0.0, err: 0.0 }; // Result;

        Q0 = gsl_sf_legendre_Q0_e( x );
        Q1 = gsl_sf_legendre_Q1_e( x );
        Qellm1 = Q0.val;
        Qell   = Q1.val;
        for ( let ell = 1; ell <= l - 1; ell++ )
        {
            Qellp1 = (x * (2 * ell + 1) * Qell - (ell) * Qellm1) / (ell + 1);
            Qellm1 = Qell;
            Qell   = Qellp1;
        }
        r.val = Qell;
        r.err = GSL_DBL_EPSILON * (l) * Math.abs( r.val );
        return r;
    }
    else
    {
        // x > 1.0
        let rat    = 0.0;
        let Qellp1 = 0.0;
        let Qell   = 0.0;
        let Qellm1 = 0.0;
        let Q0     = { val: 0.0, err: 0.0 }; // Result;
        let Q1     = { val: 0.0, err: 0.0 }; // Result;

        rat    = legendreQ_CF1_xgt1( l, 0.0, 0.0, x );
        Qellp1 = rat * GSL_SQRT_DBL_MIN;
        Qell   = GSL_SQRT_DBL_MIN;

        for ( let ell = l; ell >= 1; ell-- )
        {
            Qellm1 = (x * (2 * ell + 1) * Qell - (ell + 1) * Qellp1) / (ell);
            Qellp1 = Qell;
            Qell   = Qellm1;
        }

        if ( Math.abs( Qell ) > Math.abs( Qellp1 ) )
        {
            Q0 = gsl_sf_legendre_Q0_e( x );
            r.val = GSL_SQRT_DBL_MIN * Q0.val / Qell;
            r.err = (l) * GSL_DBL_EPSILON * Math.abs( r.val );
        }
        else
        {
            Q1 = gsl_sf_legendre_Q1_e( x );
            r.val = GSL_SQRT_DBL_MIN * Q1.val / Qellp1;
            r.err = (l) * GSL_DBL_EPSILON * Math.abs( r.val );
        }

        return r;
    }

} // gsl_sf_legendre_Ql_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

export function gsl_sf_legendre_Q0( x )
{ // gsl_sf_legendre_Q0
    return EVAL_RESULT_D( gsl_sf_legendre_Q0_e, x, "gsl_sf_legendre_Q0" );
} // gsl_sf_legendre_Q0

export function gsl_sf_legendre_Q1( x )
{ // gsl_sf_legendre_Q1
    return EVAL_RESULT_D( gsl_sf_legendre_Q1_e, x, "gsl_sf_legendre_Q1" );
} // gsl_sf_legendre_Q1

export function gsl_sf_legendre_Ql_ID( l, x )
{ // gsl_sf_legendre_Ql
    return EVAL_RESULT_ID( gsl_sf_legendre_Ql_e, { n: l, x: x }, "gsl_sf_legendre_Ql" );
} // gsl_sf_legendre_Ql

// -- ----------------------------------------------------------------------------
// EOF SF-LegendreQn.mjs