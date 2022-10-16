// SF-LegendreH3d.mjs
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
import { GSL_ROOT5_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON }   from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }      from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }       from "./SF-Machine.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { M_LNPI }                from "./SF-Math.mjs";
import { M_LN2 }                 from "./SF-Math.mjs";
import { gsl_sf_cos_e }          from "./SF-Trigonometric.mjs";
import { gsl_sf_sin_e }          from "./SF-Trigonometric.mjs";
import { gsl_sf_sin_err_e }      from "./SF-Trigonometric.mjs";
import { gsl_sf_lnsinh_e }       from "./SF-Trigonometric.mjs";
import { gsl_sf_hypot }          from "./SF-Trigonometric.mjs";
import { gsl_sf_exp_mult_err_e } from "./SF-Exponential.mjs";
import { gsl_sf_lngamma_complex_e } from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_e }         from "./SF-Gamma.mjs";
import { gsl_sf_conicalP_large_x_e } from "./SF-LegendreConical.mjs";
import { gsl_sf_conicalP_xgt1_neg_mu_largetau_e } from "./SF-LegendreConical.mjs";

import { EVAL_RESULT_DD }  from "./SF-Evaluate.mjs";
import { EVAL_RESULT_IDD } from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// See [Abbott+Schaefer, Ap.J. 308, 546 (1986)] for
// enough details to follow what is happening here.
//
//
//
// Logarithm of normalization factor, Log[N(ell,lambda)].
// N(ell,lambda) = Product[ lambda^2 + n^2, {n,0,ell} ]
//               = |Gamma(ell + 1 + I lambda)|^2  lambda sinh(Pi lambda) / Pi
// Assumes ell >= 0.
//
function legendre_H3d_lnnorm( ell, lambda )
{

    var abs_lam = 0.0;

    abs_lam = Math.abs( lambda );

    if ( abs_lam == 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( lambda > (ell + 1) / GSL_ROOT3_DBL_EPSILON )
    {
        // There is a cancellation between the sinh(Pi lambda)
        // term and the log(gamma(ell + 1 + i lambda) in the
        // result below, so we show some care and save some digits.
        // Note that the above guarantees that lambda is large,
        // since ell >= 0. We use Stirling and a simple expansion
        // of sinh.
        //
        let rat          = 0.0;
        let ln_lam2ell2  = 0.0;
        let lg_corrected = 0.0;
        let angle_terms  = 0.0;

        rat          = (ell + 1) / lambda;
        ln_lam2ell2  = 2.0 * Math.log( lambda ) + Math.log( 1.0 + rat * rat );
        lg_corrected = -2.0 * (ell + 1) + M_LNPI + ((ell) + 0.5) * ln_lam2ell2 + 1.0 / (288.0 * lambda * lambda);
        angle_terms  = lambda * 2.0 * rat * (1.0 - rat * rat / 3.0);
        return Math.log(abs_lam) + lg_corrected + angle_terms - M_LNPI;
    }
    else
    {
        let lg_r     = { val: 0.0, err: 0.0 }; // Result;
        let lg_theta = { val: 0.0, err: 0.0 }; // Result;
        let ln_sinh  = { val: 0.0, err: 0.0 }; // Result;

        gsl_sf_lngamma_complex_e( ell + 1, lambda, lg_r, lg_theta );
        ln_sinh = gsl_sf_lnsinh_e( M_PI * abs_lam );
        return Math.log(abs_lam) + ln_sinh.val + 2.0 * lg_r.val - M_LNPI;
    }

} // legendre_H3d_lnnorm
 
// ----------------------------------------------------------------------------
 
// Calculate series for small eta*lambda.
// Assumes eta > 0, lambda != 0.
//
// This is just the defining hypergeometric for the Legendre function.
//
// P^{mu}_{-1/2 + I lam}(z) = 1/Gamma(l+3/2) ((z+1)/(z-1)^(mu/2)
//                            2F1(1/2 - I lam, 1/2 + I lam; l+3/2; (1-z)/2)
// We use
//       z = cosh(eta)
// (z-1)/2 = sinh^2(eta/2)
//
// And recall
// H3d = sqrt(Pi Norm /(2 lam^2 sinh(eta))) P^{-l-1/2}_{-1/2 + I lam}(cosh(eta))
//
function legendre_H3d_series( ell, lambda, eta )
{

    const nmax      = 5000;
    var n         = 0;
    var shheta    = 0.0;
    var ln_zp1    = 0.0;
    var ln_zm1    = 0.0;
    var zeta      = 0.0;
    var term      = 0.0;
    var sum       = 0.0;
    var sum_err   = 0.0;
    var lnN       = 0.0;
    var lnpre_val = 0.0;
    var lnpre_err = 0.0;
    var lnprepow  = 0.0;
    var aR        = 0.0;
    var lg_lp32   = { val: 0.0, err: 0.0 }; // Result;
    var lnsheta   = { val: 0.0, err: 0.0 }; // Result;
    var r         = { val: 0.0, err: 0.0 }; // Result;

    shheta  = Math.sinh( 0.5 * eta );
    ln_zp1  = M_LN2 + Math.log( 1.0 + shheta * shheta );
    ln_zm1  = M_LN2 + 2.0 * Math.log( shheta );
    zeta    = -shheta * shheta;
    term    = 1.0;
    sum     = 1.0;
    sum_err = 0.0;

    lg_lp32 = gsl_sf_lngamma_e( (ell) + 3.0 / 2.0 );
    lnsheta = gsl_sf_lnsinh_e( eta );
    lnN = legendre_H3d_lnnorm( ell, lambda );
    lnprepow = 0.5 * ((ell) + 0.5) * (ln_zm1 - ln_zp1);
    lnpre_val = lnprepow + 0.5 * (lnN + M_LNPI - M_LN2 - lnsheta.val) - lg_lp32.val - Math.log( Math.abs( lambda ) );
    lnpre_err = lnsheta.err + lg_lp32.err + GSL_DBL_EPSILON * Math.abs( lnpre_val );
    lnpre_err = lnpre_err + 2.0 * GSL_DBL_EPSILON * (Math.abs( lnN ) + M_LNPI + M_LN2);
    lnpre_err = lnpre_err + 2.0 * GSL_DBL_EPSILON * (0.5 * ((ell) + 0.5) * (Math.abs( ln_zm1 ) + Math.abs( ln_zp1 )));
    n = 1;
    while ( n < nmax )
    {
        aR = (n) - 0.5;
        term = term * (aR * aR + lambda * lambda) * zeta / ((ell) + (n) + 0.5) / (n);
        sum  = sum  + term;
        sum_err = sum_err + 2.0 * GSL_DBL_EPSILON * Math.abs( term );
        if ( Math.abs( term / sum ) < 2.0 * GSL_DBL_EPSILON ) break;
        n = n + 1;
    }

    r = gsl_sf_exp_mult_err_e( lnpre_val, lnpre_err, sum, Math.abs( term ) + sum_err );
    if ( n >= nmax )
    {
        throw "SF.MaxIterationsException";
    }
    return r;

} // legendre_H3d_series

// ----------------------------------------------------------------------------

// Evaluate legendre_H3d(ell+1)/legendre_H3d(ell)
// by continued fraction. Use the Gautschi (Euler)
// equivalent series.
//
// FIXME: Maybe we have to worry about this. The a_k are
// not positive and there can be a blow-up. It happened
// for J_nu once or twice. Then we should probably use
// the method above.
//
function legendre_H3d_CF1_ser( ell, lambda, coth_eta )
{

    const maxk    = 20000;
    var pre     = 0.0;
    var tk      = 0.0;
    var sum     = 0.0;
    var rhok    = 0.0;
    var sum_err = 0.0;
    var tlk     = 0.0;
    var l1k     = 0.0;
    var ak      = 0.0;
    var k = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    pre  = gsl_sf_hypot( lambda, (ell + 1)) / ((2 * ell + 3) * coth_eta );
    tk   = 1.0;
    sum  = 1.0;
    rhok = 0.0;
    sum_err = 0.0;
 
    k = k + 1;
    while ( k < maxk )
    {
        tlk  = (2 * ell + 1 + 2 * k);
        l1k  = (ell + 1 + k);
        ak   = -(lambda * lambda + l1k * l1k) / (tlk * (tlk + 2.0) * coth_eta * coth_eta);
        rhok = -ak * (1.0 + rhok) / (1.0 + ak * (1.0 + rhok));
        tk   = tk * rhok;
        sum  = sum + tk;
        sum_err = sum_err + 2.0 * GSL_DBL_EPSILON * (k) * Math.abs( tk );
        if (Math.abs( tk / sum ) < GSL_DBL_EPSILON ) break;
        k = k + 1;
    }

    r.val = pre * sum;
    r.err = Math.abs( pre * tk );
    r.err = r.err + Math.abs( pre * sum_err );
    r.err = r.err + 4.0 * GSL_DBL_EPSILON * Math.abs( r.val );

    if ( k >= maxk )
    {
        throw "SF.MaxIterationsException";
    }
    return r;

} // legendre_H3d_CF1_ser

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_legendre_H3d_0_e( lambda, eta )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( eta < 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( eta == 0.0 || lambda == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        let lam_eta = 0.0;
        let f       = 0.0;
        let s       = { val: 0.0, err: 0.0 }; // Result;

        lam_eta = lambda * eta;
        s = gsl_sf_sin_err_e( lam_eta, 2.0 * GSL_DBL_EPSILON * Math.abs( lam_eta ) );
        if ( eta > -0.5 * GSL_LOG_DBL_EPSILON )
        {
            f = 2.0 / lambda * Math.exp( -eta );
            r.val = f * s.val;
            r.err = Math.abs( f * s.val ) * (Math.abs( eta ) + 1.0) * GSL_DBL_EPSILON;
            r.err = r.err + Math.abs( f ) * s.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val);
        }
        else
        {
            f = 1.0 / (lambda * Math.sinh( eta ));
            r.val = f * s.val;
            r.err = Math.abs( f * s.val ) * (Math.abs( eta ) + 1.0) * GSL_DBL_EPSILON;
            r.err = r.err + Math.abs( f ) * s.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        }
        return r;
    }

} // gsl_sf_legendre_H3d_0_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_H3d_1_e( lambda, eta )
{
    var xi    = 0.0;
    var lsq   = 0.0;
    var lsqp1 = 0.0;
    var r     = { val: 0.0, err: 0.0 }; // Result;

    xi    = Math.abs( eta * lambda );
    lsq   = lambda * lambda;
    lsqp1 = lsq + 1.0;

    if ( eta < 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( eta == 0.0 || lambda == 0.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else if ( xi < GSL_ROOT5_DBL_EPSILON && eta < GSL_ROOT5_DBL_EPSILON )
    {
        let etasq = 0.0;
        let xisq  = 0.0;
        let term1 = 0.0;
        let term2 = 0.0;
        let sinh_term = 0.0;
        let pre   = 0.0;

        etasq = eta * eta;
        xisq  = xi * xi;
        term1 = (etasq + xisq) / 3.0;
        term2 = -(2.0 * etasq * etasq + 5.0 * etasq * xisq + 3.0 * xisq * xisq) / 90.0;
        sinh_term = 1.0 - eta * eta / 6.0 * (1.0 - 7.0 / 60.0 * eta * eta);
        pre   = sinh_term / Math.sqrt( lsqp1 ) / eta;
        r.val = pre * (term1 + term2);
        r.err = pre * GSL_DBL_EPSILON * (Math.abs( term1 ) + Math.abs( term2 ));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else
    {
        let sin_term     = 0.0; //  Sin(xi)/xi
        let cos_term     = 0.0; //  Cos(xi)
        let coth_term    = 0.0; //  eta/Tanh(eta)
        let sinh_term    = 0.0; //  eta/Sinh(eta)
        let sin_term_err = 0.0;
        let cos_term_err = 0.0;
        let t1           = 0.0;
        let pre_val      = 0.0;
        let pre_err      = 0.0;
        let term1        = 0.0;
        let term2        = 0.0;

        if ( xi < GSL_ROOT5_DBL_EPSILON )
        {
            sin_term = 1.0 - xi * xi / 6.0 * (1.0 - xi * xi / 20.0);
            cos_term = 1.0 - 0.5 * xi * xi * (1.0 - xi * xi / 12.0);
            sin_term_err = GSL_DBL_EPSILON;
            cos_term_err = GSL_DBL_EPSILON;
        }
        else
        {
            let sin_xi_result = { val: 0.0, err: 0.0 }; // Result;
            let cos_xi_result = { val: 0.0, err: 0.0 }; // Result;

            sin_xi_result = gsl_sf_sin_e( xi );
            cos_xi_result = gsl_sf_cos_e( xi );
            sin_term = sin_xi_result.val / xi;
            cos_term = cos_xi_result.val;
            sin_term_err = sin_xi_result.err / Math.abs( xi );
            cos_term_err = cos_xi_result.err;
        }
        if ( eta < GSL_ROOT5_DBL_EPSILON )
        {
            coth_term = 1.0 + eta * eta / 3.0 * (1.0 - eta * eta / 15.0);
            sinh_term = 1.0 - eta * eta / 6.0 * (1.0 - 7.0 / 60.0 * eta * eta);
        }
        else
        {
            coth_term = eta / Math.tanh( eta );
            sinh_term = eta / Math.sinh( eta );
        }
        t1 = Math.sqrt( lsqp1 ) * eta;
        pre_val = sinh_term / t1;
        pre_err = 2.0 * GSL_DBL_EPSILON * Math.abs( pre_val );
        term1 = sin_term * coth_term;
        term2 = cos_term;
        r.val = pre_val * (term1 - term2);
        r.err = pre_err * Math.abs( term1 - term2 );
        r.err = r.err + pre_val * (sin_term_err * coth_term + cos_term_err);
        r.err = r.err + pre_val * Math.abs( term1 - term2) * (Math.abs( eta ) + 1.0) * GSL_DBL_EPSILON;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }

} // gsl_sf_legendre_H3d_1_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_H3d_e( ell, lambda, eta )
{
    var abs_lam  = 0.0;
    var lsq      = 0.0;
    var xi       = 0.0;
    var cosh_eta = 0.0;
    var r        = { val: 0.0, err: 0.0 }; // Result;

    abs_lam  = Math.abs( lambda );
    lsq      = abs_lam * abs_lam;
    xi       = abs_lam * eta;
    cosh_eta = Math.cosh( eta );

    if ( eta < 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( eta > GSL_LOG_DBL_MAX )
    {
        // cosh(eta) is too big.
        throw "SF.OverflowException";
    }
    else if ( ell == 0 )
    {
        return gsl_sf_legendre_H3d_0_e( lambda, eta );
    }
    else if ( ell == 1 )
    {
        return gsl_sf_legendre_H3d_1_e( lambda, eta );
    }
    else if ( eta == 0.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else if ( xi < 1.0 )
    {
        return legendre_H3d_series( ell, lambda, eta );
    }
    else if ( ((ell) * (ell) + lsq) / Math.sqrt( 1.0 + lsq ) / (cosh_eta * cosh_eta) < 5.0 * GSL_ROOT3_DBL_EPSILON )
    {
        // Large argument.
        //
        var lm = { Double: 0.0 };
        var P  = { val: 0.0, err: 0.0 }; // Result;

        gsl_sf_conicalP_large_x_e( -(ell) - 0.5, lambda, cosh_eta, P, lm );
        if ( P.val == 0.0 )
        {
            r.val = 0.0;
            r.err = 0.0;
            return r;
        }
        else
        {
            let lnN       = 0.0;
            let ln_abslam = 0.0;
            let lnpre_val = 0.0;
            let lnpre_err = 0.0;
            let lnsh      = { val: 0.0, err: 0.0 }; // Result;

            lnsh = gsl_sf_lnsinh_e( eta );
            lnN = legendre_H3d_lnnorm( ell, lambda );
            ln_abslam = Math.log( abs_lam );
            lnpre_val = 0.5 * (M_LNPI + lnN - M_LN2 - lnsh.val) - ln_abslam;
            lnpre_err = lnsh.err;
            lnpre_err = lnpre_err + 2.0 * GSL_DBL_EPSILON * (0.5 * (M_LNPI + M_LN2 + Math.abs( lnN )) + Math.abs( ln_abslam ));
            lnpre_err = lnpre_err + 2.0 * GSL_DBL_EPSILON * Math.abs( lnpre_val );
            r = gsl_sf_exp_mult_err_e( lnpre_val + lm.Double, lnpre_err, P.val, P.err );
            return r;
        }
    }
    else if ( abs_lam > 1000.0 * (ell) * (ell) )
    {
        // Large degree.
        //
        let lm = { Double: 0.0 };
        let r  = { val: 0.0, err: 0.0 }; // Result;
        let P  = { val: 0.0, err: 0.0 }; // Result;

        gsl_sf_conicalP_xgt1_neg_mu_largetau_e( (ell) + 0.5, lambda, cosh_eta, eta, P, lm );
        if ( P.val == 0.0 )
        {
            r.val = 0.0;
            r.err = 0.0;
            return r;
        }
        else
        {
            let lnN       = 0.0;
            let ln_abslam = 0.0;
            let lnpre_val = 0.0;
            let lnpre_err = 0.0;
            let lnsh      = { val: 0.0, err: 0.0 }; // Result;

            lnsh = gsl_sf_lnsinh_e( eta );
            lnN = legendre_H3d_lnnorm( ell, lambda );
            ln_abslam = Math.log( abs_lam );
            lnpre_val = 0.5 * (M_LNPI + lnN - M_LN2 - lnsh.val) - ln_abslam;
            lnpre_err = lnsh.err;
            lnpre_err = lnpre_err + GSL_DBL_EPSILON * (0.5 * (M_LNPI + M_LN2 + Math.abs( lnN )) + Math.abs( ln_abslam ));
            lnpre_err = lnpre_err + 2.0 * GSL_DBL_EPSILON * Math.abs( lnpre_val );
            r = gsl_sf_exp_mult_err_e( lnpre_val + lm.Double, lnpre_err, P.val, P.err );
            return r;
        }
    }
    else
    {
        // Backward recurrence.
        //
        let coth_eta      = 0.0;
        // let coth_eta_mult = 0.0;
        let coth_err_mult = 0.0;
        let Hlm1          = 0.0;
        let Hl            = 0.0;
        let Hlp1          = 0.0;
        let root_term_0   = 0.0;
        let root_term_1   = 0.0;
        let rH            = { val: 0.0, err: 0.0 }; // Result;
        let H0            = { val: 0.0, err: 0.0 }; // Result;
        let H1            = { val: 0.0, err: 0.0 }; // Result;

        coth_eta = 1.0 / Math.tanh( eta );
        coth_err_mult = Math.abs( eta ) + 1.0;
        rH = legendre_H3d_CF1_ser( ell, lambda, coth_eta );
        Hl   = GSL_SQRT_DBL_MIN;
        Hlp1 = rH.val * Hl;
        for ( let lp = ell; lp >= 1; lp-- )
        {
            root_term_0 = gsl_sf_hypot( lambda, (lp) );
            root_term_1 = gsl_sf_hypot( lambda, (lp + 1) );
            Hlm1 = ((2 * lp + 1) * coth_eta * Hl - root_term_1 * Hlp1) / root_term_0;
            Hlp1 = Hl;
            Hl   = Hlm1;
        }

        if ( Math.abs( Hl) > Math.abs( Hlp1) )
        {
            H0 = gsl_sf_legendre_H3d_0_e( lambda, eta );
            r.val = GSL_SQRT_DBL_MIN / Hl * H0.val;
            r.err = GSL_SQRT_DBL_MIN / Math.abs( Hl) * H0.err;
            r.err = r.err + Math.abs( rH.err / rH.val ) * (ell + 1) * coth_err_mult * Math.abs( r.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            return r;
        }
        else
        {
            H1 = gsl_sf_legendre_H3d_1_e( lambda, eta );
            r.val = GSL_SQRT_DBL_MIN / Hlp1 * H1.val;
            r.err = GSL_SQRT_DBL_MIN / Math.abs( Hlp1 ) * H1.err;
            r.err = r.err + Math.abs( rH.err / rH.val ) * (ell + 1) * coth_err_mult * Math.abs( r.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            return r;
        }
    }

} // gsl_sf_legendre_H3d_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_H3d_array( lmax, lambda, eta, result_array )
{

    if ( eta < 0.0 || lmax < 0 )
    {
        for ( let ell = 0; ell <= lmax; ell++ )
        {
            result_array[ell] = 0.0;
        }
        throw "SF.DomainException";
    }
    else if ( eta > GSL_LOG_DBL_MAX )
    {
        // cosh(eta) is too big.
        for ( let ell = 0; ell <= lmax; ell++ )
        {
            result_array[ell] = 0.0;
        }
        throw "SF.OverflowException";
    }
    else if ( lmax == 0 )
    {
        let H0 = { val: 0.0, err: 0.0 };
        H0 = gsl_sf_legendre_H3d_e( 0, lambda, eta );
        result_array[0] = H0.val;
        return;
    }
    else
    {
        // Not the most efficient method. But what the hell... it's simple.
        let r_Hlp1 = { val: 0.0, err: 0.0 };
        let r_Hl   = { val: 0.0, err: 0.0 };
        r_Hlp1 = gsl_sf_legendre_H3d_e( lmax,     lambda, eta );
        r_Hl   = gsl_sf_legendre_H3d_e( lmax - 1, lambda, eta );

        let coth_eta = 1.0 / Math.tanh( eta );
        let Hlp1 = r_Hlp1.val;
        let Hl   = r_Hl.val;
        let Hlm1 = 0.0;
 
        result_array[lmax]   = Hlp1;
        result_array[lmax-1] = Hl;

        for ( let ell = lmax - 1; ell > 0; ell-- )
        {
            let root_term_0 = gsl_sf_hypot( lambda, ell );
            let root_term_1 = gsl_sf_hypot( lambda, ell + 1.0 );
            Hlm1 = ((2.0 * ell + 1.0) * coth_eta*Hl - root_term_1 * Hlp1) / root_term_0;
            result_array[ell-1] = Hlm1;
            Hlp1 = Hl;
            Hl   = Hlm1;
        }

        return;
    }

} // gsl_sf_legendre_H3d_array

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_legendre_H3d_0( lambda, eta )
{ // gsl_sf_legendre_H3d_0
    return EVAL_RESULT_DD( gsl_sf_legendre_H3d_0_e, { x: lambda, y: eta }, "gsl_sf_legendre_H3d_0" );
}// gsl_sf_legendre_H3d_0;

export function gsl_sf_legendre_H3d_1( lambda, eta )
{ // gsl_sf_legendre_H3d_1
    return EVAL_RESULT_DD( gsl_sf_legendre_H3d_1_e, { x: lambda, y: eta }, "gsl_sf_legendre_H3d_1" );
}// gsl_sf_legendre_H3d_1;

export function gsl_sf_legendre_H3d( l, lambda, eta )
{ // gsl_sf_legendre_H3d
    return EVAL_RESULT_IDD( gsl_sf_legendre_H3d_e, { n: l, x: lambda, y: eta }, "gsl_sf_legendre_H3d" );
}// gsl_sf_legendre_H3d;

// -- ----------------------------------------------------------------------------
// EOF SF-LegendreH3d.mjs