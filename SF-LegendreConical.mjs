// SF-LegendreConical.mjs -- implementation of Legendre Conical functions
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

import { GSL_MODE_DEFAULT }      from "./SF-Mode.mjs";
import { GSL_DBL_MAX }           from "./SF-Machine.mjs";
import { GSL_DBL_MIN }           from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_ROOT4_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }      from "./SF-Machine.mjs";
import { GSL_SIGN }              from "./SF-Math.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { M_SQRT2 }               from "./SF-Math.mjs";
import { M_SQRTPI }              from "./SF-Math.mjs";
import { M_LNPI }                from "./SF-Math.mjs";
import { M_LN2 }                 from "./SF-Math.mjs";
import { gsl_sf_bessel_I0_scaled_e }   from "./SF-BesselI0.mjs";
import { gsl_sf_bessel_I1_scaled_e }   from "./SF-BesselI1.mjs";
import { gsl_sf_bessel_J0_e }    from "./SF-BesselJ0.mjs";
import { gsl_sf_bessel_J1_e }    from "./SF-BesselJ1.mjs";
import { gsl_sf_bessel_Jnu_e }   from "./SF-BesselJnu.mjs";
import { gsl_sf_bessel_Inu_e }   from "./SF-BesselInu.mjs";
import { gsl_sf_cos_e }          from "./SF-Trigonometric.mjs";
import { gsl_sf_cos_err_e }      from "./SF-Trigonometric.mjs";
import { gsl_sf_sin_err_e }      from "./SF-Trigonometric.mjs";
import { gsl_sf_lnsinh_e }       from "./SF-Trigonometric.mjs";
import { gsl_sf_ellint_Kcomp_e }  from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_hyperg_2F1_conj_e } from "./SF-Hypergeometric2F1.mjs";
import { gsl_sf_lngamma_e }         from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_complex_e } from "./SF-Gamma.mjs";
import { gsl_sf_exp_err_e }      from "./SF-Exponential.mjs";
import { gsl_sf_exp_mult_err_e } from "./SF-Exponential.mjs";
import { gsl_sf_exp_mult_e }     from "./SF-Exponential.mjs";
import { gsl_sf_ellint_Ecomp_e } from "./SF-EllipticIntegrals.mjs"

import { EVAL_RESULT_DD }  from "./SF-Evaluate.mjs";
import { EVAL_RESULT_IDD } from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const Root_2OverPi = 0.797884560802865355879892;
const locEPS       = 1000.0 * GSL_DBL_EPSILON;

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

const RECURSE_LARGE = 1.0e-5 * GSL_DBL_MAX;
const RECURSE_SMALL = 1.0e+5 * GSL_DBL_MIN;

// Continued fraction for f_{ell+1}/f_ell
// f_ell = P^{-mu-ell}_{-1/2 + I tau}(x),  x < 1.0
//
// Uses standard CF method from Temme's book.
//
function conicalP_negmu_xlt1_CF1( mu, ell, tau, x )
{
    const RECUR_BIG = GSL_SQRT_DBL_MAX;
    const maxiter   = 5000;
    var n         = 0;
    var xi        = 0.0;
    var Anm2      = 0.0;
    var Bnm2      = 0.0;
    var Anm1      = 0.0;
    var Bnm1      = 0.0;
    var a1        = 0.0;
    var b1        = 0.0;
    var An        = 0.0;
    var Bn        = 0.0;
    var ans       = 0.0;
    var bns       = 0.0;
    var fn        = 0.0;
    var old_fn    = 0.0;
    var del       = 0.0;
    var r         = { val: 0.0, err: 0.0 }; // Result;

    n = 1;
    xi = x / (Math.sqrt( 1.0 - x ) * Math.sqrt( 1.0 + x ));
    Anm2 = 1.0;
    Bnm2 = 0.0;
    Anm1 = 0.0;
    Bnm1 = 1.0;
    a1 = 1.0;
    b1 = 2.0 * (mu + (ell) + 1.0) * xi;
    An = b1 * Anm1 + a1 * Anm2;
    Bn = b1 * Bnm1 + a1 * Bnm2;
    fn = An / Bn;

    while ( n < maxiter )
    {
        n = n + 1;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        ans  = tau * tau + (mu - 0.5 + (ell + n)) * (mu - 0.5 + (ell + n));
        bns  = 2.0 * ((ell) + mu + (n)) * xi;
        An   = bns * Anm1 + ans * Anm2;
        Bn   = bns * Bnm1 + ans * Bnm2;
    
        if ( Math.abs( An ) > RECUR_BIG || Math.abs( Bn ) > RECUR_BIG )
        {
            An   = An   / RECUR_BIG;
            Bn   = Bn   / RECUR_BIG;
            Anm1 = Anm1 / RECUR_BIG;
            Bnm1 = Bnm1 / RECUR_BIG;
            Anm2 = Anm2 / RECUR_BIG;
            Bnm2 = Bnm2 / RECUR_BIG;
        }
    
        old_fn = fn;
        fn     = An / Bn;
        del    = old_fn / fn;
    
        if ( Math.abs( del - 1.0 ) < 2.0 * GSL_DBL_EPSILON ) break;
    }

    r.val = fn;
    r.err = 4.0 * GSL_DBL_EPSILON * (Math.sqrt( n ) + 1.0) * Math.abs( fn );

    if ( n >= maxiter )
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // conicalP_negmu_xlt1_CF1

// ----------------------------------------------------------------------------

// Continued fraction for f_{ell+1}/f_ell
// f_ell = P^{-mu-ell}_{-1/2 + I tau}(x),  x >= 1.0
//
// Uses Gautschi (Euler) equivalent series.
//
function conicalP_negmu_xgt1_CF1( mu, ell, tau, x )
{
    const maxk  = 20000;
    var gamma = 0.0;
    var pre   = 0.0;
    var tk    = 0.0;
    var sum   = 0.0;
    var rhok  = 0.0;
    var tlk   = 0.0;
    var l1k   = 0.0;
    var ak    = 0.0;
    var k     = 0;
    var r     = { val: 0.0, err: 0.0 }; // Result;

    gamma = 1.0 - 1.0 / (x * x);
    pre   = Math.sqrt( x - 1.0 ) * Math.sqrt( x + 1.0 ) / (x * (2.0 * ((ell) + mu + 1.0)));
    tk    = 1.0;
    sum   = 1.0;
    rhok  = 0.0;
 
    k = 1; 
    while ( k < maxk )
    {
        tlk  = 2.0 * ((ell) + mu + (k));
        l1k  = ((ell) + mu - 0.5 + 1.0 + (k));
        ak   = -(tau * tau + l1k * l1k) / (tlk * (tlk + 2.0)) * gamma;
        rhok = -ak * (1.0 + rhok) / (1.0 + ak * (1.0 + rhok));
        tk   = tk * rhok;
        sum  = sum + tk;
        if ( Math.abs( tk / sum ) < GSL_DBL_EPSILON ) break;
        k = k + 1;
    }
    
    r.val = pre * sum;
    r.err = Math.abs( pre * tk );
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.sqrt( k ) + 1.0) * Math.abs( pre * sum );
    
    if ( k >= maxk )
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // conicalP_negmu_xgt1_CF1

// ----------------------------------------------------------------------------

// Implementation of large negative mu asymptotic
// [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
//
//
function olver_U1( beta2, p )
{
    return (p - 1.0) / (24.0 * (1.0 + beta2)) * (3.0 + beta2 * (2.0 + 5.0 * p * (1.0 + p)));
}

// ----------------------------------------------------------------------------

function olver_U2( beta2, p )
{
    var beta4 = beta2 * beta2;
    var p2    = p*p;
    var poly1 =  4.0 * beta4 + 84.0 * beta2 - 63.0;
    var poly2 = 16.0 * beta4 + 90.0 * beta2 - 81.0;
    var poly3 = beta2 * p2 * (97.0 * beta2 - 432.0 + 77.0 * p * (beta2 - 6.0) - 385.0 * beta2 * p2 * (1.0 + p));
    return (1.0 - p) / (1152.0 * (1.0 + beta2)) * (poly1 + poly2 + poly3);
}

// ----------------------------------------------------------------------------

const U3c1 = [   -1307.0,   -1647.0,    3375.0,    3675.0 ];
const U3c2 = [   29366.0,   35835.0, -252360.0, -272630.0, 276810.0,  290499.0 ];
const U3c3 = [  -29748.0,   -8840.0, 1725295.0, 1767025.0, -7313470.0, -754778.0, 6309875.0, 6480045.0 ];
const U3c4 = [    2696.0,    -16740.0,   -524250.0,  -183975.0,
                             14670540.0,  14172939.0, -48206730.0, -48461985.0,
                             36756720.0,  37182145.0 ];
const U3c5 = [       9136.0,      22480.0,     12760.0,
                                 -252480.0,   -4662165.0,   -1705341.0,
                                92370135.0,   86244015.0, -263678415.0,
                              -260275015.0, 185910725.0,  185910725.0 ];

// ----------------------------------------------------------------------------

// Large negative mu asymptotic
// P^{-mu}_{-1/2 + I tau}, mu -> Inf
// |x| < 1
//
// [Dunster, Proc. Roy. Soc. Edinburgh 119A, 311 (1991), p. 326]
//
function gsl_sf_conicalP_xlt1_large_neg_mu_e( mu, tau, x, result, ln_multiplier )
{
    var beta  = tau / mu;
    var beta2 = beta * beta;
    var S     = beta * Math.acos( (1.0 - beta2) / (1.0 + beta2) );
    var p     = x / Math.ssqrt( beta2 * (1.0 - x * x) + 1.0 );
    var lg_mup1 = { val: 0.0, err: 0.0 };
    lg_mup1 = gsl_sf_lngamma_e( mu + 1.0 );
    var ln_pre_1 =  0.5 * mu * (S - Math.log( 1.0 + beta2 ) + Math.log( (1.0 - p) / (1.0 + p) )) - lg_mup1.val;
    var ln_pre_2 = -0.25 * Math.log( 1.0 + beta2 * (1.0 - x) );
    var ln_pre_3 = -tau * Math.atan( p * beta );
    var ln_pre = ln_pre_1 + ln_pre_2 + ln_pre_3;
    var sum   = 1.0 - olver_U1( beta2, p ) / mu + olver_U2( beta2, p ) / (mu * mu);

    if ( sum == 0.0 )
    {
        result.val = 0.0;
        result.err = 0.0;
        ln_multiplier.Double = 0.0;
        return;
    }
    else
    {
        try
        {
            let r = { val: 0.0, err: 0.0 };
            r = gsl_sf_exp_mult_e( ln_pre, sum );
            result.val = r.val;
            result.err = r.err;
            ln_multiplier.Double = 0.0;
        }
        catch ( e )
        {
            result.val = sum;
            result.err = 2.0 * GSL_DBL_EPSILON * Math.abs(sum);
            ln_multiplier.Double = ln_pre;
        }
        return;
    }
}

// ----------------------------------------------------------------------------

// Implementation of large tau asymptotic
//
// A_n^{-mu}, B_n^{-mu}  [Olver, p.465, 469]
//
//
function olver_B0_xi( mu, xi )
{
    return (1.0 - 4.0 * mu * mu) / (8.0 * xi) * (1.0 / Math.tanh( xi ) - 1.0 / xi);
} // olver_B0_xi

// ----------------------------------------------------------------------------

function olver_A1_xi( mu, xi, x )
{
    var B   = 0.0;
    var y   = 0.0;
    var s   = 0.0;
    var psi = 0.0;

    B = olver_B0_xi( mu, xi );
    if ( Math.abs( x - 1.0 ) < GSL_ROOT4_DBL_EPSILON )
    {
        y = x - 1.0;
        s = -1.0 / 3.0 + y * (2.0 / 15.0 - y * (61.0 / 945.0 - 452.0 / 14175.0 * y));
        psi = (4.0 * mu * mu - 1.0) / 16.0 * s;
    }
    else
    {
        psi = (4.0 * mu * mu - 1.0) / 16.0 * (1.0 / (x * x - 1.0) - 1.0 / (xi * xi));
    }
    return 0.5 * xi * xi * B * B + (mu + 0.5) * B - psi + mu / 6.0 * (0.25 - mu * mu);

} // olver_A1_xi

// ----------------------------------------------------------------------------

function olver_B0_th( mu, theta )
{
    return -(1.0 - 4.0 * mu * mu) / (8.0 * theta) * (1.0 / Math.tan( theta ) - 1.0 / theta);
}

// ----------------------------------------------------------------------------

function olver_A1_th( mu, theta, x )
{
    var B = olver_B0_th( mu, theta );
    var psi = 0.0;
    if ( Math.abs( x - 1.0 ) < GSL_ROOT4_DBL_EPSILON )
    {
        let y = 1.0 - x;
        let s = -1.0 / 3.0 + y * (2.0 / 15.0 - y * (61.0 / 945.0 - 452.0 / 14175.0 * y));
        psi = (4.0 * mu * mu - 1.0) / 16.0 * s;
    }
    else
    {
        psi = (4.0 * mu * mu - 1.0) / 16.0 * (1.0 / (x * x - 1.0) + 1.0 / (theta * theta));
    }
    return -0.5 * theta * theta * B * B + (mu + 0.5) * B - psi + mu / 6.0 * (0.25 - mu * mu);
}

// ----------------------------------------------------------------------------

// Large tau uniform asymptotics
// P^{-mu}_{-1/2 + I tau}
// 1 < x
// tau -> Inf 
// [Olver, p. 469]
//
export function gsl_sf_conicalP_xgt1_neg_mu_largetau_e( mu, tau, x, acosh_x, r, ln_multiplier )
{
    var xi        = acosh_x;
    var ln_xi_pre = 0.0;
    var ln_pre    = 0.0;
    var sumA      = 0.0;
    var sumB      = 0.0;
    var sum       = 0.0;
    var arg       = 0.0;
    var J_mum1    = 0.0;
    var J_mup1    = { val: 0.0, err: 0.0 }; // Result;
    var J_mu      = { val: 0.0, err: 0.0 }; // Result;
    var lnshxi    = { val: 0.0, err: 0.0 }; // Result;

    if ( xi < GSL_ROOT4_DBL_EPSILON )
    {
        ln_xi_pre = -xi * xi / 6.0; // log(1.0 - xi*xi/6.0)
    }
    else
    {
        lnshxi = gsl_sf_lnsinh_e( xi );
        ln_xi_pre = Math.log( xi ) - lnshxi.val; // log(xi/sinh(xi)
    }

    ln_pre = 0.5 * ln_xi_pre - mu * Math.log( tau );

    arg = tau * xi;

    J_mup1 = gsl_sf_bessel_Jnu_e( mu + 1.0,   arg );
    J_mu   = gsl_sf_bessel_Jnu_e( mu,         arg );
    J_mum1 = -J_mup1.val + 2.0 * mu / arg * J_mu.val; // careful of mu < 1

    sumA = 1.0 - olver_A1_xi( -mu, xi, x ) / (tau * tau);
    sumB = olver_B0_xi( -mu, xi );
    sum  = J_mu.val * sumA - xi / tau * J_mum1 * sumB;

    if ( sum == 0.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        ln_multiplier.Double = 0.0;
    }
    else
    {
        try
        {
            let r1 = { val: 0.0, err: 0.0 };
            r1 = gsl_sf_exp_mult_e( ln_pre, sum );
            r.val = sum;
            r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( sum );
            ln_multiplier.Double = ln_pre;
        }
        catch ( e )
        {
            ln_multiplier.Double = 0.0;
        }
    }

} // gsl_sf_conicalP_xgt1_neg_mu_largetau_e

// ----------------------------------------------------------------------------

// Large tau uniform asymptotics
// P^{-mu}_{-1/2 + I tau}
// -1 < x < 1
// tau -> Inf 
// [Olver, p. 473]
//
function gsl_sf_conicalP_xlt1_neg_mu_largetau_e( mu, tau, x, acos_x, result, ln_multiplier )
{
    var theta     = acos_x;
    var ln_th_pre = 0.0;
    var ln_pre    = 0.0;
    var sumA      = 0.0;
    var sumB      = 0.0;
    var sum       = 0.0;
    var sumerr    = 0.0;
    var arg       = 0.0;
    var I_mup1    = { val: 0.0, err: 0.0 };
    var I_mu      = { val: 0.0, err: 0.0 };
    var I_mum1    = 0.0;

    if ( theta < GSL_ROOT4_DBL_EPSILON )
    {
        ln_th_pre = theta * theta / 6.0;   // log(1.0 + theta*theta/6.0)
    }
    else
    {
        ln_th_pre = Math.log( theta / Math.sin( theta ) );
    }

    ln_pre = 0.5 * ln_th_pre - mu * Math.log(tau);

    arg = tau * theta;
    I_mup1 = gsl_sf_bessel_Inu_e( mu + 1.0,   arg );
    I_mu   = gsl_sf_bessel_Inu_e( mu,         arg );
    I_mum1 = I_mup1.val + 2.0 * mu / arg * I_mu.val; // careful of mu < 1

    sumA = 1.0 - olver_A1_th( -mu, theta, x ) / ( tau * tau );
    sumB = olver_B0_th( -mu, theta );
    sum  = I_mu.val * sumA - theta / tau * I_mum1 * sumB;
    sumerr  = Math.abs( I_mu.err * sumA );
    sumerr += Math.abs( I_mup1.err * theta / tau * sumB );
    sumerr += Math.abs( I_mu.err   * theta / tau * sumB * 2.0 * mu / arg );

    if(sum == 0.0)
    {
        result.val = 0.0;
        result.err = 0.0;
        ln_multiplier.Double = 0.0;
        return;
    }
    else
    {
        try
        {
            let r = { val: 0.0, err: 0.0 };
            r = gsl_sf_exp_mult_e( ln_pre, sum );
            result.val = r.val;
            result.err = r.err;
            ln_multiplier.Double = 0.0;
        }
        catch ( e )
        {
            result.val  = sum;
            result.err  = sumerr;
            result.err += GSL_DBL_EPSILON * Math.abs( sum );
            ln_multiplier.Double = ln_pre;
        }
        return;
    }
}

// ----------------------------------------------------------------------------

// Hypergeometric function which appears in the
// large x expansion below:
//
//   2F1(1/4 - mu/2 - I tau/2, 3/4 - mu/2 - I tau/2, 1 - I tau, y)
//
// Note that for the usage below y = 1/x^2;
//
function conicalP_hyperg_large_x( mu, tau, y, reF, imF )
{
    const kmax          = 1000;
    var re_a          = 0.25 - 0.5 * mu;
    var re_b          = 0.75 - 0.5 * mu;
    var re_c          = 1.0;
    var im_a          = -0.5 * tau;
    var im_b          = -0.5 * tau;
    var im_c          = -tau;
    var re_sum        = 1.0;
    var im_sum        = 0.0;
    var re_term       = 1.0;
    var im_term       = 0.0;
    var re_ak         = 0.0;
    var re_bk         = 0.0;
    var re_ck         = 0.0;
    var im_ak         = 0.0;
    var im_bk         = 0.0;
    var im_ck         = 0.0;
    var den           = 0.0;
    var re_multiplier = 0.0;
    var im_multiplier = 0.0;
    var re_tmp        = 0.0;
    var im_tmp        = 0.0;
    var asum          = 0.0;
    var k             = 0;

    k = 1;
    while ( k <= kmax )
    {
        re_ak = re_a + k - 1.0;
        re_bk = re_b + k - 1.0;
        re_ck = re_c + k - 1.0;
        im_ak = im_a;
        im_bk = im_b;
        im_ck = im_c;
        den   = re_ck * re_ck + im_ck * im_ck;
        re_multiplier = ((re_ak * re_bk - im_ak * im_bk) * re_ck + im_ck * (im_ak * re_bk + re_ak * im_bk)) / den;
        im_multiplier = ((im_ak * re_bk + re_ak * im_bk) * re_ck - im_ck * (re_ak * re_bk - im_ak * im_bk)) / den;
        re_tmp = re_multiplier * re_term - im_multiplier * im_term;
        im_tmp = im_multiplier * re_term + re_multiplier * im_term;
        asum = Math.abs( re_sum ) + Math.abs( im_sum );
        re_term = y / k * re_tmp;
        im_term = y / k * im_tmp;
        if ( Math.abs( re_term / asum ) < GSL_DBL_EPSILON && Math.abs( im_term / asum ) < GSL_DBL_EPSILON ) break;
        re_sum = re_sum + re_term;
        im_sum = im_sum + im_term;
        k = k + 1;
    }

    reF.Double = re_sum;
    imF.Double = im_sum;

    if ( k >= kmax )
    {
        throw "SF.MaxIterationsException";
    }

} // conicalP_hyperg_large_x

// ----------------------------------------------------------------------------

// P^{mu}_{-1/2 + I tau}
// x->Inf
//
export function gsl_sf_conicalP_large_x_e( mu, tau, x, r, ln_multiplier )
{
    var y     = 0.0;
    var reF   = { Double: 0.0 };
    var imF   = { Double: 0.0 };
    var angle = 0.0;
    var lnx   = 0.0;
    var lnxp1 = 0.0;
    var lnxm1 = 0.0;
    var lnpre_const = 0.0;
    var lnpre_comm  = 0.0;
    var lnpre_err   = 0.0;
    var lgr_num     = { val: 0.0, err: 0.0 }; // Result;
    var lgth_num    = { val: 0.0, err: 0.0 }; // Result;
    var lgr_den     = { val: 0.0, err: 0.0 }; // Result;
    var lgth_den    = { val: 0.0, err: 0.0 }; // Result;
    var cos_result  = { val: 0.0, err: 0.0 }; // Result;

    // 2F1 term
    //
    if ( x < 0.5 * GSL_SQRT_DBL_MAX )
    {
        y = 1.0 / (x * x);
    }
    else
    {
        y = 0.0;
    }
    conicalP_hyperg_large_x( mu, tau, y, reF, imF );

    // f = Gamma(+i tau)/Gamma(1/2 - mu + i tau)
    // FIXME: shift so it's better for tau-> 0
    //
    gsl_sf_lngamma_complex_e( 0.0, tau, lgr_num, lgth_num );
    gsl_sf_lngamma_complex_e( 0.5 - mu, tau, lgr_den, lgth_den );

    angle = lgth_num.val - lgth_den.val + Math.atan2( imF.Double, reF.Double );

    lnx   = Math.log( x );
    lnxp1 = Math.log( x + 1.0 );
    lnxm1 = Math.log( x - 1.0 );
    lnpre_const = 0.5 * M_LN2 - 0.5 * M_LNPI;
    lnpre_comm  = (mu - 0.5) * lnx - 0.5 * mu * (lnxp1 + lnxm1);
    lnpre_err   =   GSL_DBL_EPSILON * (0.5 * M_LN2 + 0.5 * M_LNPI)
                + GSL_DBL_EPSILON * Math.abs( (mu - 0.5) * lnx )
                + GSL_DBL_EPSILON * Math.abs( 0.5 * mu ) * (Math.abs( lnxp1 ) + Math.abs( lnxm1 ));

    //  result = pre*|F|*|f| * cos(angle - tau * (log(x)+M_LN2))
    //
    cos_result = gsl_sf_cos_e( angle + tau * (Math.log( x ) + M_LN2) );
    if ( cos_result.val == 0.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        return;
    }
    else
    {
        let lnFf_val  = 0.5 * Math.log( reF.Double * reF.Double + imF.Double * imF.Double ) + lgr_num.val - lgr_den.val;
        let lnFf_err  = lgr_num.err + lgr_den.err + GSL_DBL_EPSILON * Math.abs( lnFf_val );
        let lnnoc_val = lnpre_const + lnpre_comm + lnFf_val;
        let lnnoc_err = lnpre_err + lnFf_err + GSL_DBL_EPSILON * Math.abs( lnnoc_val );
        let rr = { val: 0.0, err: 0.0 };

        try
        {
            rr = gsl_sf_exp_mult_err_e( lnnoc_val, lnnoc_err, cos_result.val, cos_result.err );
            r.val = rr.val;
            r.err = rr.err;
            ln_multiplier.Double = 0.0;
        }
        catch ( e )
        {
            r.val = cos_result.val;
            r.err = cos_result.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            ln_multiplier.Double = lnnoc_val;
        }
        return;
    }

} // gsl_sf_conicalP_large_x_e

// ----------------------------------------------------------------------------

// P^{mu}_{-1/2 + I tau}  first hypergeometric representation
// -1 < x < 1
// This is more effective for |x| small, however it will work w/o
// reservation for any x < 0 because everything is positive
// definite in that case.
//
// [Kolbig,   (3)] (note typo in args of gamma functions)
// [Bateman, (22)] (correct form)
//
function conicalP_xlt1_hyperg_A( mu, tau, x )
{
    var x2      = 0.0;
    var err_amp = 0.0;
    var pre_val = 0.0;
    var pre_err = 0.0;
    var t1_val  = 0.0;
    var t1_err  = 0.0;
    var t2_val  = 0.0;
    var t2_err  = 0.0;
    var ln_g1   = { val: 0.0, err: 0.0 }; // Result;
    var ln_g2   = { val: 0.0, err: 0.0 }; // Result;
    var arg_g1  = { val: 0.0, err: 0.0 }; // Result;
    var arg_g2  = { val: 0.0, err: 0.0 }; // Result;
    var F1      = { val: 0.0, err: 0.0 }; // Result;
    var F2      = { val: 0.0, err: 0.0 }; // Result;
    var pre1    = { val: 0.0, err: 0.0 }; // Result;
    var pre2    = { val: 0.0, err: 0.0 }; // Result;
    var r       = { val: 0.0, err: 0.0 }; // Result;

    x2      = x * x;
    err_amp = 1.0 + 1.0 / (GSL_DBL_EPSILON + Math.abs( 1.0 - Math.abs( x ) ));
    pre_val = M_SQRTPI / Math.pow( (0.5 * Math.sqrt( 1.0 - x2 )), mu );
    pre_err = err_amp * GSL_DBL_EPSILON * (Math.abs( mu ) + 1.0) * Math.abs( pre_val ) ;
    
    F1 = gsl_sf_hyperg_2F1_conj_e( 0.25 - 0.5 * mu, 0.5 * tau, 0.5, x2 );
    F2 = gsl_sf_hyperg_2F1_conj_e( 0.75 - 0.5 * mu, 0.5 * tau, 1.5, x2 );
    
    gsl_sf_lngamma_complex_e( 0.75 - 0.5 * mu, -0.5 * tau, ln_g1, arg_g1 );
    gsl_sf_lngamma_complex_e( 0.25 - 0.5 * mu, -0.5 * tau, ln_g2, arg_g2 );
    
    pre1 = gsl_sf_exp_err_e( -2.0 * ln_g1.val, 2.0 * ln_g1.err );
    pre2 = gsl_sf_exp_err_e( -2.0 * ln_g2.val, 2.0 * ln_g2.err );
    pre2.val = -pre2.val * 2.0 * x;
    pre2.err =  pre2.err * 2.0 * Math.abs( x );
    pre2.err =  pre2.err + GSL_DBL_EPSILON * Math.abs( pre2.val );
    
    t1_val = pre1.val * F1.val;
    t1_err = Math.abs( pre1.val ) * F1.err + pre1.err * Math.abs( F1.val );
    t2_val = pre2.val * F2.val;
    t2_err = Math.abs( pre2.val ) * F2.err + pre2.err * Math.abs( F2.val );
    
    r.val = pre_val * (t1_val + t2_val);
    r.err = pre_val * (t1_err + t2_err);
    r.err = r.err + pre_err * Math.abs( t1_val + t2_val );
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
    
    return r;

} // conicalP_xlt1_hyperg_A

// ----------------------------------------------------------------------------

// V0, V1 from Kolbig, m = 0
//
function conicalP_0_V( u, f, tau, sgn, V0, V1 )
{
    // V0: in out { Double: 0.0 }
    // V1: in out { Double: 0.0 }

    var C = []; // ARRAY(0.. 7) OF LONG_FLOAT;
    var T = []; // ARRAY(0.. 7) OF LONG_FLOAT;
    var H = []; // ARRAY(0.. 7) OF LONG_FLOAT;
    var V = []; // ARRAY(0..11) OF LONG_FLOAT;

    T[0] = 1.0;
    H[0] = 1.0;
    V[0] = 1.0;
    for ( let i = 1; i <= 7; i++ )
    {
        T[i] = T[i-1] * u;
        H[i] = H[i-1] * (u * f);
    }
    for ( let i = 1; i <= 11; i++ )
    {
        V[i] = V[i-1] * tau;
    }

    C[0] = 1.0;
    C[1] = (H[1] - 1.0) / (8.0 * T[1]);
    C[2] = (9.0 * H[2] + 6.0 * H[1] - 15.0 - sgn * 8.0 * T[2]) / (128.0 * T[2]);
    C[3] = 5.0 * (15.0 * H[3] + 27.0 * H[2] + 21.0 * H[1] - 63.0 - sgn * T[2] * (16.0 * H[1] + 24.0)) / (1024.0 * T[3]);
    C[4] = 7.0*(525.0*H[4] + 1500.0*H[3] + 2430.0*H[2] + 1980.0*H[1] - 6435.0
              + 192.0*T[4] - sgn*T[2]*(720.0*H[2]+1600.0*H[1]+2160.0)
              ) / (32768.0*T[4]);
    C[5] = 21.0 * (2835.0 * H[5] + 11025.0 * H[4] + 24750.0 * H[3] + 38610.0 * H[2]
               + 32175.0 * H[1] - 109395.0 + T[4] * (1984.0 * H[1] + 4032.0)
               - sgn * T[2] * (4800.0 * H[3] + 15120.0 * H[2] + 26400.0 * H[1] + 34320.0)
               ) / (262144.0 * T[5]);
    C[6] = 11.0 * (218295.0 * H[6] + 1071630.0 * H[5] + 3009825.0 * H[4] + 6142500.0 * H[3]
               + 9398025.0 * H[2] + 7936110.0 * H[1] - 27776385.0
               + T[4] * (254016.0 * H[2] + 749952.0 * H[1] + 1100736.0)
               - sgn * T[2] * (441000.0 * H[4] + 1814400.0 * H[3] + 4127760.0 * H[2]
                         + 6552000.0 * H[1] + 8353800.0 + 31232.0 * T[4]
                         )
               ) / (4194304.0 * T[6]);

    V0.Double = C[0] + (-4.0 * C[3] / T[1] + C[4]) / V[4]
             + (-192.0 * C[5] / T[3] + 144.0 * C[6] / T[2]) / V[8]
             + sgn * (-C[2] / V[2]
                      + (-24.0 * C[4] / T[2] + 12.0 * C[5] / T[1] - C[6]) / V[6] 
                      + (-1920.0 * C[6] / T[4]) / V[10]
                      );
    V1.Double = C[1] / V[1] + (8.0 * (C[3] / T[2] - C[4] / T[1]) + C[5]) / V[5]
                  + (384.0 * C[5] / T[4] - 768.0 * C[6] / T[3]) / V[9]
                  + sgn * ((2.0 * C[2] / T[1] - C[3]) / V[3]
                           + (48.0 * C[4] / T[3] - 72.0 * C[5] / T[2] + 18.0 * C[6] / T[1]) / V[7]
                           + (3840.0 * C[6] / T[5]) / V[11]
                           );

} // conicalP_0_V

// ----------------------------------------------------------------------------

// V0, V1 from Kolbig, m = 1
//
function conicalP_1_V( t, f, tau, sgn, V0, V1 )
{
    // V0: in out { Double: 0.0 }
    // V1: in out { Double: 0.0 }

    var Cm1 = 0.0;
    var C   = []; // ARRAY(0.. 7) OF LONG_FLOAT = (OTHERS => 0.0);
    var S   = []; // ARRAY(0.. 7) OF LONG_FLOAT = (OTHERS => 0.0);
    var H   = []; // ARRAY(0.. 7) OF LONG_FLOAT = (OTHERS => 0.0);
    var V   = []; // ARRAY(0..11) OF LONG_FLOAT = (OTHERS => 0.0);

    S[0] = 1.0;
    H[0] = 1.0;
    V[0] = 1.0;
    for ( let i = 1; i <= 7; i++ )
    {
        S[i] = S[i-1] * t;
        H[i] = H[i-1] * (t * f);
    }
    for ( let i = 1; i <= 11;i++ )
    {
        V[i] = V[i-1] * tau;
    }

    Cm1  = -1.0;
    C[0] = 3.0 * (1.0 - H[1]) / (8.0 * S[1]);
    C[1] = (-15.0 * H[2] + 6.0 * H[1] + 9.0 + sgn * 8.0 * S[2]) / (128.0 * S[2]);
    C[2] = 3.0 * (-35.0 * H[3] - 15.0 * H[2] + 15.0 * H[1] + 35.0 + sgn * S[2] * (32.0 * H[1] + 8.0)) / (1024.0 * S[3]);
    C[3] = (-4725.0 * H[4] - 6300.0 * H[3] - 3150.0* H[2] + 3780.0 * H[1] + 10395.0
            -1216.0 * S[4] + sgn * S[2] * (6000.0 * H[2] + 5760.0 * H[1] +1680.0)) / (32768.0 * S[4]);
    C[4] = 7.0 * (-10395.0 * H[5] - 23625.0 * H[4] - 28350.0 * H[3] - 14850.0 * H[2]
         + 19305.0 * H[1] + 57915.0 -  S[4] * (6336.0 * H[1] + 6080.0)
         + sgn * S[2] * (16800.0 * H[3] + 30000.0 * H[2] + 25920.0 * H[1] + 7920.0)
         ) / (262144.0 * S[5]);
    C[5] = (-2837835.0 * H[6] - 9168390.0 * H[5] - 16372125.0 * H[4] - 18918900.0 * H[3]
         - 10135125.0 * H[2] + 13783770.0 * H[1] + 43648605.0
         - S[4] * (3044160.0 * H[2] + 5588352.0 * H[1] + 4213440.0)
         + sgn * S[2] * (5556600.0 * H[4] + 14817600.0 * H[3] + 20790000.0 * H[2]
         + 17297280.0 * H[1] + 5405400.0 + 323072.0 * S[4]
         )
         ) / (4194304.0 * S[6]);
    C[6] = 0.0;

    V0.Double = C[0] + (-4.0 * C[3] / S[1] + C[4]) / V[4]
              + (-192.0 * C[5] / S[3] + 144.0 * C[6] / S[2]) / V[8]
              + sgn * (-C[2] / V[2]
              + (-24.0 * C[4] / S[2] + 12.0 * C[5] / S[1] - C[6]) / V[6] 
              );
    V1.Double = C[1] / V[1] + (8.0 * (C[3] / S[2] - C[4] / S[1]) + C[5]) / V[5]
              + (384.0 * C[5] / S[4] - 768.0 * C[6] / S[3]) / V[9]
              + sgn * (Cm1 * V[1] + (2.0 * C[2] / S[1] - C[3]) / V[3]
              + (48.0 * C[4] / S[3] - 72.0 * C[5] / S[2] + 18.0 * C[6] / S[1]) / V[7]
              );

} // conicalP_1_V

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

// P^0_{-1/2 + I lambda}
//
export function gsl_sf_conicalP_0_e(lambda, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x <= -1.0 )
    {
        throw "SF.DomainException";
    }
    else if ( x == 1.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if ( lambda == 0.0 )
    {
        var th = 0.0;
        var s  = 0.0;
        var xi = 0.0;
        var c  = 0.0;
        var t  = 0.0;
        var K  = { val: 0.0, err: 0.0 }; // Result;

        if ( x < 1.0 )
        {
            th = Math.acos( x );
            s  = Math.sin( 0.5 * th );
            K = gsl_sf_ellint_Kcomp_e( s, GSL_MODE_DEFAULT );
            r.val = 2.0 / M_PI * K.val;
            r.err = 2.0 / M_PI * K.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            return r;
        }
        else
        {
            xi = Math.acosh( x );
            c  = Math.cosh( 0.5 * xi );
            t  = Math.tanh( 0.5 * xi );
            K = gsl_sf_ellint_Kcomp_e( t, GSL_MODE_DEFAULT );
            r.val = 2.0 / M_PI / c * K.val;
            r.err = 2.0 / M_PI / c * K.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            return r;
        }
    }
    else if ( (x <= 0.0 && lambda < 1000.0) || (x <  0.1 && lambda < 17.0) || (x <  0.2 && lambda < 5.0 ) )
    {
        return conicalP_xlt1_hyperg_A( 0.0, lambda, x );
    }
    else if ( (x <= 0.2 && lambda < 17.0) || (x <= 1.5 && lambda < 20.0) )
    {
        return gsl_sf_hyperg_2F1_conj_e( 0.5, lambda, 1.0, (1.0 - x) / 2.0 );
    }
    else if ( 1.5 < x && lambda < Math.max( x, 20.0 ) )
    {
        var P  = { val: 0.0, err: 0.0 }; // Result;
        var lm = { Double: 0.0 };

        gsl_sf_conicalP_large_x_e( 0.0, lambda, x, P, lm );
        r = gsl_sf_exp_mult_err_e( lm.Double, 2.0 * GSL_DBL_EPSILON * Math.abs( lm.Double ), P.val, P.err );
        return r;
    }
    else
    {
        let V0       = { Double: 0.0 };
        let V1       = { Double: 0.0 };
        let th       = 0.0;
        let sth      = 0.0;
        let bessterm = 0.0;
        let besserr  = 0.0;
        let arg1     = 0.0;
        let sqts     = 0.0;
        let sh       = 0.0;
        let xi       = 0.0;
        let pre_val  = 0.0;
        let pre_err  = 0.0;
        let I0       = { val: 0.0, err: 0.0 }; // Result;
        let I1       = { val: 0.0, err: 0.0 }; // Result;
        let J0       = { val: 0.0, err: 0.0 }; // Result;
        let J1       = { val: 0.0, err: 0.0 }; // Result;

        if ( x < 1.0 )
        {
            th  = Math.acos( x );
            sth = Math.sqrt( 1.0 - x * x );  // sin(th)
            I0  = gsl_sf_bessel_I0_scaled_e( th * lambda );
            I1  = gsl_sf_bessel_I1_scaled_e( th * lambda );
            conicalP_0_V( th, x / sth, lambda, -1.0, V0, V1 );
            bessterm = V0.Double * I0.val + V1.Double * I1.val;
            besserr  = Math.abs( V0.Double ) * I0.err + Math.abs( V1.Double ) * I1.err;
            arg1 = th * lambda;
            sqts = Math.sqrt( th / sth );
            r = gsl_sf_exp_mult_err_e( arg1, 4.0 * GSL_DBL_EPSILON * Math.abs( arg1 ), sqts * bessterm, sqts * besserr );
            return r;
        }
        else
        {
            sh = Math.sqrt( x - 1.0 ) * Math.sqrt( x + 1.0 );  // sinh(xi)
            xi = Math.log( x + sh );                    // xi = acosh(x)
            J0 = gsl_sf_bessel_J0_e( xi * lambda );
            J1 = gsl_sf_bessel_J1_e( xi * lambda );
            conicalP_0_V( xi, x / sh, lambda, 1.0, V0, V1 );
            bessterm = V0.Double * J0.val + V1.Double * J1.val;
            besserr  = Math.abs( V0.Double ) * J0.err + Math.abs( V1.Double ) * J1.err;
            pre_val  = Math.sqrt( xi / sh );
            pre_err  = 2.0 * Math.abs( pre_val );
            r.val = pre_val * bessterm;
            r.err = pre_val * besserr;
            r.err = r.err + pre_err * Math.abs( bessterm );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            return r;
        }
    }

} // gsl_sf_conicalP_0_e

// ----------------------------------------------------------------------------

// P^1_{-1/2 + I lambda}
//
export function gsl_sf_conicalP_1_e( lambda, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result

    if ( x <= -1.0 )
    {
        throw "SF.DomainException";
    }
    else if ( lambda == 0.0 )
    {
        let K = { val: 0.0, err: 0.0 }; // Result;
        let E = { val: 0.0, err: 0.0 }; // Result;

        if ( x == 1.0 )
        {
            r.val = 0.0;
            r.err = 0.0;
            return r;
        }
        else if ( x < 1.0 )
        {
            if ( 1.0 - x < GSL_SQRT_DBL_EPSILON )
            {
                let err_amp = Math.max( 1.0, 1.0 / (GSL_DBL_EPSILON + Math.abs( 1.0 - x )) );

                r.val = 0.25 / M_SQRT2 * Math.sqrt( 1.0 - x ) * (1.0 + 5.0 / 16.0 * (1.0 - x));
                r.err = err_amp * 3.0 * GSL_DBL_EPSILON * Math.abs( r.val );
                return r;
            }
            else
            {
                let th  = Math.acos( x );
                let s   = Math.sin( 0.5 * th );
                let c2  = 1.0 - s * s;
                let sth = Math.sin( th );
                let pre = 2.0 / (M_PI * sth);

                K = gsl_sf_ellint_Kcomp_e( s, GSL_MODE_DEFAULT );
                E = gsl_sf_ellint_Ecomp_e( s, GSL_MODE_DEFAULT );
                r.val = pre * (E.val - c2 * K.val);
                r.err = pre * (E.err + Math.abs( c2 ) * K.err);
                r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
                return r;
            }
        }
        else
        {
            if ( x - 1.0 < GSL_SQRT_DBL_EPSILON )
            {
                let err_amp = Math.max( 1.0, 1.0 / (GSL_DBL_EPSILON + Math.abs( 1.0 - x )) );

                r.val = -0.25 / M_SQRT2 * Math.sqrt( x - 1.0 ) * (1.0 - 5.0 / 16.0 * (x - 1.0));
                r.err = err_amp * 3.0 * GSL_DBL_EPSILON * Math.abs( r.val );
                return r;
            }
            else
            {
                let xi  = Math.acosh( x );
                let c   = Math.cosh( 0.5 * xi );
                let t   = Math.tanh( 0.5 * xi );
                let sxi = Math.sinh( xi );
                let pre = 2.0 / (M_PI * sxi) * c;

                K = gsl_sf_ellint_Kcomp_e( t, GSL_MODE_DEFAULT );
                E = gsl_sf_ellint_Ecomp_e( t, GSL_MODE_DEFAULT );
                r.val = pre * (E.val - K.val);
                r.err = pre * (E.err + K.err);
                r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
                return r;
            }
        }
    }
    else if ( (x <= 0.0 && lambda < 1000.0) || (x <  0.1 && lambda < 17.0) || (x <  0.2 && lambda < 5.0) )
    {
        return conicalP_xlt1_hyperg_A( 1.0, lambda, x );
    }
    else if ( (x <= 0.2 && lambda < 17.0) || (x <  1.5 && lambda < 20.0) )
    {
        let arg = 0.0;
        let sgn = 0.0;
        let pre = 0.0;
        let F   = { val: 0.0, err: 0.0 }; // Result;

        arg = Math.abs( x * x - 1.0 );
        sgn = GSL_SIGN( 1.0 - x );
        pre = 0.5 * (lambda * lambda + 0.25) * sgn * Math.sqrt( arg );
        F = gsl_sf_hyperg_2F1_conj_e( 1.5, lambda, 2.0, (1.0 - x) / 2.0 );
        r.val = pre * F.val;
        r.err = Math.abs( pre ) * F.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( 1.5 <= x && lambda < Math.max( x, 20.0 ) )
    {
        let lm = { Double: 0.0 };
        let P  = { val: 0.0, err: 0.0 }; // Result;

        gsl_sf_conicalP_large_x_e( 1.0, lambda, x, P, lm );
        r = gsl_sf_exp_mult_err_e( lm.Double, 2.0 * GSL_DBL_EPSILON * Math.abs( lm.Double ), P.val, P.err );
        return r;
    }
    else
    {
        let V0 = { Double: 0.0 };
        let V1 = { Double: 0.0};

        if ( x < 1.0 )
        {
            let sqrt_1mx = 0.0;
            let sqrt_1px = 0.0;
            let th       = 0.0;
            let sth      = 0.0;
            let bessterm = 0.0;
            let besserr  = 0.0;
            let arg1     = 0.0;
            let sqts     = 0.0;
            let I0       = { val: 0.0, err: 0.0 }; // Result;
            let I1       = { val: 0.0, err: 0.0 }; // Result;

            sqrt_1mx = Math.sqrt( 1.0 - x );
            sqrt_1px = Math.sqrt( 1.0 + x );
            th  = Math.acos( x );
            sth = sqrt_1mx * sqrt_1px;  // sin(th)
            I0 = gsl_sf_bessel_I0_scaled_e( th * lambda );
            I1 = gsl_sf_bessel_I1_scaled_e( th * lambda );
            conicalP_1_V(th, x / sth, lambda, -1.0, V0, V1);
            bessterm = V0.Double * I0.val + V1.Double * I1.val;
            besserr  =  Math.abs( V0.Double ) * I0.err + Math.abs( V1.Double ) * I1.err
                            + 2.0 * GSL_DBL_EPSILON * Math.abs( V0.Double * I0.val )
                            + 2.0 * GSL_DBL_EPSILON * Math.abs( V1.Double * I1.val );
            arg1 = th * lambda;
            sqts = Math.sqrt( th / sth );
            r = gsl_sf_exp_mult_err_e( arg1, 2.0 * GSL_DBL_EPSILON * Math.abs(arg1), sqts * bessterm, sqts * besserr );
            r.err = r.err * (1.0 / sqrt_1mx);
            return r;
        }
        else
        {
            let sqrt_xm1 = 0.0;
            let sqrt_xp1 = 0.0;
            let sh       = 0.0;
            let xi       = 0.0;
            let xi_lam   = 0.0;
            let bessterm = 0.0;
            let besserr  = 0.0;
            let pre      = 0.0;
            let J0       = { val: 0.0, err: 0.0 }; // Result;
            let J1       = { val: 0.0, err: 0.0 }; // Result;

            sqrt_xm1 = Math.sqrt( x - 1.0 );
            sqrt_xp1 = Math.sqrt( x + 1.0 );
            sh = sqrt_xm1 * sqrt_xp1;  // sinh(xi)
            xi = Math.log( x + sh );   // xi = acosh(x)
            xi_lam = xi * lambda;
            J0 = gsl_sf_bessel_J0_e( xi_lam );
            J1 = gsl_sf_bessel_J1_e( xi_lam );
            conicalP_1_V( xi, x / sh, lambda, 1.0, V0, V1 );
            bessterm = V0.Double * J0.val + V1.Double * J1.val;
            besserr  = Math.abs( V0.Double ) * J0.err + Math.abs( V1.Double ) * J1.err
                            + 512.0 * 2.0 * GSL_DBL_EPSILON * Math.abs( V0.Double * J0.val )
                            + 512.0 * 2.0 * GSL_DBL_EPSILON * Math.abs( V1.Double * J1.val )
                            + GSL_DBL_EPSILON * Math.abs( xi_lam * V0.Double * J1.val )
                            + GSL_DBL_EPSILON * Math.abs( xi_lam * V1.Double * J0.val );
            pre = Math.sqrt( xi / sh );
            r.val = pre * bessterm;
            r.err = pre * besserr * sqrt_xp1 / sqrt_xm1;
            r.err = r.err + 4.0 * GSL_DBL_EPSILON * Math.abs( r.val );
            return r;
        }
    }

} // gsl_sf_conicalP_1_e

// ----------------------------------------------------------------------------

// P^{1/2}_{-1/2 + I lambda} (x)
// [Abramowitz+Stegun 8.6.8, 8.6.12]
// checked OK [GJ] Fri May  8 12:24:36 MDT 1998 
//
export function gsl_sf_conicalP_half_e(lambda, x)
{
    var err_amp  = 0.0;
    var ac       = 0.0;
    var den      = 0.0;
    var sq_term  = 0.0;
    var ln_term  = 0.0;
    var carg_val = 0.0;
    var carg_err = 0.0;

    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= -1.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 1.0)
    {
        err_amp = 1.0 + 1.0 / (GSL_DBL_EPSILON + Math.abs(1.0 - Math.abs(x)));
        ac  = Math.acos(x);
        den = Math.sqrt(Math.sqrt(1.0 - x) * Math.sqrt(1.0 + x));
        r.val = Root_2OverPi / den * Math.cosh(ac * lambda);
        r.err = err_amp * 3.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        r.err = r.err * Math.abs(ac * lambda) + 1.0;
        return r;
    }
    else if (x == 1.0)
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        // x > 1
        err_amp = 1.0 + 1.0 / (GSL_DBL_EPSILON + Math.abs(1.0 - Math.abs(x)));
        sq_term = Math.sqrt(x - 1.0) * Math.sqrt(x + 1.0);
        ln_term = Math.log(x + sq_term);
        den = Math.sqrt(sq_term);
        carg_val = lambda * ln_term;
        carg_err = 2.0 * GSL_DBL_EPSILON * Math.abs(carg_val);
        c = gsl_sf_cos_err_e(carg_val, carg_err);
        r.val = Root_2OverPi / den * c.val;
        r.err = err_amp * Root_2OverPi / den * c.err;
        r.err = r.err + 4.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }

} // gsl_sf_conicalP_half_e

// ----------------------------------------------------------------------------

// P^{-1/2}_{-1/2 + I lambda} (x)
// [Abramowitz+Stegun 8.6.9, 8.6.14]
// checked OK [GJ] Fri May  8 12:24:43 MDT 1998 
//
export function gsl_sf_conicalP_mhalf_e(lambda, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= -1.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 1.0)
    {
        var ac  = 0.0;
        var den = 0.0;
        var arg = 0.0;
        var err_amp = 0.0;

        ac  = Math.acos(x);
        den = Math.sqrt(Math.sqrt(1.0 - x) * Math.sqrt(1.0 + x));
        arg = ac * lambda;
        err_amp = 1.0 + 1.0 / (GSL_DBL_EPSILON + Math.abs(1.0 - Math.abs(x)));
        if (Math.abs(arg) < GSL_SQRT_DBL_EPSILON)
        {
            r.val = Root_2OverPi / den * ac;
            r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            r.err = r.err * err_amp;
        }
        else
        {
            r.val = Root_2OverPi / (den * lambda) * Math.sinh(arg);
            r.err = GSL_DBL_EPSILON * (Math.abs(arg) + 1.0) * Math.abs(r.val);
            r.err = r.err * err_amp;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
        return r;
    }
    else if (x == 1.0)
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        // x > 1
        let sq_term = 0.0;
        let ln_term = 0.0;
        let den     = 0.0;
        let arg_val = 0.0;
        let arg_err = 0.0;
        let s = { val: 0.0, err: 0.0 }; // Result;

        sq_term = Math.sqrt(x - 1.0) * Math.sqrt(x + 1.0);
        ln_term = Math.log(x + sq_term);
        den     = Math.sqrt(sq_term);
        arg_val = lambda * ln_term;
        arg_err = 2.0 * GSL_DBL_EPSILON * Math.abs(arg_val);
        if (arg_val < GSL_SQRT_DBL_EPSILON)
        {
            r.val = Root_2OverPi / den * ln_term;
            r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            return r;
        }
        else
        {
            s = gsl_sf_sin_err_e(arg_val, arg_err);
            r.val = Root_2OverPi / (den * lambda) * s.val;
            r.err = Root_2OverPi / Math.abs(den * lambda) * s.err;
            r.err = r.err + 3.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            return r;
        }
    }

} // gsl_sf_conicalP_mhalf_e

// ----------------------------------------------------------------------------

export function gsl_sf_conicalP_sph_reg_e( l, lambda, x )
{

    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x <= -1.0 || l < -1 )
    {
        throw "SF.DomainException";
    }
    else if ( l == -1 )
    {
        return gsl_sf_conicalP_half_e( lambda, x );
    }
    else if ( l == 0 )
    {
        return gsl_sf_conicalP_mhalf_e( lambda, x );
    }
    else if ( x == 1.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else if ( x < 0.0 )
    {
        let c        = 0.0;
        let d        = 0.0;
        let Pellm1   = 0.0;
        let Pell     = 0.0;
        let Pellp1   = 0.0;
        let r_Pellm1 = { val: 0.0, err: 0.0 }; // Result;
        let r_Pell   = { val: 0.0, err: 0.0 }; // Result;

        c        = 1.0 / Math.sqrt( 1.0 - x * x );
        r_Pellm1 = gsl_sf_conicalP_half_e( lambda, x );  // P^( 1/2)
        r_Pell   = gsl_sf_conicalP_mhalf_e( lambda, x ); // P^(-1/2)
        Pellm1   = r_Pellm1.val;
        Pell     = r_Pell.val;
        
        for ( let ell = 0; ell <= l - 1; ell++ )
        {
            d      = (ell + 1) * (ell + 1) + lambda * lambda;
            Pellp1 = (Pellm1 - (2 * ell + 1) * c * x * Pell) / d;
            Pellm1 = Pell;
            Pell   = Pellp1;
        }
        
        r.val = Pell;
        r.err = (0.5 * (l) + 1.0) * GSL_DBL_EPSILON * Math.abs( Pell );
        r.err = r.err + GSL_DBL_EPSILON * (l) * Math.abs( r.val );
        return r;
    }
    else if ( x < 1.0 )
    {
        let xi     = 0.0;
        let d      = 0.0;
        let Pellp1 = 0.0;
        let Pell   = 0.0;
        let Pellm1 = 0.0;
        let rat    = { val: 0.0, err: 0.0 }; // Result;
        let Phf    = { val: 0.0, err: 0.0 }; // Result;

        xi     = x / (Math.sqrt( 1.0 - x ) * Math.sqrt( 1.0 + x ));
        rat    = conicalP_negmu_xlt1_CF1( 0.5, l, lambda, x );
        Phf    = gsl_sf_conicalP_half_e( lambda, x );
        Pellp1 = rat.val * GSL_SQRT_DBL_MIN;
        Pell   = GSL_SQRT_DBL_MIN;
    
        for ( let ell = l; ell >= 0; ell-- )
        {
            d      = (ell + 1) * (ell + 1) + lambda * lambda;
            Pellm1 = (2 * ell + 1) * xi * Pell + d * Pellp1;
            Pellp1 = Pell;
            Pell   = Pellm1;
        }
    
        r.val = GSL_SQRT_DBL_MIN * Phf.val / Pell;
        r.err = GSL_SQRT_DBL_MIN * Phf.err / Math.abs( Pell );
        r.err = r.err + Math.abs( rat.err / rat.val ) * (l + 1) * Math.abs( r.val );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
    
        return r;
    }
    else if ( x == 1.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        // x > 1.0
        let xi     = 0.0;
        let d      = 0.0;
        let Pellp1 = 0.0;
        let Pell   = 0.0;
        let Pellm1 = 0.0;
        let rat    = { val: 0.0, err: 0.0 }; // Result;
        let Phf    = { val: 0.0, err: 0.0 }; // Result;
        let Pmhf   = { val: 0.0, err: 0.0 }; // Result;

        xi     = x / Math.sqrt( (x - 1.0) * (x + 1.0) );
        rat    = conicalP_negmu_xgt1_CF1( 0.5, l, lambda, x );
        Pellp1 = rat.val * GSL_SQRT_DBL_MIN;
        Pell   = GSL_SQRT_DBL_MIN;

        for ( let ell = l; ell >= 0; ell-- )
        {
            d      = (ell + 1) * (ell + 1) + lambda * lambda;
            Pellm1 = (2 * ell + 1) * xi * Pell - d * Pellp1;
            Pellp1 = Pell;
            Pell   = Pellm1;
        }
        
        if ( Math.abs( Pell ) > Math.abs( Pellp1 ) )
        {
            Phf = gsl_sf_conicalP_half_e( lambda, x );
            r.val =       GSL_SQRT_DBL_MIN * Phf.val / Pell;
            r.err = 2.0 * GSL_SQRT_DBL_MIN * Phf.err / Math.abs( Pell );
            r.err = r.err + 2.0 * Math.abs( rat.err / rat.val ) * (l + 1) * Math.abs( r.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        }
        else
        {
            Pmhf = gsl_sf_conicalP_mhalf_e( lambda, x );
            r.val =       GSL_SQRT_DBL_MIN * Pmhf.val / Pellp1;
            r.err = 2.0 * GSL_SQRT_DBL_MIN * Pmhf.err / Math.abs( Pellp1 );
            r.err = r.err + 2.0 * Math.abs( rat.err / rat.val ) * (l + 1) * Math.abs( r.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        }

        return r;
    }

} // gsl_sf_conicalP_sph_reg_e

// ----------------------------------------------------------------------------

export function gsl_sf_conicalP_cyl_reg_e( m, lambda, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x <= -1.0 || m < -1 )
    {
        throw "SF.DomainException";
    }
    else if ( m == -1 )
    {
        return gsl_sf_conicalP_1_e( lambda, x );
    }
    else if ( m == 0 )
    {
        return gsl_sf_conicalP_0_e( lambda, x );
    }
    else if ( x == 1.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else if ( x < 0.0 )
    {
        let c      = 0.0;
        let d      = 0.0;
        let Pkm1   = 0.0;
        let Pk     = 0.0;
        let Pkp1   = 0.0;
        let r_Pkm1 = { val: 0.0, err: 0.0 }; // Result;
        let r_Pk   = { val: 0.0, err: 0.0 }; // Result;

        c      = 1.0 / Math.sqrt( 1.0 - x * x );
        r_Pkm1 = gsl_sf_conicalP_1_e( lambda, x ); // P^1
        r_Pk   = gsl_sf_conicalP_0_e( lambda, x ); // P^0
        Pkm1   = r_Pkm1.val;
        Pk     = r_Pk.val;

        for ( let k = 0; k <= m - 1; k++ )
        {
            d    = ((k) + 0.5) * ((k) + 0.5) + lambda * lambda;
            Pkp1 = (Pkm1 - 2.0 * (k) * c * x * Pk) / d;
            Pkm1 = Pk;
            Pk   = Pkp1;
        }

        r.val = Pk;
        r.err = (m + 2) * GSL_DBL_EPSILON * Math.abs( Pk );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );

        return r;
    }
    else if ( x < 1.0 )
    {
        let xi   = 0.0;
        let d    = 0.0;
        let Pkp1 = 0.0;
        let Pk   = 0.0;
        let Pkm1 = 0.0;
        let rat  = { val: 0.0, err: 0.0 }; // Result;
        let P0   = { val: 0.0, err: 0.0 }; // Result;

        xi   = x / (Math.sqrt( 1.0 - x ) * Math.sqrt( 1.0 + x ));
        rat  = conicalP_negmu_xlt1_CF1( 0.0, m, lambda, x );
        P0   = gsl_sf_conicalP_0_e( lambda, x );
        Pkp1 = rat.val * GSL_SQRT_DBL_MIN;
        Pk   = GSL_SQRT_DBL_MIN;

        for ( let k = m; k >= 1; k-- )
        {
            d    = ((k) + 0.5) * ((k) + 0.5) + lambda * lambda;
            Pkm1 = 2.0 * (k) * xi * Pk + d * Pkp1;
            Pkp1 = Pk;
            Pk   = Pkm1;
        }

        r.val = GSL_SQRT_DBL_MIN * P0.val / Pk;
        r.err = 2.0 * GSL_SQRT_DBL_MIN * P0.err / Math.abs( Pk );
        r.err = r.err + 2.0 * Math.abs( rat.err / rat.val ) * (m + 1) * Math.abs( r.val );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );

        return r;
    }
    else if ( x == 1.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        // x > 1.0
        let xi   = 0.0;
        let d    = 0.0;
        let Pkp1 = 0.0;
        let Pk   = 0.0;
        let Pkm1 = 0.0;
        let P0   = { val: 0.0, err: 0.0 }; // Result;
        let P1   = { val: 0.0, err: 0.0 }; // Result;
        let rat  = { val: 0.0, err: 0.0 }; // Result;

        xi   = x / Math.sqrt( (x - 1.0) * (x + 1.0) );
        rat  = conicalP_negmu_xgt1_CF1( 0.0, m, lambda, x );
        Pkp1 = rat.val * GSL_SQRT_DBL_MIN;
        Pk   = GSL_SQRT_DBL_MIN;

        for ( let k = m; k >= 0; k-- )
        {
            d    = ((k) + 0.5) * ((k) + 0.5) + lambda * lambda;
            Pkm1 = 2.0 * (k) * xi * Pk - d * Pkp1;
            Pkp1 = Pk;
            Pk   = Pkm1;
        }

        if ( Math.abs( Pk ) > Math.abs( Pkp1 ) )
        {
            P1 = gsl_sf_conicalP_1_e( lambda, x );
            r.val = GSL_SQRT_DBL_MIN * P1.val / Pk;
            r.err = 2.0 * GSL_SQRT_DBL_MIN * P1.err / Math.abs( Pk );
            r.err = r.err + 2.0 * Math.abs( rat.err / rat.val ) * (m + 2) * Math.abs( r.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        }
        else
        {
            P0 = gsl_sf_conicalP_0_e( lambda, x );
            r.val = GSL_SQRT_DBL_MIN * P0.val / Pkp1;
            r.err = 2.0 * GSL_SQRT_DBL_MIN * P0.err / Math.abs( Pkp1 );
            r.err = r.err + 2.0 * Math.abs( rat.err / rat.val ) * (m + 2) * Math.abs( r.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        }

        return r;
    }

} // gsl_sf_conicalP_cyl_reg_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

function gsl_sf_conicalP_0( lambda, x )
{ // gsl_sf_conicalP_0
    return EVAL_RESULT_DD( gsl_sf_conicalP_0_e, { x: lambda, y: x }, "gsl_sf_conicalP_0" );
} // gsl_sf_conicalP_0;

function gsl_sf_conicalP_1( lambda, x )
{ // gsl_sf_conicalP_1
    return EVAL_RESULT_DD( gsl_sf_conicalP_1_e, { x: lambda, y: x }, "gsl_sf_conicalP_1" );
} // gsl_sf_conicalP_1;

function gsl_sf_conicalP_half( lambda, x )
{ // gsl_sf_conicalP_half
    return EVAL_RESULT_DD( gsl_sf_conicalP_half_e, { x: lambda, y: x }, "gsl_sf_conicalP_half" );
} // gsl_sf_conicalP_half;

function gsl_sf_conicalP_mhalf( lambda, x )
{ // gsl_sf_conicalP_mhalf
    return EVAL_RESULT_DD( gsl_sf_conicalP_mhalf_e, { x: lambda, y:x }, "gsl_sf_conicalP_mhalf" );
} // gsl_sf_conicalP_mhalf;

function gsl_sf_conicalP_sph_reg( l, lambda, x )
{ // gsl_sf_conicalP_sph_reg
    return EVAL_RESULT_IDD( gsl_sf_conicalP_sph_reg_e, { n: l, x: lambda, y: x }, "gsl_sf_conicalP_sph_reg" );
} // gsl_sf_conicalP_sph_reg;

function gsl_sf_conicalP_cyl_reg( m, lambda, x )
{ // gsl_sf_conicalP_cyl_reg
    return EVAL_RESULT_IDD( gsl_sf_conicalP_cyl_reg_e, { n: m, x: lambda, y: x }, "gsl_sf_conicalP_cyl_reg" );
} // gsl_sf_conicalP_cyl_reg;

// ----------------------------------------------------------------------------
// EOF SF-LegendreConical.mjs
