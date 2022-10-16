// SF-Bessel.mjs
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
import { GSL_ROOT5_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }      from "./SF-Machine.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { M_SQRT2 }               from "./SF-Math.mjs";
import { gsl_sf_bessel_asymp_Mnu_e } from "./SF-BesselAF.mjs";
import { gsl_sf_bessel_asymp_thetanu_corr_e } from "./SF-BesselAF.mjs";
import { gsl_sf_poch_e }         from "./SF-Pochhammer.mjs";
import { gsl_sf_taylorcoeff_e }  from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_e }      from "./SF-Gamma.mjs";
import { gsl_sf_hypot }          from "./SF-Trigonometric.mjs";
import { gsl_sf_multiply_err_e } from "./SF-Elementary.mjs";
import { gsl_sf_exp_e }          from "./SF-Exponential.mjs";
import { gsl_sf_exp_err_e }      from "./SF-Exponential.mjs";
import { gsl_sf_bessel_Y_temme } from "./SF-BesselTemme.mjs";

// //#define CubeRoot2_  1.25992104989487316476721060728

// TYPE ALF16 IS ARRAY(0..15) OF LONG_FLOAT;

// Debye functions [Abramowitz+Stegun, 9.3.9-10]

function debye_u1(tpow)
{
    return (3.0 * tpow[1] - 5.0 * tpow[3]) / 24.0;
} // debye_u1

function debye_u2(tpow)
{
    return (81.0 * tpow[2] - 462.0 * tpow[4] + 385.0 * tpow[6]) / 1152.0;
} // debye_u2

function debye_u3(tpow)
{
    return (30375.0 * tpow[3] - 369603.0 * tpow[5] + 765765.0 * tpow[7] - 425425.0 * tpow[9]) / 414720.0;
} // debye_u3

function debye_u4(tpow)
{
    return (4465125.0 * tpow[4] - 94121676.0 * tpow[6] + 349922430.0 * tpow[8] - 
          446185740.0 * tpow[10] + 185910725.0 * tpow[12]) / 39813120.0;
} // debye_u4

function debye_u5(tpow)
{
    return (1519035525.0 * tpow[5]     - 49286948607.0 * tpow[7] + 
          284499769554.0 * tpow[9]   - 614135872350.0 * tpow[11] + 
          566098157625.0 * tpow[13]  - 188699385875.0 * tpow[15]) / 6688604160.0;
} // debye_u5

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_IJ_taylor_e(nu, x, sign, kmax, threshold)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
  
    if (nu < 0.0 || x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        if (nu == 0.0)
        {
          r.val = 1.0;
          r.err = 0.0;
        }
        else
        {
          r.val = 0.0;
          r.err = 0.0;
        }
    }
    else
    {
        const INTEGER_LAST = +2147483647;
        var k = 0;
        var N = 0;

        var y      = 0.0;
        var p      = 0.0;
        var f      = 0.0;
        var sumk   = 0.0;
        var term   = 0.0;
        var term1  = 0.0;
        var term2  = 0.0;
        var ln_pre = 0.0;
        var ln_pre_err = 0.0;

        var prefactor   = { val: 0.0, err: 0.0 }; // Result;   // (x/2)^nu / Gamma(nu+1)
        var sum         = { val: 0.0, err: 0.0 }; // Result;
        var lg          = { val: 0.0, err: 0.0 }; // Result;
        var poch_factor = { val: 0.0, err: 0.0 }; // Result;
        var tc_factor   = { val: 0.0, err: 0.0 }; // Result;

        if (nu == 0.0)
        {
            prefactor.val = 1.0;
            prefactor.err = 0.0;
        }
        else if (nu < (INTEGER_LAST - 1))
        {
            // Separate the integer part and use
            // y^nu / Gamma(nu+1) = y^N /N! y^f / (N+1)_f,
            // to control the error.
            //
            N = Math.trunc(Math.floor(nu + 0.5));
            f = nu - (N);
            poch_factor = gsl_sf_poch_e((N) + 1.0, f);
            tc_factor = gsl_sf_taylorcoeff_e(N, 0.5 * x);
            p = (0.5 * x) ** f;
            prefactor.val = tc_factor.val * p / poch_factor.val;
            prefactor.err = tc_factor.err * p / poch_factor.val;
            prefactor.err = prefactor.err + Math.abs(prefactor.val) / poch_factor.val * poch_factor.err;
            prefactor.err = prefactor.err + 2.0 * GSL_DBL_EPSILON * Math.abs(prefactor.val);
        }
        else
        {
            lg = gsl_sf_lngamma_e(nu + 1.0);
            term1  = nu * Math.log(0.5 * x);
            term2  = lg.val;
            ln_pre = term1 - term2;
            ln_pre_err = GSL_DBL_EPSILON * (Math.abs(term1) + Math.abs(term2)) + lg.err;
            prefactor = gsl_sf_exp_err_e(ln_pre, ln_pre_err);
        }
        
        // Evaluate the sum.
        // [Abramowitz+Stegun, 9.1.10]
        // [Abramowitz+Stegun, 9.6.7]
        //
        y = (sign) * 0.25 * x * x;
        sumk = 1.0;
        term = 1.0;
        
        k = 1;
        while (k <= kmax)
        {
            term = term * y / ((nu + (k)) * (k));
            sumk = sumk + term;
            if (Math.abs(term / sumk) < threshold) break;
            k = k + 1;
        }
        
        sum.val = sumk;
        sum.err = threshold * Math.abs(sumk);

        //IF (k >= kmax)
        //    stat_sum = GSL_EMAXITER;
        //ELSE
        //    stat_sum = GSL_SUCCESS;
        //END IF;
        
        r = gsl_sf_multiply_err_e(prefactor.val, prefactor.err, sum.val, sum.err);

    }

    return r;

} // gsl_sf_bessel_IJ_taylor_e

// ----------------------------------------------------------------------------

// Hankel's Asymptotic Expansion - A&S 9.2.5
//  
// x >> nu*nu+1
// error ~ O( ((nu*nu+1)/x)^4 )
//
// empirical error analysis:
//   choose  GSL_ROOT4_MACH_EPS * x > (nu*nu + 1)
//
// This is not especially useful. When the argument gets
// large enough for this to apply, the cos() and sin()
// start loosing digits. However, this seems inevitable
// for this particular method.
//
// Wed Jun 25 14:39:38 MDT 2003 [GJ]
// This function was inconsistent since the Q term did not
// go to relative order eps^2. That's why the error estimate
// originally given was screwy (it didn't make sense that the
// "empirical" error was coming out O(eps^3)).
// With Q to proper order, the error is O(eps^4).
//
// Sat Mar 15 05:16:18 GMT 2008 [BG]
// Extended to use additional terms in the series to gain
// higher accuracy.
//
export function gsl_sf_bessel_Jnu_asympx_e(nu, x)
{
    var mu  = 0.0;
    var chi = 0.0;
    var P   = 0.0;
    var Q   = 0.0;
    var k   = 0.0;
    var t   = 0.0;
    var c   = 0.0;
    var s   = 0.0;
    var pre = 0.0;
    var convP = false;
    var convQ = false;

    var r = { val: 0.0, err: 0.0 }; // Result;

    mu  = 4.0 * nu * nu;
    chi = x - (0.5 * nu + 0.25) * M_PI;
    k   = 0.0;
    t   = 1.0;
  
    while (true)
    {
        if (k != 0.0)
        {
            t = t * (-(mu - (2.0 * k - 1.0) * (2.0 * k - 1.0)) / (k * (8.0 * x)));
        }
        convP = Math.abs(t) < GSL_DBL_EPSILON * Math.abs(P);
        P = P + t;
        
        k = k + 1.0;
       
        t = t * ((mu - (2.0 * k - 1.0) * (2.0 * k - 1.0)) / (k * (8.0 * x)));
        convQ = Math.abs(t) < GSL_DBL_EPSILON * Math.abs(Q);
        Q = Q + t;
       
        // To preserve the consistency of the series we need to exit
        // when P and Q have the same number of terms
       
        if (convP && convQ && k > (nu / 2.0)) break;
       
        k = k + 1.0;
        if (k >= 1000.0) break;
    }
  
    pre = Math.sqrt(2.0 / (M_PI * x));
    c   = Math.cos(chi);
    s   = Math.sin(chi);
    
    r.val = pre * (c * P - s * Q);
    r.err = pre * GSL_DBL_EPSILON * (Math.abs(c * P) + Math.abs(s * Q) + Math.abs(t)) * (1.0 + Math.abs(x));
    // NB: final term accounts for phase error with large x
  
    return r;

} // gsl_sf_bessel_Jnu_asympx_e

// ----------------------------------------------------------------------------

// x >> nu*nu+1
//
export function gsl_sf_bessel_Ynu_asympx_e(nu, x)
{
    var ampl      = 0.0;
    var theta     = 0.0;
    var alpha     = 0.0;
    var beta      = 0.0;
    var sin_alpha = 0.0;
    var cos_alpha = 0.0;
    var sin_chi   = 0.0;
    var cos_chi   = 0.0;
    var sin_term  = 0.0;
    var sin_term_mag  = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    alpha = x;
    beta  = -0.5 * nu * M_PI;
    ampl = gsl_sf_bessel_asymp_Mnu_e(nu, x);
    theta = gsl_sf_bessel_asymp_thetanu_corr_e(nu, x);
    sin_alpha = Math.sin(alpha);
    cos_alpha = Math.cos(alpha);
    sin_chi   = Math.sin(beta + theta);
    cos_chi   = Math.cos(beta + theta);
    sin_term     = sin_alpha * cos_chi + sin_chi * cos_alpha;
    sin_term_mag = Math.abs(sin_alpha * cos_chi) + Math.abs(sin_chi * cos_alpha);

    r.val = ampl * sin_term;
    r.err = Math.abs(ampl) * GSL_DBL_EPSILON * sin_term_mag;
    r.err = r.err + Math.abs(r.val) * 2.0 * GSL_DBL_EPSILON;
  
    if (Math.abs(alpha) > 1.0/GSL_DBL_EPSILON)
    {
        r.err = r.err * 0.5 * Math.abs(alpha);
    }
    else if (Math.abs(alpha) > 1.0 / GSL_SQRT_DBL_EPSILON)
    {
        r.err = r.err * 256.0 * Math.abs(alpha) * GSL_SQRT_DBL_EPSILON;
    }
  
    return r;

} // gsl_sf_bessel_Ynu_asympx_e

// ----------------------------------------------------------------------------

// x >> nu*nu+1
//
export function gsl_sf_bessel_Inu_scaled_asympx_e( nu, x )
{
    var mu   = 4.0 * nu * nu;
    var mum1 = mu - 1.0;
    var mum9 = mu - 9.0;
    var pre  = 1.0 / Math.sqrt( 2.0 * M_PI * x );
    var r    = mu / x;
    var result = { val: 0.0, err: 0.0 };
    result.val = pre * (1.0 - mum1 / (8.0 * x) + mum1 * mum9 / (128.0 * x * x));
    result.err = 2.0 * GSL_DBL_EPSILON * Math.abs( result.val ) + pre * Math.abs( 0.1 * r * r * r );
    return result;
}

// ----------------------------------------------------------------------------

// x >> nu*nu+1
//
export function gsl_sf_bessel_Knu_scaled_asympx_e( nu, x )
{
    var mu   = 4.0 * nu * nu;
    var mum1 = mu - 1.0;
    var mum9 = mu - 9.0;
    var pre  = Math.sqrt(M_PI / (2.0 * x));
    var s    = nu / x;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    r.val = pre * (1.0 + mum1 / (8.0 * x) + mum1 * mum9 / (128.0 * x * x));
    r.err = 2.0 * GSL_DBL_EPSILON * Math.abs( r.val ) + pre * Math.abs( 0.1 * s * s * s );
    return r;

} // gsl_sf_bessel_Knu_scaled_asympx_e

// ----------------------------------------------------------------------------

// nu . Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.7]
//
// error:
//   The error has the form u_N(t)/nu^N  where  0 <= t <= 1.
//   It is not hard to show that |u_N(t)| is small for such t.
//   We have N=6 here, and |u_6(t)| < 0.025, so the error is clearly
//   bounded by 0.025/nu^6. This gives the asymptotic bound on nu
//   seen below as nu ~ 100. For general MACH_EPS it will be 
//                     nu > 0.5 / MACH_EPS^(1/6)
//   When t is small, the bound is even better because |u_N(t)| vanishes
//   as t.0. In fact u_N(t) ~ C t^N as t.0, with C ~= 0.1.
//   We write
//                     err_N <= min(0.025, C(1/(1+(x/nu)^2))^3) / nu^6
//   therefore
//                     min(0.29/nu^2, 0.5/(nu^2+x^2)) < MACH_EPS^{1/3}
//   and this is the general form.
//
// empirical error analysis, assuming 14 digit requirement:
//   choose   x > 50.000 nu   ==>  nu >   3
//   choose   x > 10.000 nu   ==>  nu >  15
//   choose   x >  2.000 nu   ==>  nu >  50
//   choose   x >  1.000 nu   ==>  nu >  75
//   choose   x >  0.500 nu   ==>  nu >  80
//   choose   x >  0.100 nu   ==>  nu >  83
//
// This makes sense. For x << nu, the error will be of the form u_N(1)/nu^N,
// since the polynomial term will be evaluated near t=1, so the bound
// on nu will become constant for small x. Furthermore, increasing x with
// nu fixed will decrease the error.
//
export function gsl_sf_bessel_Inu_scaled_asymp_unif_e(nu, x)
{
    var z         = 0.0;
    var root_term = 0.0;
    var pre       = 0.0;
    var eta       = 0.0;
    var ex_arg    = 0.0;
    var t         = 0.0;
    var sum       = 0.0;
    var tpow      = [];
    var ex_result = { val: 0.0, err: 0.0 }; // Result;
    var r         = { val: 0.0, err: 0.0 }; // Result;

    z = x / nu;
    root_term = gsl_sf_hypot(1.0, z);
    pre = 1.0 / Math.sqrt(2.0 * M_PI * nu * root_term);
    eta = root_term + Math.log(z / (1.0 + root_term));
    if (z < 1.0 / GSL_ROOT3_DBL_EPSILON)
    {
        ex_arg = nu * (-z + eta);
    }
    else
    {
        ex_arg = -0.5 * nu / z * (1.0 - 1.0 / (12.0 * z * z));
    }
    ex_result = gsl_sf_exp_e(ex_arg);
    t = 1.0 / root_term;
    tpow[0] = 1.0;
    for (let i = 1; i <= 16 - 1; i++)
    {
        tpow[i] = t * tpow[i-1];
    }
    sum = 1.0 + debye_u1(tpow) / nu + debye_u2(tpow) / (nu * nu) + debye_u3(tpow) / (nu * nu * nu)
          + debye_u4(tpow) / (nu * nu * nu * nu) + debye_u5(tpow) / (nu * nu * nu * nu * nu);
    r.val = pre * ex_result.val * sum;
    r.err = pre * ex_result.val / (nu * nu * nu * nu * nu * nu);
    r.err = r.err + pre * ex_result.err * Math.abs(sum);
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);

    return r;

} // gsl_sf_bessel_Inu_scaled_asymp_unif_e

// ----------------------------------------------------------------------------

// nu . Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.8]
//
// error:
//   identical to that above for Inu_scaled
//
export function gsl_sf_bessel_Knu_scaled_asymp_unif_e( nu, x )
{
    var z         = 0.0;
    var root_term = 0.0;
    var pre       = 0.0;
    var eta       = 0.0;
    var ex_arg    = 0.0;
    var t         = 0.0;
    var sum       = 0.0;
    var tpow      = [];
    var ex_result = { val: 0.0, err: 0.0 }; // Result;
    var r         = { val: 0.0, err: 0.0 }; // Result;

    z = x / nu;
    root_term = gsl_sf_hypot( 1.0, z ); // hypot
    pre = Math.sqrt( M_PI / (2.0 * nu * root_term) );
    eta = root_term + Math.log( z / (1.0 + root_term) );
    if ( z < 1.0 / GSL_ROOT3_DBL_EPSILON )
    {
        ex_arg = nu * (z - eta);
    }
    else
    {
        ex_arg = 0.5 * nu / z * (1.0 + 1.0 / (12.0 * z * z));
    }
    ex_result = gsl_sf_exp_e( ex_arg );
    t = 1.0 / root_term;
    tpow[0] = 1.0;
    for ( let i = 1; i <= 16 - 1; i++ )
    {
        tpow[i] = t * tpow[i-1];
    }
    sum = 1.0 - debye_u1( tpow ) / nu + debye_u2( tpow ) / (nu * nu) - debye_u3( tpow ) / (nu * nu * nu)
          + debye_u4( tpow ) / (nu * nu * nu * nu) - debye_u5( tpow ) / (nu * nu * nu * nu * nu);
    r.val = pre * ex_result.val * sum;
    r.err = pre * ex_result.err * Math.abs( sum );
    r.err = r.err + pre * ex_result.val / (nu * nu * nu * nu * nu * nu);
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );

    return r;

} // gsl_sf_bessel_Knu_scaled_asymp_unif_e

// ----------------------------------------------------------------------------

// Evaluate J_mu(x),J_{mu+1}(x) and Y_mu(x),Y_{mu+1}(x)  for |mu| < 1/2
//
export function gsl_sf_bessel_JY_mu_restricted( mu, x, Jmu, Jmup1, Ymu, Ymup1 )
{
  
    if ( x < 0.0 || Math.abs( mu ) > 0.5 )
    {
        //Jmu.val   = 0.0;
        //Jmu.err   = 0.0;
        //Jmup1.val = 0.0;
        //Jmup1.err = 0.0;
        //Ymu.val   = 0.0;
        //Ymu.err   = 0.0;
        //Ymup1.val = 0.0;
        //Ymup1.err = 0.0;
        throw "SF.DomainException";
    }
    else if ( x == 0.0 )
    {
        //IF (mu = 0.0)
        //    Jmu.val = 1.0;
        //    Jmu.err = 0.0;
        //ELSE
        //    Jmu.val = 0.0;
        //    Jmu.err = 0.0;
        //END IF;
        //Jmup1.val = 0.0;
        //Jmup1.err = 0.0;
        //Ymu.val   = 0.0;
        //Ymu.err   = 0.0;
        //Ymup1.val = 0.0;
        //Ymup1.err = 0.0;
        throw "SF.DomainException";
    }
    else
    {
        if ( x < 2.0 )
        {
            // Use Taylor series for J and the Temme series for Y.
            // The Taylor series for J requires nu > 0, so we shift
            // up one and use the recursion relation to get Jmu, in
            // case mu < 0.
            //
            var Jmup2 = { val: 0.0, err: 0.0 };
            var c = 0.0;

            Jmup1 = gsl_sf_bessel_IJ_taylor_e( mu + 1.0, x, -1, 100, GSL_DBL_EPSILON );
            Jmup2 = gsl_sf_bessel_IJ_taylor_e( mu + 2.0, x, -1, 100, GSL_DBL_EPSILON );
            c = 2.0 * (mu + 1.0) / x;
            Jmu.val = c * Jmup1.val - Jmup2.val;
            Jmu.err = c * Jmup1.err + Jmup2.err;
            Jmu.err = Jmu.err + 2.0 * GSL_DBL_EPSILON * Math.abs( Jmu.val );
            gsl_sf_bessel_Y_temme(mu, x, Ymu, Ymup1);
        }
        else if ( x < 1000.0 )
        {
            var P       = 0.0;
            var Q       = 0.0;
            var J_ratio = 0.0;
            var J_sgn   = 0.0;
            var gamma   = 0.0;
            var Jprime_J_ratio = 0.0;
            var r1 = { ratio: 0.0, sgn: 0.0 };
            var r2 = { P: 0.0, Q: 0.0 };

            r1 = gsl_sf_bessel_J_CF1( mu, x ); //, J_ratio, J_sgn );
            J_ratio = r1.ratio;
            J_sgn = r1.sgn;
            r2 = gsl_sf_bessel_JY_steed_CF2( mu, x ); //, P, Q );
            P = r2.P;
            Q = r2.Q;
            Jprime_J_ratio = mu / x - J_ratio;
            gamma = (P - Jprime_J_ratio) / Q;
            Jmu.val = J_sgn * Math.sqrt( 2.0 / (M_PI * x) / (Q + gamma * (P - Jprime_J_ratio)) );
            Jmu.err = 4.0 * GSL_DBL_EPSILON * Math.abs( Jmu.val );
            Jmup1.val = J_ratio * Jmu.val;
            Jmup1.err = Math.abs( J_ratio ) * Jmu.err;
            Ymu.val = gamma * Jmu.val;
            Ymu.err = Math.abs( gamma ) * Jmu.err;
            Ymup1.val = Ymu.val * (mu / x - P - Q / gamma);
            Ymup1.err = Ymu.err * Math.abs( mu / x - P - Q / gamma ) + 4.0 * GSL_DBL_EPSILON * Math.abs( Ymup1.val );
        }
        else
        {
            // Use asymptotics for large argument.
            //
            Jmu   = gsl_sf_bessel_Jnu_asympx_e( mu,       x );
            Jmup1 = gsl_sf_bessel_Jnu_asympx_e( mu + 1.0, x );
            Ymu   = gsl_sf_bessel_Ynu_asympx_e( mu,       x );
            Ymup1 = gsl_sf_bessel_Ynu_asympx_e( mu + 1.0, x );
        }
    }

} // gsl_sf_bessel_JY_mu_restricted;

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_J_CF1(nu, x) //, ratio, sgn)
{
    const maxiter = 10000;
    var n = 1;
    const RECUR_BIG   = GSL_SQRT_DBL_MAX;
    const RECUR_SMALL = GSL_SQRT_DBL_MIN;
    var Anm2 = 0.0;
    var Bnm2 = 0.0;
    var Anm1 = 0.0;
    var Bnm1 = 0.0;
    var a1   = 0.0;
    var An   = 0.0;
    var Bn   = 0.0;
    var ax   = 0.0;
    var fn   = 0.0;
    var dn   = 0.0;
    var s    = 0.0;
    var del  = 0.0;
    var old_fn = 0.0;
    var ratio = 0.0;
    var sgn = 0.0;
    var r = {};

    Anm2 = 1.0;
    Bnm2 = 0.0;
    Anm1 = 0.0;
    Bnm1 = 1.0;
    a1 = x / (2.0 * (nu + 1.0));
    An = Anm1 + a1 * Anm2;
    Bn = Bnm1 + a1 * Bnm2;
    fn = An / Bn;
    dn = a1;
    s  = 1.0;

    while (n < maxiter)
    {
        n = n + 1;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        ax = -x * x / (4.0 * (nu + (n) - 1.0) * (nu + (n)));
        An = Anm1 + ax * Anm2;
        Bn = Bnm1 + ax * Bnm2;
       
        if (Math.abs(An) > RECUR_BIG || Math.abs(Bn) > RECUR_BIG)
        {
            An   = An   / RECUR_BIG;
            Bn   = Bn   / RECUR_BIG;
            Anm1 = Anm1 / RECUR_BIG;
            Bnm1 = Bnm1 / RECUR_BIG;
            Anm2 = Anm2 / RECUR_BIG;
        }
        else if (Math.abs(An) < RECUR_SMALL || Math.abs(Bn) < RECUR_SMALL)
        {
            An   = An   / RECUR_SMALL;
            Bn   = Bn   / RECUR_SMALL;
            Anm1 = Anm1 / RECUR_SMALL;
            Bnm1 = Bnm1 / RECUR_SMALL;
            Anm2 = Anm2 / RECUR_SMALL;
            Bnm2 = Bnm2 / RECUR_SMALL;
        }
       
        old_fn = fn;
        fn     = An / Bn;
        del    = old_fn / fn;
       
        dn = 1.0 / (2.0 * (nu + (n)) / x - dn);
        if (dn < 0.0)
        {
            s = -s;
        }

        if (Math.abs(del - 1.0) < 2.0 * GSL_DBL_EPSILON) break;
    }
  
    // FIXME: we should return an error term here as well, because the
    // error from this recurrence affects the overall error estimate.
  
    ratio = fn;
    sgn   = s;
    r.ratio = ratio;
    r.sgn = sgn;

    if (n >= maxiter)
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // gsl_sf_bessel_J_CF1

// ----------------------------------------------------------------------------

// Evaluate the continued fraction CF1 for I_{nu+1}/I_nu
// using Gautschi (Euler) equivalent series.
//
export function gsl_sf_bessel_I_CF1_ser(nu, x)
{
    const maxk = 20000;
    var k    = 0;
    var tk   = 1.0;
    var sum  = 1.0;
    var rhok = 0.0;
    var ak   = 0.0;
  
    k = 1;  
    while (k <= maxk - 1)
    {
        ak   = 0.25 * (x / (nu + (k))) * x / (nu + (k) + 1.0);
        rhok = -ak * (1.0 + rhok) / (1.0 + ak * (1.0 + rhok));
        tk   = tk * rhok;
        sum  = sum + tk;
        if (Math.abs(tk / sum) < GSL_DBL_EPSILON) break;
        k = k + 1;
    }
  
    if (k >= maxk)
    {
        throw "SF.MaxIterationsException";
    }
  
    return x / (2.0 * (nu + 1.0)) * sum;

} // gsl_sf_bessel_I_CF1_ser

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_JY_steed_CF2( nu, x )
{

    const max_iter = 10000;
    const SMALL    = 1.0e-100;
  
    var i = 1;
  
    var x_inv = 0.0;
    var a     = 0.0;
    var p0    = 0.0;
    var q0    = 0.0;
    var br    = 0.0;
    var bi    = 0.0;
    var fact  = 0.0;
    var cr    = 0.0;
    var ci    = 0.0;
    var den   = 0.0;
    var dr    = 0.0;
    var di    = 0.0;
    var dlr   = 0.0;
    var dli   = 0.0;
    var temp  = 0.0;

    var P = 0.0;
    var Q = 0.0;

    var r = { P: 0.0, Q: 0.0 };

    x_inv = 1.0 / x;
    a = 0.25 - nu * nu;
    p0 = -0.5 * x_inv;
    q0 = 1.0;
    br = 2.0 * x;
    bi = 2.0;
    fact = a * x_inv / (p0 * p0 + q0 * q0);
    cr = br + q0 * fact;
    ci = bi + p0 * fact;
    den = br * br + bi * bi;
    dr = br / den;
    di = -bi / den;
    dlr = cr * dr - ci * di;
    dli = cr * di + ci * dr;
    temp = p0 * dlr - q0 * dli;
    q0 = p0 * dli + q0 * dlr;
    p0 = temp;
    i = 2;
    while (i <= max_iter)
    {
        a  = a + (2 * (i - 1));
        bi = bi + 2.0;
        dr = a * dr + br;
        di = a * di + bi;
        if (Math.abs(dr) + Math.abs(di) < SMALL)
        {
            dr = SMALL;
        }
        fact = a / (cr * cr + ci * ci);
        cr = br + cr * fact;
        ci = bi - ci * fact;
        if (Math.abs(cr) + Math.abs(ci) < SMALL)
        {
            cr = SMALL;
        }
        den = dr * dr + di * di;
        dr =  dr / den;
        di = -di / den;
        dlr = cr * dr - ci * di;
        dli = cr * di + ci * dr;
        temp = p0 * dlr - q0 * dli;
        q0 = p0 * dli + q0 * dlr;
        p0 = temp;
        if (Math.abs(dlr - 1.0) + Math.abs(dli) < GSL_DBL_EPSILON) break;
        i = i + 1;
    }
  
    P = p0;
    Q = q0;
  
    if (i >= max_iter)
    {
        throw "SF.MaxIterationsException";
    }

    r.P = P;
    r.Q = Q;
    return r;

} // gsl_sf_bessel_JY_steed_CF2

// ----------------------------------------------------------------------------

// Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
// to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
//
// This is unstable for small x; x > 2 is a good cutoff.
// Also requires |nu| < 1/2.
//
export function gsl_sf_bessel_K_scaled_steed_temme_CF2( nu, x )
//          K_nu: IN OUT LONG_FLOAT; K_nup1: IN OUT LONG_FLOAT; Kp_nu: IN OUT LONG_FLOAT) IS
{
    const maxiter = 10000;

    var r = { K_nu: 0.0, K_nup1: 0.0, Kp_nu: 0.0 };
  
    var i = 1;
    var bi    = 0.0;
    var di    = 0.0;
    var delhi = 0.0;
    var hi    = 0.0;
    var qi0   = 0.0;
    var qip1  = 0.0;
    var ai    = 0.0;
    var a1    = 0.0;
    var ci    = 0.0;
    var Qi    = 0.0;
    var s     = 0.0;
    var dels  = 0.0;
    var tmp   = 0.0;
  
    bi = 2.0 * (1.0 + x);
    di = 1.0 / bi;
    delhi = di;
    hi    = di;
  
    qi0  = 0.0;
    qip1 = 1.0;
  
    ai = -(0.25 - nu * nu);
    a1 = ai;
    ci = -ai;
    Qi = -ai;
  
    s = 1.0 + Qi * delhi;
    i = 2;
    while ( i <= maxiter )
    {
        ai = ai - 2.0 * (i - 1);
        ci  = -ai * ci / (i);
        tmp  = (qi0 - bi * qip1) / ai;
        qi0   = qip1;
        qip1 = tmp;
        Qi = Qi + ci * qip1;
        bi = bi + 2.0;
        di = 1.0 / (bi + ai * di);
        delhi = (bi * di - 1.0) * delhi;
        hi = hi + delhi;
        dels = Qi * delhi;
        s = s + dels;
        if ( Math.abs( dels / s ) < GSL_DBL_EPSILON ) break;
        i = i + 1;
    }
    
    hi = -hi * a1;
    
    r.K_nu   = Math.sqrt( M_PI / (2.0 * x) ) / s;
    r.K_nup1 = r.K_nu * (nu + x + 0.5 - hi) / x;
    r.Kp_nu  = -r.K_nup1 + nu / x * r.K_nu;

    if ( i >= maxiter )
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // gsl_sf_bessel_K_scaled_steed_temme_CF2;

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_cos_pi4_e(y, eps)
{
    var sy = 0.0;
    var cy = 0.0;
    var s  = 0.0;
    var d  = 0.0;
    var abs_sum = 0.0;
    var seps = 0.0;
    var ceps = 0.0;
    var e2 = 0.0;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    sy = Math.sin(y);
    cy = Math.cos(y);
    s  = sy + cy;
    d  = sy - cy;
    abs_sum = Math.abs(cy) + Math.abs(sy);

    if (Math.abs(eps) < GSL_ROOT5_DBL_EPSILON)
    {
        e2 = eps * eps;
        seps = eps * (1.0 - e2 / 6.0 * (1.0 - e2 / 20.0));
        ceps = 1.0 - e2 / 2.0 * (1.0 - e2 / 12.0);
    }
    else
    {
        seps = Math.sin(eps);
        ceps = Math.cos(eps);
    }
    r.val = (ceps * s - seps * d) / M_SQRT2;
    r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs(ceps) + Math.abs(seps)) * abs_sum / M_SQRT2;
   
   
    // Try to account for error in evaluation of sin(y), cos(y).
    // This is a little sticky because we don't really know
    // how the library routines are doing their argument reduction.
    // However, we will make a reasonable guess.
    // FIXME ?
    //
    if (y > 1.0 / GSL_DBL_EPSILON)
    {
        r.err = r.err * 0.5 * y;
    }
    else if (y > 1.0 / GSL_SQRT_DBL_EPSILON)
    {
        r.err = r.err * 256.0 * y * GSL_SQRT_DBL_EPSILON;
    }
   
    return r;

} // gsl_sf_bessel_cos_pi4_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_sin_pi4_e(y, eps)
{
    var sy = 0.0;
    var cy = 0.0;
    var s  = 0.0;
    var d  = 0.0;
    var abs_sum = 0.0;
    var seps = 0.0;
    var ceps = 0.0;
    var e2 = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    sy = Math.sin(y);
    cy = Math.cos(y);
    s = sy + cy;
    d = sy - cy;
    abs_sum = Math.abs(cy) + Math.abs(sy);

    if (Math.abs(eps) < GSL_ROOT5_DBL_EPSILON)
    {
        e2 = eps * eps;
        seps = eps * (1.0 - e2 / 6.0 * (1.0 - e2 / 20.0));
        ceps = 1.0 - e2 / 2.0 * (1.0 - e2 / 12.0);
    }
    else
    {
        seps = Math.sin(eps);
        ceps = Math.cos(eps);
    }
    r.val = (ceps * d + seps * s)/ M_SQRT2;
    r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs(ceps) + Math.abs(seps)) * abs_sum / M_SQRT2;
  
    // Try to account for error in evaluation of sin(y), cos(y).
    // See above.
    // FIXME ?
    //
    if (y > 1.0 / GSL_DBL_EPSILON)
    {
        r.err = r.err * 0.5 * y;
    }
    else if (y > 1.0 / GSL_SQRT_DBL_EPSILON)
    {
        r.err = r.err * 256.0 * y * GSL_SQRT_DBL_EPSILON;
    }
  
    return r;

} // gsl_sf_bessel_sin_pi4_e

// ----------------------------------------------------------------------------
// EOF SF-Bessel.mjs
