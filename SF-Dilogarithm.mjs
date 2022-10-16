// SF-Dilogarithm.mjs
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

import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";
import { M_PI }                 from "./SF-Math.mjs";
import { gsl_sf_hypot }         from "./SF-Trigonometric.mjs";
import { gsl_sf_clausen_e }     from "./SF-Clausen.mjs";
import { gsl_sf_complex_log_e } from "./SF-Logarithmic.mjs";

import { EVAL_RESULT_D }        from "./SF-Evaluate.mjs";

// Evaluate series for real dilog(x)
// Sum[ x^k / k^2, {k,1,Infinity}]
//
// Converges rapidly for |x| < 1/2.
//
function dilog_series_1(x)
{
    const kmax = 1000;
    var sum  = 0.0;
    var term = 0.0;
    var rk   = 0.0;

    var k = 0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    sum  = x;
    term = x;
    k = 2;
    while (k <= kmax)
    {
        rk = ((k) - 1.0) / (k);
        term = term * x;
        term = term * rk * rk;
        sum  = sum + term;
        if (Math.abs(term / sum) < GSL_DBL_EPSILON)
        {
            break;
        }
        k = k + 1;
    }
   
    r.val = sum;
    r.err = 2.0 * Math.abs(term);
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
   
    if (k >= kmax)
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // dilog_series_1

// ----------------------------------------------------------------------------

// Compute the associated series
//
//   sum_{k=1}{infty} r^k / (k^2 (k+1))
//
// This is a series which appears in the one-step accelerated
// method, which splits out one elementary function from the
// full definition of Li_2(x). See below.
//
function series_2(r)
{
    const kmax = 100;

    var ds  = 0.0;
    var rk  = 0.0;
    var sum = 0.0;

    var s = { val: 0.0, err: 0.0 }; // Result;

    rk = r;
    sum = 0.5 * r;
    for (let k = 2; k <= 10 - 1; k++)
    {
        rk = rk * r;
        ds = rk / (k * k * (k + 1));
        sum = sum + ds;
    }
    for (let k = 10; k <= kmax - 1; k++)
    {
        rk = rk * r;
        ds = rk / (k * k * (k + 1));
        sum = sum + ds;
        if (Math.abs(ds / sum) < 0.5 * GSL_DBL_EPSILON)
        {
            break;
        }
    }
   
    s.val = sum;
    s.err = 2.0 * (kmax) * GSL_DBL_EPSILON * Math.abs(sum);
   
    return s;

} // series_2

// ----------------------------------------------------------------------------

// Compute Li_2(x) using the accelerated series representation.
//
// Li_2(x) = 1 + (1-x)ln(1-x)/x + series_2(x)
//
// assumes: -1 < x < 1
//
function dilog_series_2(x)
{
    const c3 = 1.0 / 3.0;
    const c4 = 1.0 / 4.0;
    const c5 = 1.0 / 5.0;
    const c6 = 1.0 / 6.0;
    const c7 = 1.0 / 7.0;
    const c8 = 1.0 / 8.0;
    var t   = 0.0;
    var t68 = 0.0;
    var t38 = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    r = series_2(x);
    if (x > 0.01)
    {
        t = (1.0 - x) * Math.log(1.0 - x) / x;
    }
    else
    {
        t68 = c6 + x * (c7 + x * c8);
        t38 = c3 + x * (c4 + x * (c5 + x * t68));
        t = (x - 1.0) * (1.0 + x * (0.5 + x * t38));
    }
    r.val = r.val + 1.0 + t;
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(t);

    return r;

} // dilog_series_2

// ----------------------------------------------------------------------------

// Calculates Li_2(x) for real x. Assumes x >= 0.0.
function dilog_xge0(x)
{
    var t1       = 0.0;
    var t2       = 0.0;
    var t3       = 0.0;
    var log_x    = 0.0;
    var log_term = 0.0;
    var lne      = 0.0;
    var eps      = 0.0;
    var c0       = 0.0;
    var c1       = 0.0;
    var c2       = 0.0;
    var c3       = 0.0;
    var c4       = 0.0;
    var c5       = 0.0;
    var c6       = 0.0;
    var c7       = 0.0;
    var c8       = 0.0;

    var s = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x > 2.0)
    {
        s = dilog_series_2(1.0 / x);
        log_x = Math.log(x);
        t1 = M_PI * M_PI / 3.0;
        t2 = s.val;
        t3 = 0.5 * log_x * log_x;
        r.val = t1 - t2 - t3;
        r.err = GSL_DBL_EPSILON * Math.abs(log_x) + s.err;
        r.err = r.err + GSL_DBL_EPSILON * (Math.abs(t1) + Math.abs(t2) + Math.abs(t3));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x > 1.01)
    {
        s = dilog_series_2(1.0 - 1.0 / x);
        log_x    = Math.log(x);
        log_term = log_x * (Math.log(1.0 - 1.0 / x) + 0.5 * log_x);
        t1 = M_PI * M_PI / 6.0;
        t2 = s.val;
        t3 = log_term;
        r.val = t1 + t2 - t3;
        r.err = GSL_DBL_EPSILON * Math.abs(log_x) + s.err;
        r.err = r.err + GSL_DBL_EPSILON * (Math.abs(t1) + Math.abs(t2) + Math.abs(t3));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x > 1.0)
    {
        // series around x = 1.0
        eps = x - 1.0;
        lne = Math.log(eps);
        c0 = M_PI * M_PI / 6.0;
        c1 =   1.0 - lne;
        c2 = -(1.0 - 2.0 * lne) / 4.0;
        c3 =  (1.0 - 3.0 * lne) / 9.0;
        c4 = -(1.0 - 4.0 * lne) / 16.0;
        c5 =  (1.0 - 5.0 * lne) / 25.0;
        c6 = -(1.0 - 6.0 * lne) / 36.0;
        c7 =  (1.0 - 7.0 * lne) / 49.0;
        c8 = -(1.0 - 8.0 * lne) / 64.0;
        r.val = c0 + eps * (c1 + eps * (c2 + eps * (c3 + eps * (c4 + eps * (c5 + eps * (c6 + eps * (c7 + eps * c8)))))));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x == 1.0)
    {
        r.val = M_PI * M_PI / 6.0;
        r.err = 2.0 * GSL_DBL_EPSILON * M_PI * M_PI / 6.0;
    }
    else if (x > 0.5)
    {
        s = dilog_series_2(1.0 - x);
        log_x = Math.log(x);
        t1 = M_PI * M_PI / 6.0;
        t2 = s.val;
        t3 = log_x * Math.log(1.0 - x);
        r.val = t1 - t2 - t3;
        r.err = GSL_DBL_EPSILON * Math.abs(log_x) + s.err;
        r.err = r.err + GSL_DBL_EPSILON * (Math.abs(t1) + Math.abs(t2) + Math.abs(t3));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x > 0.25)
    {
        r = dilog_series_2(x);
    }
    else if (x > 0.0)
    {
        r = dilog_series_1(x);
    }
    else
    {
        // x == 0.0
        r.val = 0.0;
        r.err = 0.0;
    }

    return r;

} // dilog_xge0

// ----------------------------------------------------------------------------

// Evaluate the series representation for Li2(z):
//
//   Li2(z) = Sum[ |z|^k / k^2 Exp[i k arg(z)], {k,1,Infinity}]
//   |z|    = r
//   arg(z) = theta
//   
// Assumes 0 < r < 1.
// It is used only for small r.
//
function dilogc_series_1( r, x, y, real_result, imag_result )
{
    var kmax = 0;

    var cos_theta = 0.0;
    var sin_theta = 0.0;
    var alpha     = 0.0;
    var beta      = 0.0;
    var ck        = 0.0;
    var sk        = 0.0;
    var rk        = 0.0;
    var real_sum  = 0.0;
    var imag_sum  = 0.0;
    var dr        = 0.0;
    var di        = 0.0;
    var ck_tmp    = 0.0;

    cos_theta = x / r;
    sin_theta = y / r;
    alpha = 1.0 - cos_theta;
    beta  = sin_theta;
    ck = cos_theta;
    sk = sin_theta;
    rk = r;
    real_sum = r * ck;
    imag_sum = r * sk;
    kmax = 50 + Math.trunc( 22.0 / (-Math.log( r )) ); // tuned for double-precision

    for ( let k = 2; k <= kmax - 1; k++ )
    {
        ck_tmp = ck;
        ck = ck - (alpha * ck + beta * sk);
        sk = sk - (alpha * sk - beta * ck_tmp);
        rk = rk * r;
        dr = rk / (k * k) * ck;
        di = rk / (k * k) * sk;
        real_sum = real_sum + dr;
        imag_sum = imag_sum + di;
        if ( Math.abs( (dr * dr + di * di) / (real_sum * real_sum + imag_sum * imag_sum) ) < GSL_DBL_EPSILON * GSL_DBL_EPSILON ) break;
    }
   
    real_result.val = real_sum;
    real_result.err = 2.0 * (kmax) * GSL_DBL_EPSILON * Math.abs( real_sum );
    imag_result.val = imag_sum;
    imag_result.err = 2.0 * (kmax) * GSL_DBL_EPSILON * Math.abs( imag_sum );

} // dilogc_series_1

// ----------------------------------------------------------------------------

// Compute
//
//   sum_{k=1}{infty} z^k / (k^2 (k+1))
//
// This is a series which appears in the one-step accelerated
// method, which splits out one elementary function from the
// full definition of Li_2.
//
function series_2_c( r, x, y, sum_re, sum_im )
{
    var kmax = 0;

    var cos_theta = 0.0;
    var sin_theta = 0.0;
    var alpha     = 0.0;
    var beta      = 0.0;
    var ck        = 0.0;
    var sk        = 0.0;
    var rk        = 0.0;
    var real_sum  = 0.0;
    var imag_sum  = 0.0;
    var dr        = 0.0;
    var di        = 0.0;
    var ck_tmp    = 0.0;

    cos_theta = x / r;
    sin_theta = y / r;
    alpha = 1.0 - cos_theta;
    beta  = sin_theta;
    ck = cos_theta;
    sk = sin_theta;
    rk = r;
    real_sum = 0.5 * r * ck;
    imag_sum = 0.5 * r * sk;
    kmax = 30 + Math.trunc( 18.0 / (-Math.log( r )) ); // tuned for double-precision

    for ( let k = 2; k <= kmax - 1; k++ )
    {
        ck_tmp = ck;
        ck = ck - (alpha * ck + beta * sk);
        sk = sk - (alpha * sk - beta * ck_tmp);
        rk = rk * r;
        dr = rk / (k * k * (k + 1)) * ck;
        di = rk / (k * k * (k + 1)) * sk;
        real_sum = real_sum + dr;
        imag_sum = imag_sum + di;
        if ( Math.abs( (dr * dr + di * di) / (real_sum * real_sum + imag_sum * imag_sum) ) < GSL_DBL_EPSILON * GSL_DBL_EPSILON ) break;
    }
   
    sum_re.val = real_sum;
    sum_re.err = 2.0 * (kmax) * GSL_DBL_EPSILON * Math.abs( real_sum );
    sum_im.val = imag_sum;
    sum_im.err = 2.0 * (kmax) * GSL_DBL_EPSILON * Math.abs( imag_sum );

} // series_2_c

// ----------------------------------------------------------------------------

// Compute Li_2(z) using the one-step accelerated series.
//
// Li_2(z) = 1 + (1-z)ln(1-z)/z + series_2_c(z)
//
// z = r exp(i theta)
// assumes: r < 1
// assumes: r > epsilon, so that we take no special care with log(1-z)
//
function dilogc_series_2( r, x, y, real_dl, imag_dl )
{
    var t_x = 0.0;
    var t_y = 0.0;
    var r_x = 0.0;
    var r_y = 0.0;

    var sum_re       = { val: 0.0, err: 0.0 }; // Result;
    var sum_im       = { val: 0.0, err: 0.0 }; // Result;
    var ln_omz_r     = { val: 0.0, err: 0.0 }; // Result;
    var ln_omz_theta = { val: 0.0, err: 0.0 }; // Result;

    if ( r == 0.0 )
    {
        real_dl.val = 0.0;
        imag_dl.val = 0.0;
        real_dl.err = 0.0;
        imag_dl.err = 0.0;
    }
    else
    {
        series_2_c( r, x, y, sum_re, sum_im );
       
        // t = ln(1-z)/z
        gsl_sf_complex_log_e( 1.0 - x, -y, ln_omz_r, ln_omz_theta );
        t_x = ( ln_omz_r.val * x + ln_omz_theta.val * y) / (r * r);
        t_y = (-ln_omz_r.val * y + ln_omz_theta.val * x) / (r * r);
       
        // r = (1-z) ln(1-z)/z
        r_x = (1.0 - x) * t_x + y * t_y;
        r_y = (1.0 - x) * t_y - y * t_x;
       
        real_dl.val = sum_re.val + r_x + 1.0;
        imag_dl.val = sum_im.val + r_y;
        real_dl.err = sum_re.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( real_dl.val ) + Math.abs( r_x ));
        imag_dl.err = sum_im.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( imag_dl.val ) + Math.abs( r_y ));
    }

} // dilogc_series_2

// ----------------------------------------------------------------------------

// Evaluate a series for Li_2(z) when |z| is near 1.
// This is uniformly good away from z=1.
//
//   Li_2(z) = Sum[ a^n/n! H_n(theta), {n, 0, Infinity}]
//
// where
//   H_n(theta) = Sum[ e^(i m theta) m^n / m^2, {m, 1, Infinity}]
//   a = ln(r)
//
//  H_0(t) = Gl_2(t) + i Cl_2(t)
//  H_1(t) = 1/2 ln(2(1-c)) + I atan2(-s, 1-c)
//  H_2(t) = -1/2 + I/2 s/(1-c)
//  H_3(t) = -1/2 /(1-c)
//  H_4(t) = -I/2 s/(1-c)^2
//  H_5(t) = 1/2 (2 + c)/(1-c)^2
//  H_6(t) = I/2 s/(1-c)^5 (8(1-c) - s^2 (3 + c))
//
function dilogc_series_3( r, x, y, real_result, imag_result )
{
    var t         = 0.0;
    var theta     = 0.0;
    var cos_theta = 0.0;
    var sin_theta = 0.0;
    var a         = 0.0;
    var omc       = 0.0;
    var omc2      = 0.0;
    var an        = 0.0;
    var nfact     = 0.0;
    var sum_re    = 0.0;
    var sum_im    = 0.0;

    var H_re = []; //[7];
    var H_im = []; //[7];

    var Him0 = { val: 0.0, err: 0.0 }; // Result;

    theta = Math.atan2( y, x ); // atan2
    cos_theta = x / r;
    sin_theta = y / r;
    a = Math.log( r );
    omc = 1.0 - cos_theta;
    omc2 = omc * omc;

    H_re[0] = M_PI * M_PI / 6.0 + 0.25 * (theta * theta - 2.0 * M_PI * Math.abs( theta ));
    Him0 = gsl_sf_clausen_e( theta );
    H_im[0] = Him0.val;
   
    H_re[1] = -0.5 * Math.log( 2.0 * omc );
    H_im[1] = -Math.atan2( -sin_theta, omc ); // atan2
   
    H_re[2] = -0.5;
    H_im[2] = 0.5 * sin_theta / omc;
   
    H_re[3] = -0.5 / omc;
    H_im[3] = 0.0;
   
    H_re[4] = 0.0;
    H_im[4] = -0.5 * sin_theta / omc2;
   
    H_re[5] = 0.5 * (2.0 + cos_theta) / omc2;
    H_im[5] = 0.0;
   
    H_re[6] = 0.0;
    H_im[6] = 0.5 * sin_theta / (omc2 * omc2 * omc) * (8.0 * omc - sin_theta * sin_theta * (3.0 + cos_theta));
   
    sum_re = H_re[0];
    sum_im = H_im[0];
    an = 1.0;
    nfact = 1.0;
    for ( let n = 1; n <= 6; n++ )
    {
        an = an * a;
        nfact = nfact * (n);
        t = an / nfact;
        sum_re = sum_re + t * H_re[n];
        sum_im = sum_im + t * H_im[n];
    }
   
    real_result.val = sum_re;
    real_result.err = 2.0 * 6.0 * GSL_DBL_EPSILON * Math.abs( sum_re ) + Math.abs( an / nfact );
    imag_result.val = sum_im;
    imag_result.err = 2.0 * 6.0 * GSL_DBL_EPSILON * Math.abs( sum_im ) + Him0.err + Math.abs( an / nfact );

} // dilogc_series_3

// ----------------------------------------------------------------------------

// Calculate complex dilogarithm Li_2(z) in the fundamental region,
// which we take to be the intersection of the unit disk with the
// half-space x < MAGIC_SPLIT_VALUE. It turns out that 0.732 is a
// nice choice for MAGIC_SPLIT_VALUE since then points mapped out
// of the x > MAGIC_SPLIT_VALUE region and into another part of the
// unit disk are bounded in radius by MAGIC_SPLIT_VALUE itself.
//
// If |z| < 0.98 we use a direct series summation. Otherwise z is very
// near the unit circle, and the series_2 expansion is used; see above.
// Because the fundamental region is bounded away from z = 1, this
// works well.
//
function dilogc_fundamental( r, x, y, real_dl, imag_dl )
{

    if ( r > 0.98 )
    {  
        dilogc_series_3( r, x, y, real_dl, imag_dl );
    }
    else if ( r > 0.25 )
    {
        dilogc_series_2( r, x, y, real_dl, imag_dl );
    }
    else
    {
        dilogc_series_1( r, x, y, real_dl, imag_dl );
    }

} // dilogc_fundamental

// ----------------------------------------------------------------------------

// Compute Li_2(z) for z in the unit disk, |z| < 1. If z is outside
// the fundamental region, which means that it is too close to z = 1,
// then it is reflected into the fundamental region using the identity
//
//   Li2(z) = -Li2(1-z) + zeta(2) - ln(z) ln(1-z).
//
function dilogc_unitdisk( x, y, real_dl, imag_dl )
{
    const MAGIC_SPLIT_VALUE = 0.732;
    var zeta2  = M_PI * M_PI / 6.0;
    var r      = 0.0;
    var x_tmp  = 0.0;
    var y_tmp  = 0.0;
    var r_tmp  = 0.0;
    var lnz    = 0.0;
    var lnomz  = 0.0;
    var argz   = 0.0;
    var argomz = 0.0;

    var result_re_tmp = { val: 0.0, err: 0.0 }; // Result;
    var result_im_tmp = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_hypot( x, y );
    //r = hypot(x, y);

    if ( x > MAGIC_SPLIT_VALUE )
    {
        // Reflect away from z = 1 if we are too close. The magic value
        // insures that the reflected value of the radius satisfies the
        // related inequality r_tmp < MAGIC_SPLIT_VALUE.
        //
        x_tmp = 1.0 - x;
        y_tmp =     - y;
        r_tmp = gsl_sf_hypot( x_tmp, y_tmp );
        //r_tmp = hypot(x_tmp, y_tmp);
        // const double cos_theta_tmp = x_tmp/r_tmp;
        // const double sin_theta_tmp = y_tmp/r_tmp;
       
        dilogc_fundamental( r_tmp, x_tmp, y_tmp, result_re_tmp, result_im_tmp );
       
        lnz    =  Math.log( r );              //  log(|z|)
        lnomz  =  Math.log( r_tmp );          //  log(|1-z|)
        argz   =  Math.atan2( y, x );         //  arg(z) assuming principal branch
        argomz =  Math.atan2( y_tmp, x_tmp ); //  arg(1-z)
        real_dl.val = -result_re_tmp.val + zeta2 - lnz * lnomz + argz * argomz;
        real_dl.err =  result_re_tmp.err;
        real_dl.err =  real_dl.err + 2.0 * GSL_DBL_EPSILON * (zeta2 + Math.abs( lnz * lnomz ) + Math.abs( argz * argomz ));
        imag_dl.val = -result_im_tmp.val - argz * lnomz - argomz * lnz;
        imag_dl.err =  result_im_tmp.err;
        imag_dl.err =  imag_dl.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( argz * lnomz ) + Math.abs( argomz * lnz ));
    }
    else
    {
        dilogc_fundamental( r, x, y, real_dl, imag_dl );
    }

} // dilogc_unitdisk

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_dilog_e(x)
{
    var d1 = { val: 0.0, err: 0.0 }; // Result;
    var d2 = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (x >= 0.0)
    {
        r = dilog_xge0(x);
    }
    else
    {
        d1 = dilog_xge0(   -x);
        d2 = dilog_xge0(x * x);
        r.val = -d1.val + 0.5 * d2.val;
        r.err =  d1.err + 0.5 * d2.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_dilog_e

// ----------------------------------------------------------------------------

export function gsl_sf_complex_dilog_xy_e( x, y, real_dl, imag_dl )
{
   
    const zeta2 = M_PI * M_PI / 6.0;
    const r2    = x * x + y * y;
    var theta = 0.0;
    var term1 = 0.0;
    var term2 = 0.0;
    var r     = 0.0;
    var x_tmp = 0.0;
    var y_tmp = 0.0;
    var theta_abs = 0.0;
    var theta_sgn = 0.0;
    var ln_minusz_re = 0.0;
    var ln_minusz_im = 0.0;
    var lmz2_re = 0.0;
    var lmz2_im = 0.0;

    var result_re_tmp = { val: 0.0, err: 0.0 }; // Result;
    var result_im_tmp = { val: 0.0, err: 0.0 }; // Result;

    var t = { val: 0.0, err: 0.0 };

    if ( y == 0.0 )
    {
        if ( x >= 1.0 )
        {
            imag_dl.val = -M_PI * Math.log( x );
            imag_dl.err = 2.0 * GSL_DBL_EPSILON * Math.abs( imag_dl.val );
        }
        else
        {
            imag_dl.val = 0.0;
            imag_dl.err = 0.0;
        }
        t = gsl_sf_dilog_e( x );
        real_dl.val = t.val;
        real_dl.err = t.err;
    }
    else if ( Math.abs( r2 - 1.0 ) < GSL_DBL_EPSILON )
    {
        // Lewin A.2.4.1 and A.2.4.2
        theta = Math.atan2( y, x );
        term1 = theta * theta / 4.0;
        term2 = M_PI * Math.abs( theta ) / 2.0;
        real_dl.val = zeta2 + term1 - term2;
        real_dl.err = 2.0 * GSL_DBL_EPSILON * (zeta2 + term1 + term2);
        imag_dl = gsl_sf_clausen_e( theta );
    }
    else if ( r2 < 1.0 )
    {
        dilogc_unitdisk( x, y, real_dl, imag_dl );
    }
    else
    {
        // Reduce argument to unit disk.
        r = Math.sqrt( r2 );
        x_tmp =  x / r2;
        y_tmp = -y / r2;
        // const double r_tmp = 1.0/r;
       
        dilogc_unitdisk( x_tmp, y_tmp, result_re_tmp, result_im_tmp );
       
        // Unwind the inversion.
        //
        //  Li_2(z) + Li_2(1/z) = -zeta(2) - 1/2 ln(-z)^2
        //
        theta = Math.atan2( y, x );
        theta_abs = Math.abs( theta );
        if ( theta < 0.0 )
        {
            theta_sgn = -1.0;
        }
        else
        {
            theta_sgn =  1.0;
        }
        ln_minusz_re = Math.log( r );
        ln_minusz_im = theta_sgn * (theta_abs - M_PI);
        lmz2_re = ln_minusz_re * ln_minusz_re - ln_minusz_im * ln_minusz_im;
        lmz2_im = 2.0 * ln_minusz_re * ln_minusz_im;
        real_dl.val = -result_re_tmp.val - 0.5 * lmz2_re - zeta2;
        real_dl.err =  result_re_tmp.err + 2.0 * GSL_DBL_EPSILON * (0.5 * Math.abs( lmz2_re ) + zeta2);
        imag_dl.val = -result_im_tmp.val - 0.5 * lmz2_im;
        imag_dl.err =  result_im_tmp.err + 2.0 * GSL_DBL_EPSILON * Math.abs( lmz2_im );
    }

} // gsl_sf_complex_dilog_xy_e

// ----------------------------------------------------------------------------

export function gsl_sf_complex_dilog_e( r, theta, real_dl, imag_dl )
{

    var cos_theta = 0.0;
    var sin_theta = 0.0;
    var x = 0.0;
    var y = 0.0;

    cos_theta = Math.cos( theta );
    sin_theta = Math.sin( theta );
    x = r * cos_theta;
    y = r * sin_theta;
    gsl_sf_complex_dilog_xy_e( x, y, real_dl, imag_dl );

} // gsl_sf_complex_dilog_e

// ----------------------------------------------------------------------------

export function gsl_sf_complex_spence_xy_e( x, y, real_sp, imag_sp )
{
    var oms_x = 0.0;
    var oms_y = 0.0;

    oms_x = 1.0 - x;
    oms_y =     - y;
    gsl_sf_complex_dilog_xy_e( oms_x, oms_y, real_sp, imag_sp );

} // gsl_sf_complex_spence_xy_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_dilog( x )
{ // gsl_sf_dilog
    return EVAL_RESULT_D( gsl_sf_dilog_e, x, "gsl_sf_dilog" );
} // gsl_sf_dilog

// ----------------------------------------------------------------------------
// EOF SF-Dilogarithm.mjs
