// SF-Hypergeometric1F1.mjs
// These routines compute the confluent hypergeometric function 1F1
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

import { M_PI }                       from "./SF-Math.mjs";
import { M_E }                        from "./SF-Math.mjs";
import { M_LN10 }                     from "./SF-Math.mjs";
import { GSL_DBL_MAX }                from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }            from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }            from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON }        from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }           from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }       from "./SF-Machine.mjs";
import { gsl_sf_lngamma_sgn_e }       from "./SF-Gamma.mjs";
import { gsl_sf_beta_e }              from "./SF-Beta.mjs";
import { gsl_sf_exp_mult_err_e }      from "./SF-Exponential.mjs";
import { gsl_sf_exp_mult_err_e10_e }  from "./SF-Exponential.mjs";
import { gsl_sf_exp_e }               from "./SF-Exponential.mjs";
import { gsl_sf_expm1_e }             from "./SF-Exponential.mjs";
import { gsl_sf_exp_err_e }           from "./SF-Exponential.mjs";
import { gsl_sf_exprel_e }            from "./SF-Exponential.mjs";
import { gsl_sf_exprel_2_e }          from "./SF-Exponential.mjs";
import { gsl_sf_exprel_n_e }          from "./SF-Exponential.mjs";
import { gsl_sf_bessel_I1_scaled_e }  from "./SF-BesselI1.mjs";
import { gsl_sf_bessel_Inu_scaled_e } from "./SF-BesselInu.mjs";
import { gsl_sf_bessel_In_scaled }    from "./SF-BesselIn.mjs";
import { gsl_sf_bessel_J1_e }         from "./SF-BesselJ1.mjs";
import { gsl_sf_lngamma_e }           from "./SF-Gamma.mjs";
import { gsl_sf_lnfact_e }            from "./SF-Gamma.mjs";
import { gsl_sf_lnbeta_e }            from "./SF-Beta.mjs";
import { gsl_sf_log_1plusx }          from "./SF-Logarithmic.mjs";
import { gsl_sf_hyperg_2F0_series_e } from "./SF-Hypergeometric.mjs";
import { gsl_sf_hyperg_1F1_series_e } from "./SF-Hypergeometric.mjs";
import { gsl_sf_hyperg_U_e10_e }      from "./SF-HypergeometricU.mjs";
import { gsl_sf_laguerre_n_e }        from "./SF-Laguerre.mjs";
import { gsl_sf_multiply_err_e }      from "./SF-Elementary.mjs";

import { EVAL_RESULT_3D }         from "./SF-Evaluate.mjs";
import { EVAL_RESULT_IID }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const H1F1_INT_THRESHOLD = 100.0 * GSL_DBL_EPSILON;


// Asymptotic result for 1F1(a, b, x)  x . -Infinity.
// Assumes b-a != neg integer and b != neg integer.
//
function hyperg_1F1_asymp_negx( a, b, x )
{
    var lg_b        = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var lg_bma      = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var F           = { val: 0.0, err: 0.0 }; // Result;
    var r           = { val: 0.0, err: 0.0 }; // Result;
    var sgn_b       = 0.0;
    var sgn_bma     = 0.0;
    var ln_term_val = 0.0;
    var ln_term_err = 0.0;
    var ln_pre_val  = 0.0;
    var ln_pre_err  = 0.0;

    lg_b   = gsl_sf_lngamma_sgn_e( b     ); //,     lg_b,   sgn_b );
    lg_bma = gsl_sf_lngamma_sgn_e( b - a ); //, lg_bma, sgn_bma );
    sgn_b   = lg_b.sign;
    sgn_bma = lg_bma.sign;
  
    F = gsl_sf_hyperg_2F0_series_e( a, 1.0 + a - b, -1.0 / x, -1 );
    if ( F.val != 0.0 )
    {
        ln_term_val = a * Math.log( -x );
        ln_term_err = 2.0 * GSL_DBL_EPSILON * (Math.abs( a ) + Math.abs( ln_term_val ));
        ln_pre_val = lg_b.val - lg_bma.val - ln_term_val;
        ln_pre_err = lg_b.err + lg_bma.err + ln_term_err;
        r = gsl_sf_exp_mult_err_e( ln_pre_val, ln_pre_err, sgn_bma * sgn_b * F.val, F.err );
    }
    else
    {
        r.val = 0.0;
        r.err = 0.0;
    }

    return r;

} // hyperg_1F1_asymp_negx

// ----------------------------------------------------------------------------

// Asymptotic result for 1F1(a, b, x)  x . +Infinity
// Assumes b != neg integer and a != neg integer
//
function hyperg_1F1_asymp_posx( a, b, x )
{
    var lg_b = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var lg_a = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var F    = { val: 0.0, err: 0.0 }; // Result;
    var r    = { val: 0.0, err: 0.0 }; // Result;
    var sgn_b       = 0.0;
    var sgn_a       = 0.0;
    var lnx         = 0.0;
    var ln_term_val = 0.0;
    var ln_term_err = 0.0;
    var ln_pre_val  = 0.0;
    var ln_pre_err  = 0.0;

    lg_b = gsl_sf_lngamma_sgn_e( b ); //, lg_b, sgn_b );
    lg_a = gsl_sf_lngamma_sgn_e( a ); //, lg_a, sgn_a );
    sgn_b = lg_b.sign;
    sgn_a = lg_a.sign;
  
    F = gsl_sf_hyperg_2F0_series_e( b - a, 1.0 - a, 1.0 / x, -1 );
    if ( F.val != 0.0 )
    {
        lnx = Math.log( x );
        ln_term_val = (a - b) * lnx;
        ln_term_err = 2.0 * GSL_DBL_EPSILON * (Math.abs( a ) + Math.abs( b )) * Math.abs( lnx )
                    + 2.0 * GSL_DBL_EPSILON * Math.abs( a - b );
        ln_pre_val = lg_b.val - lg_a.val + ln_term_val + x;
        ln_pre_err = lg_b.err + lg_a.err + ln_term_err + 2.0 * GSL_DBL_EPSILON * Math.abs( x );
        r = gsl_sf_exp_mult_err_e( ln_pre_val, ln_pre_err, sgn_a * sgn_b * F.val, F.err );
    }
    else
    {
        r.val = 0.0;
        r.err = 0.0;
    }

    return r;

} // hyperg_1F1_asymp_posx

// ----------------------------------------------------------------------------

// Asymptotic result from Slater 4.3.7 
// 
// To get the general series, write M(a,b,x) as
//
//  M(a,b,x)=sum ((a)_n/(b)_n) (x^n / n!)
//
// and expand (b)_n in inverse powers of b as follows
//
// -log(1/(b)_n) = sum_(k=0)^(n-1) log(b+k)
//             = n log(b) + sum_(k=0)^(n-1) log(1+k/b)
//
// Do a taylor expansion of the log in 1/b and sum the resulting terms
// using the standard algebraic formulas for finite sums of powers of
// k.  This should then give
//
// M(a,b,x) = sum_(n=0)^(inf) (a_n/n!) (x/b)^n * (1 - n(n-1)/(2b) 
//                          + (n-1)n(n+1)(3n-2)/(24b^2) + ...
//
// which can be summed explicitly. The trick for summing it is to take
// derivatives of sum_(i=0)^(inf) a_n*y^n/n! = (1-y)^(-a);
//
// [BJG 16/01/2007]
function hyperg_1F1_largebx( a, b, x )
{
    var y  = x / b;
    var f  = Math.exp( -a * gsl_sf_log_1plusx( -y ) ); //log1p(-y));
    var t1 = -((a * (a + 1.0)) / (2.0 * b)) * Math.pow((y / (1.0 - y)), 2.0);
    var t2 = (1.0 / (24.0 * b * b)) * ((a * (a + 1.0) * y * y) / Math.pow((1.0 - y), 4.0)) * (12.0 + 8.0 * (2.0 * a + 1.0) * y + (3.0 * a * a - a - 2.0) * y * y);
    var t3 = (-1.0 / (48.0 * b * b * b * Math.pow((1.0 - y), 6.0))) * a * ((a + 1.0) * ((y * ((a + 1.0) * (a * (y * (y * ((y * (a - 2.0) + 16.0) * (a - 1.0)) + 72.0)) + 96.0)) + 24.0) * (y ** 2)));

    var r  = { val: 0.0, err: 0.0 }; // Result;

    r.val = f * (1.0 + t1 + t2 + t3);
    r.err = 2.0 * Math.abs( f * t3 ) + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
    return r;

} // hyperg_1F1_largebx

// ----------------------------------------------------------------------------
 
// Asymptotic result for x < 2b-4a, 2b-4a large.
// [Abramowitz+Stegun, 13.5.21]
//
// assumes 0 <= x/(2b-4a) <= 1
//
function hyperg_1F1_large2bm4a( a, b, x )
{
    var eta     = 0.0;
    var cos2th  = 0.0;
    var sin2th  = 0.0;
    var th      = 0.0;
    var pre_h   = 0.0;
    var t1      = 0.0;
    var t2      = 0.0;
    var lnpre_val = 0.0;
    var lnpre_err = 0.0;
    var s1      = 0.0;
    var s2      = 0.0;
    var ser_val = 0.0;
    var ser_err = 0.0;
    var lg_b    = { val: 0.0, err: 0.0 }; // Result;
    var r       = { val: 0.0, err: 0.0 }; // Result;

    eta    = 2.0 * b - 4.0 * a;
    cos2th = x / eta;
    sin2th = 1.0 - cos2th;
    th = Math.acos( Math.sqrt( cos2th ) );
    pre_h = 0.25 * M_PI * M_PI * eta * eta * cos2th * sin2th;
    lg_b = gsl_sf_lngamma_e( b );
    t1 = 0.5 * (1.0 - b) * Math.log( 0.25 * x * eta );
    t2 = 0.25 * Math.log( pre_h );
    lnpre_val = lg_b.val + 0.5 * x + t1 - t2;
    lnpre_err = lg_b.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( 0.5 * x ) + Math.abs( t1 ) + Math.abs( t2 ));
    //#if SMALL_ANGLE
    //  const double eps = asin(sqrt(cos2th));  /* theta = pi/2 - eps */
    //  double s1 = (fmod(a, 1.0) = 0.0) ? 0.0 : sin(a*M_PI);
    //  double eta_reduc = (fmod(eta + 1, 4.0) = 0.0) ? 0.0 : fmod(eta + 1, 8.0);
    //  double phi1 = 0.25*eta_reduc*M_PI;
    //  double phi2 = 0.25*eta*(2*eps + sin(2.0*eps));
    //  double s2 = sin(phi1 - phi2);
    //#else
    s1 = Math.sin( a * M_PI );
    s2 = Math.sin( 0.25 * eta * (2.0 * th - Math.sin( 2.0 * th )) + 0.25 * M_PI );
    //#endif
    ser_val = s1 + s2;
    ser_err = 2.0 * GSL_DBL_EPSILON * (Math.abs( s1 ) + Math.abs( s2 ));
    r = gsl_sf_exp_mult_err_e( lnpre_val, lnpre_err, ser_val, ser_err );
    return r;

} // hyperg_1F1_large2bm4a

// ----------------------------------------------------------------------------

// Luke's rational approximation.
// See [Luke, Algorithms for the Computation of Mathematical Functions, p.182]
//
// Like the case of the 2F1 rational approximations, these are
// probably guaranteed to converge for x < 0, barring gross
// numerical instability in the pre-asymptotic regime.
//
function hyperg_1F1_luke( a, c, xin )
{
    const RECUR_BIG = 1.0e+50;
    const nmax = 5000;
    var n  = 0;
    var x  = 0.0;
    var x3 = 0.0;
    var t0 = 0.0;
    var t1 = 0.0;
    var t2 = 0.0;
    var F  = 0.0;
    var prec = 0.0;
    var Bnm3 = 0.0;
    var Bnm2 = 0.0;
    var Bnm1 = 0.0;
    var Anm3 = 0.0;
    var Anm2 = 0.0;
    var Anm1 = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    n = 3;
    x  = -xin;
    x3 = x * x * x;
    t0 = a / c;
    t1 = (a + 1.0) / (2.0 * c);
    t2 = (a + 2.0) / (2.0 * (c + 1.0));
    F  = 1.0;
  
    Bnm3 = 1.0;                                  // B0
    Bnm2 = 1.0 + t1 * x;                         // B1
    Bnm1 = 1.0 + t2 * x * (1.0 + t1 / 3.0 * x);  // B2
   
    Anm3 = 1.0;                                                                // A0
    Anm2 = Bnm2 - t0 * x;                                                      // A1
    Anm1 = Bnm1 - t0 * (1.0 + t2 * x) * x + t0 * t1 * (c / (c + 1.0)) * x * x; // A2
  
    while ( true )
    {
        var npam1 = (n) + a - 1.0;
        var npcm1 = (n) + c - 1.0;
        var npam2 = (n) + a - 2.0;
        var npcm2 = (n) + c - 2.0;
        var tnm1  = (2 * n - 1);
        var tnm3  = (2 * n - 3);
        var tnm5  = (2 * n - 5);
        var F1 =  ((n) - a - 2.0) / (2.0 * tnm3 * npcm1);
        var F2 =  ((n) + a) * npam1 / (4.0 * tnm1 * tnm3 * npcm2 * npcm1);
        var F3 = -npam2 * npam1 * ((n) - a - 2.0) / (8.0 * tnm3 * tnm3 * tnm5 * ((n) + c - 3.0) * npcm2 * npcm1);
        var E  = -npam1 * ((n) - c - 1.0) / (2.0 * tnm3 * npcm2 * npcm1);
        
        var An = (1.0 + F1 * x) * Anm1 + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
        var Bn = (1.0 + F1 * x) * Bnm1 + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
        var r1 = An / Bn;

        prec = Math.abs( (F - r1) / F );
        F = r1;
        
        if ( prec < GSL_DBL_EPSILON || n > nmax ) break;
        
        if ( Math.abs( An ) > RECUR_BIG || Math.abs( Bn ) > RECUR_BIG )
        {
            An   = An   / RECUR_BIG;
            Bn   = Bn   / RECUR_BIG;
            Anm1 = Anm1 / RECUR_BIG;
            Bnm1 = Bnm1 / RECUR_BIG;
            Anm2 = Anm2 / RECUR_BIG;
            Bnm2 = Bnm2 / RECUR_BIG;
            Anm3 = Anm3 / RECUR_BIG;
            Bnm3 = Bnm3 / RECUR_BIG;
        }
        else if ( Math.abs( An ) < 1.0 / RECUR_BIG || Math.abs( Bn ) < 1.0 / RECUR_BIG )
        {
            An   = An   * RECUR_BIG;
            Bn   = Bn   * RECUR_BIG;
            Anm1 = Anm1 * RECUR_BIG;
            Bnm1 = Bnm1 * RECUR_BIG;
            Anm2 = Anm2 * RECUR_BIG;
            Bnm2 = Bnm2 * RECUR_BIG;
            Anm3 = Anm3 * RECUR_BIG;
            Bnm3 = Bnm3 * RECUR_BIG;
        }
        
        n = n + 1;
        Bnm3 = Bnm2;
        Bnm2 = Bnm1;
        Bnm1 = Bn;
        Anm3 = Anm2;
        Anm2 = Anm1;
        Anm1 = An;
    }
  
    r.val = F;
    r.err = 2.0 * Math.abs( F * prec );
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * (n - 1) * Math.abs( F );
  
    return r;

} // hyperg_1F1_luke

// ----------------------------------------------------------------------------

// Series for 1F1(1,b,x)
// b > 0
//
function hyperg_1F1_1_series( b, x )
{
    var sum_val = 0.0;
    var sum_err = 0.0;
    var term    = 0.0;
    var n       = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    sum_val = 1.0;
    sum_err = 0.0;
    term = 1.0;
    n    = 1.0;
    while ( Math.abs( term / sum_val ) > 0.25 * GSL_DBL_EPSILON )
    {
        term = term * (x / (b + n - 1.0));
        sum_val = sum_val + term;
        sum_err = sum_err + 8.0 * GSL_DBL_EPSILON * Math.abs( term ) + GSL_DBL_EPSILON * Math.abs( sum_val );
        n = n + 1.0;
    }
    r.val = sum_val;
    r.err = sum_err;
    r.err = r.err + 2.0 *  Math.abs( term );
    return r;

} // hyperg_1F1_1_series

// ----------------------------------------------------------------------------

// 1F1(1,b,x)
// b >= 1, b integer
//
function hyperg_1F1_1_int( b, x )
{

    if ( b < 1 )
    {
        throw "SF.DomainException";
    }
    else if ( b == 1 )
    {
        return gsl_sf_exp_e( x );
    }
    else if ( b == 2 )
    {
        return gsl_sf_exprel_e( x );
    }
    else if ( b == 3 )
    {
        return gsl_sf_exprel_2_e( x );
    }
    else
    {
        return gsl_sf_exprel_n_e( b - 1, x );
    }

} // hyperg_1F1_1_int

// ----------------------------------------------------------------------------

// 1F1(1,b,x)
// b >=1, b real
//
// checked OK: [GJ] Thu Oct  1 16:46:35 MDT 1998
//
function hyperg_1F1_1( b, x )
{
    var ax = 0.0;
    var ib = 0.0;
    const INTEGER_LAST = +2147483647;

    var r  = { val: 0.0, err: 0.0 }; // Result;

    ax = Math.abs( x );
    ib = Math.floor( b + 0.1 );
  
    if ( b < 1.0 )
    {
        throw "SF.DomainException";
    }
    else if ( b == 1.0 )
    {
        return gsl_sf_exp_e( x );
    }
    else if ( b >= 1.4 * ax )
    {
        return hyperg_1F1_1_series( b, x );
    }
    else if ( Math.abs( b - ib ) < H1F1_INT_THRESHOLD && ib < (INTEGER_LAST) )
    {
        return hyperg_1F1_1_int( Math.trunc( ib ), x );
    }
    else if ( x > 0.0 )
    {
        if ( x > 100.0 && b < 0.75 * x )
        {
            return hyperg_1F1_asymp_posx( 1.0, b, x );
        }
        else if ( b < 1.0e+05 )
        {
            // Recurse backward on b, from a
            // chosen offset point. For x > 0,
            // which holds here, this should
            // be a stable direction.
            //
            var off     = 0.0;
            var bp      = 0.0;
            var err_rat = 0.0;
            var M       = { val: 0.0, err: 0.0 }; // Result;

            off = Math.ceil( 1.4 * x - b ) + 1.0;
            bp  = b + off;
            M = hyperg_1F1_1_series( bp, x );
            err_rat = M.err / Math.abs( M.val );
            while ( bp > b + 0.1 )
            {
                // M(1,b-1) = x/(b-1) M(1,b) + 1
                bp = bp - 1.0;
                M.val = 1.0 + x / bp * M.val;
            }
            r.val = M.val;
            r.err = err_rat * Math.abs( M.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( off ) + 1.0) * Math.abs( M.val );
            return r;
        }
        else if ( Math.abs( x ) < Math.abs( b ) && Math.abs( x ) < Math.sqrt( Math.abs( b ) ) * Math.abs( b - x ) )
        {
            return hyperg_1F1_largebx( 1.0, b, x );
        }
        else if ( Math.abs( x ) > Math.abs( b ) )
        {
            return hyperg_1F1_1_series( b, x );
        }
        else
        {
            return hyperg_1F1_large2bm4a( 1.0, b, x );
        }
    }
    else
    {
        // x <= 0 and b not large compared to |x|
        //
        if ( ax < 10.0 && b < 10.0 )
        {
            return hyperg_1F1_1_series( b, x );
        }
        else if ( ax >= 100.0 && Math.max( Math.abs( 2.0 - b ), 1.0 ) < 0.99 * ax )
        {
            return hyperg_1F1_asymp_negx( 1.0, b, x );
        }
        else
        {
            return hyperg_1F1_luke( 1.0, b, x );
        }
    }

} // hyperg_1F1_1

// ----------------------------------------------------------------------------

// 1F1(a,b,x)/Gamma(b) for b.0
// [limit of Abramowitz+Stegun 13.3.7]
//
function hyperg_1F1_renorm_b0( a, x )
{
    var eta = a * x;
    var r   = { val: 0.0, err: 0.0 }; // Result;

    if ( eta > 0.0 )
    {
        var root_eta = Math.sqrt( eta );
        var I1_scaled = { val: 0.0, err: 0.0 }; // Result;

        I1_scaled = gsl_sf_bessel_I1_scaled_e( 2.0 * root_eta );
        if ( I1_scaled.val <= 0.0 )
        {
            r.val = 0.0;
            r.err = 0.0;
            //code = GSL_ERROR_SELECT_2(stat_I, GSL_EDOM);
            throw "SF.DomainException";
        }
        else
        {
            // Note that 13.3.7 contains higher terms which are zeroth order
            // in b.  These make a non-negligible contribution to the sum.
            // With the first correction term, the I1 above is replaced by
            // I1 + (2/3)*a*(x/(4a))**(3/2)*I2(2*root_eta).  We will add
            // this as part of the result and error estimate.
            var corr1   = (2.0 / 3.0) * a * (x / Math.pow( 4.0 * a, 1.5) ) * gsl_sf_bessel_In_scaled( 2, 2.0 * root_eta );
            var lnr_val =  0.5 * x + 0.5 * Math.log( eta ) + Math.abs( 2.0 * root_eta ) + Math.log( I1_scaled.val + corr1 );
            var lnr_err =  GSL_DBL_EPSILON * (1.5 * Math.abs( x ) + 1.0) + Math.abs( (I1_scaled.err + corr1) / I1_scaled.val );

            return gsl_sf_exp_err_e( lnr_val, lnr_err );
        }
    }
    else if ( eta == 0.0 )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        // eta < 0
        let root_eta = Math.sqrt( -eta );
        let J1 = { val: 0.0, err: 0.0 }; // Result;

        J1 = gsl_sf_bessel_J1_e( 2.0 * root_eta );
        if ( J1.val <= 0.0 )
        {
            r.val = 0.0;
            r.err = 0.0;
            //code = GSL_ERROR_SELECT_2(stat_J, GSL_EDOM);
            throw "SF.DomainException";
        }
        else
        {
            let t1 = 0.5 * x;
            let t2 = 0.5 * Math.log( -eta );
            let t3 = Math.abs( x );
            let t4 = Math.log( J1.val );
            let lnr_val = t1 + t2 + t3 + t4;
            let lnr_err = GSL_DBL_EPSILON * (1.5 * Math.abs( x ) + 1.0) + Math.abs( J1.err / J1.val );
            let ex = { val: 0.0, err: 0.0 }; // Result;

            ex = gsl_sf_exp_err_e( lnr_val, lnr_err );
            r.val = -ex.val;
            r.err =  ex.err;
            return r;
        }
    }

} // hyperg_1F1_renorm_b0

// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

// 1F1'(a,b,x)/1F1(a,b,x)
// Uses Gautschi's series transformation of the
// continued fraction. This is apparently the best
// method for getting this ratio in the stable region.
// The convergence is monotone and supergeometric
// when b > x.
// Assumes a >= -1.
//
function hyperg_1F1_CF1_p_ser( a, b, x )
{
    var r = { val: 0.0, err: 0.0 }; // LONG_FLOAT = 0.0;

    if ( a == 0.0 )
    {
        r = 0.0;
        return r;
    }
    else
    {
        const maxiter = 5000;
        let sum  = 1.0;
        let pk   = 1.0;
        let rhok = 0.0;
        let ak   = 0.0;
        let k    = 0;

        k = 1;
        while ( k <= maxiter - 1 )
        {
            ak = (a + (k)) * x / ((b - x + (k - 1)) * (b - x + (k)));
            rhok = -ak * (1.0 + rhok) / (1.0 + ak * (1.0 + rhok));
            pk  = pk * rhok;
            sum = sum + pk;
            if ( Math.abs( pk / sum ) < 2.0 * GSL_DBL_EPSILON ) break;
            k = k + 1;
        }
        r = a / (b - x) * sum;
        if ( k >= maxiter )
        {
            throw "SF.MaxIterationsException";
        }
        else
        {
            return r;
        }
    }

} // hyperg_1F1_CF1_p_ser

// ----------------------------------------------------------------------------

// 1F1(a,b,x)
// |a| <= 1, b > 0
//
function hyperg_1F1_small_a_bgt0( a, b, x )
{
    var bma       = b - a;
    var oma       = 1.0 - a;
    var ap1mb     = 1.0 + a - b;
    var abs_bma   = Math.abs( bma );
    var abs_oma   = Math.abs( oma );
    var abs_ap1mb = Math.abs( ap1mb );
    var ax        = Math.abs( x );

    var r = { val: 0.0, err: 0.0 }; // Result;
  
    if ( a == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if ( a == 1.0 && b >= 1.0 )
    {
        return hyperg_1F1_1( b, x );
    }
    else if ( a == -1.0 )
    {
        r.val = 1.0 + a / b * x;
        r.err = GSL_DBL_EPSILON * (1.0 + Math.abs( a / b * x ));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( b >= 1.4 * ax )
    {
        return gsl_sf_hyperg_1F1_series_e( a, b, x );
    }
    else if ( x > 0.0 )
    {
        if ( x > 100.0 && abs_bma * abs_oma < 0.5 * x )
        {
            return hyperg_1F1_asymp_posx( a, b, x );
        }
        else if ( b < 5.0e+06 )
        {
            // Recurse backward on b from
            // a suitably high point.
            //
            let b_del   = 0.0;
            let Mbp1    = 0.0;
            let Mb      = 0.0;
            let Mbm1    = 0.0;
            let bp      = 0.0;
            let err_rat = 0.0;
            let r_Mbp1  = { val: 0.0, err: 0.0 }; // Result;
            let r_Mb    = { val: 0.0, err: 0.0 }; // Result;

            b_del = Math.ceil( 1.4 * x - b ) + 1.0;
            bp = b + b_del;
            r_Mbp1 = gsl_sf_hyperg_1F1_series_e( a, bp + 1.0, x );
            r_Mb   = gsl_sf_hyperg_1F1_series_e( a, bp,       x );
            err_rat = Math.abs( r_Mbp1.err / r_Mbp1.val ) + Math.abs( r_Mb.err / r_Mb.val );
            Mbp1 = r_Mbp1.val;
            Mb   = r_Mb.val;
            while ( bp > b + 0.1 )
            {
                // Do backward recursion.
                Mbm1 = ((x + bp - 1.0) * Mb - x * (bp - a) / bp * Mbp1) / (bp - 1.0);
                bp   = bp - 1.0;
                Mbp1 = Mb;
                Mb   = Mbm1;
            }
            r.val = Mb;
            r.err = err_rat * (Math.abs( b_del ) + 1.0) * Math.abs( Mb );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( Mb );
            return r;
        }
        else if ( Math.abs( x ) < Math.abs( b ) && Math.abs( a * x ) < Math.sqrt( Math.abs( b ) ) * Math.abs( b - x ))
        {
            return hyperg_1F1_largebx( a, b, x );
        }
        else
        {
            return hyperg_1F1_large2bm4a( a, b, x );
        }
    }
    else
    {
        // x < 0 and b not large compared to |x|
        //
        if ( ax < 10.0 && b < 10.0 )
        {
            return gsl_sf_hyperg_1F1_series_e( a, b, x );
        }
        else if ( ax >= 100.0 && Math.max( abs_ap1mb, 1.0 ) < 0.99 * ax )
        {
            return hyperg_1F1_asymp_negx( a, b, x );
        }
        else
        {
            return hyperg_1F1_luke( a, b, x );
        }
    }

} // hyperg_1F1_small_a_bgt0

// ----------------------------------------------------------------------------

// 1F1(b+eps,b,x)
// |eps|<=1, b > 0
//
function hyperg_1F1_beps_bgt0( eps, b, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( b > Math.abs( x ) && Math.abs( eps ) < GSL_SQRT_DBL_EPSILON )
    {
        // If b-a is very small and x/b is not too large we can
        // use this explicit approximation.
        //
        // 1F1(b+eps,b,x) = exp(ax/b) (1 - eps x^2 (v2 + v3 x + ...) + ...)
        //
        //   v2 = a/(2b^2(b+1))
        //   v3 = a(b-2a)/(3b^3(b+1)(b+2))
        //   ...
        //
        // See [Luke, Mathematical Functions and Their Approximations, p.292]
        //
        // This cannot be used for b near a negative integer or zero.
        // Also, if x/b is large the deviation from exp(x) behaviour grows.
        //
        let a  = 0.0;
        let v2 = 0.0;
        let v3 = 0.0;
        let v  = 0.0;
        let f  = 0.0;
        let exab = { val: 0.0, err: 0.0 }; // Result;

        a = b + eps;
        exab = gsl_sf_exp_e( a * x / b );
        v2 = a / (2.0 * b * b * (b + 1.0));
        v3 = a * (b - 2.0 * a) / (3.0 * b * b * b * (b + 1.0) * (b + 2.0));
        v  = v2 + v3 * x;
        f  = (1.0 - eps * x * x * v);
        r.val = exab.val * f;
        r.err = exab.err * Math.abs( f );
        r.err = r.err + Math.abs( exab.val ) * GSL_DBL_EPSILON * (1.0 + Math.abs( eps * x * x * v ));
        r.err = r.err + 4.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else
    {
        // Otherwise use a Kummer transformation to reduce
        // it to the small a case.
        //
        let Kummer_1F1 = { val: 0.0, err: 0.0 }; // Result;

        Kummer_1F1 = hyperg_1F1_small_a_bgt0( -eps, b, -x );
        if ( Kummer_1F1.val != 0.0 )
        {
            r = gsl_sf_exp_mult_err_e( x, 2.0 * GSL_DBL_EPSILON * Math.abs( x ), Kummer_1F1.val, Kummer_1F1.err );
        }
        else
        {
            r.val = 0.0;
            r.err = 0.0;
        }
        return r;
    }

} // hyperg_1F1_beps_bgt0

// ----------------------------------------------------------------------------

// 1F1(a,2a,x) = Gamma(a + 1/2) E(x) (|x|/4)^(-a+1/2) scaled_I(a-1/2,|x|/2)
//
// E(x) = exp(x) x > 0
//      = 1      x < 0
//
// a >= 1/2
//
function hyperg_1F1_beq2a_pos( a, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        let ln_term   = 0.0;
        let lnpre_val = 0.0;
        let lnpre_err = 0.0;
        let I  = { val: 0.0, err: 0.0 }; // Result;
        let lg = { val: 0.0, err: 0.0 }; // Result;

        I = gsl_sf_bessel_Inu_scaled_e( a - 0.5, 0.5 * Math.abs( x ) );
        lg = gsl_sf_lngamma_e( a + 0.5 );
        ln_term   = (0.5 - a) * Math.log( 0.25 * Math.abs( x ) );
        lnpre_val = lg.val + Math.max( x, 0.0 ) + ln_term;
        lnpre_err = lg.err + GSL_DBL_EPSILON * (Math.abs( ln_term ) + Math.abs( x ));
        r = gsl_sf_exp_mult_err_e( lnpre_val, lnpre_err, I.val, I.err );
        return r;
    }

} // hyperg_1F1_beq2a_pos

// ----------------------------------------------------------------------------

// Handle the case of a and b both positive integers.
// Assumes a > 0 and b > 0.
//
function hyperg_1F1_ab_posint( a, b, x )
{
    var ax = Math.abs( x );
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if ( a == b )
    {
        return gsl_sf_exp_e( x );             // 1F1(a,a,x)
    }
    else if ( a == 1 )
    {
        return gsl_sf_exprel_n_e( b - 1, x ); // 1F1(1,b,x)
    }
    else if ( b == a + 1 )
    {
        let K = { val: 0.0, err: 0.0 }; // Result;

        K = gsl_sf_exprel_n_e( a, -x );  // 1F1(1,1+a,-x)
        return gsl_sf_exp_mult_err_e( x, 2.0 * GSL_DBL_EPSILON * Math.abs( x ), K.val, K.err );
    }
    else if ( a == b + 1 )
    {
        let ex = { val: 0.0, err: 0.0 }; // Result;

        ex = gsl_sf_exp_e( x );
        r.val = ex.val * (1.0 + x / (b));
        r.err = ex.err * (1.0 + x / (b));
        r.err = r.err + ex.val * GSL_DBL_EPSILON * (1.0 + Math.abs( x / (b) ));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( a == b + 2 )
    {
        let ex   = { val: 0.0, err: 0.0 }; // Result;
        let poly = 0.0;

        ex = gsl_sf_exp_e( x );
        poly = (1.0 + x / (b) * (2.0 + x / (b + 1)));
        r.val = ex.val * poly;
        r.err = ex.err * Math.abs( poly );
        r.err = r.err + ex.val * GSL_DBL_EPSILON * (1.0 + Math.abs( x / (b) ) * (2.0 + Math.abs( x / (b + 1) )));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( b == 2 * a )
    {
        return hyperg_1F1_beq2a_pos( (a), x );  // 1F1(a,2a,x)
    }
    else if ( (b < 10 && a < 10 && ax < 5.0) || ((b) > (a) * ax) || (b > a && ax < 5.0) )
    {
        return gsl_sf_hyperg_1F1_series_e( (a), (b), x );
    }
    else if ( b > a && (b) >= (2 * a) + x )
    {
        // Use the Gautschi CF series, then
        // recurse backward to a=0 for normalization.
        // This will work for either sign of x.
        //
        let rap  = 0.0;
        let ra   = 0.0;
        let Ma   = 0.0;
        let Map1 = 0.0;
        let Mnp1 = 0.0;
        let Mn   = 0.0;
        let Mnm1 = 0.0;

        rap = hyperg_1F1_CF1_p_ser( (a), (b), x );
        ra = 1.0 + x / (a) * rap;
        Ma   = GSL_SQRT_DBL_MIN;
        Map1 = ra * Ma;
        Mnp1 = Map1;
        Mn   = Ma;
        for ( let n = a; n >= 1; n-- )
        {
            Mnm1 = ((n) * Mnp1 - ((2 * n - b) + x) * Mn) / (b - n);
            Mnp1 = Mn;
            Mn   = Mnm1;
        }
        r.val = Ma / Mn;
        r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( a ) + 1) * Math.abs( Ma / Mn );
        return r;
    }
    else if ( b > a && (b) < (2 * a) + x && (b) > x )
    {
        // Use the Gautschi series representation of
        // the continued fraction. Then recurse forward
        // to the a=b line for normalization. This will
        // work for either sign of x, although we do need
        // to check for b > x, for when x is positive.
        //
        let rap  = 0.0;
        let ra   = 0.0;
        let Ma   = 0.0;
        let Map1 = 0.0;
        let Mnm1 = 0.0;
        let Mn   = 0.0;
        let Mnp1 = 0.0;
        let ex   = { val: 0.0, err: 0.0 }; // Result;

        rap = hyperg_1F1_CF1_p_ser( (a), (b), x) ;
        ra = 1.0 + x / (a) * rap;
        
        Ma   = GSL_SQRT_DBL_MIN;
        Map1 = ra * Ma;
        Mnm1 = Ma;
        Mn   = Map1;
        for( let n = a + 1; n <= b - 1; n++ )
        {
            Mnp1 = ((b - n) * Mnm1 + ((2 * n - b) + x) * Mn) / (n);
            Mnm1 = Mn;
            Mn   = Mnp1;
        }
        
        ex = gsl_sf_exp_e( x );  // 1F1(b,b,x)
        r.val = ex.val * Ma / Mn;
        r.err = ex.err * Math.abs( Ma / Mn );
        r.err = r.err + 4.0 * GSL_DBL_EPSILON * (Math.abs( b - a ) + 1) * Math.abs( r.val );
        return r;
    }
    else if ( x >= 0.0 )
    {
        if ( b < a )
        {
            // The point b,b is below the b=2a+x line.
            // Forward recursion on a from b,b+1 is possible.
            // Note that a > b + 1 as well, since we already tried a = b + 1.
            //
            if ( x + Math.log( Math.abs( x / (b) ) ) < GSL_LOG_DBL_MAX - 2.0 )
            {
                let ex   = Math.exp( x );
                let Mnm1 = ex;                 // 1F1(b,b,x)
                let Mn   = ex * (1.0 + x / (b)); // 1F1(b+1,b,x)
                let Mnp1 = 0.0;

                for ( let n = b + 1; n <= a - 1; n++ )
                {
                    Mnp1 = ((b - n) * Mnm1 + ((2 * n - b) + x) * Mn) / (n);
                    Mnm1 = Mn;
                    Mn   = Mnp1;
                }
                r.val = Mn;
                r.err = (x + 1.0) * GSL_DBL_EPSILON * Math.abs( Mn );
                r.err = r.err * (Math.abs( a - b ) + 1);
                return r;
            }
            else
            {
                throw "SF.OverflowException";
            }
        }
        else
        {
            // b > a
            // b < 2a + x 
            // b <= x (otherwise we would have finished above)
            //
            // Gautschi anomalous convergence region. However, we can
            // recurse forward all the way from a=0,1 because we are
            // always underneath the b=2a+x line.
            //
            let r_Mn = { val: 0.0, err: 0.0 }; // Result;
            let Mnm1 = 1.0; // 1F1(0,b,x)
            let Mn   = 0.0; // 1F1(1,b,x)
            let Mnp1 = 0.0;

            r_Mn = gsl_sf_exprel_n_e( b - 1, x );
            Mn = r_Mn.val;
            for ( let n = 1; n <= a - 1; n++ )
            {
                Mnp1 = ((b - n) * Mnm1 + ((2 * n - b) + x) * Mn) / (n);
                Mnm1 = Mn;
                Mn   = Mnp1;
            }
            r.val = Mn;
            r.err = Math.abs( Mn ) * (1 + Math.abs( a )) * Math.abs( r_Mn.err / r_Mn.val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( Mn );
            return r;
        }
    }
    else
    {
        // x < 0
        // b < a (otherwise we would have tripped one of the above)
        //
       
        if ( (a) <= 0.5 * ((b) - x) || (a) >= -x )
        {
            // Gautschi continued fraction is in the anomalous region,
            // so we must find another way. We recurse down in b,
            // from the a=b line.
            //
            let ex    = Math.exp(x);
            let Manp1 = ex;
            let Man   = ex * (1.0 + x / (a - 1));
            let Manm1 = 0.0;

            for ( let n = a - 1; n >= b + 1; n-- )
            {
                Manm1 = (-(n) * ((1 - n) - x) * Man - x * (n - a) * Manp1) / (n * (n - 1));
                Manp1 = Man;
                Man   = Manm1;
            }
            r.val = Man;
            r.err = (Math.abs( x ) + 1.0) * GSL_DBL_EPSILON * Math.abs( Man );
            r.err = r.err * (Math.abs( b - a ) + 1);
            return r;
        }
        else
        {
            // Pick a0 such that b ~= 2a0 + x, then
            // recurse down in b from a0,a0 to determine
            // the values near the line b=2a+x. Then recurse
            // forward on a from a0.
            //
            let a0 = Math.trunc( Math.ceil( 0.5 * ((b) - x) ) );
            let Ma0b   = 0.0;  // M(a0,b)
            let Ma0bp1 = 0.0;  // M(a0,b+1)
            let Ma0p1b = 0.0;  // M(a0+1,b)
            let Mnm1   = 0.0;
            let Mn     = 0.0;
            let Mnp1   = 0.0;
            let ex     = 0.0;
            let Ma0np1 = 0.0;
            let Ma0n   = 0.0;
            let Ma0nm1 = 0.0;

            ex = Math.exp( x );
            Ma0np1 = ex;
            Ma0n   = ex * (1.0 + x / (a0 - 1));
            for ( let n = a0 - 1; n >= b + 1; n-- )
            {
                Ma0nm1 = (-(n) * ((1 - n) - x) * Ma0n - x * (n - a0) * Ma0np1) / (n * (n - 1));
                Ma0np1 = Ma0n;
                Ma0n = Ma0nm1;
            }
            Ma0bp1 = Ma0np1;
            Ma0b   = Ma0n;
            Ma0p1b = ((b) * ((a0) + x) * Ma0b + x * (a0 - b) * Ma0bp1) / (a0 * b);
            
            // Initialise the recurrence correctly BJG
            
            if ( a0 >= a )
            {
                Mn = Ma0b;
            }
            else if ( a0 + 1 >= a )
            {
                Mn = Ma0p1b;
            }
            else
            {
                Mnm1 = Ma0b;
                Mn   = Ma0p1b;
            
                for ( let n = a0 + 1; n <= a - 1; n++ )
                {
                    Mnp1 = ((b - n) * Mnm1 + ((2 * n - b) + x) * Mn) / (n);
                    Mnm1 = Mn;
                    Mn   = Mnp1;
                }
            }
            
            r.val = Mn;
            r.err = (Math.abs( x ) + 1.0) * GSL_DBL_EPSILON *  Math.abs( Mn );
            r.err = r.err * (Math.abs( b - a ) + 1);
            return r;
        }
    }

} // hyperg_1F1_ab_posint

// ----------------------------------------------------------------------------

// Evaluate a <= 0, a integer, cases directly. (Polynomial; Horner)
// When the terms are all positive, this
// must work. We will assume this here.
//
function hyperg_1F1_a_negint_poly( a, b, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( a == 0 )
    {
        r.val = 1.0;
        r.err = 1.0;
        return r;
    }
    else
    {
        var N    = -a;
        var poly = 1.0;
        var t    = 0.0;
        var s    = 0.0;

        for ( let k = N - 1; k >= 0; k-- )
        {
            t = (a + k) / (b + k) * (x / (k + 1));
            s = t + 1.0 / poly;
            if ( s > 0.9 * GSL_DBL_MAX / poly )
            {
                throw "SF.OverflowException";
            }
            else
            {
                poly = poly * s;  // P_n = 1 + t_n P_{n-1}
            }
        }
        r.val = poly;
        r.err = 2.0 * (Math.sqrt( (N) ) + 1.0) * GSL_DBL_EPSILON * Math.abs( poly );
        return r;
    }

} // hyperg_1F1_a_negint_poly

// ----------------------------------------------------------------------------

// Evaluate negative integer a case by relation
// to Laguerre polynomials. This is more general than
// the direct polynomial evaluation, but is safe
// for all values of x.
//
// 1F1(-n,b,x) = n!/(b)_n Laguerre[n,b-1,x]
//             = n B(b,n) Laguerre[n,b-1,x]
//
// assumes b is not a negative integer
//
function hyperg_1F1_a_negint_lag( a, b, x )
{
    var n   = -a;
    var lag = { val: 0.0, err: 0.0 }; // Result;
    var r   = { val: 0.0, err: 0.0 }; // Result;

    lag = gsl_sf_laguerre_n_e( n, b - 1.0, x );
    if ( b < 0.0 )
    {
        let lnfact    = { val: 0.0, err: 0.0 }; // Result;
        let lng1      = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let lng2      = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let s1        = 0.0;
        let s2        = 0.0;
        let lnpre_val = 0.0;
        let lnpre_err = 0.0;

        lnfact = gsl_sf_lnfact_e( n );
        lng1 = gsl_sf_lngamma_sgn_e( b + (n) ); //, lng1, s1 );
        lng2 = gsl_sf_lngamma_sgn_e( b ); //, lng2, s2 );
        s1 = lng1.sign;
        s2 = lng2.sign;
        lnpre_val = lnfact.val - (lng1.val - lng2.val);
        lnpre_err = lnfact.err + lng1.err + lng2.err + 2.0 * GSL_DBL_EPSILON * Math.abs( lnpre_val );
        r = gsl_sf_exp_mult_err_e( lnpre_val, lnpre_err, s1 * s2 * lag.val, lag.err );
        return r;
    }
    else
    {
        let lnbeta = { val: 0.0, err: 0.0 }; // Result;

        lnbeta = gsl_sf_lnbeta_e( b, (n) );
        if ( Math.abs(lnbeta.val) < 0.1 )
        {
            // As we have noted, when B(x,y) is near 1,
            // evaluating log(B(x,y)) is not accurate.
            // Instead we evaluate B(x,y) directly.
            //
            let ln_term_val = Math.log( 1.25 * (n) );
            let ln_term_err = 2.0 * GSL_DBL_EPSILON * ln_term_val;
            let beta = { val: 0.0, err: 0.0 }; // Result;

            beta = gsl_sf_beta_e( b, (n) );
            r = gsl_sf_exp_mult_err_e( ln_term_val, ln_term_err, lag.val, lag.err );
            r.val = r.val * beta.val / 1.25;
            r.err = r.err * beta.val / 1.25;
            return r;
        }
        else
        {
            // B(x,y) was not near 1, so it is safe to use
            // the logarithmic values.
            //
            let ln_n = Math.log( (n) );
            let ln_term_val = lnbeta.val + ln_n;
            let ln_term_err = lnbeta.err + 2.0 * GSL_DBL_EPSILON * Math.abs( ln_n );

            r = gsl_sf_exp_mult_err_e( ln_term_val, ln_term_err, lag.val, lag.err );
            return r;
        }
    }

} // hyperg_1F1_a_negint_lag

// ----------------------------------------------------------------------------


// Assumes a <= -1,  b <= -1, and b <= a.
//
function hyperg_1F1_ab_negint( a, b, x )
{
    var K = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if ( x > 0.0 )
    {
        return hyperg_1F1_a_negint_poly( a, b, x );
    }
    else
    {
        // Apply a Kummer transformation to make x > 0 so
        // we can evaluate the polynomial safely. Of course,
        // this assumes b <= a, which must be true for
        // a<0 and b<0, since otherwise the thing is undefined.
        //
        K = hyperg_1F1_a_negint_poly( b - a, b, -x );
        return gsl_sf_exp_mult_err_e( x, 2.0 * GSL_DBL_EPSILON * Math.abs( x ), K.val, K.err );
    }

} // hyperg_1F1_ab_negint

// ----------------------------------------------------------------------------

// [Abramowitz+Stegun, 13.1.3]
//
// M(a,b,x) = Gamma(1+a-b)/Gamma(2-b) x^(1-b) *
//            { Gamma(b)/Gamma(a) M(1+a-b,2-b,x) - (b-1) U(1+a-b,2-b,x) }
//
// b not an integer >= 2
// a-b not a negative integer
//
function hyperg_1F1_U( a, b, x )
{
    var bp         = 0.0;
    var ap         = 0.0;
    var sg_ap      = 0.0;
    var t1         = 0.0;
    var lnpre_val  = 0.0;
    var lnpre_err  = 0.0;
    var lnc1_val   = 0.0;
    var lnc1_err   = 0.0;
    var sg_2mbp    = 0.0;
    var sg_1papmbp = 0.0;
    var ombp       = 0.0;
    var Uee_val    = 0.0;
    var Uee_err    = 0.0;
    var Mee_val    = 0.0;
    var Mee_err    = 0.0;
    var lg_ap      = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var lg_bp      = { val: 0.0, err: 0.0 }; // Result;
    var lg_2mbp    = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var lg_1papmbp = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var M          = { val: 0.0, err: 0.0 }; // Result;
    var U          = { val: 0.0, err: 0.0 }; // Result;
    var term_M     = { val: 0.0, err: 0.0 }; // Result;
    var r          = { val: 0.0, err: 0.0 }; // Result;

    bp = 2.0 - b;
    ap = a - b + 1.0;
  
    
    lg_ap = gsl_sf_lngamma_sgn_e( ap ); //, lg_ap, sg_ap );
    sg_ap = lg_ap.sign;
    lg_bp = gsl_sf_lngamma_e( bp );
    t1 = (bp - 1.0) * Math.log( x );
    lnpre_val = lg_ap.val - lg_bp.val + t1;
    lnpre_err = lg_ap.err + lg_bp.err + 2.0 * GSL_DBL_EPSILON * Math.abs( t1 );
  
    lg_2mbp = gsl_sf_lngamma_sgn_e( 2.0 - bp ); //,      lg_2mbp,    sg_2mbp );
    sg_2mbp = lg_2mbp.sign;
    lg_1papmbp = gsl_sf_lngamma_sgn_e( 1.0 + ap - bp ); //, lg_1papmbp, sg_1papmbp );
    sg_1papmbp = lg_1papmbp.sign;
    lnc1_val = lg_2mbp.val - lg_1papmbp.val;
    lnc1_err = lg_2mbp.err + lg_1papmbp.err + GSL_DBL_EPSILON * (Math.abs( lg_2mbp.val ) + Math.abs( lg_1papmbp.val ));
  
    M = gsl_sf_hyperg_1F1_e( ap, bp, x );
    U = gsl_sf_hyperg_U_e10_e( ap, bp, x );
  
    term_M = gsl_sf_exp_mult_err_e10_e( lnc1_val, lnc1_err, sg_2mbp * sg_1papmbp * M.val, M.err );
  
    ombp = 1.0 - bp;
    Uee_val = (U.e10) * M_LN10;
    Uee_err = 2.0 * GSL_DBL_EPSILON * Math.abs( Uee_val );
    Mee_val = (term_M.e10) * M_LN10;
    Mee_err = 2.0 * GSL_DBL_EPSILON * Math.abs( Mee_val );
  
    // Do a little dance with the exponential prefactors
    // to avoid overflows in intermediate results.
    //
    if ( Uee_val > Mee_val )
    {
        let factorM_val = Math.exp( Mee_val - Uee_val );
        let factorM_err = 2.0 * GSL_DBL_EPSILON * (Math.abs( Mee_val - Uee_val ) + 1.0) * factorM_val;
        let inner_val   = term_M.val * factorM_val - ombp * U.val;
        let inner_err   = term_M.err * factorM_val + Math.abs( ombp ) * U.err
                                + Math.abs( term_M.val ) * factorM_err
                                + GSL_DBL_EPSILON * (Math.abs( term_M.val * factorM_val ) + Math.abs( ombp * U.val ));

        r = gsl_sf_exp_mult_err_e( lnpre_val + Uee_val, lnpre_err + Uee_err, sg_ap * inner_val, inner_err );
    }
    else
    {
        let factorU_val = Math.exp( Uee_val - Mee_val );
        let factorU_err = 2.0 * GSL_DBL_EPSILON * (Math.abs( Mee_val - Uee_val ) + 1.0) * factorU_val;
        let inner_val   = term_M.val - ombp * factorU_val * U.val;
        let inner_err   = term_M.err + Math.abs( ombp * factorU_val * U.err )
                                + Math.abs( ombp * factorU_err * U.val )
                                + GSL_DBL_EPSILON * (Math.abs( term_M.val ) + Math.abs( ombp * factorU_val * U.val ));

        r = gsl_sf_exp_mult_err_e( lnpre_val + Mee_val, lnpre_err + Mee_err, sg_ap * inner_val, inner_err );
    }
  
    return r;

} // hyperg_1F1_U

// ----------------------------------------------------------------------------

// Handle case of generic positive a, b.
// Assumes b-a is not a negative integer.
//
function hyperg_1F1_ab_pos( a, b, x )
{
    var ax = Math.abs( x );
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if ( (b < 10.0 && a < 10.0 && ax < 5.0) || (b > a * ax) || (b > a && ax < 5.0) )
    {
        return gsl_sf_hyperg_1F1_series_e( a, b, x );
    }
    else if ( x < -100.0 && Math.max( Math.abs( a ), 1.0 ) * Math.max( Math.abs( 1.0 + a - b ), 1.0 ) < 0.7 * Math.abs( x ) )
    {
        // Large negative x asymptotic.
        //
        return hyperg_1F1_asymp_negx( a, b, x );
    }
    else if ( x > 100.0 && Math.max( Math.abs( b - a ), 1.0 ) * Math.max( Math.abs( 1.0 - a ), 1.0 ) < 0.7 * Math.abs( x ) )
    {
        // Large positive x asymptotic.
        //
        return hyperg_1F1_asymp_posx( a, b, x );
    }
    else if ( Math.abs( b - a ) <= 1.0 )
    {
        // Directly handle b near a.
        //
        return hyperg_1F1_beps_bgt0( a - b, b, x );  // a = b + eps
    }
    else if ( b > a && b >= 2.0 * a + x )
    {
        // Use the Gautschi CF series, then
        // recurse backward to a near 0 for normalization.
        // This will work for either sign of x.
        //
        let n    = 0.0;
        let rap  = 0.0;
        let ra   = 0.0;
        let Ma   = 0.0;
        let Mn   = 0.0;
        let Map1 = 0.0;
        let Mnp1 = 0.0;
        let Mnm1 = 0.0;
        let Mn_true = { val: 0.0, err: 0.0 }; // Result;

        rap = hyperg_1F1_CF1_p_ser( a, b, x );
        ra = 1.0 + x / a * rap;
    
        Ma   = GSL_SQRT_DBL_MIN;
        Map1 = ra * Ma;
        Mnp1 = Map1;
        Mn   = Ma;
        n = a;
        while ( n > 0.5 )
        {
            Mnm1 = (n * Mnp1 - (2.0 * n - b + x) * Mn) / (b - n);
            Mnp1 = Mn;
            Mn   = Mnm1;
            n    = n - 1.0;
        }
    
        Mn_true = hyperg_1F1_small_a_bgt0( n, b, x );
    
        r.val = (Ma / Mn) * Mn_true.val;
        r.err = Math.abs( Ma / Mn ) * Mn_true.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( a ) + 1.0) * Math.abs( r.val );
        return r;
    }
    else if ( b > a && b < 2.0 * a + x && b > x )
    {
        // Use the Gautschi series representation of
        // the continued fraction. Then recurse forward
        // to near the a=b line for normalization. This will
        // work for either sign of x, although we do need
        // to check for b > x, which is relevant when x is positive.
        //
        let Mn_true = { val: 0.0, err: 0.0 }; // Result;
        let rap  = 0.0;
        let ra   = 0.0;
        let Ma   = 0.0;
        let Mn   = 0.0;
        let Mnm1 = 0.0;
        let Mnp1 = 0.0;
        let n    = 0.0;

        rap = hyperg_1F1_CF1_p_ser( a, b, x );
        ra   = 1.0 + x / a * rap;
        Ma   = GSL_SQRT_DBL_MIN;
        Mnm1 = Ma;
        Mn   = ra * Mnm1;
        n = a + 1.0;
        while ( n < b - 0.5 )
        {
            Mnp1 = ((b - n) * Mnm1 + (2.0 * n - b + x) * Mn) / n;
            Mnm1 = Mn;
            Mn   = Mnp1;
            n    = n + 1.0;
        }
        Mn_true = hyperg_1F1_beps_bgt0( n - b, b, x );
        r.val = Ma / Mn * Mn_true.val;
        r.err = Math.abs( Ma / Mn ) * Mn_true.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( b - a ) + 1.0) * Math.abs( r.val );
        return r;
    }
    else if ( x >= 0.0 )
    {
        if ( b < a )
        {
            // Forward recursion on a from a=b+eps-1,b+eps.
            //
            let N    = Math.floor( a - b );
            let eps  = a - b - N;
            let M0   = 0.0;
            let M1   = 0.0;
            let Mam1 = 0.0;
            let Ma   = 0.0;
            let Map1 = 0.0;
            let ap   = 0.0;
            let start_pair = 0.0;
            let minim_pair = 0.0;
            let pair_ratio = 0.0;
            let rat_0 = 0.0;
            let rat_1 = 0.0;
            let r_M0 = { val: 0.0, err: 0.0 }; // Result;
            let r_M1 = { val: 0.0, err: 0.0 }; // Result;

            r_M0 = hyperg_1F1_beps_bgt0( eps - 1.0, b, x );
            r_M1 = hyperg_1F1_beps_bgt0( eps,       b, x );
            M0 = r_M0.val;
            M1 = r_M1.val;
            
            Mam1 = M0;
            Ma   = M1;
            start_pair = Math.abs( M0 ) + Math.abs( M1 );
            minim_pair = GSL_DBL_MAX;
            rat_0 = Math.abs( r_M0.err / r_M0.val );
            rat_1 = Math.abs( r_M1.err / r_M1.val );
            ap = b + eps;
            while ( ap < a - 0.1 )
            {
                Map1 = ((b - ap) * Mam1 + (2.0 * ap - b + x) * Ma) / ap;
                Mam1 = Ma;
                Ma   = Map1;
                minim_pair = Math.min( Math.abs( Mam1 ) + Math.abs( Ma ), minim_pair );
                ap = ap + 1.0;
            }
            pair_ratio = start_pair / minim_pair;
            r.val = Ma;
            r.err = 2.0 * (rat_0 + rat_1 + GSL_DBL_EPSILON) * (Math.abs( b - a ) + 1.0) * Math.abs( Ma );
            r.err = r.err + 2.0 * (rat_0 + rat_1) * pair_ratio * pair_ratio * Math.abs( Ma );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( Ma );
            return r;
        }
        else
        {
            // b > a
            // b < 2a + x 
            // b <= x
            //
            // Recurse forward on a from a=eps,eps+1.
            //
            let n    = 0.0;
            let eps  = 0.0;
            let M0   = 0.0;
            let M1   = 0.0;
            let Mnm1 = 0.0;
            let Mn   = 0.0;
            let Mnp1 = 0.0;
            let ap   = 0.0;
            let start_pair = 0.0;
            let minim_pair = 0.0;
            let pair_ratio = 0.0;
            let rat_0 = 0.0;
            let rat_1 = 0.0;
            let r_Mnm1 = { val: 0.0, err: 0.0 }; // Result;
            let r_Mn   = { val: 0.0, err: 0.0 }; // Result;

            eps    = a - Math.floor( a );
            r_Mnm1 = hyperg_1F1_small_a_bgt0( eps,       b, x );
            r_Mn   = hyperg_1F1_small_a_bgt0( eps + 1.0, b, x );
            Mnm1   = r_Mnm1.val;
            Mn     = r_Mn.val;
            
            start_pair = Math.abs( Mn ) + Math.abs( Mnm1 );
            minim_pair = GSL_DBL_MAX;
            rat_0 = Math.abs( r_Mnm1.err / r_Mnm1.val );
            rat_1 = Math.abs( r_Mn.err / r_Mn.val );
            n = eps + 1.0;
            while ( n < a - 0.1 )
            {
                Mnp1 = ((b - n) * Mnm1 + (2.0 * n - b + x) * Mn) / n;
                Mnm1 = Mn;
                Mn   = Mnp1;
                minim_pair = Math.min( Math.abs( Mn ) + Math.abs( Mnm1 ), minim_pair );
                n = n + 1.0;
            }
            pair_ratio = start_pair / minim_pair;
            r.val = Mn;
            r.err = 2.0 * (rat_0 + rat_1 + GSL_DBL_EPSILON) * (Math.abs( a ) + 1.0) * Math.abs( Mn );
            r.err = r.err + 2.0 * (rat_0 + rat_1) * pair_ratio * pair_ratio * Math.abs( Mn );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( Mn );
            return r;
        }
    }
    else
    {
        // x < 0
        // b < a
        //
    
        if ( a <= 0.5 * (b - x) || a >= -x )
        {
            // Recurse down in b, from near the a=b line, b=a+eps,a+eps-1.
            //
            let n       = 0.0;
            let eps     = 1.0 + Math.floor( a - b ) - a + b;
            let Manp1   = 0.0;
            let Man     = 0.0;
            let Manm1   = 0.0;
            let start_pair = 0.0;
            let minim_pair = 0.0;
            let pair_ratio = 0.0;
            let rat_0   = 0.0;
            let rat_1   = 0.0;
            let r_Manp1 = { val: 0.0, err: 0.0 }; // Result;
            let r_Man   = { val: 0.0, err: 0.0 }; // Result;

            r_Manp1 = hyperg_1F1_beps_bgt0( -eps,      a + eps,       x );
            r_Man   = hyperg_1F1_beps_bgt0( 1.0 - eps, a + eps - 1.0, x );
            Manp1 = r_Manp1.val;
            Man   = r_Man.val;
            start_pair = Math.abs( Manp1 ) + Math.abs( Man );
            minim_pair = GSL_DBL_MAX;
            rat_0 = Math.abs( r_Manp1.err / r_Manp1.val );
            rat_1 = Math.abs( r_Man.err / r_Man.val );
            n = a + eps - 1.0;
            while ( n > b + 0.1 )
            {
                Manm1 = (-n * (1.0 - n - x) * Man - x * (n - a) * Manp1) / (n * (n - 1.0));
                Manp1 = Man;
                Man = Manm1;
                minim_pair = Math.min( Math.abs( Manp1 ) + Math.abs( Man ), minim_pair );
                n = n - 1.0;
            }
            
            // FIXME: this is a nasty little hack; there is some
            // (transient?) instability in this recurrence for some
            // values. I can tell when it happens, which is when
            // this pair_ratio is large. But I do not know how to
            // measure the error in terms of it. I guessed quadratic
            // below, but it is probably worse than that.
            //
            pair_ratio = start_pair / minim_pair;
            r.val = Man;
            r.err = 2.0 * (rat_0 + rat_1 + GSL_DBL_EPSILON) * (Math.abs( b - a ) + 1.0) * Math.abs( Man );
            r.err = r.err * (pair_ratio * pair_ratio + 1.0);
            return r;
        }
        else
        {
            // Pick a0 such that b ~= 2a0 + x, then
            // recurse down in b from a0,a0 to determine
            // the values near the line b=2a+x. Then recurse
            // forward on a from a0.
            //
            let epsa = a - Math.floor( a );
            let a0   = Math.floor( 0.5 * (b - x) ) + epsa;
            let epsb = 1.0 + Math.floor( a0 - b ) - a0 + b;
            let Ma0b = 0.0;
            let Ma0bp1    = 0.0;
            let Ma0p1b    = 0.0;
            let Mnm1      = 0.0;
            let Mn        = 0.0;
            let Mnp1      = 0.0;
            let n         = 0.0;
            let err_rat   = 0.0;
            let Ma0np1    = 0.0;
            let Ma0n      = 0.0;
            let Ma0nm1    = 0.0;
            let r_Ma0np1  = { val: 0.0, err: 0.0 }; // Result;
            let r_Ma0n    = { val: 0.0, err: 0.0 }; // Result;
 
            r_Ma0np1 = hyperg_1F1_beps_bgt0( -epsb,      a0 + epsb,       x );
            r_Ma0n   = hyperg_1F1_beps_bgt0( 1.0 - epsb, a0 + epsb - 1.0, x );
            Ma0np1 = r_Ma0np1.val;
            Ma0n   = r_Ma0n.val;
            
            err_rat = Math.abs( r_Ma0np1.err / r_Ma0np1.val ) + Math.abs( r_Ma0n.err / r_Ma0n.val );
            n = a0 + epsb - 1.0;
            while ( n > b + 0.1 )
            {
                Ma0nm1 = (-n * (1.0 - n - x) * Ma0n - x * (n - a0) * Ma0np1) / (n * (n - 1.0));
                Ma0np1 = Ma0n;
                Ma0n   = Ma0nm1;
                n      = n - 1.0;
            }
            Ma0bp1 = Ma0np1;
            Ma0b   = Ma0n;
            Ma0p1b = (b * (a0 + x) * Ma0b + x * (a0 - b) * Ma0bp1) / (a0 * b); // right-down hook
                
            // Initialise the recurrence correctly BJG
            
            if ( a0 >= a - 0.1 )
            {
                Mn = Ma0b;
            }
            else if ( a0 + 1.0 >= a - 0.1 )
            {
                Mn = Ma0p1b;
            }
            else
            {
                Mnm1 = Ma0b;
                Mn   = Ma0p1b;
            
                n = a0 + 1.0;
                while ( n < a - 0.1 )
                {
                    Mnp1 = ((b - n) * Mnm1 + (2.0 * n - b + x) * Mn) / n;
                    Mnm1 = Mn;
                    Mn   = Mnp1;
                    n    = n + 1.0;
                }
            }
            
            r.val = Mn;
            r.err = (err_rat + GSL_DBL_EPSILON) * (Math.abs( b - a ) + 1.0) * Math.abs( Mn );
            return r;
        }
    }

} // hyperg_1F1_ab_pos

// ----------------------------------------------------------------------------

// Assumes b != integer
// Assumes a != integer when x > 0
// Assumes b-a != neg integer when x < 0
//
function hyperg_1F1_ab_neg( a, b, x )
{
    var bma    = b - a;
    var abs_x  = Math.abs( x );
    var abs_a  = Math.abs( a );
    var abs_b  = Math.abs( b );
    var size_a = Math.max( abs_a, 1.0 );
    var size_b = Math.max( abs_b, 1.0 );
    var bma_integer = (bma - Math.floor( bma + 0.5 ) < H1F1_INT_THRESHOLD);

    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( (abs_a < 10.0 && abs_b < 10.0 && abs_x < 5.0) || (b > 0.8 * Math.max( Math.abs( a ), 1.0 ) * Math.abs( x )) )
    {
        return gsl_sf_hyperg_1F1_series_e( a, b, x );
    }
    else if ( x > 0.0 && size_b > size_a && size_a * Math.log( M_E * x / size_b ) < GSL_LOG_DBL_EPSILON + 7.0 )
    {
        // Series terms are positive definite up until
        // there is a sign change. But by then the
        // terms are small due to the last condition.
        //
        return gsl_sf_hyperg_1F1_series_e( a, b, x );
    }
    else if ( (abs_x < 5.0 && Math.abs( bma ) < 10.0 && abs_b < 10.0) || (b > 0.8 * Math.max( Math.abs( bma ), 1.0 ) * abs_x) )
    {
        // Use Kummer transformation to render series safe.
        //
        var Kummer_1F1 = { val: 0.0, err: 0.0 }; // Result;

        Kummer_1F1 = gsl_sf_hyperg_1F1_series_e( bma, b, -x );
        r = gsl_sf_exp_mult_err_e( x, GSL_DBL_EPSILON * Math.abs( x ), Kummer_1F1.val, Kummer_1F1.err );
        return r;
    }
    else if ( x < -30.0 && Math.max( Math.abs( a ), 1.0 ) * Math.max( Math.abs( 1.0 + a - b ), 1.0 ) < 0.99 * Math.abs( x ) )
    {
        // Large negative x asymptotic.
        // Note that we do not check if b-a is a negative integer.
        //
        return hyperg_1F1_asymp_negx( a, b, x );
    }
    else if ( x > 100.0 && Math.max( Math.abs( bma ), 1.0 ) * Math.max( Math.abs( 1.0 - a ), 1.0 ) < 0.99 * Math.abs( x ) )
    {
        // Large positive x asymptotic.
        // Note that we do not check if a is a negative integer.
        //
        return hyperg_1F1_asymp_posx( a, b, x );
    }
    else if ( x > 0.0 && ! (bma_integer && bma > 0.0) )
    {
        return hyperg_1F1_U( a, b, x );
    }
    else
    {
        // FIXME:  if all else fails, try the series... BJG
        if ( x < 0.0 )
        {
            // Apply Kummer Transformation
            var K_factor = 0.0;

            r = gsl_sf_hyperg_1F1_series_e( b - a, b, -x );
            K_factor = Math.exp( x );
            r.val = r.val * K_factor;
            r.err = r.err * K_factor;
            return r;
        }
        else
        {
            return gsl_sf_hyperg_1F1_series_e( a, b, x );
        }
       
        // Sadness...
        // result.val = 0.0;
        // result.err = 0.0;
        // GSL_ERROR ("error", GSL_EUNIMPL);
    }

} // hyperg_1F1_ab_neg

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hyperg_1F1_int_e( a, b, x )
{
    var Kummer_1F1 = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if ( a == b )
    {
        return gsl_sf_exp_e( x );
    }
    else if ( b == 0 )
    {
        throw "SF.DomainException";
    }
    else if ( a == 0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if ( b < 0 && (a < b || a > 0) )
    {
        // Standard domain error due to singularity.
        throw "SF.DomainException";
    }
    else if ( x > 100.0 && Math.max( 1.0, (Math.abs( b - a )) ) * Math.max( 1.0, (Math.abs( 1 - a )) ) < 0.5 * x )
    {
        // x . +Inf asymptotic
        return hyperg_1F1_asymp_posx( a, b, x );
    }
    else if ( x < -100.0 && Math.max( 1.0, (Math.abs( a )) ) * Math.max( 1.0, (Math.abs( 1 + a - b )) ) < 0.5 * Math.abs( x ) )
    {
        // x . -Inf asymptotic
        return hyperg_1F1_asymp_negx( a, b, x );
    }
    else if ( a < 0 && b < 0 )
    {
        return hyperg_1F1_ab_negint( a, b, x );
    }
    else if ( a < 0 && b > 0 )
    {
        // Use Kummer to reduce it to the positive integer case.
        // Note that b > a, strictly, since we already trapped b = a.
        //
        Kummer_1F1 = hyperg_1F1_ab_posint( b - a, b, -x );
        return gsl_sf_exp_mult_err_e( x, GSL_DBL_EPSILON * Math.abs( x ), Kummer_1F1.val, Kummer_1F1.err );
    } 
    else
    {
        // a > 0 and b > 0
        return hyperg_1F1_ab_posint( a, b, x );
    }

} // gsl_sf_hyperg_1F1_int_e

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_1F1_e( a, b, x )
{
    var bma     = 0.0;
    var rinta   = 0.0;
    var rintb   = 0.0;
    var rintbma = 0.0;

    var a_integer       = false;
    var b_integer       = false;
    var bma_integer     = false;
    var b_neg_integer   = false;
    var a_neg_integer   = false;
    var bma_neg_integer = false;

    const INTEGER_FIRST = -2147483648;
    const INTEGER_LAST  = +2147483647;

    var r = { val: 0.0, err: 0.0 }; // Result;

    bma     = b - a;
    rinta   = Math.floor( a + 0.5 );
    rintb   = Math.floor( b + 0.5 );
    rintbma = Math.floor( bma + 0.5 );
    a_integer   = ( Math.abs(a-rinta) < H1F1_INT_THRESHOLD && rinta > (INTEGER_FIRST) && rinta < (INTEGER_LAST) );
    b_integer   = ( Math.abs(b-rintb) < H1F1_INT_THRESHOLD && rintb > (INTEGER_FIRST) && rintb < (INTEGER_LAST) );
    bma_integer = ( Math.abs(bma-rintbma) < H1F1_INT_THRESHOLD && rintbma > (INTEGER_FIRST) && rintbma < (INTEGER_LAST) );
    b_neg_integer   = ( b < -0.1 && b_integer );
    a_neg_integer   = ( a < -0.1 && a_integer );
    bma_neg_integer = ( bma < -0.1 &&  bma_integer );

    if ( x == 0.0 )
    {
        // Testing for this before testing a and b
        // is somewhat arbitrary. The result is that
        // we have 1F1(a,0,0) = 1.
        //
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if ( b == 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( a == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if ( a == b )
    {
        // case: a=b; exp(x)
        // It's good to test exact equality now.
        // We also test approximate equality later.
        //
        return gsl_sf_exp_e( x );
    }
    else if ( Math.abs( b ) < H1F1_INT_THRESHOLD && Math.abs( a ) < H1F1_INT_THRESHOLD )
    {
        // a and b near zero: 1 + a/b (exp(x)-1)
        //
       
        // Note that neither a nor b is zero, since
        // we eliminated that with the above tests.
        //
        let exm1 = { val: 0.0, err: 0.0 }; // Result;
        let hx   = { val: 0.0, err: 0.0 }; // Result;
        let sa   = 0.0;
        let sb   = 0.0;
        let lnab = 0.0;

        exm1 = gsl_sf_expm1_e( x );
        sa = ( a > 0.0 ) ? 1.0: -1.0;
        sb = ( b > 0.0 ) ? 1.0: -1.0;
        lnab = Math.log( Math.abs( a / b ) ); // safe
        hx = gsl_sf_exp_mult_err_e( lnab, GSL_DBL_EPSILON * Math.abs( lnab ), sa * sb * exm1.val, exm1.err );
        if ( hx.val == GSL_DBL_MAX ) // FIXME: excessive paranoia ? what is DBL_MAX+1 ?
        {
            r.val = hx.val;
        }
        else
        {
            r.val = 1.0 + hx.val;
        }
        r.err = hx.err;
        return r;
    }
    else if ( Math.abs( b ) < H1F1_INT_THRESHOLD && Math.abs( x * a ) < 1.0 )
    {
        // b near zero and a not near zero
        //
        let m_arg    = 1.0 / (0.5 * b);
        let F_renorm = { val: 0.0, err: 0.0 }; // Result;

        F_renorm = hyperg_1F1_renorm_b0( a, x );
        r = gsl_sf_multiply_err_e( m_arg, 2.0 * GSL_DBL_EPSILON * m_arg, 0.5 * F_renorm.val, 0.5 * F_renorm.err );
        return r;
    }
    else if ( a_integer && b_integer )
    {
        // Check for reduction to the integer case.
        // Relies on the arbitrary "near an integer" test.
        //
        return gsl_sf_hyperg_1F1_int_e( rinta, rintb, x );
    }
    else if ( b_neg_integer && ! (a_neg_integer && a > b) )
    {
        // Standard domain error due to
        // uncancelled singularity.
        //
        throw "SF.DomainException";
    }
    else if ( a_neg_integer )
    {
        return hyperg_1F1_a_negint_lag( rinta, b, x );
    }
    else if ( b > 0.0 )
    {
        if ( -1.0 <= a && a <= 1.0 )
        {
            // Handle small a explicitly.
            //
            return hyperg_1F1_small_a_bgt0( a, b, x );
        }
        else if ( bma_neg_integer )
        {
            // Catch this now, to avoid problems in the
            // generic evaluation code.
            //
            let Kummer_1F1 = { val: 0.0, err: 0.0 }; // Result;

            Kummer_1F1 = hyperg_1F1_a_negint_lag( rintbma, b, -x );
            r = gsl_sf_exp_mult_err_e( x, GSL_DBL_EPSILON * Math.abs( x ), Kummer_1F1.val, Kummer_1F1.err );
            return r;
        }
        else if ( a < 0.0 && Math.abs( x ) < 100.0 )
        {
            // Use Kummer to reduce it to the generic positive case.
            // Note that b > a, strictly, since we already trapped b = a.
            // Also b-(b-a)=a, and a is not a negative integer here,
            // so the generic evaluation is safe.
            //
            let Kummer_1F1 = { val: 0.0, err: 0.0 }; // Result;

            Kummer_1F1 = hyperg_1F1_ab_pos( b - a, b, -x );
            r = gsl_sf_exp_mult_err_e( x, GSL_DBL_EPSILON * Math.abs( x ), Kummer_1F1.val, Kummer_1F1.err );
            return r;
        }
        else if ( a > 0.0 )
        {
            // a > 0.0
            return hyperg_1F1_ab_pos( a, b, x );
        }
        else
        {
            return gsl_sf_hyperg_1F1_series_e( a, b, x );
        }
    }
    else
    {
        // b < 0.0
  
        if ( bma_neg_integer && x < 0.0 )
        {
            // Handle this now to prevent problems
            // in the generic evaluation.
            //
            let K = { val: 0.0, err: 0.0 }; // Result;

            if ( a < 0.0 )
            {
                // Kummer transformed version of safe polynomial.
                // The condition a < 0 is equivalent to b < b-a,
                // which is the condition required for the series
                // to be positive definite here.
                //
                K = hyperg_1F1_a_negint_poly( rintbma, b, -x );
            }
            else
            {
                // Generic eval for negative integer a.
                K = hyperg_1F1_a_negint_lag( rintbma, b, -x );
            }
            r = gsl_sf_exp_mult_err_e( x, GSL_DBL_EPSILON * Math.abs( x ), K.val, K.err );
            return r;
        }
        else if ( a > 0.0 )
        {
            // Use Kummer to reduce it to the generic negative case.
            //
            let K = { val: 0.0, err: 0.0 }; // Result;

            K = hyperg_1F1_ab_neg( b - a, b, -x );
            r = gsl_sf_exp_mult_err_e( x, GSL_DBL_EPSILON * Math.abs( x ), K.val, K.err );
            return r;
        }
        else
        {
            return hyperg_1F1_ab_neg( a, b, x );
        }
    }

} // gsl_sf_hyperg_1F1_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hyperg_1F1_int( m, n, x )
{ // gsl_sf_hyperg_1F1_int
    return EVAL_RESULT_IID( gsl_sf_hyperg_1F1_int_e, { n: m, m: n, x: x }, "gsl_sf_hyperg_1F1_int" );
} // gsl_sf_hyperg_1F1_int

export function gsl_sf_hyperg_1F1( a, b, x )
{ // gsl_sf_hyperg_1F1
    return EVAL_RESULT_3D( gsl_sf_hyperg_1F1_e, { x: a, y: b, z: x}, "gsl_sf_hyperg_1F1" );
} // gsl_sf_hyperg_1F1

// ----------------------------------------------------------------------------
// EOF SF-Hypergeometric1F1.mjs
