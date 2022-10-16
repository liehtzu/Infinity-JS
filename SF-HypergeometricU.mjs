// SF-HypergeometricU.mjs
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

import { M_PI }                   from "./SF-Math.mjs";
import { M_LN10 }                 from "./SF-Math.mjs";
import { M_LNPI }                 from "./SF-Math.mjs";
import { M_SQRT2 }                from "./SF-Math.mjs";
import { GSL_IS_ODD }             from "./SF-Math.mjs";
import { GSL_SIGN }               from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }        from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MAX }       from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }       from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }   from "./SF-Machine.mjs";
import { gsl_sf_gamma_e }         from "./SF-Gamma.mjs";
import { gsl_sf_gammainv_e }      from "./SF-Gamma.mjs";
import { gsl_sf_lnfact_e }        from "./SF-Gamma.mjs";
import { gsl_sf_pow_int_e }       from "./SF-Power.mjs";
import { gsl_sf_poch_e }          from "./SF-Pochhammer.mjs";
import { gsl_sf_pochrel_e }       from "./SF-Pochhammer.mjs";
import { gsl_sf_exprel_e }        from "./SF-Exponential.mjs";
import { gsl_sf_exp_e }           from "./SF-Exponential.mjs";
import { gsl_sf_exp_err_e }       from "./SF-Exponential.mjs";
import { gsl_sf_exp_e10_e }       from "./SF-Exponential.mjs";
import { gsl_sf_result_smash_e }  from "./SF-Results.mjs";
import { gsl_sf_exp_mult_err_e10_e } from "./SF-Exponential.mjs";
import { gsl_sf_exp_err_e10_e }   from "./SF-Exponential.mjs";
import { gsl_sf_bessel_lnKnu_e}   from "./SF-BesselKnu.mjs";
import { gsl_sf_laguerre_n_e }    from "./SF-Laguerre.mjs";
import { gsl_sf_hyperg_U_large_b_e } from "./SF-Hypergeometric.mjs"

// ----------------------------------------------------------------------------

const INT_THRESHOLD = 1000.0 * GSL_DBL_EPSILON;
const INTEGER_LAST  = 2147483647;

function SERIES_EVAL_OK( a, b, x )
{ // SERIES_EVAL_OK
    return ((Math.abs( a ) < 5 && b < 5 && x < 2.0) || (Math.abs(a) <  10 && b < 10 && x < 1.0));
} // SERIES_EVAL_OK

function ASYMP_EVAL_OK( a, b, x )
{ // ASYMP_EVAL_OK
    return (Math.max( Math.abs( (a) ), 1.0 ) * Math.max( Math.abs( (1 + a - b) ), 1.0 ) < 0.99 * Math.abs( x ));
} // ASYMP_EVAL_OK

// Log[U(a,2a,x)]
// [Abramowitz+stegun, 13.6.21]
// Assumes x > 0, a > 1/2.
//
function hyperg_lnU_beq2a( a, x )
{
    var lx    = 0.0;
    var nu    = 0.0;
    var lnpre = 0.0;
    var lnK   = { val: 0.0, err: 0.0 }; // Result;
    var r     = { val: 0.0, err: 0.0 }; // Result;

    lx = Math.log( x );
    nu = a - 0.5;
    lnpre = 0.5 * (x - M_LNPI) - nu * lx;
    lnK = gsl_sf_bessel_lnKnu_e( nu, 0.5 * x );
    r.val = lnpre + lnK.val;
    r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( 0.5 * x ) + 0.5 * M_LNPI + Math.abs( nu * lx ));
    r.err = r.err + lnK.err;
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
    return r;

} // hyperg_lnU_beq2a

// ----------------------------------------------------------------------------

// Evaluate u_{N+1}/u_N by Steed's continued fraction method.
//
// u_N = Gamma[a+N]/Gamma[a] U(a + N, b, x)
//
// u_{N+1}/u_N = (a+N) U(a+N+1,b,x)/U(a+N,b,x)
//
// PROCEDURE hyperg_U_CF1(a: LONG_FLOAT; b: LONG_FLOAT; N: INTEGER; x: LONG_FLOAT; r: IN OUT LONG_FLOAT; count: IN OUT INTEGER) IS
function hyperg_U_CF1( a, b, N, x ) // r: IN OUT LONG_FLOAT; count: IN OUT INTEGER) IS
{
    const RECUR_BIG = GSL_SQRT_DBL_MAX;
    const maxiter   = 20000;
    var n0        = 1;
    var Anm2      = 0.0;
    var Bnm2      = 0.0;
    var Anm1      = 0.0;
    var Bnm1      = 0.0;
    var a1        = 0.0;
    var b1        = 0.0;
    var An        = 0.0;
    var Bn        = 0.0;
    var ann       = 0.0;
    var bnn       = 0.0;
    var fn        = 0.0;
    var old_fn    = 0.0;
    var del       = 0.0;

    var r = { val: 0.0, err: 0.0, count: 0 };

    n0 = 1;
    Anm2 = 1.0;
    Bnm2 = 0.0;
    Anm1 = 0.0;
    Bnm1 = 1.0;
    a1 = -(a + (N));
    b1 =  (b - 2.0 * a - x - 2.0 * (N + 1));
    An = b1 * Anm1 + a1 * Anm2;
    Bn = b1 * Bnm1 + a1 * Bnm2;
    fn = An / Bn;
  
    while ( n0 < maxiter )
    {
        n0   = n0 + 1;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        ann  = -(a + (N + n0) - b) * (a + (N + n0 - 1));
        bnn  =  (b - 2.0 * a - x - 2.0 * (N + n0));
        An   = bnn * Anm1 + ann * Anm2;
        Bn   = bnn * Bnm1 + ann * Bnm2;
        
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
        
        if ( Math.abs( del - 1.0 ) < 10.0 * GSL_DBL_EPSILON ) break;
    }
    
    r.val   = fn;
    r.count = n0;
  
    if ( n0 >= maxiter )
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // hyperg_U_CF1

// ----------------------------------------------------------------------------

// Large x asymptotic for  x^a U(a,b,x)
// Based on SLATEC D9CHU() [W. Fullerton]
//
// Uses a rational approximation due to Luke.
// See [Luke, Algorithms for the Computation of Special Functions, p. 252]
//     [Luke, Utilitas Math. (1977)]
//
// z^a U(a,b,z) ~ 2F0(a,1+a-b,-1/z)
//
// This assumes that a is not a negative integer and
// that 1+a-b is not a negative integer. If one of them
// is, then the 2F0 actually terminates, the above
// relation is an equality, and the sum should be
// evaluated directly [see below].
//
function d9chu( a, b, x )
{
    const EPS     = 8.0 * GSL_DBL_EPSILON;  // EPS = 4.0D0*D1MACH(4)
    const maxiter = 500;
    var aa      = []; //: ARRAY (0..3) OF LONG_FLOAT;
    var bb      = []; //: ARRAY (0..3) OF LONG_FLOAT;
    var i       = 0;
  
    var bp   = 1.0 + a - b;
    var ab   = a * bp;
    var ct2  = 2.0 * (x - ab);
    var sab  = a + bp;
    var ct3  = sab + 1.0 + ab;
    var anbn = ct3 + sab + 3.0;
    var ct1  = 1.0 + 2.0 * x / anbn;

    var c2   = 0.0;
    var g1   = 0.0;
    var g2   = 0.0;
    var g3   = 0.0;
    var d1z  = 0.0;
    var x2i1 = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    bb[0] = 1.0;
    aa[0] = 1.0;
  
    bb[1] = 1.0 + 2.0 * x / ct3;
    aa[1] = 1.0 + ct2 / ct3;
    
    bb[2] = 1.0 + 6.0 * ct1 * x / ct3;
    aa[2] = 1.0 + 6.0 * ab / anbn + 3.0 * ct1 * ct2 / ct3;
  
    i = 4;
    while ( i <= maxiter - 1 )
    {
        x2i1 = (2 * i - 3);
        ct1  = x2i1 / (x2i1 - 2.0);
        anbn = anbn + x2i1 + sab;
        ct2  = (x2i1 - 1.0) / anbn;
        c2   = x2i1 * ct2 - 1.0;
        d1z  = 2.0 * x2i1 * x / anbn;
        
        ct3 = sab * ct2;
        g1  = d1z + ct1 * (c2 + ct3);
        g2  = d1z - c2;
        g3  = ct1 * (1.0 - ct3 - 2.0 * ct2);
        
        bb[3] = g1 * bb[2] + g2 * bb[1] + g3 * bb[0];
        aa[3] = g1 * aa[2] + g2 * aa[1] + g3 * aa[0];
        
        if ( Math.abs( aa[3] * bb[0] - aa[0] * bb[3] ) < EPS * Math.abs( bb[3] * bb[0] ) ) break;
        
        for ( let j = 0; j <= 3 - 1; j++ )
        {
            aa[j] = aa[j+1];
            bb[j] = bb[j+1];
        }
        i = i + 1;
    }
    
    r.val = aa[3] / bb[3];
    r.err = 8.0 * GSL_DBL_EPSILON * Math.abs( r.val );
    
    if ( i >= maxiter )
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // d9chu

// ----------------------------------------------------------------------------

// Evaluate asymptotic for z^a U(a,b,z) ~ 2F0(a,1+a-b,-1/z)
// We check for termination of the 2F0 as a special case.
// Assumes x > 0.
// Also assumes a,b are not too large compared to x.
//
function hyperg_zaU_asymp( a, b, x )
{
    var ap = a;
    var bp = 1.0 + a - b;
    var rintap = Math.floor( ap + 0.5 );
    var rintbp = Math.floor( bp + 0.5 );
    var ap_neg_int = ( ap < 0.0 && Math.abs( ap - rintap ) < INT_THRESHOLD );
    var bp_neg_int = ( bp < 0.0 && Math.abs( bp - rintbp ) < INT_THRESHOLD );

    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( ap_neg_int || bp_neg_int )
    {
        // Evaluate 2F0 polynomial.
        //
        let mxi  = -1.0 / x;
        let nmax = -Math.trunc( Math.min( ap, bp ) - 0.1 );
        let tn   = 1.0;
        let sum  = 1.0;
        let n    = 1.0;
        let sum_err = 0.0;
        let apn  = 0.0;
        let bpn  = 0.0;

        while ( n <= nmax )
        {
            apn = (ap + n - 1.0);
            bpn = (bp + n - 1.0);
            tn  = tn * ((apn / n) * mxi) * bpn;
            sum = sum + tn;
            sum_err = sum_err + 2.0 * GSL_DBL_EPSILON * Math.abs(tn);
            n = n + 1.0;
        }
        r.val = sum;
        r.err = sum_err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( nmax ) + 1.0) * Math.abs( sum );
        return r;
    }
    else
    {
        return d9chu( a, b, x );
    }

} // hyperg_zaU_asymp

// -- ----------------------------------------------------------------------------

// Evaluate finite sum which appears below.
//
function hyperg_U_finite_sum( N, a, b, x, xeps )
{
    var sum_val = 0.0;
    var sum_err = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( N <= 0 )
    {
        let t_val = 1.0;
        let t_err = 0.0;
        let xi1   = 0.0;
        let mult  = 0.0;
        let poch  = { val: 0.0, err: 0.0 }; // Result;

        sum_val = 1.0;
        sum_err = 0.0;
        for ( let i = 1; i <= -N; i++ )
        {
            xi1  = (i - 1);
            mult = (a + xi1) * x / ((b + xi1) * (xi1 + 1.0));
            t_val = t_val * mult;
            t_err = t_err + Math.abs( mult ) * t_err + Math.abs( t_val ) * 8.0 * 2.0 * GSL_DBL_EPSILON;
            sum_val = sum_val + t_val;
            sum_err = sum_err + t_err;
        }
        
        poch = gsl_sf_poch_e( 1.0 + a - b, -a );
        
        r.val = sum_val * poch.val;
        r.err = Math.abs( sum_val ) * poch.err + sum_err * Math.abs( poch.val );
        r.err = r.err + Math.abs( poch.val ) * (Math.abs( N ) + 2) * GSL_DBL_EPSILON * Math.abs( sum_val );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        r.err = r.err * 2.0; // FIXME: fudge factor... why is the error estimate too small?
        return r;
    }
    else
    {
        var M = N - 2;

        if ( M < 0 )
        {
            r.val = 0.0;
            r.err = 0.0;
            return r;
        }
        else
        {
            let gbm1  = { val: 0.0, err: 0.0 }; // Result;
            let gamr  = { val: 0.0, err: 0.0 }; // Result;
            let t_val = 1.0;
            let t_err = 0.0;
            let mult  = 0.0;

            sum_val = 1.0;
            sum_err = 0.0;
            for ( let i = 1; i <= M; i++ )
            {
                mult = (a - b + (i)) * x / ((1.0 - b + (i)) * (i));
                t_val = t_val * mult;
                t_err = t_err + t_err * Math.abs( mult ) + Math.abs( t_val ) * 8.0 * 2.0 * GSL_DBL_EPSILON;
                sum_val = sum_val + t_val;
                sum_err = sum_err + t_err;
            }
            
            gbm1 = gsl_sf_gamma_e( b - 1.0 );
            gamr = gsl_sf_gammainv_e( a );
            
            let powx1N = { val: 0.0, err: 0.0 }; // Result;
            let pe_val = 0.0;
            let pe_err = 0.0;
            let coeff_val = 0.0;
            let coeff_err = 0.0;

            powx1N = gsl_sf_pow_int_e( x, 1 - N );
            pe_val = powx1N.val * xeps;
            pe_err = powx1N.err * Math.abs( xeps ) + 2.0 * GSL_DBL_EPSILON * Math.abs( pe_val );
            coeff_val = gbm1.val * gamr.val * pe_val;
            coeff_err = gbm1.err * Math.abs( gamr.val * pe_val )
                        + gamr.err * Math.abs( gbm1.val * pe_val )
                        + Math.abs( gbm1.val * gamr.val ) * pe_err
                        + 2.0 * GSL_DBL_EPSILON * Math.abs( coeff_val );
            
            r.val = sum_val * coeff_val;
            r.err = Math.abs( sum_val ) * coeff_err + sum_err * Math.abs( coeff_val );
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * ( M + 2 ) * Math.abs( r.val );
            r.err = r.err * 2.0; // FIXME: fudge factor... why is the error estimate too small?
            return r;
        }
    }

} // hyperg_U_finite_sum

// ----------------------------------------------------------------------------

// Based on SLATEC DCHU() [W. Fullerton]
// Assumes x > 0.
// This is just a series summation method, and
// it is not good for large a.
//
// I patched up the window for 1+a-b near zero. [GJ]
//
function hyperg_U_series( a, b, x )
{
    const EPS      = 2.0 * GSL_DBL_EPSILON;  // EPS = D1MACH(3)
    const SQRT_EPS = M_SQRT2 * GSL_SQRT_DBL_EPSILON;

    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( Math.abs( 1.0 + a - b ) < SQRT_EPS )
    {
        // Original Comment: ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X
        //
        // We can however do the following:
        // U(a,b,x) = U(a,a+1,x) when 1+a-b=0
        // and U(a,a+1,x) = x^(-a).
        //
        var lnr = -a * Math.log( x );

        r = gsl_sf_exp_e( lnr );
        r.err = r.err + 2.0 * SQRT_EPS * Math.abs( r.val );
        return r;
    }
    else
    {
        let aintb       = 0.0;
        let beps        = 0.0;
        let lnx         = 0.0;
        let xeps        = 0.0;
        let xi          = 0.0;
        let sarg        = 0.0;
        let sfact       = 0.0;
        let factor_val  = 0.0;
        let factor_err  = 0.0;
        let b0_val      = 0.0;
        let b0_err      = 0.0;
        let N           = 0;
        let istrt       = 0;
        let sum         = { val: 0.0, err: 0.0 }; // Result;
        let gamr        = { val: 0.0, err: 0.0 }; // Result;
        let powx        = { val: 0.0, err: 0.0 }; // Result;
        let pochai      = { val: 0.0, err: 0.0 }; // Result;
        let gamri1      = { val: 0.0, err: 0.0 }; // Result;
        let gamrni      = { val: 0.0, err: 0.0 }; // Result;
        let pochaxibeps = { val: 0.0, err: 0.0 }; // Result;
        let gamrxi1beps = { val: 0.0, err: 0.0 }; // Result;

        if ( b < 0.0 )
        {
            aintb =  Math.ceiling( b - 0.5 );
        }
        else
        {
            aintb =  Math.floor( b + 0.5 );
        }
        beps = b - aintb;
        N = Math.trunc( aintb );

        lnx  = Math.log( x );
        xeps = Math.exp( -beps * lnx );

        // Evaluate finite sum.
        //
        sum = hyperg_U_finite_sum( N, a, b, x, xeps );

        // Evaluate infinite sum.
        //
        if ( N < 1 )
        {
            istrt = 1 - N;
        }
        else
        {
            istrt = 0;
        }
        xi = (istrt);

        gamr = gsl_sf_gammainv_e( 1.0 + a - b );
        powx = gsl_sf_pow_int_e( x, istrt );
        sarg = beps * M_PI;
        if ( sarg != 0.0 )
        {
            sfact = sarg / Math.sin( sarg );
        }
        else
        {
            sfact = 1.0;
        }
        factor_val = sfact * gamr.val * powx.val;
        if ( GSL_IS_ODD( N ) )
        {
            factor_val = -factor_val;
        }
        factor_err = Math.abs( gamr.val ) * powx.err + Math.abs( powx.val ) * gamr.err
                    + 2.0 * GSL_DBL_EPSILON * Math.abs( factor_val );

        pochai = gsl_sf_poch_e( a, xi );
        gamri1 = gsl_sf_gammainv_e( xi + 1.0 );
        gamrni = gsl_sf_gammainv_e( aintb + xi );

        pochaxibeps = gsl_sf_poch_e( a, xi - beps );
        gamrxi1beps = gsl_sf_gammainv_e( xi + 1.0 - beps );

        b0_val = factor_val * pochaxibeps.val * gamrni.val * gamrxi1beps.val;
        b0_err =  Math.abs( factor_val * pochaxibeps.val * gamrni.val ) * gamrxi1beps.err
                    + Math.abs( factor_val * pochaxibeps.val * gamrxi1beps.val ) * gamrni.err
                    + Math.abs( factor_val * gamrni.val * gamrxi1beps.val ) * pochaxibeps.err
                    + Math.abs( pochaxibeps.val * gamrni.val * gamrxi1beps.val ) * factor_err
                    + 2.0 * GSL_DBL_EPSILON * Math.abs( b0_val );

        if ( Math.abs( xeps - 1.0 ) < 0.5 )
        {
            //
            // X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE
            // CAREFUL IN EVALUATING THE DIFFERENCES.
            //
            let i = 0;
            let pch1ai       = { val: 0.0, err: 0.0 }; // Result;
            let pch1i        = { val: 0.0, err: 0.0 }; // Result;
            let poch1bxibeps = { val: 0.0, err: 0.0 }; // Result;
            let dexprl       = { val: 0.0, err: 0.0 }; // Result;
            let t_val        = 0.0;
            let t_err        = 0.0;
            let xn           = 0.0;
            let c0_t1_val    = 0.0;
            let c0_t1_err    = 0.0;
            let c0_t2_val    = 0.0;
            let c0_t2_err    = 0.0;
            let c0_val       = 0.0;
            let c0_err       = 0.0;
            let xeps1_val    = 0.0;
            let xeps1_err    = 0.0;
            let dchu_val     = 0.0;
            let dchu_err     = 0.0;

            pch1ai = gsl_sf_pochrel_e( a + xi, -beps );
            pch1i = gsl_sf_pochrel_e( xi + 1.0 - beps, beps );
            poch1bxibeps = gsl_sf_pochrel_e( b + xi, -beps );
            c0_t1_val = beps * pch1ai.val * pch1i.val;
            c0_t1_err = Math.abs( beps ) * Math.abs( pch1ai.val ) * pch1i.err
                        + Math.abs( beps ) * Math.abs( pch1i.val )  * pch1ai.err
                        + 2.0 * GSL_DBL_EPSILON * Math.abs( c0_t1_val );
            c0_t2_val = -poch1bxibeps.val + pch1ai.val - pch1i.val + c0_t1_val;
            c0_t2_err =  poch1bxibeps.err + pch1ai.err + pch1i.err + c0_t1_err
                        + 2.0 * GSL_DBL_EPSILON * Math.abs( c0_t2_val );
            c0_val = factor_val * pochai.val * gamrni.val * gamri1.val * c0_t2_val;
            c0_err =  Math.abs( factor_val * pochai.val * gamrni.val * gamri1.val ) * c0_t2_err
                    + Math.abs( factor_val * pochai.val * gamrni.val * c0_t2_val ) * gamri1.err
                    + Math.abs( factor_val * pochai.val * gamri1.val * c0_t2_val ) * gamrni.err
                    + Math.abs( factor_val * gamrni.val * gamri1.val * c0_t2_val ) * pochai.err
                    + Math.abs( pochai.val * gamrni.val * gamri1.val * c0_t2_val ) * factor_err
                    + 2.0 * GSL_DBL_EPSILON * Math.abs( c0_val );
            //
            //  XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
            //
            dexprl = gsl_sf_exprel_e( -beps * lnx );
            xeps1_val = lnx * dexprl.val;
            xeps1_err = 2.0 * GSL_DBL_EPSILON * (1.0 + Math.abs( beps * lnx )) * Math.abs( dexprl.val )
                        + Math.abs( lnx ) * dexprl.err
                        + 2.0 * GSL_DBL_EPSILON * Math.abs( xeps1_val );
            dchu_val = sum.val + c0_val + xeps1_val*b0_val;
            dchu_err = sum.err + c0_err
                        + Math.abs( xeps1_val ) * b0_err + xeps1_err * Math.abs( b0_val )
                        + Math.abs( b0_val * lnx ) * dexprl.err
                        + 2.0 * GSL_DBL_EPSILON * (Math.abs( sum.val ) + Math.abs( c0_val ) + Math.abs( xeps1_val * b0_val ));
            xn = (N);

            i = 1;
            while ( i <= 2000 - 1 )
            {
                let xi              = 0.0;
                let xi1             = 0.0;
                let tmp             = 0.0;
                let b0_multiplier   = 0.0;
                let c0_multiplier_1 = 0.0;
                let c0_multiplier_2 = 0.0;

                xi  = (istrt + i);
                xi1 = (istrt + i - 1);
                tmp = (a - 1.0) * (xn + 2.0 * xi - 1.0) + xi * (xi - beps);
                b0_multiplier = (a + xi1 - beps) * x / ((xn + xi1) * (xi - beps));
                c0_multiplier_1 = (a + xi1) * x / ((b + xi1) * xi);
                c0_multiplier_2 = tmp / (xi * (b + xi1) * (a + xi1 - beps));
                b0_val = b0_val * b0_multiplier;
                b0_err = b0_err + Math.abs( b0_multiplier ) * b0_err + Math.abs(b0_val) * 8.0 * 2.0 * GSL_DBL_EPSILON;
                c0_val = c0_multiplier_1 * c0_val - c0_multiplier_2 * b0_val;
                c0_err =  Math.abs( c0_multiplier_1 ) * c0_err
                            + Math.abs( c0_multiplier_2 ) * b0_err
                            + Math.abs( c0_val ) * 8.0 * 2.0 * GSL_DBL_EPSILON
                            + Math.abs( b0_val * c0_multiplier_2 ) * 16.0 * 2.0 * GSL_DBL_EPSILON;
                t_val = c0_val + xeps1_val * b0_val;
                t_err = c0_err + Math.abs( xeps1_val ) * b0_err;
                t_err = t_err + Math.abs( b0_val * lnx ) * dexprl.err;
                t_err = t_err + Math.abs( b0_val ) * xeps1_err;
                dchu_val = dchu_val + t_val;
                dchu_err = dchu_err + t_err;
                if ( Math.abs( t_val ) < EPS * Math.abs( dchu_val ) ) break;
                i = i + 1;
            }

            r.val = dchu_val;
            r.err = 2.0 * dchu_err;
            r.err = r.err + 2.0 * Math.abs( t_val );
            r.err = r.err + 4.0 * GSL_DBL_EPSILON * (i + 2) * Math.abs( dchu_val );
            r.err = r.err * 2.0; // FIXME: fudge factor

            if ( i >= 2000 )
            {
                throw "SF.MaxIterationsException";
            }

            return r;
        }
        else
        {
            //
            // X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE
            // STRAIGHTFORWARD FORMULATION IS STABLE.
            //
            let i = 0;
            let dchu_val = 0.0;
            let dchu_err = 0.0;
            let t_val    = 0.0;
            let t_err    = 0.0;
            let a0_val   = 0.0;
            let a0_err   = 0.0;
            let dgamrbxi = { val: 0.0, err: 0.0 }; // Result;

            dgamrbxi = gsl_sf_gammainv_e( b + xi );
            a0_val = factor_val * pochai.val * dgamrbxi.val * gamri1.val / beps;
            a0_err =  Math.abs( factor_val * pochai.val * dgamrbxi.val / beps ) * gamri1.err
                    + Math.abs( factor_val * pochai.val * gamri1.val / beps ) * dgamrbxi.err
                    + Math.abs( factor_val * dgamrbxi.val * gamri1.val / beps ) * pochai.err
                    + Math.abs( pochai.val * dgamrbxi.val * gamri1.val / beps ) * factor_err
                    + 2.0 * GSL_DBL_EPSILON * Math.abs( a0_val );

            b0_val = xeps * b0_val / beps;
            b0_err = Math.abs( xeps / beps ) * b0_err + 4.0 * GSL_DBL_EPSILON * Math.abs( b0_val );
            dchu_val = sum.val + a0_val - b0_val;
            dchu_err = sum.err + a0_err + b0_err
                        + 2.0 * GSL_DBL_EPSILON * (Math.abs( sum.val ) + Math.abs( a0_val ) + Math.abs( b0_val ));

            i = 1;
            while ( i <= 2000 - 1 )
            {
                let xi  = 0.0;
                let xi1 = 0.0;
                let a0_multiplier = 0.0;
                let b0_multiplier = 0.0;

                xi = (istrt + i);
                xi1 = (istrt + i - 1);
                a0_multiplier = (a + xi1) * x / ((b + xi1) * xi);
                b0_multiplier = (a + xi1 - beps) * x / ((aintb + xi1) * (xi - beps));
                a0_val = a0_val * a0_multiplier;
                a0_err = a0_err + Math.abs( a0_multiplier ) * a0_err;
                b0_val = b0_val * b0_multiplier;
                b0_err = b0_err + Math.abs( b0_multiplier ) * b0_err;
                t_val = a0_val - b0_val;
                t_err = a0_err + b0_err;
                dchu_val = dchu_val + t_val;
                dchu_err = dchu_err + t_err;
                if ( Math.abs( t_val ) < EPS * Math.abs( dchu_val ) ) break;
                i = i + 1;
            }

            r.val = dchu_val;
            r.err = 2.0 * dchu_err;
            r.err = r.err + 2.0 * Math.abs( t_val );
            r.err = r.err + 4.0 * GSL_DBL_EPSILON * (i + 2) * Math.abs( dchu_val );
            r.err = r.err * 2.0; // FIXME: fudge factor

            if ( i >= 2000 )
            {
                throw "SF.MaxIterationsException";
            }

            return r;
        }
    }

} // hyperg_U_series

// ----------------------------------------------------------------------------

// Assumes b > 0 and x > 0.
//
function hyperg_U_small_ab( a, b, x )
{
    var p     = 0.0;
    var asymp = { val: 0.0, err: 0.0 }; // Result;
    var r     = { val: 0.0, err: 0.0 }; // Result;

    if ( a == -1.0 )
    {
        // U(-1,c+1,x) = Laguerre[c,0,x] = -b + x
        //
        r.val = -b + x;
        r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( b ) + Math.abs( x ));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else if ( a == 0.0 )
    {
        // U(0,c+1,x) = Laguerre[c,0,x] = 1
        //
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else if ( ASYMP_EVAL_OK( Math.trunc( a ), Math.trunc( b ), x ) )
    {
        p = Math.pow( x, -a );
        asymp = hyperg_zaU_asymp( a, b, x );
        r.val = asymp.val * p;
        r.err = asymp.err * p;
        r.err = r.err +Math.abs( asymp.val ) * GSL_DBL_EPSILON * Math.abs( a ) * p;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else
    {
        return hyperg_U_series( a, b, x );
    }

} // hyperg_U_small_ab

// ----------------------------------------------------------------------------

// Assumes b > 0 and x > 0.
//
//function hyperg_U_small_a_bgt0( a: LONG_FLOAT; b: LONG_FLOAT; x: LONG_FLOAT; r: IN OUT Result; ln_multiplier: IN OUT LONG_FLOAT) IS
function hyperg_U_small_a_bgt0( a, b, x ) //, r, ln_multiplier )
{
    var r = { val: 0.0, err: 0.0, lnm: 0.0 };

    if ( a == 0.0 )
    {
        r.val = 1.0;
        r.err = 1.0;
        r.lnm = 0.0;
        //ln_multiplier = 0.0;
        return r;
    }
    else if ( (b > 5000.0 && x < 0.90 * Math.abs( b )) || (b >  500.0 && x < 0.50 * Math.abs( b )) )
    {
        return gsl_sf_hyperg_U_large_b_e( a, b, x ); //, r, ln_multiplier );
    }
    else if ( b > 15.0 )
    {
        // Recurse up from b near 1.
        //
        var eps    = b - Math.floor( b );
        var b0     = 1.0 + eps;
        var r_Ubm1 = { val: 0.0, err: 0.0 }; // Result;
        var r_Ub   = { val: 0.0, err: 0.0 }; // Result;
        var Ubm1   = 0.0;
        var Ub     = 0.0;
        var Ubp1   = 0.0;
        var bp     = 0.0;

        r_Ubm1 = hyperg_U_small_ab( a, b0,       x );
        r_Ub   = hyperg_U_small_ab( a, b0 + 1.0, x );
        Ubm1   = r_Ubm1.val;
        Ub     = r_Ub.val;

        bp = b0 + 1.0;
        while ( bp < b - 0.1 )
        {
            Ubp1 = ((1.0 + a - bp) * Ubm1 + (bp + x - 1.0) * Ub) / x;
            Ubm1 = Ub;
            Ub   = Ubp1;
            bp   = bp + 1.0;
        }
        r.val = Ub;
        r.err = (Math.abs( r_Ubm1.err / r_Ubm1.val ) + Math.abs( r_Ub.err / r_Ub.val )) * Math.abs( Ub );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( b - b0 ) + 1.0) * Math.abs( Ub );
        r.lnm = 0.0;
        //ln_multiplier = 0.0;
        return r;
    }
    else
    {
        //ln_multiplier = 0.0;
        r = hyperg_U_small_ab( a, b, x );
        r.lnm = 0.0;
        return r;
    }

} // hyperg_U_small_a_bgt0

// ----------------------------------------------------------------------------

// We use this to keep track of large
// dynamic ranges in the recursions.
// This can be important because sometimes
// we want to calculate a very large and
// a very small number and the answer is
// the product, of order 1. This happens,
// for instance, when we apply a Kummer
// transform to make b positive and
// both x and b are large.
//
function RESCALE_2( u0, u1, factor, count )
{
    var au0 = Math.abs( u0 );
    var r = { u0: 0.0, u1: 0.0, count: 0 };

    if ( au0 > factor )
    {
        r.u0 = u0 / factor;
        r.u1 = u1 / factor;
        r.count = count + 1;
    }
    else if ( au0 < 1.0 / factor )
    {
        r.u0 = u0 * factor;
        r.u1 = u1 * factor;
        r.count = count - 1;
    }
    else
    {
        r.u0 = u0;
        r.u1 = u1;
        r.count = count;
    }
    return r;

} // RESCALE_2

// ----------------------------------------------------------------------------

// Specialization to b >= 1, for integer parameters.
// Assumes x > 0.
//
export function hyperg_U_int_bge1( a, b, x )
{
    var r = { val: 0.0, err: 0.0, e10: 0 }; // Result;

    if ( a == 0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        r.e10 = 0;
        return r;
    }
    else if ( a == -1 )
    {
        r.val = -(b) + x;
        r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( b ) + Math.abs( x ));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        r.e10 = 0;
        return r;
    }
    else if ( b == a + 1 )
    {
        // U(a,a+1,x) = x^(-a)
        //
        return gsl_sf_exp_e10_e( -(a) * Math.log( x ) );
    }
    else if ( ASYMP_EVAL_OK( a, b, x ) )
    {
        var ln_pre_val = -(a) * Math.log( x );
        var ln_pre_err = 2.0 * GSL_DBL_EPSILON * Math.abs( ln_pre_val );
        var asymp = { val: 0.0, err: 0.0 }; // Result;

        asymp = hyperg_zaU_asymp( (a), (b), x );
        r = gsl_sf_exp_mult_err_e10_e( ln_pre_val, ln_pre_err, asymp.val, asymp.err );
        return r;
    }
    else if ( SERIES_EVAL_OK( a, b, x ) )
    {
        let ser = { val: 0.0, err: 0.0 }; // Result;

        ser = hyperg_U_series( (a), (b), x );
        r.val = ser.val;
        r.err = ser.err;
        r.e10 = 0;
        return r;
    }
    else if ( a < 0 )
    {
        // Recurse backward from a = -1,0.
        //
        let scale_count  = 0;
        let scale_factor = GSL_SQRT_DBL_MAX;
        let lnm     = { val: 0.0, err: 0.0 }; // Result;
        let y       = { val: 0.0, err: 0.0 }; // Result;
        let lnscale = 0.0;
        let Uap1    = 1.0;     // U(0,b,x)
        let Ua      = -(b) + x;  // U(-1,b,x)
        let Uam1    = 0.0;
        let rr = { u0: 0.0, u1: 0.0, count: 0 };

        for ( let ap = -1; ap > a; ap-- )
        {
            Uam1 = (ap) * (b - ap - 1) * Uap1 + (x + (2 * ap - b)) * Ua;
            Uap1 = Ua;
            Ua   = Uam1;
            rr = RESCALE_2( Ua, Uap1, scale_factor, scale_count );
            Ua = rr.u0;
            Uap1 = rr.u1;
            scale_count = rr.count;
        }
    
        lnscale = Math.log( scale_factor );
        lnm.val = (scale_count) * lnscale;
        lnm.err = 2.0 * GSL_DBL_EPSILON * Math.abs( lnm.val );
        y.val = Ua;
        y.err = 4.0 * GSL_DBL_EPSILON * (Math.abs( (a) ) + 1.0) * Math.abs( Ua );
        r = gsl_sf_exp_mult_err_e10_e( lnm.val, lnm.err, y.val, y.err );
        return r;
    }
    else if ( (b) >= 2.0 * (a) + x )
    {
        // Recurse forward from a = 0,1.
        //
        let scale_count  = 0;
        let scale_factor = GSL_SQRT_DBL_MAX;
        let lnscale      = 0.0;
        let lm           = 0.0;
        let Uam1         = 0.0;
        let Ua           = 0.0;
        let Uap1         = 0.0;
        let r_Ua         = { val: 0.0, err: 0.0 }; // Result;
        let lnm          = { val: 0.0, err: 0.0 }; // Result;
        let y            = { val: 0.0, err: 0.0 }; // Result;
        let rlm          = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
        let rr = { u0: 0.0, u1: 0.0, count: 0 };

        rlm = hyperg_U_small_a_bgt0( 1.0, (b), x ); //, r_Ua, lm );  // U(1,b,x)
        r_Ua.val = rlm.val;
        r_Ua.err = rlm.err;
        lm = rlm.lnm;
        Uam1 = 1.0;  // U(0,b,x)
        Ua   = r_Ua.val;
        Uam1 = Uam1 * Math.exp( -lm );
    
        for ( let ap = 1; ap <= a - 1; ap++ )
        {
            Uap1 = -(Uam1 + ((b - 2 * ap) - x) * Ua) / ((ap) * (1 + ap - b));
            Uam1 = Ua;
            Ua   = Uap1;
            rr = RESCALE_2( Ua, Uam1, scale_factor, scale_count );
            Ua = rr.u0;
            Uam1 = rr.u1;
            scale_count = rr.count;
        }
    
        lnscale = Math.log( scale_factor );
        lnm.val = lm + (scale_count) * lnscale;
        lnm.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( lm ) + Math.abs( (scale_count) * lnscale) );
        y.val = Ua;
        y.err = Math.abs( r_Ua.err / r_Ua.val ) * Math.abs( Ua );
        y.err = y.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( a ) + 1) * Math.abs( Ua );
        r = gsl_sf_exp_mult_err_e10_e( lnm.val, lnm.err, y.val, y.err );
        return r;
    }
    else
    {
        if ( (b) <= x )
        {
            // Recurse backward either to the b=a+1 line
            // or to a=0, whichever we hit.
            //
            let scale_factor = GSL_SQRT_DBL_MAX;
            let scale_count  = 0;
            let ru           = 0.0;
            let CF1_count    = 0;
            let a_target     = 0;
            let lnU_target   = 0.0;
            let Ua           = 0.0;
            let Uap1         = 0.0;
            let Uam1         = 0.0;

            let rur = { val: 0.0, err: 0.0, count: 0 };
            let rr = { u0: 0.0, u1: 0.0, count: 0 };

            if ( b < a + 1 )
            {
                a_target = b - 1;
                lnU_target = -(a_target) * Math.log( x );
            }
            else
            {
                a_target = 0;
                lnU_target = 0.0;
            }
               
            rur = hyperg_U_CF1( (a), (b), 0, x ); //, ru, CF1_count );
            ru = rur.val;
            CF1_count = rur.count;
            
            Ua   = 1.0;
            Uap1 = ru / (a) * Ua;
            for ( let ap = a; ap >= a_target + 1; ap-- )
            {
                Uam1 = -(((b - 2 * ap) - x) * Ua + (ap) * (1 + ap - b) * Uap1);
                Uap1 = Ua;
                Ua   = Uam1;
                rr = RESCALE_2( Ua, Uap1, scale_factor, scale_count );
                Ua = rr.u0;
                Uap1 = rr.u1;
                scale_count = rr.count;
            }
            
            if ( Ua == 0.0 )
            {
                //result.val = 0.0;
                //result.err = 0.0;
                //result.e10 = 0;
                throw "SF.DivideByZeroException";
            }
            else
            {
                let lnscl = -(scale_count) * Math.log( scale_factor );
                let lnpre_val = lnU_target + lnscl;
                let lnpre_err = 2.0 * GSL_DBL_EPSILON * (Math.abs( lnU_target ) + Math.abs( lnscl ));
                let oUa_err   = 2.0 * (Math.abs( a_target - a ) + CF1_count + 1) * GSL_DBL_EPSILON * Math.abs( 1.0 / Ua );

                r = gsl_sf_exp_mult_err_e10_e( lnpre_val, lnpre_err, 1.0 / Ua, oUa_err );
                return r;
            }
        }
        else
        {
            // Recurse backward to near the b=2a+x line, then
            // determine normalization by either direct evaluation
            // or by a forward recursion. The direct evaluation
            // is needed when x is small (which is precisely
            // when it is easy to do).
            //
            let scale_factor = GSL_SQRT_DBL_MAX;
            let scale_count_for = 0;
            let scale_count_bck = 0;
            let a0 = 1;
            let a1 = a0 + Math.trunc( Math.ceil( 0.5 * ((b) - x) - (a0) ) );
            let Ua1_bck_val = 0.0;
            let Ua1_bck_err = 0.0;
            let Ua1_for_val = 0.0;
            let Ua1_for_err = 0.0;
            let lm_for = { val: 0.0, err: 0.0 }; // Result;

            // Recurse back to determine U(a1,b), sans normalization.
            //
            let CF1_count = 0;
            let ru   = 0.0;
            let Ua   = 0.0;
            let Uap1 = 0.0;
            let Uam1 = 0.0;

            let rur = { val: 0.0, err: 0.0, count: 0 };
            let rr = { u0: 0.0, u1: 0.0, count: 0 };

            rur = hyperg_U_CF1( (a), (b), 0, x ); //, ru, CF1_count );
            ru = rur.val;
            CF1_count = rur.count;
            Ua   = 1.0;
            Uap1 = ru / (a) * Ua;
            for ( let ap = a; ap >= a1 + 1; ap-- )
            {
                Uam1 = -(((b - 2 * ap) - x) * Ua + (ap) * (1 + ap - b) * Uap1);
                Uap1 = Ua;
                Ua   = Uam1;
                rr = RESCALE_2( Ua, Uap1, scale_factor, scale_count_bck );
                Ua = rr.u0;
                Uap1 = rr.u1;
                scale_count_bck = rr.count;
            }
            Ua1_bck_val = Ua;
            Ua1_bck_err = 2.0 * GSL_DBL_EPSILON * (Math.abs( a1 - a ) + CF1_count + 1) * Math.abs( Ua );

            if ( b == 2 * a1 && a1 > 1 )
            {
                // This can happen when x is small, which is
                // precisely when we need to be careful with
                // this evaluation.
                //
                lm_for = hyperg_lnU_beq2a( (a1), x );
                Ua1_for_val = 1.0;
                Ua1_for_err = 0.0;
            }
            else if ( b == 2 * a1 - 1 && a1 > 1 )
            {
                // Similar to the above. Happens when x is small.
                // Use
                //   U(a,2a-1) = (x U(a,2a) - U(a-1,2(a-1))) / (2a - 2)
                //
                let lnU00 = { val: 0.0, err: 0.0 }; // Result;
                let lnU12 = { val: 0.0, err: 0.0 }; // Result;
                let U00   = { val: 0.0, err: 0.0 }; // Result;
                let U12   = { val: 0.0, err: 0.0 }; // Result;

                lnU00 = hyperg_lnU_beq2a( (a1) - 1.0, x );
                lnU12 = hyperg_lnU_beq2a( (a1),       x );
                if ( lnU00.val > lnU12.val )
                {
                    lm_for.val = lnU00.val;
                    lm_for.err = lnU00.err;
                    U00.val = 1.0;
                    U00.err = 0.0;
                    U12 = gsl_sf_exp_err_e( lnU12.val - lm_for.val, lnU12.err + lm_for.err );
                }
                else
                {
                    lm_for.val = lnU12.val;
                    lm_for.err = lnU12.err;
                    U12.val = 1.0;
                    U12.err = 0.0;
                    U00 = gsl_sf_exp_err_e( lnU00.val - lm_for.val, lnU00.err + lm_for.err );
                }
                Ua1_for_val = (x * U12.val - U00.val) / (2 * a1 - 2);
                Ua1_for_err = (Math.abs( x ) * U12.err + U00.err) / (Math.abs( 2 * a1 - 2 ));
                Ua1_for_err = Ua1_for_err + 2.0 * GSL_DBL_EPSILON * Math.abs( Ua1_for_val );
            }
            else
            {
                // Recurse forward to determine U(a1,b) with
                // absolute normalization.
                //
                let r_Ua = { val: 0.0, err: 0.0 }; // Result;
                let Uam1 = 1.0;  // U(a0-1,b,x) = U(0,b,x)
                let Ua   = 0.0;
                let Uap1 = 0.0;
                let ap   = 0;
                let lm_for_local = 0.0;
                let rr = { u0: 0.0, u1: 0.0, count: 0 };
                let rlm = { val: 0.0, err: 0.0, lnm: 0.0 };

                rlm = hyperg_U_small_a_bgt0( (a0), (b), x ); //, r_Ua, lm_for_local ); // U(1,b,x)
                r_Ua.val = rlm.val;
                r_Ua.err = rlm.err;
                lm_for_local = rlm.lnm;
                Ua = r_Ua.val;
                Uam1 = Uam1 * Math.exp( -lm_for_local );
                lm_for.val = lm_for_local;
                lm_for.err = 0.0;
                
                for ( let ap = a0; ap <= a1 - 1; ap++ )
                {
                    Uap1 = -(Uam1 + ((b - 2 * ap) - x) * Ua) / (ap * (1 + ap - b));
                    Uam1 = Ua;
                    Ua   = Uap1;
                    rr = RESCALE_2( Ua, Uam1, scale_factor, scale_count_for );
                    Ua = rr.u0;
                    Uam1 = rr.u1;
                    scale_count_for = rr.count;
                }
                Ua1_for_val = Ua;
                Ua1_for_err = Math.abs( Ua ) * Math.abs( r_Ua.err / r_Ua.val );
                Ua1_for_err = Ua1_for_err + 2.0 * GSL_DBL_EPSILON * (Math.abs( a1 - a0 ) + 1) * Math.abs( Ua1_for_val );
            }

            // Now do the matching to produce the final result.
            //
            if ( Ua1_bck_val == 0.0 )
            {
                //result.val = 0.0;
                //result.err = 0.0;
                //result.e10 = 0;
                throw "SF.DivideByZeroException";
            }
            else if ( Ua1_for_val == 0.0 )
            {
                // Should never happen.
                throw "SF.UnderflowException";
            }
            else
            {
                let lns = (scale_count_for - scale_count_bck) * Math.log( scale_factor );
                let ln_for_val = Math.log( Math.abs( Ua1_for_val ) );
                let ln_for_err = GSL_DBL_EPSILON + Math.abs( Ua1_for_err / Ua1_for_val );
                let ln_bck_val = Math.log( Math.abs( Ua1_bck_val ) );
                let ln_bck_err = GSL_DBL_EPSILON + Math.abs( Ua1_bck_err / Ua1_bck_val );
                let lnr_val = lm_for.val + ln_for_val - ln_bck_val + lns;
                let lnr_err = lm_for.err + ln_for_err + ln_bck_err
                            + 2.0 * GSL_DBL_EPSILON * (Math.abs( lm_for.val ) + Math.abs( ln_for_val ) + Math.abs( ln_bck_val ) + Math.abs( lns ));
                let sgn = GSL_SIGN( Ua1_for_val ) * GSL_SIGN( Ua1_bck_val );

                r = gsl_sf_exp_err_e10_e( lnr_val, lnr_err );
                r.val = r.val * sgn;
                return r;
            }
        }
    }

} // hyperg_U_int_bge1

// ----------------------------------------------------------------------------

// Handle b >= 1 for generic a,b values.
//
function hyperg_U_bge1( a, b, x )
{
    const rinta = Math.floor( a + 0.5 );
    var a_neg_integer = (a < 0.0 && Math.abs( a - rinta ) < INT_THRESHOLD);

    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( a == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        r.e10 = 0;
        return r;
    }
    else if ( a_neg_integer && Math.abs( rinta ) < INTEGER_LAST )
    {
        // U(-n,b,x) = (-1)^n n! Laguerre[n,b-1,x]
        //
        let n      = -Math.trunc( rinta );
        let sgn    = 0.0;
        let lnfact = { val: 0.0, err: 0.0 }; // Result;
        let L      = { val: 0.0, err: 0.0 }; // Result;

        if ( GSL_IS_ODD( n ) )
        {
            sgn = -1.0;
        }
        else
        {
            sgn = 1.0;
        }
        L = gsl_sf_laguerre_n_e( n, b - 1.0, x );
        lnfact = gsl_sf_lnfact_e( n );
        r = gsl_sf_exp_mult_err_e10_e( lnfact.val, lnfact.err, sgn * L.val, L.err );
        return r;
    }
    else if ( ASYMP_EVAL_OK( Math.trunc( a ), Math.trunc( b ), x ) )
    {
        let ln_pre_val = -a * Math.log( x );
        let ln_pre_err = 2.0 * GSL_DBL_EPSILON * Math.abs( ln_pre_val );
        let asymp      = { val: 0.0, err: 0.0 }; // Result;
        let stat_asymp = 0;

        asymp = hyperg_zaU_asymp( a, b, x );
        r = gsl_sf_exp_mult_err_e10_e( ln_pre_val, ln_pre_err, asymp.val, asymp.err );
        return r;
    }
    else if ( Math.abs( a ) <= 1.0 )
    {
        let rU = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
        let ln_multiplier = 0.0;

        rU = hyperg_U_small_a_bgt0( a, b, x ); //, rU, ln_multiplier );
        ln_multiplier = rU.lnm;
        r = gsl_sf_exp_mult_err_e10_e( ln_multiplier, 2.0 * GSL_DBL_EPSILON * Math.abs(ln_multiplier), rU.val, rU.err );
        return r;
    }
    else if ( SERIES_EVAL_OK( Math.trunc( a ), Math.trunc( b ), x ) )
    {
        let ser = { val: 0.0, err: 0.0 }; // Result;

        ser = hyperg_U_series( a, b, x );
        r.val = ser.val;
        r.err = ser.err;
        r.e10 = 0;
        return r;
    }
    else if ( a < 0.0 )
    {
        // Recurse backward on a and then upward on b.
        //
        const scale_factor = GSL_SQRT_DBL_MAX;
        let scale_count  = 0;
        let a0     = a - Math.floor( a ) - 1.0;
        let b0     = b - Math.floor( b ) + 1.0;
        let lm_0   = 0.0;
        let lm_1   = 0.0;
        let lm_max = 0.0;
        let r_Uap1 = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
        let r_Ua   = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
        let Uap1   = 0.0;
        let Ua     = 0.0;
        let Uam1   = 0.0;
        let ap     = 0.0;
        let rr = { u0: 0.0, u1: 0.0, count: 0 };

        r_Uap1 = hyperg_U_small_a_bgt0( a0 + 1.0, b0, x ); //, r_Uap1, lm_0 );
        r_Ua   = hyperg_U_small_a_bgt0( a0,       b0, x ); //, r_Ua,   lm_1 );
        lm_0 = r_Uap1.lnm;
        lm_1 = r_Ua.lnm;
        Uap1 = r_Uap1.val;
        Ua   = r_Ua.val;
        lm_max = Math.max( lm_0, lm_1 );
        Uap1 = Uap1 * Math.exp( lm_0 - lm_max );
        Ua   = Ua   * Math.exp( lm_1 - lm_max );

        // Downward recursion on a.
        //
        ap = a0;
        while ( ap > a + 0.1 )
        {
            Uam1 = ap * (b0 - ap - 1.0) * Uap1 + (x + 2.0 * ap - b0) * Ua;
            Uap1 = Ua;
            Ua   = Uam1;
            rr = RESCALE_2(Ua, Uap1, scale_factor, scale_count);
            Ua = rr.u0;
            Uap1 = rr.u1;
            scale_count = rr.count;
            ap = ap - 1.0;
        }

        if ( b < 2.0 )
        {
            // b = b0, so no recursion necessary
            //
            let lnscale = Math.log( scale_factor );
            let lnm     = { val: 0.0, err: 0.0 }; // Result;
            let y       = { val: 0.0, err: 0.0 }; // Result;

            lnm.val = lm_max + (scale_count) * lnscale;
            lnm.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( lm_max ) + (scale_count) * Math.abs( lnscale ));
            y.val = Ua;
            y.err = Math.abs( r_Uap1.err / r_Uap1.val ) * Math.abs( Ua );
            y.err = y.err + Math.abs( r_Ua.err / r_Ua.val ) * Math.abs( Ua );
            y.err = y.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( a - a0 ) + 1.0) * Math.abs( Ua );
            y.err = y.err * (Math.abs( lm_0 - lm_max ) + 1.0);
            y.err = y.err * (Math.abs( lm_1 - lm_max ) + 1.0);
            r = gsl_sf_exp_mult_err_e10_e( lnm.val, lnm.err, y.val, y.err );
        }
        else
        {
            // Upward recursion on b.
            //
            let err_mult = Math.abs( b - b0 ) + Math.abs( a - a0 ) + 1.0;
            let lnscale  = Math.log( scale_factor );
            let lnm      = { val: 0.0, err: 0.0 }; // Result;
            let y        = { val: 0.0, err: 0.0 }; // Result;
            let Ubm1     = Ua;                                 // U(a,b0)
            let Ub       = (a * (b0 - a - 1.0) * Uap1 + (a + x) * Ua) / x; // U(a,b0+1)
            let Ubp1     = 0.0;
            let bp       = 0.0;

            bp = b0 + 1.0;
            while ( bp < b - 0.1 )
            {
                Ubp1 = ((1.0 + a - bp) * Ubm1 + (bp + x - 1.0) * Ub) / x;
                Ubm1 = Ub;
                Ub   = Ubp1;
                rr = RESCALE_2(Ub, Ubm1, scale_factor, scale_count);
                Ub = rr.u0;
                Ubm1 = rr.u1;
                scale_count = rr.count;
                bp = bp + 1.0;
            }
            
            lnm.val = lm_max + (scale_count) * lnscale;
            lnm.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( lm_max ) + Math.abs( (scale_count) * lnscale) );
            y.val = Ub;
            y.err = 2.0 * err_mult * Math.abs( r_Uap1.err / r_Uap1.val ) * Math.abs( Ub );
            y.err = y.err + 2.0 * err_mult * Math.abs( r_Ua.err / r_Ua.val ) * Math.abs( Ub );
            y.err = y.err + 2.0 * GSL_DBL_EPSILON * err_mult * Math.abs( Ub );
            y.err = y.err * (Math.abs( lm_0 - lm_max ) + 1.0);
            y.err = y.err * (Math.abs( lm_1 - lm_max ) + 1.0);
            r = gsl_sf_exp_mult_err_e10_e( lnm.val, lnm.err, y.val, y.err );
        }
        return r;
    }
    else if ( b >= 2.0 * a + x )
    {
        // Recurse forward from a near zero.
        // Note that we cannot cross the singularity at
        // the line b=a+1, because the only way we could
        // be in that little wedge is if a < 1. But we
        // have already dealt with the small a case.
        //
        let scale_count  = 0;
        let scale_factor = GSL_SQRT_DBL_MAX;
        let a0      = a - Math.floor( a );
        let lnscale = 0.0;
        let lm_0    = 0.0;
        let lm_1    = 0.0;
        let lm_max  = 0.0;
        let r_Uam1  = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
        let r_Ua    = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
        let lnm     = { val: 0.0, err: 0.0 }; // Result;
        let y       = { val: 0.0, err: 0.0 }; // Result;
        let Uam1    = 0.0;
        let Ua      = 0.0;
        let Uap1    = 0.0;
        let ap      = 0.0;
        let rr = { u0: 0.0, u1: 0.0, count: 0 };

        r_Uam1 = hyperg_U_small_a_bgt0( a0 - 1.0, b, x ); //, r_Uam1, lm_0);
        r_Ua   = hyperg_U_small_a_bgt0( a0,       b, x ); //, r_Ua,   lm_1);
        lm_0 = r_Uam1.lnm;
        lm_1 = r_Ua.lnm;
        Uam1 = r_Uam1.val;
        Ua   = r_Ua.val;
        lm_max = Math.max( lm_0, lm_1 );
        Uam1 = Uam1 * Math.exp( lm_0 - lm_max );
        Ua   = Ua   * Math.exp( lm_1 - lm_max );

        ap = a0;
        while ( ap < a - 0.1 )
        {
            Uap1 = -(Uam1 + (b - 2.0 * ap - x) * Ua) / (ap * (1.0 + ap - b));
            Uam1 = Ua;
            Ua   = Uap1;
            rr = RESCALE_2(Ua, Uam1, scale_factor, scale_count);
            Ua = rr.u0;
            Uam1 = rr.u1;
            scale_count = rr.count;
            ap = ap + 1.0;
        }

        lnscale = Math.log( scale_factor );
        lnm.val = lm_max + (scale_count) * lnscale;
        lnm.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( lm_max ) + Math.abs( (scale_count) * lnscale) );
        y.val = Ua;
        y.err = Math.abs( r_Uam1.err / r_Uam1.val ) * Math.abs( Ua );
        y.err = y.err + Math.abs( r_Ua.err / r_Ua.val ) * Math.abs( Ua );
        y.err = y.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( a - a0 ) + 1.0) * Math.abs( Ua );
        y.err = y.err * (Math.abs( lm_0 - lm_max ) + 1.0);
        y.err = y.err * (Math.abs( lm_1 - lm_max ) + 1.0);
        r = gsl_sf_exp_mult_err_e10_e( lnm.val, lnm.err, y.val, y.err );
        return r;
    }
    else
    {
        if ( b <= x )
        {
            // Recurse backward to a near zero.
            //
            let scale_factor = GSL_SQRT_DBL_MAX;
            let scale_count  = 0;
            let a0        = a - Math.floor( a );
            let lnm       = { val: 0.0, err: 0.0 }; // Result;
            let y         = { val: 0.0, err: 0.0 }; // Result;
            let U0        = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
            let lnscale   = 0.0;
            let lm_0      = 0.0;
            let Uap1      = 0.0;
            let Ua        = 0.0;
            let Uam1      = 0.0;
            let ap        = 0.0;
            let ru        = 0.0;
            let s         = 0.0;
            let CF1_count = 0;
            let rr = { u0: 0.0, u1: 0.0, count: 0 };
            let rur = { val: 0.0, err: 0.0, count: 0 };

            rur = hyperg_U_CF1(a, b, 0, x ); //, ru, CF1_count); // ***
            ru = rur.val;
            CF1_count = rur.count;
            s = ru / a;
            Ua   = GSL_SQRT_DBL_MIN;
            Uap1 = s * Ua;
            ap = a;
            while ( ap > a0 + 0.1 )
            {
                Uam1 = -((b - 2.0 * ap - x) * Ua + ap * (1.0 + ap - b) * Uap1);
                Uap1 = Ua;
                Ua   = Uam1;
                rr = RESCALE_2( Ua, Uap1, scale_factor, scale_count );
                Ua = rr.u0;
                Uap1 = rr.u1;
                scale_count = rr.count;
                ap = ap - 1.0;
            }
        
            U0 = hyperg_U_small_a_bgt0( a0, b, x ); //, U0, lm_0 );
            lm_0 = U0.lnm;
        
            lnscale = Math.log( scale_factor );
            lnm.val = lm_0 - (scale_count) * lnscale;
            lnm.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( lm_0 ) + Math.abs( (scale_count) * lnscale) );
            y.val = GSL_SQRT_DBL_MIN * (U0.val / Ua);
            y.err = GSL_SQRT_DBL_MIN * (U0.err / Math.abs( Ua ));
            y.err = y.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( a0 - a ) + (CF1_count) + 1.0) * Math.abs( y.val );
            r = gsl_sf_exp_mult_err_e10_e( lnm.val, lnm.err, y.val, y.err );
            return r;
        }
        else
        {
            // Recurse backward to near the b=2a+x line, then
            // forward from a near zero to get the normalization.
            //
            let scale_count_for = 0;
            let scale_count_bck = 0;
            let scale_factor    = GSL_SQRT_DBL_MAX;
            let eps       = a - Math.floor( a );
            let a0        = 0.0;
            let a1        = 0.0;
            let lnm       = { val: 0.0, err: 0.0 }; // Result;
            let y         = { val: 0.0, err: 0.0 }; // Result;
            let lm_for    = 0.0;
            let lnscale   = 0.0;
            let Ua1_bck   = 0.0;
            let Ua1_for   = 0.0;
            let CF1_count = 0;

            if ( eps == 0.0 )
            {
                a0 = 1.0;
            }
            else
            {
                a0 = eps;
            }
            a1 = a0 + Math.ceil( 0.5 * (b - x) - a0 );

            // Recurse back to determine U(a1,b), sans normalization.
            //
            let Uap1 = 0.0;
            let Ua   = 0.0;
            let Uam1 = 0.0;
            let ap   = 0.0;
            let ru   = 0.0;
            let r    = 0.0;

            let rr = { u0: 0.0, u1: 0.0, count: 0 };
            let rur = { val: 0.0, err: 0.0, count: 0 };

            rur = hyperg_U_CF1( a, b, 0, x ); //, ru, CF1_count ); // ***
            ru = rur.val;
            CF1_count = rur.count;
            r = ru / a;
            Ua   = GSL_SQRT_DBL_MIN;
            Uap1 = r * Ua;
            ap = a;
            while ( ap > a1 + 0.1 )
            {
                Uam1 = -((b - 2.0 * ap - x) * Ua + ap * (1.0 + ap - b) * Uap1);
                Uap1 = Ua;
                Ua   = Uam1;
                rr = RESCALE_2( Ua, Uap1, scale_factor, scale_count_bck );
                Ua = rr.u0;
                Uap1 = rr.u1;
                scale_count_bck = rr.count;
                ap = ap - 1.0;
            }
            Ua1_bck = Ua;

            // Recurse forward to determine U(a1,b) with
            // absolute normalization.
            //
            let r_Uam1 = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
            let r_Ua   = { val: 0.0, err: 0.0, lnm: 0.0 }; // Result;
            let lm_0   = 0.0;
            let lm_1   = 0.0;

            r_Uam1 = hyperg_U_small_a_bgt0( a0 - 1.0, b, x ); //, r_Uam1, lm_0);
            r_Ua   = hyperg_U_small_a_bgt0( a0,       b, x ); //, r_Ua,   lm_1);
            lm_0  = r_Uam1.lnm;
            lm_1  = r_Ua.lnm;
            Uam1 = r_Uam1.val;
            Ua   = r_Ua.val;
            
            lm_for = Math.max( lm_0, lm_1 );
            Uam1 = Uam1 * Math.exp( lm_0 - lm_for );
            Ua   = Ua   * Math.exp( lm_1 - lm_for );
            ap = a0;
            while ( ap < a1 - 0.1 )
            {
                Uap1 = -(Uam1 + (b - 2.0 * ap - x) * Ua) / (ap * (1.0 + ap - b));
                Uam1 = Ua;
                Ua   = Uap1;
                rr = RESCALE_2(Ua, Uam1, scale_factor, scale_count_for);
                Ua = rr.u0;
                Uam1 = rr.u1;
                scale_count_for = rr.count;
                ap = ap + 1.0;
            }
            Ua1_for = Ua;

            lnscale = Math.log( scale_factor );
            lnm.val = lm_for + (scale_count_for - scale_count_bck) * lnscale;
            lnm.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( lm_for ) + Math.abs( (scale_count_for - scale_count_bck) ) * Math.abs( lnscale ));
            y.val = GSL_SQRT_DBL_MIN * Ua1_for / Ua1_bck;
            y.err = 2.0 * GSL_DBL_EPSILON * (Math.abs( a - a0 ) + (CF1_count) + 1.0) * Math.abs( y.val );
            r = gsl_sf_exp_mult_err_e10_e( lnm.val, lnm.err, y.val, y.err );
            return r;
        }
    }

} // hyperg_U_bge1

// --*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hyperg_U_int_e10_e( a, b, x )
{

    if ( x <= 0.0 )
    {
        throw "SF.DomainException";
    }
    else
    {
        if ( b >= 1 )
        {
            return hyperg_U_int_bge1( a, b, x );
        }
        else
        {
            // Use the reflection formula
            // U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
            //
            var U = { val: 0.0, err: 0.0 }; // Result;
            var r = { val: 0.0, err: 0.0 }; // Result;
            var ln_x = Math.log(x);
            var ap = 1 + a - b;
            var bp = 2 - b;
            var ln_pre_val = 0.0;
            var ln_pre_err = 0.0;

            U = hyperg_U_int_bge1( ap, bp, x );
            ln_pre_val = (1.0 - (b)) * ln_x;
            ln_pre_err = 2.0 * GSL_DBL_EPSILON * (Math.abs( b ) + 1) * Math.abs( ln_x );
            ln_pre_err = ln_pre_err + 2.0 * GSL_DBL_EPSILON * Math.abs( 1.0 - (b) ); // error in log(x)
            r = gsl_sf_exp_mult_err_e10_e( ln_pre_val + (U.e10) * M_LN10, ln_pre_err, U.val, U.err );
            return r;
        }
    }

} // gsl_sf_hyperg_U_int_e10_e

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_U_e10_e( a, b, x )
{
    var rinta = Math.floor( a + 0.5 );
    var rintb = Math.floor( b + 0.5 );
    var a_integer = (Math.abs( a - rinta ) < INT_THRESHOLD);
    var b_integer = (Math.abs( b - rintb ) < INT_THRESHOLD);

    var r = { val: 0.0, err: 0.0, e10: 0 }; // Result;

    if ( x <= 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( a == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        r.e10 = 0;
        return r;
    }
    else if ( a_integer && b_integer )
    {
        return gsl_sf_hyperg_U_int_e10_e( Math.trunc( rinta ), Math.trunc( rintb ), x );
    }
    else
    {
        if ( b >= 1.0 )
        {
            // Use b >= 1 function.
            //
            return hyperg_U_bge1( a, b, x );
        }
        else
        {
            // Use the reflection formula
            // U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)
            //
            var lnx = Math.log(x);
            var ln_pre_val = (1.0 - b) * lnx;
            var ln_pre_err = Math.abs(lnx) * 2.0 * GSL_DBL_EPSILON * (1.0 + Math.abs(b));
            var ap = 1.0 + a - b;
            var bp = 2.0 - b;
            var U = { val: 0.0, err: 0.0 }; // Result;

            U = hyperg_U_bge1( ap, bp, x );
            r = gsl_sf_exp_mult_err_e10_e( ln_pre_val + (U.e10) * M_LN10, ln_pre_err, U.val, U.err);
            return r;
        }
    }

} // gsl_sf_hyperg_U_e10_e

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_U_int_e( a, b, x )
{
    var re = { val: 0.0, err: 0.0, e10: 0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    re = gsl_sf_hyperg_U_int_e10_e( a, b, x );
    r  = gsl_sf_result_smash_e( re );
    return r;

} // gsl_sf_hyperg_U_int_e

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_U_e( a, b, x )
{
    var re = { val: 0.0, err: 0.0, e10: 0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    re = gsl_sf_hyperg_U_e10_e(a, b, x);
    r  = gsl_sf_result_smash_e(re);
    return r;

} // gsl_sf_hyperg_U_e


// --*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

// function gsl_sf_hyperg_U_int(a: INTEGER; b: INTEGER; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_hyperg_U_int
//     return EVAL_RESULT(gsl_sf_hyperg_U_int_e'Access, (a, b, x), "gsl_sf_hyperg_U_int");
// END gsl_sf_hyperg_U_int;

// function gsl_sf_hyperg_U(a: LONG_FLOAT; b: LONG_FLOAT; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_hyperg_U
//     return EVAL_RESULT(gsl_sf_hyperg_U_e'Access, (a, b, x), "gsl_sf_hyperg_U");
// END gsl_sf_hyperg_U;

// ----------------------------------------------------------------------------
// EOF SF-HypergeometricU.mjs
