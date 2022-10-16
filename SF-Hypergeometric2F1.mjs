// SF-Hypergeometric2F1.mjs
// These routines compute the confluent hypergeometric function 2F1
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

import { M_EULER }                from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }        from "./SF-Machine.mjs";
import { GSL_ROOT5_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }        from "./SF-Machine.mjs";
import { GSL_IS_ODD }             from "./SF-Math.mjs";
import { gsl_sf_lngamma_e }       from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_sgn_e }   from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_complex_e }   from "./SF-Gamma.mjs";
import { gsl_sf_exp_mult_err_e }  from "./SF-Exponential.mjs";
import { gsl_sf_hyperg_1F1_e }    from "./SF-Hypergeometric1F1.mjs";
import { gsl_sf_exp_err_e }       from "./SF-Exponential.mjs";
import { gsl_sf_psi_e }           from "./SF-Psi.mjs";

import { EVAL_RESULT_4D }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const locEPS = 1000.0 * GSL_DBL_EPSILON;


// Assumes c != negative integer.
//
function hyperg_2F1_series( a, b, c, x )
{
    var sum_pos = 1.0;
    var sum_neg = 0.0;
    var del_pos = 1.0;
    var del_neg = 0.0;
    var del     = 1.0;
    var k       = 0.0;
    var i       = 0;
    var r       = { val: 0.0, err: 0.0 }; // Result;

    if ( Math.abs( c ) < GSL_DBL_EPSILON )
    {
        //r.val = 0.0; // FIXME: ??
        //r.err = 1.0;
        throw "SF.DomainException";
    }

    i = 0;
    while ( true )
    {
        i = i + 1;
        if ( i > 30000 )
        {
            //r.val = sum_pos - sum_neg;
            //r.err = del_pos + del_neg;
            //r.err = r.err + 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
            //r.err = r.err + 2.0 * GSL_DBL_EPSILON * (2.0 * Sqrt(k) + 1.0) * Math.abs(r.val);
            throw "SF.MaxIterationsException";
        }
        del = del * ((a + k) * (b + k) * x / ((c + k) * (k + 1.0))); // Gauss series

        if ( del > 0.0 )
        {
            del_pos = del;
            sum_pos = sum_pos + del;
        }
        else if ( del == 0.0 )
        {
            // Exact termination (a or b was a negative integer).
            del_pos = 0.0;
            del_neg = 0.0;
            break;
        }
        else
        {
            del_neg = -del;
            sum_neg = sum_neg - del;
        }

        k = k + 1.0;
        if ( ! (Math.abs((del_pos + del_neg) / (sum_pos - sum_neg)) > GSL_DBL_EPSILON) ) break;
    }

    r.val = sum_pos - sum_neg;
    r.err = del_pos + del_neg;
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * (2.0 * Math.sqrt( k ) + 1.0) * Math.abs( r.val );

    return r;

} // hyperg_2F1_series

// ----------------------------------------------------------------------------

// a = aR + i aI, b = aR - i aI
function hyperg_2F1_conj_series( aR, aI, c, x )
{

    if ( c == 0.0 )
    {
        //r.val = 0.0; // FIXME: should be Inf
        //r.err = 0.0;
        throw "SF.DomainException";
    }
    else
    {
        var sum_pos = 0.0;
        var sum_neg = 0.0;
        var del_pos = 0.0;
        var del_neg = 0.0;
        var del     = 0.0;
        var k       = 0.0;
        var r       = { val: 0.0, err: 0.0 }; // Result;

        sum_pos = 1.0;
        sum_neg = 0.0;
        del_pos = 1.0;
        del_neg = 0.0;
        del = 1.0;
        k = 0.0;
        while ( true )
        {
            del = del * (((aR + k) * (aR + k) + aI * aI) / ((k + 1.0) * (c + k)) * x);

            if ( del >= 0.0 )
            {
                del_pos = del;
                sum_pos = sum_pos + del;
            }
            else
            {
                del_neg = -del;
                sum_neg = sum_neg - del;
            }

            if ( k > 30000 )
            {
                //r.val = sum_pos - sum_neg;
                //r.err = del_pos + del_neg;
                //r.err = r.err + 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
                //r.err = r.err + 2.0 * GSL_DBL_EPSILON * (2.0 * Sqrt(k) + 1.0) * Math.abs(r.val);
                throw "SF.MaxIterationsException";
            }

            k = k + 1.0;
            if ( ! (Math.abs( (del_pos + del_neg) / (sum_pos - sum_neg) ) > GSL_DBL_EPSILON) ) break;
        }
        
        r.val = sum_pos - sum_neg;
        r.err = del_pos + del_neg;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (2.0 * Math.sqrt( k ) + 1.0) * Math.abs( r.val );
        
        return r;
    }

} // hyperg_2F1_conj_series

// ----------------------------------------------------------------------------

// Luke's rational approximation. The most accesible
// discussion is in [Kolbig, CPC 23, 51 (1981)].
// The convergence is supposedly guaranteed for x < 0.
// You have to read Luke's books to see this and other
// results. Unfortunately, the stability is not so
// clear to me, although it seems very efficient when
// it works.
//
function hyperg_2F1_luke( a, b, c, xin )
{
    const RECUR_BIG = 1.0e+50;
    const nmax = 20000;
    var n    = 0;
    var x    = 0.0;
    var x3   = 0.0;
    var t0   = 0.0;
    var t1   = 0.0;
    var t2   = 0.0;
    var F    = 0.0;
    var prec = 0.0;
    var Bnm3 = 0.0;
    var Bnm2 = 0.0;
    var Bnm1 = 0.0;
    var Anm3 = 0.0;
    var Anm2 = 0.0;
    var Anm1 = 0.0;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    n  = 3;
    x  = -xin;
    x3 = x * x * x;
    t0 = a * b / c;
    t1 = (a + 1.0) * (b + 1.0) / (2.0 * c);
    t2 = (a + 2.0) * (b + 2.0) / (2.0 * (c + 1.0));
    F  = 1.0;
  
    Bnm3 = 1.0;                                  // B0
    Bnm2 = 1.0 + t1 * x;                         // B1
    Bnm1 = 1.0 + t2 * x * (1.0 + t1 / 3.0 * x);  // B2
   
    Anm3 = 1.0;                                                                // A0
    Anm2 = Bnm2 - t0 * x;                                                      // A1
    Anm1 = Bnm1 - t0 * (1.0 + t2 * x) * x + t0 * t1 * (c / (c + 1.0)) * x * x; // A2

    while ( true )
    {
        var npam1 = 0.0;
        var npbm1 = 0.0;
        var npcm1 = 0.0;
        var npam2 = 0.0;
        var npbm2 = 0.0;
        var npcm2 = 0.0;
        var tnm1  = 0.0;
        var tnm3  = 0.0;
        var tnm5  = 0.0;
        var n2    = 0.0;
        var F1    = 0.0;
        var F2    = 0.0;
        var F3    = 0.0;
        var E     = 0.0;
        var An    = 0.0;
        var Bn    = 0.0;
        var s     = 0.0;

        npam1 = (n) + a - 1.0;
        npbm1 = (n) + b - 1.0;
        npcm1 = (n) + c - 1.0;
        npam2 = (n) + a - 2.0;
        npbm2 = (n) + b - 2.0;
        npcm2 = (n) + c - 2.0;
        tnm1  = (2 * n - 1);
        tnm3  = (2 * n - 3);
        tnm5  = (2 * n - 5);
        n2 = (n) * (n);
        F1 =  (3.0 * n2 + (a + b - 6.0) * (n) + 2.0 - a * b - 2.0 * (a + b)) / (2.0 * tnm3 * npcm1);
        F2 = -(3.0 * n2 - (a + b + 6.0) * (n) + 2.0 - a * b) * npam1 * npbm1 / (4.0 * tnm1 * tnm3 * npcm2 * npcm1);
        F3 = (npam2 * npam1 * npbm2 * npbm1 * ((n) - a - 2.0) * ((n) - b - 2.0)) / (8.0 * tnm3 * tnm3 * tnm5 * ((n) + c - 3.0) * npcm2 * npcm1);
        E  = -npam1 * npbm1 * ((n) - c - 1.0) / (2.0 * tnm3 * npcm2 * npcm1);
    
        An = (1.0 + F1 * x) * Anm1 + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
        Bn = (1.0 + F1 * x) * Bnm1 + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
        s  = An / Bn;
    
        prec = Math.abs( (F - s) / F );
        F    = s;

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
    r.err = 2.0 * Math.abs( prec * F );
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * (n + 1) * Math.abs( F );

    // FIXME: just a hack: there's a lot of shit going on here
    r.err = r.err * 8.0 * (Math.abs( a ) + Math.abs( b ) + 1.0);

    //stat_iter = (n >= nmax ? GSL_EMAXITER : GSL_SUCCESS );
    if ( n >= nmax )
    {
        throw "SF.MaxIterationsException";
    }

    return r;

} // hyperg_2F1_luke

// ----------------------------------------------------------------------------

// Luke's rational approximation for the
// case a = aR + i aI, b = aR - i aI.
//
function hyperg_2F1_conj_luke( aR, aI, c, xin )
{
    const RECUR_BIG = 1.0e+50;
    const nmax      = 10000;
    var n         = 0;
    var x         = 0.0;
    var x3        = 0.0;
    var atimesb   = 0.0;
    var apb       = 0.0;
    var t0        = 0.0;
    var t1        = 0.0;
    var t2        = 0.0;
    var F         = 0.0;
    var prec      = 0.0;
    var Bnm3      = 0.0;
    var Bnm2      = 0.0;
    var Bnm1      = 0.0;
    var Anm3      = 0.0;
    var Anm2      = 0.0;
    var Anm1      = 0.0;
    var r         = { val: 0.0, err: 0.0 }; // Result;

    x       = -xin;
    x3      = x * x * x;
    atimesb = aR * aR + aI * aI;
    apb     = 2.0 * aR;
    t0      = atimesb / c;
    t1      = (atimesb +       apb + 1.0) / (2.0 * c);
    t2      = (atimesb + 2.0 * apb + 4.0) / (2.0 * (c + 1.0));
    F       = 1.0;
    
    Bnm3 = 1.0;                                  // B0
    Bnm2 = 1.0 + t1 * x;                         // B1
    Bnm1 = 1.0 + t2 * x * (1.0 + t1 / 3.0 * x);  // B2
    
    Anm3 = 1.0;                                                                // A0
    Anm2 = Bnm2 - t0 * x;                                                      // A1
    Anm1 = Bnm1 - t0 * (1.0 + t2 * x) * x + t0 * t1 * (c / (c + 1.0)) * x * x; // A2

    n = 3;
    while ( true )
    {
        let nm1         = 0.0;
        let nm2         = 0.0;
        let npam1_npbm1 = 0.0;
        let npam2_npbm2 = 0.0;
        let npcm1       = 0.0;
        let npcm2       = 0.0;
        let tnm1        = 0.0;
        let tnm3        = 0.0;
        let tnm5        = 0.0;
        let n2          = 0.0;
        let F1          = 0.0;
        let F2          = 0.0;
        let F3          = 0.0;
        let E           = 0.0;
        let An          = 0.0;
        let Bn          = 0.0;
        let r1          = 0.0;

        nm1 = (n - 1);
        nm2 = (n - 2);
        npam1_npbm1 = atimesb + nm1 * apb + nm1 * nm1;
        npam2_npbm2 = atimesb + nm2 * apb + nm2 * nm2;
        npcm1 = nm1 + c;
        npcm2 = nm2 + c;
        tnm1  = (2 * n - 1);
        tnm3  = (2 * n - 3);
        tnm5  = (2 * n - 5);
        n2 = (n) * (n);
        F1 =  (3.0 * n2 + (apb - 6.0) * n + 2 - atimesb - 2.0 * apb) / (2.0 * tnm3 * npcm1);
        F2 = -(3.0 * n2 - (apb + 6.0) * n + 2 - atimesb) * npam1_npbm1 / (4.0 * tnm1 * tnm3 * npcm2 * npcm1);
        F3 = (npam2_npbm2 * npam1_npbm1 * (nm2 * nm2 - nm2 * apb + atimesb)) / (8.0 * tnm3 * tnm3 * tnm5 * ((n) + c - 3.0) * npcm2 * npcm1);
        E  = -npam1_npbm1 * ((n) - c - 1.0) / (2.0 * tnm3 * npcm2 * npcm1);
        
        An = (1.0 + F1 * x) * Anm1 + (E + F2 * x) * x * Anm2 + F3 * x3 * Anm3;
        Bn = (1.0 + F1 * x) * Bnm1 + (E + F2 * x) * x * Bnm2 + F3 * x3 * Bnm3;
        r1 = An / Bn;
        
        prec = Math.abs( F - r1 ) / Math.abs( F );
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
    r.err = 2.0 * Math.abs( prec * F );
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * (n + 1) * Math.abs( F );
    
    // FIXME: see above
    r.err = r.err * 8.0 * (Math.abs( aR ) + Math.abs( aI ) + 1.0);

    if ( n >= nmax )
    {
        throw "SF.MaxIterationsException";
    }
    
    return r;

} // hyperg_2F1_conj_luke

// ----------------------------------------------------------------------------

// Do the reflection described in [Moshier, p. 334].
// Assumes a,b,c != neg integer.
//
function hyperg_2F1_reflect( a, b, c, x )
{
    var d         = 0.0;
    var intd      = 0;
    var d_integer = false;
    var r         = { val: 0.0, err: 0.0 }; // Result;

    d = c - a - b;
    intd  = Math.trunc( Math.floor( d + 0.5 ) );
    d_integer = (Math.abs( d - intd ) < locEPS);

    if ( d_integer )
    {
        var ln_omx  = 0.0;
        var ad      = 0.0;
        var sgn_2   = 0.0;
        var d1      = 0.0;
        var d2      = 0.0;
        var F1      = { val: 0.0, err: 0.0 }; // Result;
        var F2      = { val: 0.0, err: 0.0 }; // Result;
        var lng_c   = { val: 0.0, err: 0.0 }; // Result;
        var lng_ad2 = { val: 0.0, err: 0.0 }; // Result;
        var lng_bd2 = { val: 0.0, err: 0.0 }; // Result;
        var ok      = false;

        ln_omx = Math.log( 1.0 - x );
        ad = Math.abs( d );

        if ( d >= 0.0 )
        {
            d1 = d;
            d2 = 0.0;
        }
        else
        {
            d1 = 0.0;
            d2 = d;
        }

        ok = true;
        try
        {
            lng_ad2 = gsl_sf_lngamma_e( a + d2 );
            lng_bd2 = gsl_sf_lngamma_e( b + d2 );
            lng_c   = gsl_sf_lngamma_e( c );
        }
        catch ( e )
        {
            ok = false;
        }
    
        // Evaluate F1.
        //
        if ( ad < GSL_DBL_EPSILON )
        {
            // d = 0
            F1.val = 0.0;
            F1.err = 0.0;
        }
        else
        {
            var lng_ad  = { val: 0.0, err: 0.0 }; // Result;
            var lng_ad1 = { val: 0.0, err: 0.0 }; // Result;
            var lng_bd1 = { val: 0.0, err: 0.0 }; // Result;
            var sum1        = 0.0;
            var term        = 0.0;
            var ln_pre1_val = 0.0;
            var ln_pre1_err = 0.0;
            var j = 0;
            try
            {
                lng_ad  = gsl_sf_lngamma_e( ad );
                lng_ad1 = gsl_sf_lngamma_e( a + d1 );
                lng_bd1 = gsl_sf_lngamma_e( b + d1 );

                // Gamma functions in the denominator are ok.
                // Proceed with evaluation.
                //
                sum1 = 1.0;
                term = 1.0;
                ln_pre1_val = lng_ad.val + lng_c.val + d2 * ln_omx - lng_ad1.val - lng_bd1.val;
                ln_pre1_err = lng_ad.err + lng_c.err + lng_ad1.err + lng_bd1.err + GSL_DBL_EPSILON * Math.abs( ln_pre1_val );

                // Do F1 sum.
                //
                for (let i = 1; i <= Math.trunc(ad) - 1; i++ )
                {
                    j    = i - 1;
                    term = term * (a + d2 + (j)) * (b + d2 + (j)) / (1.0 + d2 + (j)) / (i) * (1.0 - x);
                    sum1 = sum1 + term;
                }
    
                F1 = gsl_sf_exp_mult_err_e( ln_pre1_val, ln_pre1_err, sum1, GSL_DBL_EPSILON * Math.abs( sum1 ) );
            }
            catch ( e )
            { 
                // Gamma functions in the denominator were not ok.
                // So the F1 term is zero.
                //
                F1.val = 0.0;
                F1.err = 0.0;
            }
        } // end F1 evaluation

        // Evaluate F2.
        //
        if ( ok )
        {
            // Gamma functions in the denominator are ok.
            // Proceed with evaluation.
            //
            const maxiter     = 2000;
            let psi_1       = -M_EULER;
            let psi_1pd     = { val: 0.0, err: 0.0 }; // Result;
            let psi_apd1    = { val: 0.0, err: 0.0 }; // Result;
            let psi_bpd1    = { val: 0.0, err: 0.0 }; // Result;
            let psi_val     = 0.0;
            let psi_err     = 0.0;
            let fact        = 0.0;
            let sum2_val    = 0.0;
            let sum2_err    = 0.0;
            let ln_pre2_val = 0.0;
            let ln_pre2_err = 0.0;
            let delta1      = 0.0;
            let term1       = 0.0;
            let term2       = 0.0;
            let j           = 0;

            psi_1pd  = gsl_sf_psi_e( 1.0 + ad );
            psi_apd1 = gsl_sf_psi_e( a + d1 );
            psi_bpd1 = gsl_sf_psi_e( b + d1 );

            psi_val = psi_1 + psi_1pd.val - psi_apd1.val - psi_bpd1.val - ln_omx;
            psi_err = psi_1pd.err + psi_apd1.err + psi_bpd1.err + GSL_DBL_EPSILON * Math.abs( psi_val );
            fact = 1.0;
            sum2_val = psi_val;
            sum2_err = psi_err;
            ln_pre2_val = lng_c.val + d1*ln_omx - lng_ad2.val - lng_bd2.val;
            ln_pre2_err = lng_c.err + lng_ad2.err + lng_bd2.err + GSL_DBL_EPSILON * Math.abs( ln_pre2_val );

            // Do F2 sum.
            //
            j = 1;
            while ( j < maxiter )
            {
                // values for psi functions use recurrence; Abramowitz+Stegun 6.3.5
                term1    = 1.0 / (j)  + 1.0 / (ad + (j));
                term2    = 1.0 / (a + d1 + (j) - 1.0) + 1.0 / (b + d1 + (j) - 1.0);
                delta1   = 0.0;
                psi_val  = psi_val + term1 - term2;
                psi_err  = psi_err + GSL_DBL_EPSILON * (Math.abs( term1 ) + Math.abs( term2 ));
                fact     = fact * ((a + d1 + (j) - 1.0) * (b + d1 + (j) - 1.0) / ((ad + (j)) * (j)) * (1.0 - x));
                delta1   = fact * psi_val;
                sum2_val = sum2_val + delta1;
                sum2_err = sum2_err + Math.abs( fact * psi_err ) + GSL_DBL_EPSILON * Math.abs( delta1 );
                if ( Math.abs( delta1 ) < GSL_DBL_EPSILON * Math.abs(sum2_val) ) break;
                j = j + 1;
            }

            if ( j >= maxiter )
            {
                throw "SF.MaxIterationsException";
            }

            if ( sum2_val == 0.0 )
            {
                F2.val = 0.0;
                F2.err = 0.0;
            }
            else
            {
                F2 = gsl_sf_exp_mult_err_e( ln_pre2_val, ln_pre2_err, sum2_val, sum2_err );
            }
        }
        else
        {
            // Gamma functions in the denominator not ok.
            // So the F2 term is zero.
            //
            F2.val = 0.0;
            F2.err = 0.0;
        } // end F2 evaluation

        if ( GSL_IS_ODD( intd ) )
        {
            sgn_2 = -1.0;
        }
        else
        {
            sgn_2 = 1.0;
        }
        r.val = F1.val + sgn_2 * F2.val;
        r.err = F1.err + F2. err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( F1.val ) + Math.abs( F2.val ));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    else
    {
        // d not an integer
        let sgn1     = 0.0;
        let sgn2     = 0.0;
        let sgn_g1ca = 0.0;
        let sgn_g1cb = 0.0;
        let sgn_g2a  = 0.0;
        let sgn_g2b  = 0.0;
        let sgn_gc   = 0.0;
        let sgn_gd   = 0.0;
        let sgn_gmd  = 0.0;
        let ln_pre1_val = 0.0;
        let ln_pre2_val = 0.0;
        let ln_pre1_err = 0.0;
        let ln_pre2_err = 0.0;
        let pre1     = { val: 0.0, err: 0.0 }; // Result;
        let pre2     = { val: 0.0, err: 0.0 }; // Result;
        let F1       = { val: 0.0, err: 0.0 }; // Result;
        let F2       = { val: 0.0, err: 0.0 }; // Result;
        let ln_g1ca  = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let ln_g1cb  = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let ln_g2a   = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let ln_g2b   = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let ln_gc    = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let ln_gd    = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let ln_gmd   = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        let ok1      = false;
        let ok2      = false;

        // These gamma functions appear in the denominator, so we
        // catch their harmless domain errors and set the terms to zero.
        //
        ok1 = true;
        try
        {
            ln_g1ca = gsl_sf_lngamma_sgn_e( c - a ); //, ln_g1ca, sgn_g1ca );
            ln_g1cb = gsl_sf_lngamma_sgn_e( c - b ); //, ln_g1cb, sgn_g1cb );
            sgn_g1ca = ln_g1ca.sign;
            sgn_g1cb = ln_g1cb.sign;
        }
        catch ( e )
        {
            ok1 = false;
        }
        ok2 = true;
        try
        {
            ln_g2a = gsl_sf_lngamma_sgn_e( a ); //, ln_g2a, sgn_g2a );
            ln_g2b = gsl_sf_lngamma_sgn_e( b ); //, ln_g2b, sgn_g2b );
            sgn_g2a = ln_g2a.sign;
            sgn_g2b = ln_g2b.sign;
        }
        catch ( e )
        {
            ok2 = false;
        }

        ln_gc  = gsl_sf_lngamma_sgn_e(  c ); //, ln_gc,  sgn_gc );
        ln_gd  = gsl_sf_lngamma_sgn_e(  d ); //, ln_gd,  sgn_gd );
        ln_gmd = gsl_sf_lngamma_sgn_e( -d ); //, ln_gmd, sgn_gmd );
        sgn_gc = ln_gc.sign;
        sgn_gd = ln_gd.sign;
        sgn_gmd = ln_gmd.sign;

        sgn1 = sgn_gc * sgn_gd  * sgn_g1ca * sgn_g1cb;
        sgn2 = sgn_gc * sgn_gmd * sgn_g2a  * sgn_g2b;

        if ( ok1 && ok2 )
        {
            ln_pre1_val = ln_gc.val + ln_gd.val  - ln_g1ca.val - ln_g1cb.val;
            ln_pre2_val = ln_gc.val + ln_gmd.val - ln_g2a.val  - ln_g2b.val + d * Math.log( 1.0 - x );
            ln_pre1_err = ln_gc.err + ln_gd.err + ln_g1ca.err + ln_g1cb.err;
            ln_pre2_err = ln_gc.err + ln_gmd.err + ln_g2a.err  + ln_g2b.err;
            if ( ln_pre1_val < GSL_LOG_DBL_MAX && ln_pre2_val < GSL_LOG_DBL_MAX)
            {
                pre1 = gsl_sf_exp_err_e( ln_pre1_val, ln_pre1_err );
                pre2 = gsl_sf_exp_err_e( ln_pre2_val, ln_pre2_err );
                pre1.val = pre1.val * sgn1;
                pre2.val = pre2.val * sgn2;
            }
            else
            {
                throw "SF.OverflowException";
            }
        }
        else if ( ok1 && ! ok2 )
        {
            ln_pre1_val = ln_gc.val + ln_gd.val - ln_g1ca.val - ln_g1cb.val;
            ln_pre1_err = ln_gc.err + ln_gd.err + ln_g1ca.err + ln_g1cb.err;
            if ( ln_pre1_val < GSL_LOG_DBL_MAX )
            {
                pre1 = gsl_sf_exp_err_e( ln_pre1_val, ln_pre1_err );
                pre1.val = pre1.val * sgn1;
                pre2.val = 0.0;
                pre2.err = 0.0;
            }
            else
            {
                throw "SF.OverflowException";
            }
        }
        else if ( ! ok1 && ok2 )
        {
            ln_pre2_val = ln_gc.val + ln_gmd.val - ln_g2a.val - ln_g2b.val + d * Math.log( 1.0 - x );
            ln_pre2_err = ln_gc.err + ln_gmd.err + ln_g2a.err + ln_g2b.err;
            if ( ln_pre2_val < GSL_LOG_DBL_MAX )
            {
                pre1.val = 0.0;
                pre1.err = 0.0;
                pre2 = gsl_sf_exp_err_e( ln_pre2_val, ln_pre2_err );
                pre2.val = pre2.val * sgn2;
            }
            else
            {
                throw "SF.OverflowException";
            }
        }
        else
        {
            //pre1.val = 0.0;
            //pre2.val = 0.0;
            throw "SF.UnderflowException";
        }

        F1 = hyperg_2F1_series(     a,     b, 1.0 - d, 1.0 - x );
        F2 = hyperg_2F1_series( c - a, c - b, 1.0 + d, 1.0 - x );

        r.val = pre1.val * F1.val + pre2.val * F2.val;
        r.err = Math.abs( pre1.val * F1.err ) + Math.abs( pre2.val * F2.err );
        r.err = r.err + Math.abs( pre1.err * F1.val ) + Math.abs( pre2.err * F2.val );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs( pre1.val * F1.val ) + Math.abs( pre2.val * F2.val ));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs( r.val );

        return r;
    }

} // hyperg_2F1_reflect

// ----------------------------------------------------------------------------

function pow_omx( x, p )
{
    var ln_omx    = 0.0;
    var ln_result = 0.0;

    if ( Math.abs( x ) < GSL_ROOT5_DBL_EPSILON )
    {
        ln_omx = -x * (1.0 + x * (1.0 / 2.0 + x * (1.0 / 3.0 + x / 4.0 + x * x / 5.0)));
    }
    else
    {
        ln_omx = Math.log( 1.0 - x );
    }
    ln_result = p * ln_omx;
    return gsl_sf_exp_err_e( ln_result, GSL_DBL_EPSILON * Math.abs( ln_result ) );

} // pow_omx

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

export function gsl_sf_hyperg_2F1_e( a, b, c, x )
{
    var d     = c - a - b;
    var rinta = Math.floor( a + 0.5 );
    var rintb = Math.floor( b + 0.5 );
    var rintc = Math.floor( c + 0.5 );
    var a_neg_integer = ( a < 0.0 && Math.abs( a - rinta ) < locEPS );
    var b_neg_integer = ( b < 0.0 && Math.abs( b - rintb ) < locEPS );
    var c_neg_integer = ( c < 0.0 && Math.abs( c - rintc ) < locEPS );
    var r = { val: 0.0, err: 0.0 }; // Result;
    var ap = 0.0;
    var bp = 0.0;

    r.val = 0.0;
    r.err = 0.0;

    // Handle x == 1.0 RJM

    if ( Math.abs( x - 1.0 ) < locEPS && (c - a - b) > 0.0 && c != 0.0 && ! c_neg_integer )
    {
        try
        {
            let lngamc   = { val: 0.0, err: 0.0, sign: 0 }; // Result;
            let lngamcab = { val: 0.0, err: 0.0 }; // Result;
            let lngamca  = { val: 0.0, err: 0.0, sign: 0 }; // Result;
            let lngamcb  = { val: 0.0, err: 0.0, sign: 0 }; // Result;
            let lngamc_sgn  = 0.0;
            let lngamca_sgn = 0.0;
            let lngamcb_sgn = 0.0;

            lngamc = gsl_sf_lngamma_sgn_e( c ); //, lngamc, lngamc_sgn );
            lngamcab = gsl_sf_lngamma_e( c - a - b );
            lngamca = gsl_sf_lngamma_sgn_e( c - a ); //, lngamca, lngamca_sgn );
            lngamcb = gsl_sf_lngamma_sgn_e( c - b ); //, lngamcb, lngamcb_sgn );
    
            r = gsl_sf_exp_err_e( lngamc.val + lngamcab.val - lngamca.val - lngamcb.val, lngamc.err + lngamcab.err + lngamca.err + lngamcb.err );
    
            r.val = r.val * lngamc.sign / (lngamca.sign * lngamcb.sign);
            return r;
        }
        catch ( e )
        {
            throw "SF.DomainException";
        }
    }
  
    if ( x < -1.0 || 1.0 <= x )
    {
        throw "SF.DomainException";
    }

    if ( c_neg_integer )
    {
        if ( ! (a_neg_integer && a > c + 0.1) )
        {
            throw "SF.DomainException";
        }
        if ( ! (b_neg_integer && b > c + 0.1) )
        {
            throw "SF.DomainException";
        }
    }

    if ( Math.abs( c - b ) < locEPS || Math.abs( c - a ) < locEPS )
    {
        return pow_omx( x, d );  // (1-x)^(c-a-b)
    }

    if ( a >= 0.0 && b >= 0.0 && c >=0.0 && x >= 0.0 && x < 0.995 )
    {
        // Series has all positive definite
        // terms and x is not close to 1.
        //
        return hyperg_2F1_series( a, b, c, x );
    }

    if ( Math.abs( a ) < 10.0 && Math.abs( b ) < 10.0 )
    {
        // a and b are not too large, so we attempt
        // variations on the series summation.

        if ( a_neg_integer )
        {
            return hyperg_2F1_series( rinta, b, c, x );
        }
        if ( b_neg_integer )
        {
            return hyperg_2F1_series( a, rintb, c, x );
        }

        if ( x < -0.25 )
        {
            return hyperg_2F1_luke( a, b, c, x );
        }
        else if ( x < 0.5 )
        {
            return hyperg_2F1_series( a, b, c, x );
        }
        else
        {
            if ( Math.abs( c ) > 10.0 )
            {
                return hyperg_2F1_series( a, b, c, x );
            }
            else
            {
                return hyperg_2F1_reflect( a, b, c, x );
            }
        }
    }
    else
    {
        // Either a or b or both large.
        // Introduce some new variables ap,bp so that bp is
        // the larger in magnitude.

        if ( Math.abs( a ) > Math.abs( b ) )
        {
            bp = a;
            ap = b;
        }
        else
        {
            bp = b;
            ap = a;
        }

        if ( x < 0.0 )
        {
            // What the hell, maybe Luke will converge.
            return hyperg_2F1_luke( a, b, c, x );
        }

        if ( Math.max( Math.abs( a ), 1.0 ) * Math.abs( bp ) * Math.abs( x ) < 2.0 * Math.abs( c ) )
        {
            // If c is large enough or x is small enough,
            // we can attempt the series anyway.
            return hyperg_2F1_series( a, b, c, x );
        }

        if ( Math.abs( bp * bp * x * x ) < 0.001 * Math.abs( bp ) && Math.abs( a ) < 10.0 )
        {
            // The famous but nearly worthless "large b" asymptotic.
            r = gsl_sf_hyperg_1F1_e( a, c, bp * x );
            r.err = 0.001 * Math.abs( r.val );
            return r;
        }

        // We give up.
        r.val = 0.0;
        r.err = 0.0;
        throw "SF.NotImplementedException";
    }

} // gsl_sf_hyperg_2F1_e

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_2F1_conj_e( aR, aI, c, x )
{
    var ax            = 0.0;
    var rintc         = 0.0;
    var c_neg_integer = false;
    var r             = { val: 0.0, err: 0.0 }; // Result;

    ax = Math.abs( x );
    rintc = Math.floor( c + 0.5 );
    c_neg_integer = ( c < 0.0 && Math.abs( c - rintc ) < locEPS );
    
    r.val = 0.0;
    r.err = 0.0;
    
    if ( ax >= 1.0 || c_neg_integer || c == 0.0 )
    {
        throw "SF.DomainException";
    }

    if ( (ax < 0.25 && Math.abs( aR ) < 20.0 && Math.abs( aI ) < 20.0) || (c > 0.0 && x > 0.0) )
    {
        return hyperg_2F1_conj_series( aR, aI, c, x );
    }
    else if ( Math.abs( aR ) < 10.0 && Math.abs( aI ) < 10.0 )
    {
        if ( x < -0.25 )
        {
            return hyperg_2F1_conj_luke( aR, aI, c, x );
        }
        else
        {
            return hyperg_2F1_conj_series( aR, aI, c, x );
        }
    }
    else
    {
        if ( x < 0.0 )
        {
            // What the hell, maybe Luke will converge.
            //
            return hyperg_2F1_conj_luke( aR, aI, c, x );
        }

        // Give up.
        //r.val = 0.0;
        //r.err = 0.0;
        throw "SF.NotImplementedException";
    }

} // gsl_sf_hyperg_2F1_conj_e

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_2F1_renorm_e( a, b, c, x )
{
    var rinta = 0.0;
    var rintb = 0.0;
    var rintc = 0.0;
    var a_neg_integer = false;
    var b_neg_integer = false;
    var c_neg_integer = false;
    var r = { val: 0.0, err: 0.0 }; // Result;

    rinta = Math.floor( a + 0.5 );
    rintb = Math.floor( b + 0.5 );
    rintc = Math.floor( c + 0.5 );
    a_neg_integer = ( a < 0.0 && Math.abs( a - rinta ) < locEPS );
    b_neg_integer = ( b < 0.0 && Math.abs( b - rintb ) < locEPS );
    c_neg_integer = ( c < 0.0 && Math.abs( c - rintc ) < locEPS );
  
    if ( c_neg_integer )
    {
        if ( (a_neg_integer && a > c + 0.1) || (b_neg_integer && b > c + 0.1) )
        {
            // 2F1 terminates early
            r.val = 0.0;
            r.err = 0.0;
            return r;
        }
        else
        {
            // 2F1 does not terminate early enough, so something survives
            // [Abramowitz+Stegun, 15.1.2]
            try
            {
                let g1 = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
                let g2 = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
                let g3 = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
                let g4 = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
                let g5 = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
                let sg = 0.0;
                let F  = { val: 0.0, err: 0.0 }; // Result;
                let ln_pre_val = 0.0;
                let ln_pre_err = 0.0;

                g1 = gsl_sf_lngamma_sgn_e( a - c + 1.0 ); //, g1, s1 );
                g2 = gsl_sf_lngamma_sgn_e( b - c + 1.0 ); //, g2, s2 );
                g3 = gsl_sf_lngamma_sgn_e( a ); //, g3, s3 );
                g4 = gsl_sf_lngamma_sgn_e( b ); //, g4, s4 );
                g5 = gsl_sf_lngamma_sgn_e( -c + 2.0 ); //, g5, s5 );

                F = gsl_sf_hyperg_2F1_e( a - c + 1.0, b - c + 1.0, -c + 2.0, x );
                ln_pre_val = g1.val + g2.val - g3.val - g4.val - g5.val;
                ln_pre_err = g1.err + g2.err + g3.err + g4.err + g5.err;
                sg = g1.sign * g2.sign * g3.sign * g4.sign * g5.sign;
                r = gsl_sf_exp_mult_err_e( ln_pre_val, ln_pre_err, sg * F.val, F.err );
                return r;
            }
            catch ( e )
            {
                throw "SF.DomainException";
            }
        }
    }
    else
    {
        // generic c
        let F   = { val: 0.0, err: 0.0 }; // Result;
        let lng  = { val: 0.0, err: 0.0, sign: 0.0 };

        lng = gsl_sf_lngamma_sgn_e( c ); //, lng, sgn );
        F = gsl_sf_hyperg_2F1_e( a, b, c, x );
        r = gsl_sf_exp_mult_err_e( -lng.val, lng.err, lng.sign * F.val, F.err ) ;
        return r;
    }

} // gsl_sf_hyperg_2F1_renorm_e

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_2F1_conj_renorm_e( aR, aI, c, x )
{
    var rintc = 0.0;
    var rinta = 0.0;
    var a_neg_integer = false;
    var c_neg_integer = false;
    var r = { val: 0.0, err: 0.0 }; // Result;

    rintc = Math.floor( c  + 0.5 );
    rinta = Math.floor( aR + 0.5 );
    a_neg_integer = ( aR < 0.0 && Math.abs( aR - rinta ) < locEPS && aI  == 0.0 );
    c_neg_integer = (  c < 0.0 && Math.abs( c  - rintc ) < locEPS );
    
    if ( c_neg_integer )
    {
        if ( a_neg_integer && aR > c + 0.1 )
        {
            // 2F1 terminates early
            r.val = 0.0;
            r.err = 0.0;
            return r;
        }
        else
        {
            // 2F1 does not terminate early enough, so something survives
            // [Abramowitz+Stegun, 15.1.2]
            try
            {
                let g1 = { val: 0.0, err: 0.0 }; // Result;
                let g2 = { val: 0.0, err: 0.0 }; // Result;
                let g3 = { val: 0.0, err: 0.0 }; // Result;
                let a1 = { val: 0.0, err: 0.0 }; // Result;
                let a2 = { val: 0.0, err: 0.0 }; // Result;
                let F  = { val: 0.0, err: 0.0 }; // Result;
                let ln_pre_val = 0.0;
                let ln_pre_err = 0.0;

                gsl_sf_lngamma_complex_e( aR - c + 1.0, aI, g1, a1 );
                gsl_sf_lngamma_complex_e( aR, aI, g2, a2 );
                g3 = gsl_sf_lngamma_e( -c + 2.0 );
                F = gsl_sf_hyperg_2F1_conj_e( aR - c + 1.0, aI, -c + 2.0, x );
                ln_pre_val = 2.0 * (g1.val - g2.val) - g3.val;
                ln_pre_err = 2.0 * (g1.err + g2.err) + g3.err;
                r = gsl_sf_exp_mult_err_e( ln_pre_val, ln_pre_err, F.val, F.err );
                return r;
            }
            catch ( e )
            {
                throw "SF.DomainException";
            }
        }
    }
    else
    {
        // generic c
        let F   = { val: 0.0, err: 0.0 }; // Result;
        let lng = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;

        lng = gsl_sf_lngamma_sgn_e( c ); //, lng, sgn );
        F = gsl_sf_hyperg_2F1_conj_e( aR, aI, c, x );
        r = gsl_sf_exp_mult_err_e( -lng.val, lng.err, lng.sign * F.val, F.err );
        return r;
    }

} // gsl_sf_hyperg_2F1_conj_renorm_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hyperg_2F1( a, b, c, x )
{ // gsl_sf_hyperg_2F1
    return EVAL_RESULT_4D( gsl_sf_hyperg_2F1_e, { x: a, dx: b, y: c, dy: x }, "gsl_sf_hyperg_2F1" );
} // gsl_sf_hyperg_2F1;

export function gsl_sf_hyperg_2F1_conj(aR, aI, c, x )
{ // gsl_sf_hyperg_2F1_conj
    return EVAL_RESULT_4D( gsl_sf_hyperg_2F1_conj_e, { x: aR, dx: aI, y: c, dy: x }, "gsl_sf_hyperg_2F1_conj" );
} // gsl_sf_hyperg_2F1_conj;

export function gsl_sf_hyperg_2F1_renorm(a, b, c, x )
{ // gsl_sf_hyperg_2F1_renorm
    return EVAL_RESULT_4D( gsl_sf_hyperg_2F1_renorm_e, { x: a, dx: b, y: c, dy: x }, "gsl_sf_hyperg_2F1_renorm" );
} // gsl_sf_hyperg_2F1_renorm;

export function gsl_sf_hyperg_2F1_conj_renorm(aR, aI, c, x )
{ // gsl_sf_hyperg_2F1_conj_renorm
    return EVAL_RESULT_4D( gsl_sf_hyperg_2F1_conj_renorm_e, { x: aR, dx: aI, y: c, dy: x }, "gsl_sf_hyperg_2F1_conj_renorm" );
} // gsl_sf_hyperg_2F1_conj_renorm;

// ----------------------------------------------------------------------------
// EOF SF-Hypergeometric2F1.mjs
