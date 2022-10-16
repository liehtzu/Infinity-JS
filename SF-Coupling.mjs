// SF-Coupling.mjs
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
// Translation to Ada: Igor Izvarin

import { GSL_IS_ODD }      from "./SF-Math.mjs";
import { gsl_sf_choose_e } from "./SF-Gamma.mjs";
import { gsl_sf_fact_e }   from "./SF-Gamma.mjs";
import { GSL_DBL_EPSILON } from "./SF-Machine.mjs";

function locMax3 (a, b, c )
{ // locMax3
    return Math.max( Math.max( a, b ), c );
} // locMax3

// ----------------------------------------------------------------------------

function locMin3( a, b, c )
{ // locMin3
    return Math.min( Math.min( a, b ), c );
} // locMin3

// ----------------------------------------------------------------------------

function locMin5( a, b, c, d, e )
{
    var f = Math.min( a, b );
    var g = Math.min( c, d );
    var h = Math.min( f, g );
    return Math.min( e, h );
} // locMin5

// ----------------------------------------------------------------------------

// See: [Thompson, Atlas for Computing Mathematical Functions]

function delta1( ta, tb, tc )
{
    var f1 = { val: 0.0, err: 0.0 };
    var f2 = { val: 0.0, err: 0.0 };
    var f3 = { val: 0.0, err: 0.0 };
    var f4 = { val: 0.0, err: 0.0 };
    var r  = { val: 0.0, err: 0.0 };

    try
    {
        f1 = gsl_sf_fact_e( (ta + tb - tc) / 2 );
        f2 = gsl_sf_fact_e( (ta + tc - tb) / 2 );
        f3 = gsl_sf_fact_e( (tb + tc - ta) / 2 );
        f4 = gsl_sf_fact_e( (ta + tb + tc) / 2 + 1 );
        r.val = f1.val * f2.val * f3.val / f4.val;
        r.err = 4.0 * GSL_DBL_EPSILON * Math.abs( r.val );
        return r;
    }
    catch ( e )
    {
        throw "SF.OverflowException";
    }

} // delta1

// ----------------------------------------------------------------------------

function triangle_selection_fails( two_ja, two_jb, two_jc )
{ // triangle_selection_fails
    return ((two_jb < Math.abs( two_ja - two_jc )) || (two_jb > two_ja + two_jc));
} // triangle_selection_fails

// ----------------------------------------------------------------------------

function m_selection_fails( two_ja, two_jb, two_jc, two_ma, two_mb, two_mc )
{ // m_selection_fails
    return (
           Math.abs( two_ma ) > two_ja 
        || Math.abs( two_mb ) > two_jb
        || Math.abs( two_mc ) > two_jc
        || GSL_IS_ODD( two_ja + two_ma )
        || GSL_IS_ODD( two_jb + two_mb )
        || GSL_IS_ODD( two_jc + two_mc )
        || (two_ma + two_mb + two_mc) != 0
            );
} // m_selection_fails

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_coupling_3j_e( two_ja,
                              two_jb,
                              two_jc,
                              two_ma,
                              two_mb,
                              two_mc )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( two_ja < 0 || two_jb < 0 || two_jc < 0 )
    {
        throw "SF.DomainException";
    }
    else if ( triangle_selection_fails( two_ja, two_jb, two_jc )
        || m_selection_fails( two_ja, two_jb, two_jc, two_ma, two_mb, two_mc )
        )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        var jca  = (-two_ja + two_jb + two_jc) / 2;
        var jcb  = ( two_ja - two_jb + two_jc) / 2;
        var jcc  = ( two_ja + two_jb - two_jc) / 2;
        var jmma = ( two_ja - two_ma) / 2;
        var jmmb = ( two_jb - two_mb) / 2;
        var jmmc = ( two_jc - two_mc) / 2;
        var jpma = ( two_ja + two_ma) / 2;
        var jpmb = ( two_jb + two_mb) / 2;
        var jpmc = ( two_jc + two_mc) / 2;
        var jsum = ( two_ja + two_jb + two_jc) / 2;
        var kmin = locMax3( 0, jpmb - jmmc, jmma - jpmc );
        var kmax = locMin3( jcc, jmma, jpmb );
        var sign = 0;
        var sum_pos = 0.0;
        var sum_neg = 0.0;
        var norm    = 0.0;
        var term    = 0.0;
        var bc1  = { val: 0.0, err: 0.0 }; // Result;
        var bc2  = { val: 0.0, err: 0.0 }; // Result;
        var bc3  = { val: 0.0, err: 0.0 }; // Result;
        var bcn1 = { val: 0.0, err: 0.0 }; // Result;
        var bcn2 = { val: 0.0, err: 0.0 }; // Result;
        var bcd1 = { val: 0.0, err: 0.0 }; // Result;
        var bcd2 = { val: 0.0, err: 0.0 }; // Result;
        var bcd3 = { val: 0.0, err: 0.0 }; // Result;
        var bcd4 = { val: 0.0, err: 0.0 }; // Result;

        if ( GSL_IS_ODD( kmin - jpma + jmmb ) )
        {
            sign = -1;
        }
        else
        {
            sign = 1;
        }

        try
        {
            bcn1 = gsl_sf_choose_e( two_ja,   jcc  );
            bcn2 = gsl_sf_choose_e( two_jb,   jcc  );
            bcd1 = gsl_sf_choose_e( jsum + 1, jcc  );
            bcd2 = gsl_sf_choose_e( two_ja,   jmma );
            bcd3 = gsl_sf_choose_e( two_jb,   jmmb );
            bcd4 = gsl_sf_choose_e( two_jc,   jpmc );
        }
        catch ( e )
        {
            throw "SF.OverflowException";
        }
        
        norm = Math.sqrt( bcn1.val * bcn2.val )
                / Math.sqrt( bcd1.val * bcd2.val * bcd3.val * bcd4.val * ((two_jc) + 1.0) );
        
        for ( let k = kmin; k <= kmax; k++ )
        {
            try
            {
                bc1 = gsl_sf_choose_e( jcc, k );
                bc2 = gsl_sf_choose_e( jcb, jmma - k );
                bc3 = gsl_sf_choose_e( jca, jpmb - k );
            }
            catch ( e )
            { 
                throw "SF.OverflowException";
            }
            
            term = bc1.val * bc2.val * bc3.val;
            
            if ( sign < 0 )
            {
                sum_neg = sum_neg + norm * term;
            }
            else
            {
                sum_pos = sum_pos + norm * term;
            }
            
            sign = -sign;
        }
        
        r.val = sum_pos - sum_neg;
        r.err = 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (kmax - kmin) * Math.abs( r.val );
        
        return r;
    }

} // gsl_sf_coupling_3j_e

// ----------------------------------------------------------------------------

//int
//gsl_sf_coupling_6j_INCORRECT_e(int two_ja, int two_jb, int two_jc,
//                               int two_jd, int two_je, int two_jf,
//                               gsl_sf_result * result)
//{
//  return gsl_sf_coupling_6j_e(two_ja, two_jb, two_je, two_jd, two_jc, two_jf, result);
//}

// ----------------------------------------------------------------------------

export function gsl_sf_coupling_6j_e( two_ja,
                              two_jb,
                              two_jc,
                              two_jd,
                              two_je,
                              two_jf )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( two_ja < 0 || two_jb < 0 || two_jc < 0
       || two_jd < 0 || two_je < 0 || two_jf < 0
       )
    {
        throw "SF.DomainException";
    }
    else if ( triangle_selection_fails( two_ja, two_jb, two_jc )
          || triangle_selection_fails( two_ja, two_je, two_jf )
          || triangle_selection_fails( two_jb, two_jd, two_jf )
          || triangle_selection_fails( two_je, two_jd, two_jc )
          )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        var n1        = { val: 0.0, err: 0.0 }; // Result;
        var d1        = { val: 0.0, err: 0.0 }; // Result;
        var d2        = { val: 0.0, err: 0.0 }; // Result;
        var d3        = { val: 0.0, err: 0.0 }; // Result;
        var d4        = { val: 0.0, err: 0.0 }; // Result;
        var d5        = { val: 0.0, err: 0.0 }; // Result;
        var d6        = { val: 0.0, err: 0.0 }; // Result;
        var norm      = 0.0;
        var tk        = 0;
        var tkmin     = 0;
        var tkmax     = 0;
        var phase     = 0.0;
        var sum_pos   = 0.0;
        var sum_neg   = 0.0;
        var sumsq_err = 0.0;
        var term      = 0.0;
        var term_err  = 0.0;
        var den_1     = { val: 0.0, err: 0.0 }; // Result;
        var den_2     = { val: 0.0, err: 0.0 }; // Result;
        var d1_a      = { val: 0.0, err: 0.0 }; // Result;
        var d1_b      = { val: 0.0, err: 0.0 }; // Result;

        try
        {
            d1 = delta1( two_ja, two_jb, two_jc );
            d2 = delta1( two_ja, two_je, two_jf );
            d3 = delta1( two_jb, two_jd, two_jf );
            d4 = delta1( two_je, two_jd, two_jc );
        }
        catch ( e )
        {
            throw "SF.OverflowException";
        }
        norm = Math.sqrt( d1.val ) * Math.sqrt( d2.val ) * Math.sqrt( d3.val ) * Math.sqrt( d4.val );
        
        tkmin = locMax3( 0,
                        two_ja + two_jd - two_jc - two_jf,
                        two_jb + two_je - two_jc - two_jf );
        
        tkmax = locMin5( two_ja + two_jb + two_je + two_jd + 2,
                        two_ja + two_jb - two_jc,
                        two_je + two_jd - two_jc,
                        two_ja + two_je - two_jf,
                        two_jb + two_jd - two_jf );
        
        if ( GSL_IS_ODD( (two_ja + two_jb + two_je + two_jd + tkmin) / 2 ) )
        {
            phase = -1.0;
        }
        else
        {
            phase =  1.0;
        }
        
        tk = tkmin;
        while ( tk <= tkmax )
        {
            try
            {
                n1 = gsl_sf_fact_e( (two_ja + two_jb + two_je + two_jd - tk) / 2 + 1 );
                d1_a = gsl_sf_fact_e( tk / 2 );
                d1_b = gsl_sf_fact_e( (two_jc + two_jf - two_ja - two_jd + tk) / 2 );
                d2 = gsl_sf_fact_e( (two_jc + two_jf - two_jb - two_je + tk) / 2 );
                d3 = gsl_sf_fact_e( (two_ja + two_jb - two_jc - tk) / 2 );
                d4 = gsl_sf_fact_e( (two_je + two_jd - two_jc - tk) / 2 );
                d5 = gsl_sf_fact_e( (two_ja + two_je - two_jf - tk) / 2 );
                d6 = gsl_sf_fact_e( (two_jb + two_jd - two_jf - tk) / 2 );
            }
            catch ( e )
            { 
                throw "SF.OverflowException";
            }
            
            d1.val = d1_a.val * d1_b.val;
            d1.err = d1_a.err * Math.abs( d1_b.val ) + Math.abs( d1_a.val ) * d1_b.err;
            
            den_1.val = d1.val * d2.val * d3.val;
            den_1.err = d1.err * Math.abs( d2.val * d3.val );
            den_1.err = den_1.err + d2.err * Math.abs( d1.val * d3.val );
            den_1.err = den_1.err + d3.err * Math.abs( d1.val * d2.val );
            
            den_2.val = d4.val * d5.val * d6.val;
            den_2.err = d4.err * Math.abs( d5.val * d6.val );
            den_2.err = den_2.err + d5.err * Math.abs( d4.val * d6.val );
            den_2.err = den_2.err + d6.err * Math.abs( d4.val * d5.val );
            
            term  = phase * n1.val / den_1.val / den_2.val;
            phase = -phase;
            term_err = n1.err / Math.abs( den_1.val ) / Math.abs( den_2.val );
            term_err = term_err + Math.abs( term / den_1.val ) * den_1.err;
            term_err = term_err + Math.abs( term / den_2.val ) * den_2.err;
            
            if ( term >= 0.0 )
            {
                sum_pos = sum_pos + norm * term;
            }
            else
            {
                sum_neg = sum_neg - norm * term;
            }
            
            sumsq_err = sumsq_err + norm * norm * term_err * term_err;
            tk = tk + 2;
        }
        
        r.val = sum_pos - sum_neg;
        r.err = 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
        r.err = r.err + Math.sqrt( sumsq_err / (0.5 * (tkmax - tkmin) + 1.0) );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (tkmax - tkmin + 2) * Math.abs( r.val );
        
        return r;
    }

} // gsl_sf_coupling_6j_e

// ----------------------------------------------------------------------------

//int
//gsl_sf_coupling_RacahW_e(int two_ja, int two_jb, int two_jc,
//                         int two_jd, int two_je, int two_jf,
//                         gsl_sf_result * result)
//{
//  int status = gsl_sf_coupling_6j_e(two_ja, two_jb, two_je, two_jd, two_jc, two_jf, result);
//  int phase_sum = (two_ja + two_jb + two_jc + two_jd)/2;
//  result.val *= ( GSL_IS_ODD(phase_sum) ? -1.0 : 1.0 );
//  return status;
//}

// ----------------------------------------------------------------------------

export function gsl_sf_coupling_9j_e( two_ja,
                              two_jb,
                              two_jc,
                              two_jd,
                              two_je,
                              two_jf,
                              two_jg,
                              two_jh,
                              two_ji )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( two_ja < 0 || two_jb < 0 || two_jc < 0
       || two_jd < 0 || two_je < 0 || two_jf < 0
       || two_jg < 0 || two_jh < 0 || two_ji < 0
       )
    {
        throw "SF.DomainException";
    }
    else if ( triangle_selection_fails( two_ja, two_jb, two_jc )
            || triangle_selection_fails( two_jd, two_je, two_jf )
            || triangle_selection_fails( two_jg, two_jh, two_ji )
            || triangle_selection_fails( two_ja, two_jd, two_jg )
            || triangle_selection_fails( two_jb, two_je, two_jh )
            || triangle_selection_fails( two_jc, two_jf, two_ji )
       )
    {
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        var tk        = 0;
        var tkmin     = 0;
        var tkmax     = 0;
        var sum_pos   = 0.0;
        var sum_neg   = 0.0;
        var sumsq_err = 0.0;
        var phase     = 0.0;
        var term      = 0.0;
        var term_err  = 0.0;
        var s1        = { val: 0.0, err: 0.0 }; // Result;
        var s2        = { val: 0.0, err: 0.0 }; // Result;
        var s3        = { val: 0.0, err: 0.0 }; // Result;

        tkmin = locMax3( Math.abs( two_ja - two_ji ), Math.abs( two_jh - two_jd ), Math.abs( two_jb - two_jf ) );
        tkmax = locMin3( two_ja + two_ji, two_jh + two_jd, two_jb + two_jf );
        tk = tkmin;
        while ( tk <= tkmax )
        {
            try
            {
                s1 = gsl_sf_coupling_6j_e( two_ja, two_ji, tk,  two_jh, two_jd, two_jg );
                s2 = gsl_sf_coupling_6j_e( two_jb, two_jf, tk,  two_jd, two_jh, two_je );
                s3 = gsl_sf_coupling_6j_e( two_ja, two_ji, tk,  two_jf, two_jb, two_jc );
            }
            catch ( e )
            {
                throw "SF.OverflowException";
            }
            term = s1.val * s2.val * s3.val;
            term_err = s1.err * Math.abs( s2.val * s3.val );
            term_err = term_err + s2.err * Math.abs( s1.val * s3.val );
            term_err = term_err + s3.err * Math.abs( s1.val * s2.val );
            
            if ( term >= 0.0 )
            {
                sum_pos = sum_pos + (tk + 1) * term;
            }
            else
            {
                sum_neg = sum_neg - (tk + 1) * term;
            }
            
            sumsq_err = sumsq_err + ((tk + 1) * term_err) * ((tk + 1) * term_err);
            tk = tk + 2;
        }
        
        if ( GSL_IS_ODD( tkmin ) )
        {
            phase = -1.0;
        }
        else
        {
            phase =  1.0;
        }
        
        r.val = phase * (sum_pos - sum_neg);
        r.err = 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
        r.err = r.err + Math.sqrt( sumsq_err / (0.5 * (tkmax - tkmin) + 1.0) );
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (tkmax - tkmin + 2) * Math.abs( r.val );
        
        return r;
    }

} // 

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

// export function gsl_sf_coupling_3j( two_ja, two_jb, two_jc, two_ma, two_mb, two_mc )
// { // gsl_sf_coupling_3j
//     return EVAL_RESULT( gsl_sf_coupling_3j_e, (two_ja, two_jb, two_jc, two_ma, two_mb, two_mc), "gsl_sf_coupling_3j" );
// } // gsl_sf_coupling_3j


//double gsl_sf_coupling_6j_INCORRECT(int two_ja, int two_jb, int two_jc,
//                                    int two_jd, int two_je, int two_jf)
//{
//  EVAL_RESULT(gsl_sf_coupling_6j_INCORRECT_e(two_ja, two_jb, two_jc,
//                                             two_jd, two_je, two_jf,
//                                             &result));
//}


// export function gsl_sf_coupling_6j( two_ja, two_jb, two_jc, two_jd, two_je, two_jf )
// { // gsl_sf_coupling_6j
//     return EVAL_RESULT( gsl_sf_coupling_6j_e, (two_ja, two_jb, two_jc, two_jd, two_je, two_jf), "gsl_sf_coupling_6j" );
// } // gsl_sf_coupling_6j


//double gsl_sf_coupling_RacahW(int two_ja, int two_jb, int two_jc,
//                          int two_jd, int two_je, int two_jf)
//{
//  EVAL_RESULT(gsl_sf_coupling_RacahW_e(two_ja, two_jb, two_jc,
//                                      two_jd, two_je, two_jf,
//                                      &result));
//}


// export function gsl_sf_coupling_9j( two_ja, two_jb, two_jc,
//                            two_jd, two_je, two_jf,
//                            two_jg, two_jh, two_ji )
// { // gsl_sf_coupling_9j
//     return EVAL_RESULT( gsl_sf_coupling_9j_e, (two_ja, two_jb, two_jc,
//                                                      two_jd, two_je, two_jf,
//                                                      two_jg, two_jh, two_ji),
//                                                      "gsl_sf_coupling_9j" );
// } // gsl_sf_coupling_9j

// ----------------------------------------------------------------------------
// EOF SF-Coupling.mjs
