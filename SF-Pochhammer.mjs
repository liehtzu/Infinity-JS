// SF-Pochhammer.mjs
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

import { GSL_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX } from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN } from "./SF-Machine.mjs";
import { M_PI } from "./SF-Math.mjs";
import { M_LN2 } from "./SF-Math.mjs";
import { M_SQRT2 } from "./SF-Math.mjs";
import { M_SQRT3 } from "./SF-Math.mjs";
import { GSL_SIGN } from "./SF-Math.mjs";
import { gsl_sf_expm1_e } from "./SF-Exponential.mjs";
import { gsl_sf_exp_err_e } from "./SF-Exponential.mjs";
import { gsl_sf_log_1plusx_e } from "./SF-Logarithmic.mjs";
import { GSL_SF_GAMMA_XMAX } from "./SF-Gamma.mjs";
import { gsl_sf_gammainv_e } from "./SF-Gamma.mjs";
import { gsl_sf_psi_e } from "./SF-Psi.mjs";

const bern =//: CONSTANT ARRAY(0..20) OF LONG_FLOAT = --[21]
    [
    0.0,   // no element 0
   +0.833333333333333333333333333333333e-01,
   -0.138888888888888888888888888888888e-02,
   +0.330687830687830687830687830687830e-04,
   -0.826719576719576719576719576719576e-06,
   +0.208767569878680989792100903212014e-07,
   -0.528419013868749318484768220217955e-09,
   +0.133825365306846788328269809751291e-10,
   -0.338968029632258286683019539124944e-12,
   +0.858606205627784456413590545042562e-14,
   -0.217486869855806187304151642386591e-15,
   +0.550900282836022951520265260890225e-17,
   -0.139544646858125233407076862640635e-18,
   +0.353470703962946747169322997780379e-20,
   -0.895351742703754685040261131811274e-22,
   +0.226795245233768306031095073886816e-23,
   -0.574472439520264523834847971943400e-24,
   +0.145517247561486490186626486727132e-26,
   -0.368599494066531017818178247990866e-28,
   +0.933673425709504467203255515278562e-30,
   -0.236502241570062993455963519636983e-31
    ];


// ((a)_x - 1)/x in the "small x" region where
// cancellation must be controlled.
//
// Based on SLATEC DPOCH1().
//
//
// When Math.abs(X) is so small that substantial cancellation will occur if
// the straightforward formula is used, we use an expansion due
// to Fields and discussed by Y. L. Luke, The Special Functions and Their
// Approximations, Vol. 1, Academic Press, 1969, page 34.
//
// The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
//        (A+(X-1)/2)**X-- polynomial in (A+(X-1)/2)**(-2) .
// In order to maintain significance in POCH1, we write for positive a
//        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
//                       = 1.0 + Q*EXPREL(Q) .
// Likewise the polynomial is written
//        POLY = 1.0 + X*POLY1(A,X) .
// Thus,
//        POCH1(A,X) = (POCH(A,X) - 1) / X
//                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
//
//
function pochrel_smallx(a, x)
{
    // SQTBIG = 1.0D0/SQRT(24.0D0*D1MACH(1))
    // ALNEPS = LOG(D1MACH(3))
    const SQTBIG = 1.0 / (2.0 * M_SQRT2 * M_SQRT3 * GSL_SQRT_DBL_MIN);
    const ALNEPS = GSL_LOG_DBL_EPSILON - M_LN2;

    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0)
    {
        r = gsl_sf_psi_e(a);
    }
    else
    {
        var incr   = 0;
        var bp     = 0.0;
        var b      = 0.0;
        var dpoch1 = 0.0;
        var var1   = 0.0;
        var alnvar = 0.0;
        var q      = 0.0;
        var poly1  = 0.0;
        var binv   = 0.0;
        var dexprl = { val: 0.0, err: 0.0 }; // Result;

        if (a < -0.5)
        {
            bp = 1.0 - a - x;
        }
        else
        {
            bp = a;
        }
        if (bp < 10.0)
        {
            incr = Math.trunc(11.0 - bp);
        }
        else
        {
            incr = 0;
        }
        b      = bp + (incr);
        var1   = b + 0.5 * (x - 1.0);
        alnvar = Math.log(var1);
        q      = x * alnvar;
        poly1  = 0.0;

        if (var1 < SQTBIG)
        {
            var nterms = 0;
            var var2   = 0.0;
            var rho    = 0.0;
            var term   = 0.0;
            var gbk    = 0.0;
            var gbern  = []; //: ARRAY (0..23) OF LONG_FLOAT; --[24];

            nterms = Math.trunc(-0.5 * ALNEPS / alnvar + 1.0);
            var2   = (1.0 / var1) / var1;
            rho    = 0.5 * (x + 1.0);
            term   = var2;
            
            gbern[1] = 1.0;
            gbern[2] = -rho / 12.0;
            poly1 = gbern[2] * term;
            
            if (nterms > 20)
            {
                // NTERMS IS TOO BIG, MAYBE D1MACH(3) IS BAD
                // nterms = 20;
                //result.val = 0.0;
                //result.err = 0.0;
                //GSL_ERROR("error", GSL_ESANITY);
                throw "SF.SanityException";
            }
            
            for (let k = 2; k <= nterms; k++)
            {
                gbk = 0.0;
                for (let j = 1; j <= k; j++)
                {
                    gbk = gbk + bern[k-j+1] * gbern[j];
                }
                gbern[k+1] = -rho * gbk / (k);
                
                term  = term * ((2 * k - 2) - x) * ((2 * k - 1) - x) * var2;
                poly1 = poly1 + gbern[k+1] * term;
            }
        }
        
        dexprl = gsl_sf_expm1_e(q);
        dexprl.val = dexprl.val / q;
        poly1 = poly1 * (x - 1.0);
        dpoch1 = dexprl.val * (alnvar + q * poly1) + poly1;
        
        for (let i = incr - 1; i >= 0; i--)
        {
            // WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
            // TO OBTAIN DPOCH1(BP,X).
            binv   = 1.0 / (bp + (i));
            dpoch1 = (dpoch1 - binv) / (1.0 + x * binv);
        }
        
        if (bp == a)
        {
            r.val = dpoch1;
            r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs((incr)) + 1.0) * Math.abs(r.val);
        }
        else
        {
            // WE HAVE DPOCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A
            // REFLECTION FORMULA TO OBTAIN DPOCH1(A,X).
            var sinpxx = Math.sin(M_PI * x) / x;
            var sinpx2 = Math.sin(0.5 * M_PI * x);
            var t1     = sinpxx / Math.tan(M_PI * b);
            var t2     = 2.0 * sinpx2 * (sinpx2 / x);
            var trig   = t1 - t2;

            r.val = dpoch1 * (1.0 + x * trig) + trig;
            r.err = (Math.abs(dpoch1 * x) + 1.0) * GSL_DBL_EPSILON * (Math.abs(t1) + Math.abs(t2));
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs((incr)) + 1.0) * Math.abs(r.val);
        }
    }

    return r;

} // pochrel_smallx

// ----------------------------------------------------------------------------

// Assumes a>0 and a+x>0.

function lnpoch_pos(a, x)
{
    var absx = Math.abs(x);
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ((absx > 0.1 * a) || (absx * Math.log(Math.max(a, 2.0)) > 0.1))
    {
        if ((a < GSL_SF_GAMMA_XMAX) && (a + x < GSL_SF_GAMMA_XMAX))
        {
            // If we can do it by calculating the gamma functions
            // directly, then that will be more accurate than
            // doing the subtraction of the logs.
            //
            var g1 = { val: 0.0, err: 0.0 }; // Result;
            var g2 = { val: 0.0, err: 0.0 }; // Result;

            g1 = gsl_sf_gammainv_e(a);
            g2 = gsl_sf_gammainv_e(a + x);
            r.val = -Math.log(g2.val / g1.val);
            r.err = g1.err / Math.abs(g1.val) + g2.err / Math.abs(g2.val);
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
        else
        {
            // Otherwise we must do the subtraction.
            //
            var lg1 = { val: 0.0, err: 0.0 }; // Result;
            var lg2 = { val: 0.0, err: 0.0 }; // Result;

            lg1 = gsl_sf_lngamma_e(a);
            lg2 = gsl_sf_lngamma_e(a + x);
            r.val = lg2.val - lg1.val;
            r.err = lg2.err + lg1.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
    }
    else if ((absx < 0.1 * a) && (a > 15.0))
    {
        // Be careful about the implied subtraction.
        // Note that both a+x and and a must be
        // large here since a is not small
        // and x is not relatively large.
        // So we calculate using Stirling for Log[Gamma(z)].
        //
        //   Log[Gamma(a+x)/Gamma(a)] = x(Log[a]-1) + (x+a-1/2)Log[1+x/a]
        //                              + (1/(1+eps)   - 1) / (12 a)
        //                              - (1/(1+eps)^3 - 1) / (360 a^3)
        //                              + (1/(1+eps)^5 - 1) / (1260 a^5)
        //                              - (1/(1+eps)^7 - 1) / (1680 a^7)
        //                              + ...
        //
        const eps   = x / a;
        const den   = 1.0 + eps;
        const d3    = den * den * den;
        const d5    = d3 * den * den;
        const d7    = d5 * den * den;
        const c1    = -eps / den;
        const c3    = -eps * (3.0 + eps * (3.0 + eps)) / d3;
        const c5    = -eps * (5.0 + eps * (10.0 + eps * (10.0 + eps * (5.0 + eps)))) / d5;
        const c7    = -eps * (7.0 + eps * (21.0 + eps * (35.0 + eps * (35.0 + eps * (21.0 + eps * (7.0 + eps)))))) / d7;
        const p8    = Math.pow(1.0 + eps, 8); //gsl_sf_pow_int(1.0 + eps, 8);
        const c8    = 1.0 / p8                 - 1.0;  // these need not
        const c9    = 1.0 / (p8 * (1.0 + eps)) - 1.0;  // be very accurate
        const a4    = a * a * a * a;
        const a6    = a4 * a * a;
        const ser_1 = c1 + c3 / (30.0 * a * a) + c5 / (105.0 * a4) + c7 / (140.0 * a6);
        const ser_2 = c8 / (99.0 * a6 * a * a) - 691.0 / 360360.0 * c9 / (a6 * a4);
        const ser   = (ser_1 + ser_2) / (12.0 * a);
        const term1 = x * Math.log(a / M_E);
        var term2 = 0.0;
        var ln_1peps = { val: 0.0, err: 0.0 }; // Result;

        ln_1peps = gsl_sf_log_1plusx_e(eps);  // log(1 + x/a)
        term2 = (x + a - 0.5) * ln_1peps.val;
    
        r.val = term1 + term2 + ser;
        r.err = GSL_DBL_EPSILON * Math.abs(term1);
        r.err = r.err + Math.abs((x + a - 0.5) * ln_1peps.err);
        r.err = r.err + Math.abs(ln_1peps.val) * GSL_DBL_EPSILON * (Math.abs(x) + Math.abs(a) + 0.5);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var poch_rel = { val: 0.0, err: 0.0 }; // Result;
        var eps      = 0.0;

        poch_rel = pochrel_smallx(a, x);
        eps = x * poch_rel.val;
        r = gsl_sf_log_1plusx_e(eps);
        r.err = 2.0 * Math.abs(x * poch_rel.err / (1.0 + eps));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // lnpoch_pos


// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_lnpoch_e(a, x)
{

    if ((a <= 0.0) || (a + x <= 0.0))
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        return { val: 0.0, err: 0.0 };//, 0);
    }
    else
    {
        return lnpoch_pos(a, x);
    }

} // gsl_sf_lnpoch_e

// ----------------------------------------------------------------------------

export function gsl_sf_lnpoch_sgn_e(a, x) //, r, sgn)
{
    var sin_1  = 0.0;
    var sin_2  = 0.0;
    var lnterm = 0.0;
    //var s_apn  = 0.0;
    //var s_a    = 0.0;
    var sgn    = 0.0;

    var lnp_pos = { val: 0.0, err: 0.0 }; // Result;
    var r       = { val: 0.0, err: 0.0 }; // Result;

    var lg_apn  = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var lg_a    = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var rs      = { val: 0.0, err: 0.0, sign: 0.0 };

    if ((a == 0.0) || (a + x == 0.0))
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        sgn = 1.0;
        r.val = 0.0;
        r.err = 0.0;
    }
    else if ((a > 0.0) && (a + x > 0.0))
    {
        sgn = 1.0;
        r = lnpoch_pos(a, x);
    }
    else if ((a < 0.0) && (a + x < 0.0))
    {
        // Reduce to positive case using reflection.
        //
        sin_1 = Math.sin(M_PI * (1.0 - a));
        sin_2 = Math.sin(M_PI * (1.0 - a - x));
        if ((sin_1 == 0.0) || (sin_2 == 0.0))
        {
            throw "SF.DomainException";
        }
        else
        {
            lnp_pos = lnpoch_pos(1.0 - a, -x);
            lnterm = Math.log(Math.abs(sin_1 / sin_2));
            r.val = lnterm - lnp_pos.val;
            r.err = lnp_pos.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs(1.0 - a) + Math.abs(1.0 - a - x)) * Math.abs(lnterm);
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            sgn = GSL_SIGN(sin_1 * sin_2);
        }
    }
    else
    {
        // Evaluate gamma ratio directly.
        //
        lgn_apn = gsl_sf_lngamma_sgn_e(a + x); //, lg_apn, s_apn);
        lgn_a = gsl_sf_lngamma_sgn_e(a); //, lg_a, s_a);
        r.val = lg_apn.val - lg_a.val;
        r.err = lg_apn.err + lg_a.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        sgn = s_a * s_apn;
    }

    rs.val = r.val;
    rs.err = r.err;
    rs.sign = sgn;
    return rs;

} // gsl_sf_lnpoch_sgn_e

// ----------------------------------------------------------------------------

export function gsl_sf_poch_e(a, x)
{
    var sgn = 0.0;

    var lnpoch = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var r      = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0)
    {
        r.val = 1.0;
        r.err = 0.0;
    }
    else
    {
        lnpoch = gsl_sf_lnpoch_sgn_e(a, x); //, lnpoch, sgn);
        r = gsl_sf_exp_err_e(lnpoch.val, lnpoch.err);
        r.val = r.val * lnpoch.sign;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_poch_e

// ----------------------------------------------------------------------------

export function gsl_sf_pochrel_e(a, x)
{
    var absx = Math.abs(x);
    var absa = Math.abs(a);
    var sgn  = 0.0;
    var el   = 0.0;

    var lnpoch = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var r      = { val: 0.0, err: 0.0 }; // Result;
  
    if ((absx > 0.1 * absa) || (absx * Math.log(Math.max(absa, 2.0)) > 0.1))
    {
        lnpoch = gsl_sf_lnpoch_sgn_e(a, x); //, lnpoch, sgn);
        if (lnpoch.val > GSL_LOG_DBL_MAX)
        {
            throw "SF.OverflowException";
        }
        else
        {
            el = Math.exp(lnpoch.val);
            r.val = (lnpoch.sign * el - 1.0) / x;
            r.err = Math.abs(r.val) * (lnpoch.err + 2.0 * GSL_DBL_EPSILON);
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs(lnpoch.sign * el) + 1.0) / Math.abs(x);
        }
    }
    else
    {
        r = pochrel_smallx(a, x);
    }

    return r;

} // gsl_sf_pochrel_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

// function gsl_sf_lnpoch(a: LONG_FLOAT; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_lnpoch
//     return EVAL_RESULT(gsl_sf_lnpoch_e'Access, (a, x), "gsl_sf_lnpoch");
// END gsl_sf_lnpoch;

// function gsl_sf_poch(a: LONG_FLOAT; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_poch
//     return EVAL_RESULT(gsl_sf_poch_e'Access, (a, x), "gsl_sf_poch");
// END gsl_sf_poch;

// function gsl_sf_pochrel(a: LONG_FLOAT; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_pochrel
//     return EVAL_RESULT(gsl_sf_pochrel_e'Access, (a, x), "gsl_sf_pochrel");
// END gsl_sf_pochrel;

// ----------------------------------------------------------------------------
// EOF SF-Pochhammer.mjs
