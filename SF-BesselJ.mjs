// SF-BesselJ.adb
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

import { M_E }                             from "./SF-Math.mjs";
import { M_PI }                            from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }                 from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }            from "./SF-Machine.mjs";
import { GSL_ROOT4_DBL_EPSILON }           from "./SF-Machine.mjs";
import { GSL_ROOT6_DBL_EPSILON }           from "./SF-Machine.mjs";
import { GSL_DBL_MIN }                     from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MIN }                from "./SF-Machine.mjs";
import { gsl_sf_sin_e }                    from "./SF-Trigonometric.mjs";
import { gsl_sf_cos_e }                    from "./SF-Trigonometric.mjs";
import { gsl_sf_hypot }                    from "./SF-Trigonometric.mjs";
import { gsl_sf_bessel_IJ_taylor_e }       from "./SF-Bessel.mjs";
import { gsl_sf_bessel_J_CF1 }             from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Jnu_asympx_e }      from "./SF-Bessel.mjs";
import { gsl_sf_bessel_Jnu_asymp_Olver_e } from "./SF-BesselOlver.mjs";

import { EVAL_RESULT_D }                   from "./SF-Evaluate.mjs";
import { EVAL_RESULT_ID }                  from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_j0s_e(x)
{
    var ax = Math.abs(x);
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (ax < 0.5)
    {
        const y  = x * x;
        const c1 = -1.0 / 6.0;
        const c2 =  1.0 / 120.0;
        const c3 = -1.0 / 5040.0;
        const c4 =  1.0 / 362880.0;
        const c5 = -1.0 / 39916800.0;
        const c6 =  1.0 / 6227020800.0;

        r.val = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * c6)))));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var sin_result = { val: 0.0, err: 0.0 }; // Result;

        sin_result = gsl_sf_sin_e(x);
        r.val = sin_result.val / x;
        r.err = Math.abs(sin_result.err / x);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_j0s_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_j1s_e(x)
{
    var ax = Math.abs(x);
    var r  = { val: 0.0, err: 0.0 }; // Result;
   
    if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (ax < 3.1 * GSL_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (ax < 0.25)
    {
        const y  = x * x;
        const c1 = -1.0 / 10.0;
        const c2 =  1.0 / 280.0;
        const c3 = -1.0 / 15120.0;
        const c4 =  1.0 / 1330560.0;
        const c5 = -1.0 / 172972800.0;
        const sum = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * c5))));

        r.val = x / 3.0 * sum;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var cos_result = { val: 0.0, err: 0.0 }; // Result;
        var sin_result = { val: 0.0, err: 0.0 }; // Result;
        var cos_x      = 0.0;
        var sin_x      = 0.0;

        cos_result = gsl_sf_cos_e(x);
        sin_result = gsl_sf_sin_e(x);
        cos_x = cos_result.val;
        sin_x = sin_result.val;
        r.val = (sin_x / x - cos_x) / x;
        r.err = (Math.abs(sin_result.err / x) + Math.abs(cos_result.err)) / Math.abs(x);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs(sin_x / (x * x)) + Math.abs(cos_x / x));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_j1s_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_j2s_e(x)
{
    var ax = Math.abs(x);
    var r  = { val: 0.0, err: 0.0 }; // Result;
    
    if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (ax < 4.0 * GSL_SQRT_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (ax < 1.3)
    {
        const y  = x * x;
        const c1 = -1.0 / 14.0;
        const c2 =  1.0 / 504.0;
        const c3 = -1.0 / 33264.0;
        const c4 =  1.0 / 3459456.0;
        const c5 = -1.0 / 518918400.0;
        const c6 =  1.0 / 105859353600.0;
        const c7 = -1.0 / 28158588057600.0;
        const c8 =  1.0 / 9461285587353600.0;
        const c9 = -1.0 / 3916972233164390400.0;
        const sum = 1.0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * (c7 + y * (c8 + y * c9))))))));

        r.val = y / 15.0 * sum;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        var cos_result = { val: 0.0, err: 0.0 }; // Result;
        var sin_result = { val: 0.0, err: 0.0 }; // Result;
        var cos_x      = 0.0;
        var sin_x      = 0.0;
        var f          = 0.0;

        cos_result = gsl_sf_cos_e(x);
        sin_result = gsl_sf_sin_e(x );
        cos_x = cos_result.val;
        sin_x = sin_result.val;
        f = (3.0 / (x * x) - 1.0);
        r.val = (f * sin_x - 3.0 * cos_x / x) / x;
        r.err = Math.abs(f * sin_result.err / x) + Math.abs((3.0 * cos_result.err / x) / x);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs(f * sin_x / x) + 3.0 * Math.abs(cos_x / (x * x)));
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_bessel_j2s_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_jl_e(l, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (l < 0 || x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        if (l > 0)
        {
            r.val = 0.0;
        }
        else
        {
            r.val = 1.0;
        }
        r.err = 0.0;
    }
    else if (l == 0)
    {
        r = gsl_sf_bessel_j0s_e(x);
    }
    else if (l == 1)
    {
        r = gsl_sf_bessel_j1s_e(x);
    }
    else if (l == 2)
    {
        r = gsl_sf_bessel_j2s_e(x);
    }
    else if (x * x < 10.0 * ((l) + 0.5) / M_E)
    {
        let b   = { val: 0.0, err: 0.0 }; // Result;
        let pre = 0.0;

        b = gsl_sf_bessel_IJ_taylor_e((l) + 0.5, x, -1, 50, GSL_DBL_EPSILON);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = pre * b.val;
        r.err = pre * b.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (GSL_ROOT4_DBL_EPSILON * x > (l * l + l + 1))
    {
        let b   = { val: 0.0, err: 0.0 }; // Result;
        let pre = 0.0;

        b = gsl_sf_bessel_Jnu_asympx_e((l) + 0.5, x);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = pre * b.val;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val) + pre * b.err;
    }
    else if ((l) > 1.0 / GSL_ROOT6_DBL_EPSILON)
    {
        let b   = { val: 0.0, err: 0.0 }; // Result;
        let pre = 0.0;

        b = gsl_sf_bessel_Jnu_asymp_Olver_e((l) + 0.5, x);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = pre * b.val;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val) + pre * b.err;
    }
    else if (x > 1000.0 && x > (l * l))
    {
        // We need this path to avoid feeding large x to CF1 below;
        let b   = { val: 0.0, err: 0.0 }; // Result;
        let pre = 0.0;

        b = gsl_sf_bessel_Jnu_asympx_e((l) + 0.5, x);
        pre = Math.sqrt((0.5 * M_PI) / x);
        r.val = pre * b.val;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val) + pre * b.err;
    }
    else
    {
        let pre    = 0.0;
        let sgn    = 0.0;
        let ratio  = 0.0;
        let jellp1 = 0.0;
        let jell   = 0.0;
        let jellm1 = 0.0;
        let j0_result = { val: 0.0, err: 0.0 }; // Result;
        let j1_result = { val: 0.0, err: 0.0 }; // Result;
        let rs = { val: 0.0, err: 0.0, sign: 0.0 };

        // The CF1 call will hit 10000 iterations for x > 10000 + l
        rs = gsl_sf_bessel_J_CF1((l) + 0.5, x); //, ratio, sgn);
        ratio = rs.ratio;
        sgn = rs.sign;
        jellp1 = GSL_SQRT_DBL_EPSILON * ratio;
        jell   = GSL_SQRT_DBL_EPSILON;
        for (let ell = l; ell >= 1; ell--)
        {
            jellm1 = -jellp1 + (2 * ell + 1) / x * jell;
            jellp1 = jell;
            jell   = jellm1;
        }
        
        if (Math.abs(jell) > Math.abs(jellp1))
        {
            j0_result = gsl_sf_bessel_j0s_e(x);
            pre = GSL_SQRT_DBL_EPSILON / jell;
            r.val = j0_result.val * pre;
            r.err = j0_result.err * Math.abs(pre);
            r.err = r.err + 4.0 * GSL_DBL_EPSILON * (0.5 * (l) + 1.0) * Math.abs(r.val);
        }
        else
        {
            j1_result = gsl_sf_bessel_j1s_e(x);
            pre = GSL_SQRT_DBL_EPSILON / jellp1;
            r.val = j1_result.val * pre;
            r.err = j1_result.err * Math.abs(pre);
            r.err = r.err + 4.0 * GSL_DBL_EPSILON * (0.5 * (l) + 1.0) * Math.abs(r.val);
        }
    }

    return r;

} //  gsl_sf_bessel_jl_e

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_jl_array( lmax, x, result_array )
{
  
    if ( lmax < 0 || x < 0.0 )
    {
        //FOR j IN 0 .. lmax LOOP
        //    result_array(j) = 0.0;
        //END LOOP;
        throw "SF.DomainException";
    }
    else if ( x == 0.0 )
    {
        for ( let j = 1; j <= lmax; j++ )
        {
            result_array[j] = 0.0;
        }
        result_array[0] = 1.0;
    }
    else
    {
        var r_jellp1 = { val: 0.0, err: 0.0 }; // Result;
        var r_jell   = { val: 0.0, err: 0.0 }; // Result;
        var jellp1   = 0.0;
        var jell     = 0.0;
        var jellm1   = 0.0;

        r_jellp1 = gsl_sf_bessel_jl_e( lmax + 1, x );
        r_jell   = gsl_sf_bessel_jl_e( lmax,     x );
        jellp1 = r_jellp1.val;
        jell   = r_jell.val;
        result_array[lmax] = jell;
        for ( let ell = lmax; ell >= 1; ell-- )
        {
            jellm1 = -jellp1 + (2 * ell + 1) / x * jell;
            jellp1 = jell;
            jell   = jellm1;
            result_array[ell-1] = jellm1;
        }
    }

} // gsl_sf_bessel_jl_array

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_jl_steed_array( lmax, x, jl_x )
{
  
    if ( lmax < 0 || x < 0.0 )
    {
        //FOR j IN 0 .. lmax LOOP
        //    jl_x(j) = 0.0;
        //END LOOP;
        throw "SF.DomainException";
    }
    else if ( x == 0.0 )
    {
        for ( let j = 1; j <= lmax; j++ )
        {
            jl_x[j] = 0.0;
        }
        jl_x[0] = 1.0;
    }
    else if ( x < 2.0 * GSL_ROOT4_DBL_EPSILON )
    {
        // first two terms of Taylor series
        var inv_fact = 1.0;  // 1/(1 3 5 ... (2l+1))
        var x_l      = 1.0;  // x^l

        for ( let l = 0; l <= lmax; l++ )
        {
            jl_x[l]  = x_l * inv_fact;
            jl_x[l]  = jl_x[l] * (1.0 - 0.5 * x * x / (2 * l + 3));
            inv_fact = inv_fact / (2 * l + 3);
            x_l      = x_l      * x;
        }
    }
    else
    {
        // Steed/Barnett algorithm [Comp. Phys. Comm. 21, 297 (1981)]
        var x_inv = 1.0 / x;
        var W     = 2.0 * x_inv;
        var F     = 1.0;
        var FP    = (lmax + 1) * x_inv;
        var B     = 2.0 * FP + x_inv;
        var end0  = B + 20000.0 * W;
        var D     = 1.0 / B;
        var del   = -D;

        var XP2   = 0.0;
        var PL    = 0.0;
        var L     = 0;

        FP = FP + del;
        
        // continued fraction
        while ( true )
        {
            B = B + W;
            D = 1.0 / (B - D);
            del = del * (B * D - 1.0);
            FP = FP + del;
            if ( D < 0.0 )
            {
                F = -F;
            }
            if ( B > end0 )
            {
                throw "SF.MaxIterationsException";
            }
            if ( ! (Math.abs( del ) >= Math.abs( FP ) * GSL_DBL_EPSILON) ) break;
        }
        
        FP = FP * F;
        
        if ( lmax > 0 )
        {
            // downward recursion
            XP2 = FP;
            PL = (lmax) * x_inv;
            L  = lmax;
            jl_x[lmax] = F;
            for ( let LP = 1; LP <= lmax; LP++ )
            {
                jl_x[L-1] = PL * jl_x[L] + XP2;
                FP = PL * jl_x[L-1] - jl_x[L];
                XP2 = FP;
                PL = PL - x_inv;
                L = L - 1;
            }
            F = jl_x[0];
        }
        
        // normalization
        W = x_inv / gsl_sf_hypot( FP, F );
        jl_x[0] = W * F;
        if ( lmax > 0 )
        {
            for ( let L = 1; L <= lmax; L++ )
            {
                jl_x[L] = jl_x[L] * W;
            }
        }
    }

} // gsl_sf_bessel_jl_steed_array

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_bessel_j0s( x )
{ // gsl_sf_bessel_j0s
    return EVAL_RESULT_D( gsl_sf_bessel_j0s_e, x, "gsl_sf_bessel_j0s" );
} // gsl_sf_bessel_j0s;

export function gsl_sf_bessel_j1s( x )
{ // gsl_sf_bessel_j1s
    return EVAL_RESULT_D( gsl_sf_bessel_j1s_e, x, "gsl_sf_bessel_j1s" );
} // gsl_sf_bessel_j1s;

export function gsl_sf_bessel_j2s( x )
{ // gsl_sf_bessel_j2s
    return EVAL_RESULT_D( gsl_sf_bessel_j2s_e, x, "gsl_sf_bessel_j2s" );
} // gsl_sf_bessel_j2s;

export function gsl_sf_bessel_jl( l, x )
{ // gsl_sf_bessel_jl
    return EVAL_RESULT_ID( gsl_sf_bessel_jl_e, { i: l, x: x }, "gsl_sf_bessel_jl" );
} // gsl_sf_bessel_jl;

// ----------------------------------------------------------------------------
// EOF SF-BesselJ.mjs
