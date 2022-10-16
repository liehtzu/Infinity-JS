// SF-LegendrePolynomials.mjs
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

import { GSL_LOG_DBL_MIN }      from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_IS_ODD }           from "./SF-Math.mjs";
import { M_PI }                 from "./SF-Math.mjs";
import { M_LNPI }               from "./SF-Math.mjs";
import { gsl_sf_log_1plusx_e }  from "./SF-Logarithmic.mjs";
import { gsl_sf_bessel_J0_e }   from "./SF-BesselJ0.mjs";
import { gsl_sf_bessel_Jn_e }   from "./SF-BesselJn.mjs";
import { gsl_sf_lnpoch_e }      from "./SF-Pochhammer.mjs";

// Calculate P_m^m(x) from the analytic result:
//   P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
//            = 1 , m = 0
//
function legendre_Pmm(m, x)
{

    if (m == 0)
    {
        return 1.0;
    }
    else
    {
        var p_mm        = 0.0;
        var root_factor = 0.0;
        var fact_coeff  = 0.0;

        p_mm = 1.0;
        root_factor = Math.sqrt(1.0 - x) * Math.sqrt(1.0 + x);
        fact_coeff = 1.0;
        for (let i = 1; i <= m; i++)
        {
            p_mm = -p_mm * fact_coeff * root_factor;
            fact_coeff = fact_coeff + 2.0;
        }
        return p_mm;
    }

} // legendre_Pmm

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_legendre_P1_e(x)
{
    return { val: x, err: 0.0, e10: 0 };

} // gsl_sf_legendre_P1_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_P2_e(x)
{

    return { val: 0.5 * (3.0 * x * x - 1.0), err: GSL_DBL_EPSILON * (Math.abs(3.0 * x * x) + 1.0), e10: 0 };

} // gsl_sf_legendre_P2_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_P3_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r.val = 0.5 * x * (5.0 * x * x - 3.0);
    r.err = GSL_DBL_EPSILON * (Math.abs(r.val) + 0.5 * Math.abs(x) * (Math.abs(5.0 * x * x) + 3.0));
    return r;

} // gsl_sf_legendre_P3_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_Pl_e(l, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (l < 0 || x < -1.0 || x > 1.0)
    {
        throw "SF.DomainException";
    }
    else if (l == 0)
    {
        r.val = 1.0;
        r.err = 0.0;
    }
    else if (l == 1)
    {
        r.val = x;
        r.err = 0.0;
    }
    else if (l == 2)
    {
        r.val = 0.5 * (3.0 * x * x - 1.0);
        r.err = GSL_DBL_EPSILON * (Math.abs(3.0 * x * x) + 1.0);
        //result.err = 3.0 * GSL_DBL_EPSILON * Math.abs(result.val);
        // removed this old bogus estimate [GJ)
    }
    else if (x == 1.0)
    {
        r.val = 1.0;
        r.err = 0.0;
    }
    else if (x == -1.0)
    {
        if (GSL_IS_ODD(l))
        {
            r.val = -1.0;
        }
        else
        {
            r.val = 1.0;
        }
        r.err = 0.0;
    }
    else if (l < 100000)
    {
        // upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2}
        var p_ellm2 = 1.0;    // P_0(x)
        var p_ellm1 = x;      // P_1(x)
        var p_ell   = p_ellm1;
        var e_ellm2 = GSL_DBL_EPSILON;
        var e_ellm1 = Math.abs(x) * GSL_DBL_EPSILON;
        var e_ell   = e_ellm1;

        for (let ell = 2; ell <= l; ell++)
        {
            p_ell   = (x * (2 * ell - 1) * p_ellm1 - (ell - 1) * p_ellm2) / (ell);
            p_ellm2 = p_ellm1;
            p_ellm1 = p_ell;
            
            e_ell   = 0.5 * (Math.abs(x) * (2 * ell - 1) * e_ellm1 + (ell - 1) * e_ellm2) / (ell);
            e_ellm2 = e_ellm1;
            e_ellm1 = e_ell;
        }
        
        r.val = p_ell;
        r.err = e_ell + (l) * Math.abs(p_ell) * GSL_DBL_EPSILON;
    }
    else
    {
        // Asymptotic expansion.
        // [Olver, p. 473)
        //
        var pre = 0.0;
        var B00 = 0.0;
        var c1  = 0.0;
        var u   = 0.0;
        var th  = 0.0;
        var J0  = { val: 0.0, err: 0.0 }; // Result;
        var Jm1 = { val: 0.0, err: 0.0 }; // Result;
        var sin_th = 0.0;
        var cot_th = 0.0;

        u  = (l) + 0.5;
        th = Math.acos(x);
        J0  = gsl_sf_bessel_J0_e(u * th);
        Jm1 = gsl_sf_bessel_Jn_e(-1, u * th);
        // B00 = 1/8 (1 - th cot(th) / th^2
        // pre = sqrt(th/sin(th))
        //
        if (th < GSL_ROOT4_DBL_EPSILON)
        {
            B00 = (1.0 + th * th / 15.0) / 24.0;
            pre = 1.0 + th * th / 12.0;
        }
        else
        {
            sin_th = Math.sqrt(1.0 - x * x);
            cot_th = x / sin_th;
            B00 = 1.0 / 8.0 * (1.0 - th * cot_th) / (th * th);
            pre = Math.sqrt(th / sin_th);
        }
        
        c1 = th / u * B00;
        
        r.val = pre * (J0.val + c1 * Jm1.val);
        r.err = pre * (J0.err + Math.abs(c1) * Jm1.err);
        r.err = r.err + GSL_SQRT_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_legendre_Pl_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_Pl_array(lmax, x,  result_array)
{

    if (lmax < 0 || x < -1.0 || x > 1.0)
    {
        throw "SF.DomainException";
    }
    else if (lmax == 0)
    {
        result_array[0] = 1.0;
    }
    else if (lmax == 1)
    {
        result_array[0] = 1.0;
        result_array[1] = x;
    }
    else
    {
        // upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2}
        var p_ellm2 = 0.0;
        var p_ellm1 = 0.0;
        var p_ell   = 0.0;

        p_ellm2 = 1.0;  // P_0(x)
        p_ellm1 = x;    // P_1(x)
        p_ell   = p_ellm1;

        result_array[0] = 1.0;
        result_array[1] = x;

        for (let ell = 2; ell <= lmax; ell++)
        {
            p_ell = (x * (2 * ell - 1) * p_ellm1 - (ell - 1) * p_ellm2) / (ell);
            p_ellm2 = p_ellm1;
            p_ellm1 = p_ell;
            result_array[ell] = p_ell;
        }
    }

} // gsl_sf_legendre_Pl_array

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_Pl_deriv_array(lmax, x, result_array, result_deriv_array)
{
    var pre    = 0.0;
    var sgn    = 0.0;
    var diff_a = 0.0;
    var diff_b = 0.0;

    gsl_sf_legendre_Pl_array(lmax, x, result_array);

    if (lmax >= 0)
    {
        result_deriv_array[0] = 0.0;
    }
    if (lmax >= 1)
    {
        result_deriv_array[1] = 1.0;
    }

    if (Math.abs(x - 1.0) * (lmax + 1) * (lmax + 1) <  GSL_SQRT_DBL_EPSILON)
    {
        // x is near 1
        for (let ell = 2; ell <= lmax; ell++)
        {
            pre = 0.5 * (ell) * (ell + 1);
            result_deriv_array[ell] = pre * (1.0 - 0.25 * (1.0 - x) * (ell + 2) * (ell - 1));
        }
    }
    else if (Math.abs(x + 1.0) * (lmax + 1) * (lmax + 1) <  GSL_SQRT_DBL_EPSILON)
    {
        // x is near -1
        for (let ell = 2; ell <= lmax; ell++)
        {
            // derivative is odd in x for even ell
            if (GSL_IS_ODD(ell))
            {
                sgn = 1.0;
            }
            else
            {
                sgn = -1.0;
            }
            pre = sgn * 0.5 * (ell) * (ell + 1);
            result_deriv_array[ell] = pre * (1.0 - 0.25 * (1.0 + x) * (ell + 2) * (ell - 1));
        }
    }
    else
    {
        diff_a = 1.0 + x;
        diff_b = 1.0 - x;
        for (let ell = 2; ell <= lmax; ell++)
        {
            result_deriv_array[ell] = -(ell) * (x * result_array[ell] - result_array[ell-1]) / (diff_a * diff_b);
        }
    }

} // gsl_sf_legendre_Pl_deriv_array

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_Plm_e(l, m, x)
{
    var dif = 0.0;
    var sum = 0.0;
    var t_d = 0.0;
    var t_s = 0.0;
    var exp_check = 0.0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    // If l is large and m is large, then we have to worry
    // about overflow. Calculate an approximate exponent which
    // measures the normalization of this thing.
    //
    dif = (l - m);
    sum = (l + m);
    if (dif == 0.0)
    {
        t_d = 0.0;
    }
    else
    {
        t_d = 0.5 * dif * (Math.log(dif) - 1.0);
    }
    if (dif == 0.0)
    {
        t_s = 0.0;
    }
    else
    {
        t_s = 0.5 * sum * (Math.log(sum) - 1.0);
    }
    exp_check = 0.5 * Math.log(2.0 * (l) + 1.0) + t_d - t_s;


    if (m < 0 || l < m || x < -1.0 || x > 1.0)
    {
        throw "SF.DomainException";
    }
    else if (exp_check < GSL_LOG_DBL_MIN + 10.0)
    {
        // Bail out.
        throw "SF.OverflowException";
    }
    else
    {
        // Account for the error due to the
        // representation of 1-x.
        //
        var err_amp = 0.0;
        var p_mm    = 0.0;
        var p_mmp1  = 0.0;

        err_amp = 1.0 / (GSL_DBL_EPSILON + Math.abs(1.0 - Math.abs(x)));

        // P_m^m(x) and P_{m+1}^m(x)
        p_mm   = legendre_Pmm(m, x);
        p_mmp1 = x * (2 * m + 1) * p_mm;

        if (l == m)
        {
            r.val = p_mm;
            r.err = err_amp * 2.0 * GSL_DBL_EPSILON * Math.abs(p_mm);
        }
        else if (l == m + 1)
        {
            r.val = p_mmp1;
            r.err = err_amp * 2.0 * GSL_DBL_EPSILON * Math.abs(p_mmp1);
        }
        else
        {
            // upward recurrence: (l-m) P(l,m) = (2l-1) z P(l-1,m) - (l+m-1) P(l-2,m)
            // start at P(m,m), P(m+1,m)
            //
            var p_ellm2 = 0.0;
            var p_ellm1 = 0.0;
            var p_ell   = 0.0;

            p_ellm2 = p_mm;
            p_ellm1 = p_mmp1;
            p_ell   = 0.0;

            for (let ell = m + 2; ell <= l; ell++)
            {
                p_ell = (x * (2 * ell - 1) * p_ellm1 - (ell + m - 1) * p_ellm2) / (ell - m);
                p_ellm2 = p_ellm1;
                p_ellm1 = p_ell;
            }

            r.val = p_ell;
            r.err = err_amp * (0.5 * (l - m) + 1.0) * GSL_DBL_EPSILON * Math.abs(p_ell);
        }
    }

    return r;

} // gsl_sf_legendre_Plm_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_Plm_array(lmax, m, x, result_array)
{
    var dif = 0.0;
    var sum = 0.0;
    var t_d = 0.0;
    var t_s = 0.0;
    var exp_check = 0.0;

    // If l is large and m is large, then we have to worry
    // about overflow. Calculate an approximate exponent which
    // measures the normalization of this thing.
    //
    dif = (lmax - m);
    sum = (lmax + m);
    if (dif == 0.0)
    {
        t_d = 0.0;
        t_s = 0.0;
    }
    else
    {
        t_d = 0.5 * dif * (Math.log(dif) - 1.0);
        t_s = 0.5 * sum * (Math.log(sum) - 1.0);
    }
    exp_check = 0.5 * Math.log(2.0 * (lmax) + 1.0) + t_d - t_s;

    if (m < 0 || lmax < m || x < -1.0 || x > 1.0)
    {
        throw "SF.DomainException";
    }
    else if (m > 0 && (x == 1.0 || x == -1.0))
    {
        for (let ell = m; ell <= lmax; ell++)
        {
            result_array[ell-m] = 0.0;
        }
        return;
    }
    else if (exp_check < GSL_LOG_DBL_MIN + 10.0)
    {
        // Bail out.
        throw "SF.OverflowException";
    }
    else
    {
        var p_mm    = 0.0;
        var p_mmp1  = 0.0;
        var p_ellm2 = 0.0;
        var p_ellm1 = 0.0;
        var p_ell   = 0.0;

        p_mm   = legendre_Pmm(m, x);
        p_mmp1 = x * (2 * m + 1) * p_mm;
    
        if (lmax == m)
        {
            result_array[0] = p_mm;
            return;
        }
        else if (lmax == m + 1)
        {
            result_array[0] = p_mm;
            result_array[1] = p_mmp1;
            return;
        }
        else
        {
            p_ellm2 = p_mm;
            p_ellm1 = p_mmp1;
            p_ell = 0.0;
    
            result_array[0] = p_mm;
            result_array[1] = p_mmp1;
    
            for (let ell = m + 2; ell <= lmax; ell++)
            {
                p_ell = (x * (2 * ell - 1) * p_ellm1 - (ell + m - 1) * p_ellm2) / (ell - m);
                p_ellm2 = p_ellm1;
                p_ellm1 = p_ell;
                result_array[ell-m] = p_ell;
            }
    
            return;
        }
    }

} // gsl_sf_legendre_Plm_array

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_Plm_deriv_array(lmax, m, x, result_array, result_deriv_array)
{
    var sgn    = 0.0;
    var diff_a = 0.0;
    var diff_b = 0.0;

    if (m < 0 || m > lmax)
    {
        throw "SF.DomainException"; // WITH "m < 0 or m > lmax";
    }
    else if (m == 0)
    {
        // It is better to do m=0 this way, so we can more easily
        // trap the divergent case which can occur when m = 1.
        //
        gsl_sf_legendre_Pl_deriv_array(lmax, x, result_array, result_deriv_array);
    }
    else
    {
        gsl_sf_legendre_Plm_array(lmax, m, x, result_array);

        if (m == 1 && (1.0 - Math.abs(x) < GSL_DBL_EPSILON))
        {
            // This divergence is real and comes from the cusp-like
            // behaviour for m = 1. For example, P(1,1) = - Sqrt(1-x^2).
            //
            throw "SF.OverflowException"; // WITH "divergence near |x| = 1.0 since m = 1";
        }
        else if (m == 2 && (1.0 - Math.abs(x) < GSL_DBL_EPSILON))
        {
            // m = 2 gives a finite nonzero result for |x| near 1
            if (Math.abs(x - 1.0) < GSL_DBL_EPSILON)
            {
                for (let ell = m; ell <= lmax; ell++)
                {
                    result_deriv_array[ell-m] = -0.25 * x * (ell - 1) * (ell) * (ell + 1) * (ell + 2);
                }
            }
            else if (Math.abs(x + 1.0) < GSL_DBL_EPSILON)
            {
                for (let ell = m; ell <= lmax; ell++)
                {
                    if (GSL_IS_ODD(ell))
                    {
                        sgn = 1.0;
                    }
                    else
                    {
                        sgn = -1.0;
                    }
                    result_deriv_array[ell-m] = -0.25 * sgn * x * (ell - 1) * (ell) * (ell + 1) * (ell + 2);
                }
            }
            return;
        }
        else 
        {
            // m > 2 is easier to deal with since the endpoints always vanish
            if (1.0 - Math.abs(x) < GSL_DBL_EPSILON)
            {
                for (let ell = m; ell <= lmax; ell++)
                {
                    result_deriv_array[ell-m] = 0.0;
                }
                return;
            }
            else
            {
                diff_a = 1.0 + x;
                diff_b = 1.0 - x;
                result_deriv_array[0] = -(m) * x / (diff_a * diff_b) * result_array[0];
                if (lmax - m >= 1)
                {
                    result_deriv_array[1] = (2 * m + 1) * (x * result_deriv_array[0] + result_array[0]);
                }
                for (let ell = m + 2; ell <= lmax; ell++)
                {
                    result_deriv_array[ell-m] = - ((ell) * x * result_array[ell-m] - (ell + m) * result_array[ell-1-m]) / (diff_a * diff_b);
                }
                return;
            }
        }
    }

} // gsl_sf_legendre_Plm_deriv_array

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_sphPlm_e(l, m, x)
{
    var pre = 0.0;
    var P = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (m < 0 || l < m || x < -1.0 || x > 1.0)
    {
        throw "SF.DomainException";
    }
    else if (m == 0)
    {
        P = gsl_sf_legendre_Pl_e(l, x);
        pre = Math.sqrt((2 * l + 1) / (4.0 * M_PI));
        r.val = pre * P.val;
        r.err = pre * P.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }
    else if (x == 1.0 || x == -1.0)
    {
        // m > 0 here
        r.val = 0.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        // m > 0 and |x| < 1 here
        var lncirc = { val: 0.0, err: 0.0 }; // Result;
        var lnpoch = { val: 0.0, err: 0.0 }; // Result;
        var ex_pre = { val: 0.0, err: 0.0 }; // Result;
        var lnpre_val  = 0.0;
        var lnpre_err  = 0.0;
        var sr         = 0.0;
        var sgn        = 0.0;
        var y_mmp1_factor = 0.0;
        var y_mm       = 0.0;
        var y_mm_err   = 0.0;
        var y_mmp1     = 0.0;
        var y_mmp1_err = 0.0;

        // Starting value for recursion.
        // Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
        //
        if (GSL_IS_ODD(m))
        {
            sgn = -1.0;
        }
        else
        {
            sgn = 1.0;
        }
        y_mmp1_factor = x * Math.sqrt(2.0 * (m) + 3.0);
        lncirc = gsl_sf_log_1plusx_e(-x * x);
        lnpoch = gsl_sf_lnpoch_e((m), 0.5);  // Gamma(m+1/2)/Gamma(m)
        lnpre_val = -0.25 * M_LNPI + 0.5 * (lnpoch.val + (m) * lncirc.val);
        lnpre_err = 0.25 * M_LNPI * GSL_DBL_EPSILON + 0.5 * (lnpoch.err + Math.abs((m)) * lncirc.err);
        // Compute exp(ln_pre) with error term, avoiding call to gsl_sf_exp_err BJG
        ex_pre.val = Math.exp(lnpre_val);
        ex_pre.err = 2.0 * (Math.sinh(lnpre_err) + GSL_DBL_EPSILON) * ex_pre.val;
        sr     = Math.sqrt((2.0 + 1.0 / (m)) / (4.0 * M_PI));
        y_mm   = sgn * sr * ex_pre.val;
        y_mm_err = 2.0 * GSL_DBL_EPSILON * Math.abs(y_mm) + sr * ex_pre.err;
        y_mm_err = y_mm_err * (1.0 + 1.0 / (GSL_DBL_EPSILON + Math.abs(1.0 - x)));
        y_mmp1 = y_mmp1_factor * y_mm;
        y_mmp1_err = Math.abs(y_mmp1_factor) * y_mm_err;

        if (l == m)
        {
            r.val = y_mm;
            r.err = y_mm_err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(y_mm);
            return r;
        }
        else if (l == m + 1)
        {
            r.val = y_mmp1;
            r.err = y_mmp1_err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(y_mmp1);
            return r;
        }
        else
        {
            var y_ell     = 0.0;
            var y_ell_err = 0.0;
            var rat1      = 0.0;
            var rat2      = 0.0;
            var factor1   = 0.0;
            var factor2   = 0.0;

            // Compute Y_l^m, l > m+1, upward recursion on l.
            for (let ell = m + 2; ell <= l; ell++)
            {
                rat1    = (ell - m) / (ell + m);
                rat2    = (ell - m - 1) / (ell + m - 1);
                factor1 = Math.sqrt(rat1 * (2 * ell + 1) * (2 * ell - 1));
                factor2 = Math.sqrt(rat1 * rat2 * (2 * ell + 1) / (2 * ell - 3));
                y_ell   = (x * y_mmp1 * factor1 - (ell + m - 1) * y_mm * factor2) / (ell - m);
                y_mm    = y_mmp1;
                y_mmp1  = y_ell;
                
                y_ell_err  = 0.5 * (Math.abs(x * factor1) * y_mmp1_err + Math.abs((ell + m - 1) * factor2) * y_mm_err) / Math.abs((ell - m));
                y_mm_err   = y_mmp1_err;
                y_mmp1_err = y_ell_err;
            }
            r.val = y_ell;
            r.err = y_ell_err + (0.5 * (l - m) + 1.0) * GSL_DBL_EPSILON * Math.abs(y_ell);
            return r;
        }
    }

} // gsl_sf_legendre_sphPlm_e

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_sphPlm_array(lmax, m, x, result_array)
{

    if (m < 0 || lmax < m || x < -1.0 || x > 1.0)
    {
        throw "SF.DomainException";
    }
    else if (m > 0 && (x == 1.0 || x == -1.0))
    {
        for (let ell = m; ell <= lmax; ell++)
        {
            result_array[ell-m] = 0.0;
        }
        return;
    }
    else
    {
        var y_mm   = 0.0;
        var y_mmp1 = 0.0;

        if (m == 0)
        {
            y_mm   = 0.5 / M_SQRTPI;          // Y00 = 1/sqrt(4pi)
            y_mmp1 = x * M_SQRT3 * y_mm;
        }
        else
        {
            // |x| < 1 here
            var lncirc = { val: 0.0, err: 0.0 }; // Result;
            var lnpoch = { val: 0.0, err: 0.0 }; // Result;
            var lnpre  = 0.0;
            var sgn    = 0.0;

            if (GSL_IS_ODD(m))
            {
                sgn = -1.0;
            }
            else
            {
                sgn = 1.0;
            }
            lncirc = gsl_sf_log_1plusx_e(-x * x);
            lnpoch = gsl_sf_lnpoch_e((m), 0.5);  // Gamma(m+1/2)/Gamma(m)
            lnpre  = -0.25 * M_LNPI + 0.5 * (lnpoch.val + (m) * lncirc.val);
            y_mm   = Math.sqrt((2.0 + 1.0 / (m)) / (4.0 * M_PI)) * sgn * Math.exp(lnpre);
            y_mmp1 = x * Math.sqrt(2.0 * (m) + 3.0) * y_mm;
        }

        if (lmax == m)
        {
            result_array[0] = y_mm;
            return;
        }
        else if (lmax == m + 1)
        {
            result_array[0] = y_mm;
            result_array[1] = y_mmp1;
            return;
        }
        else
        {
            var y_ell   = 0.0;
            var rat1    = 0.0;
            var rat2    = 0.0;
            var factor1 = 0.0;
            var factor2 = 0.0;

            result_array[0] = y_mm;
            result_array[1] = y_mmp1;
            
            // Compute Y_l^m, l > m+1, upward recursion on l.
            for (let ell = m + 2; ell <= lmax; ell++)
            {
                rat1    = (ell - m) / (ell + m);
                rat2    = (ell - m - 1) / (ell + m - 1);
                factor1 = Math.sqrt(rat1 * (2 * ell + 1) * (2 * ell - 1));
                factor2 = Math.sqrt(rat1 * rat2 * (2 * ell + 1) / (2 * ell - 3));
                y_ell   = (x * y_mmp1 * factor1 - (ell + m - 1) * y_mm * factor2) / (ell - m);
                y_mm    = y_mmp1;
                y_mmp1  = y_ell;
                result_array[ell-m] = y_ell;
            }
        }

        return;
    }

} // gsl_sf_legendre_sphPlm_array

// ----------------------------------------------------------------------------

export function gsl_sf_legendre_sphPlm_deriv_array(lmax, m, x, result_array, result_deriv_array)
{
    var diff_a    = 0.0;
    var diff_b    = 0.0;
    var c1        = 0.0;
    var prefactor = 0.0;

    if (m < 0 || lmax < m || x < -1.0 || x > 1.0)
    {
        throw "SF.DomainException";
    }
    else if (m == 0)
    {
        // m = 0 is easy to trap
        gsl_sf_legendre_Pl_deriv_array(lmax, x, result_array, result_deriv_array);
        for (let ell = 0; ell <= lmax; ell++)
        {
            prefactor = Math.sqrt((2 * ell + 1) / (4.0 * M_PI));
            result_array[ell] = result_array[ell] * prefactor;
            result_deriv_array[ell] = result_deriv_array[ell] * prefactor;
        }
        return;
    }
    else if (m == 1)
    {
        // Trapping m = 1 is necessary because of the possible divergence.
        // Recall that this divergence is handled properly in ..._Plm_deriv_array(),
        // and the scaling factor is not large for small m, so we just scale.
        //
        gsl_sf_legendre_Plm_deriv_array(lmax, m, x, result_array, result_deriv_array);
        for (let ell = 1; ell <= lmax; ell++)
        {
            prefactor = Math.sqrt((2 * ell + 1) / (ell + 1) / (4.0 * M_PI * (ell)));
            result_array[ell-1] = result_array[ell-1] * prefactor;
            result_deriv_array[ell-1] = result_deriv_array[ell-1] * prefactor;
        }
        return;
    }
    else
    {
        // as for the derivative of P_lm, everything is regular for m >= 2

        gsl_sf_legendre_sphPlm_array(lmax, m, x, result_array);

        if (1.0 - Math.abs(x) < GSL_DBL_EPSILON)
        {
            for (let ell = m; ell <= lmax; ell++)
            {
                result_deriv_array[ell-m] = 0.0;
            }
            return;
        }
        else
        {
            diff_a = 1.0 + x;
            diff_b = 1.0 - x;
            result_deriv_array[0] = -(m) * x / (diff_a * diff_b) * result_array[0];
            if (lmax - m >= 1)
            {
                result_deriv_array[1] = Math.sqrt(2.0 * (m) + 3.0) * (x * result_deriv_array[0] + result_array[0]);
            }
            for (let ell = m + 2; ell <= lmax; ell++)
            {
                c1 = Math.sqrt(((2 * ell + 1) / (2 * ell - 1)) * ((ell - m) / (ell + m)));
                result_deriv_array[ell-m] = -((ell) * x * result_array[ell-m] - c1 * (ell + m) * result_array[ell-1-m]) / (diff_a * diff_b);
            }
            return;
        }
    }

} // gsl_sf_legendre_sphPlm_deriv_array

//int
//gsl_sf_legendre_array_size(const int lmax, const int m)
//{
//  return lmax-m+1;
//}

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

// function gsl_sf_legendre_P1(x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_legendre_P1
//     return EVAL_RESULT(gsl_sf_legendre_P1_e'Access, (x), "gsl_sf_legendre_P1");
// END gsl_sf_legendre_P1;

// function gsl_sf_legendre_P2(x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_legendre_P2
//     return EVAL_RESULT(gsl_sf_legendre_P2_e'Access, (x), "gsl_sf_legendre_P2");
// END gsl_sf_legendre_P2;

// function gsl_sf_legendre_P3(x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_legendre_P3
//     return EVAL_RESULT(gsl_sf_legendre_P3_e'Access, (x), "gsl_sf_legendre_P3");
// END gsl_sf_legendre_P3;

// function gsl_sf_legendre_Pl(l: INTEGER; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_legendre_Pl
//     return EVAL_RESULT(gsl_sf_legendre_Pl_e'Access, (l, x), "gsl_sf_legendre_Pl");
// END gsl_sf_legendre_Pl;

// function gsl_sf_legendre_Plm(l: INTEGER; m: INTEGER; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_legendre_Plm
//     return EVAL_RESULT(gsl_sf_legendre_Plm_e'Access, (l, m, x), "gsl_sf_legendre_Plm");
// END gsl_sf_legendre_Plm;

//function gsl_sf_legendre_sphPlm(l: INTEGER; m: INTEGER; x: LONG_FLOAT) return LONG_FLOAT IS
//BEGIN -- gsl_sf_legendre_sphPlm
//    return EVAL_RESULT(gsl_sf_legendre_sphPlm_e'Access, (l, m, x), "gsl_sf_legendre_sphPlm");
//END gsl_sf_legendre_sphPlm;

// ----------------------------------------------------------------------------
// EOF SF-LegendrePolynomials.mjs
