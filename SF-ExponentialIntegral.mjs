// SF-ExponentialIntegral.mjs
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

import { GSL_LOG_DBL_MIN }       from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON }   from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { gsl_sf_gamma_inc_e }    from "./SF-GammaIncomplete.mjs";
import { cheb_eval_e }           from "./SF-Chebyshev.mjs";

import { EVAL_RESULT_D }         from "./SF-Evaluate.mjs";
import { EVAL_RESULT_ID }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

//
// Chebyshev expansions: based on SLATEC e1.f, W. Fullerton
// 
// Series for AE11       on the interval -1.00000D-01 to  0.
//                                        with weighted error   1.76E-17
//                                         log weighted error  16.75
//                               significant figures required  15.70
//                                    decimal places required  17.55
//
//
// Series for AE12       on the interval -2.50000D-01 to -1.00000D-01
//                                        with weighted error   5.83E-17
//                                         log weighted error  16.23
//                               significant figures required  15.76
//                                    decimal places required  16.93
//
//
// Series for E11        on the interval -4.00000D+00 to -1.00000D+00
//                                        with weighted error   1.08E-18
//                                         log weighted error  17.97
//                               significant figures required  19.02
//                                    decimal places required  18.61
//
//
// Series for E12        on the interval -1.00000D+00 to  1.00000D+00
//                                        with weighted error   3.15E-18
//                                         log weighted error  17.50
//                        approx significant figures required  15.8
//                                    decimal places required  18.10
//
//
// Series for AE13       on the interval  2.50000D-01 to  1.00000D+00
//                                        with weighted error   2.34E-17
//                                         log weighted error  16.63
//                               significant figures required  16.14
//                                    decimal places required  17.33
//
//
// Series for AE14       on the interval  0.          to  2.50000D-01
//                                        with weighted error   5.41E-17
//                                         log weighted error  16.27
//                               significant figures required  15.38
//                                    decimal places required  16.97
//

const AE11_data =
    [
    0.121503239716065790,
   -0.065088778513550150,
    0.004897651357459670,
   -0.000649237843027216,
    0.000093840434587471,
    0.000000420236380882,
   -0.000008113374735904,
    0.000002804247688663,
    0.000000056487164441,
   -0.000000344809174450,
    0.000000058209273578,
    0.000000038711426349,
   -0.000000012453235014,
   -0.000000005118504888,
    0.000000002148771527,
    0.000000000868459898,
   -0.000000000343650105,
   -0.000000000179796603,
    0.000000000047442060,
    0.000000000040423282,
   -0.000000000003543928,
   -0.000000000008853444,
   -0.000000000000960151,
    0.000000000001692921,
    0.000000000000607990,
   -0.000000000000224338,
   -0.000000000000200327,
   -0.000000000000006246,
    0.000000000000045571,
    0.000000000000016383,
   -0.000000000000005561,
   -0.000000000000006074,
   -0.000000000000000862,
    0.000000000000001223,
    0.000000000000000716,
   -0.000000000000000024,
   -0.000000000000000201,
   -0.000000000000000082,
    0.000000000000000017
    ];
const AE11_cs = { length: 38, c: AE11_data, order: 38, a: -1.0, b: 1.0, order_sp: 20 };

const AE12_data =
    [
    0.582417495134726740,
   -0.158348850905782750,
   -0.006764275590323141,
    0.005125843950185725,
    0.000435232492169391,
   -0.000143613366305483,
   -0.000041801320556301,
   -0.000002713395758640,
    0.000001151381913647,
    0.000000420650022012,
    0.000000066581901391,
    0.000000000662143777,
   -0.000000002844104870,
   -0.000000000940724197,
   -0.000000000177476602,
   -0.000000000015830222,
    0.000000000002905732,
    0.000000000001769356,
    0.000000000000492735,
    0.000000000000093709,
    0.000000000000010707,
   -0.000000000000000537,
   -0.000000000000000716,
   -0.000000000000000244,
   -0.000000000000000058
    ];
const AE12_cs = { length: 24, c: AE12_data, order: 24, a: -1.0, b: 1.0, order_sp: 15 };

const E11_data =
    [
  -16.11346165557149402600,
    7.79407277874268027690,
   -1.95540581886314195070,
    0.37337293866277945612,
   -0.05692503191092901938,
    0.00721107776966009185,
   -0.00078104901449841593,
    0.00007388093356262168,
   -0.00000620286187580820,
    0.00000046816002303176,
   -0.00000003209288853329,
    0.00000000201519974874,
   -0.00000000011673686816,
    0.00000000000627627066,
   -0.00000000000031481541,
    0.00000000000001479904,
   -0.00000000000000065457,
    0.00000000000000002733,
   -0.00000000000000000108
    ];
const E11_cs = { length: 18, c: E11_data, order: 18, a: -1.0, b: 1.0, order_sp: 13 };

const E12_data =
    [
   -0.03739021479220279500,
    0.04272398606220957700,
   -0.13031820798497005440,
    0.01441912402469889073,
   -0.00134617078051068022,
    0.00010731029253063780,
   -0.00000742999951611943,
    0.00000045377325690753,
   -0.00000002476417211390,
    0.00000000122076581374,
   -0.00000000005485141480,
    0.00000000000226362142,
   -0.00000000000008635897,
    0.00000000000000306291,
   -0.00000000000000010148,
    0.00000000000000000315
    ];
const E12_cs = { length: 15, c: E12_data, order: 15, a: -1.0, b: 1.0, order_sp: 10 };

const AE13_data =
    [
   -0.605773246640603460,
   -0.112535243483660900,
    0.013432266247902779,
   -0.001926845187381145,
    0.000309118337720603,
   -0.000053564132129618,
    0.000009827812880247,
   -0.000001885368984916,
    0.000000374943193568,
   -0.000000076823455870,
    0.000000016143270567,
   -0.000000003466802211,
    0.000000000758754209,
   -0.000000000168864333,
    0.000000000038145706,
   -0.000000000008733026,
    0.000000000002023672,
   -0.000000000000474132,
    0.000000000000112211,
   -0.000000000000026804,
    0.000000000000006457,
   -0.000000000000001568,
    0.000000000000000383,
   -0.000000000000000094,
    0.000000000000000023
    ];
const AE13_cs = { length: 24, c: AE13_data, order: 24, a: -1.0, b: 1.0, order_sp: 15 };

const AE14_data =
    [
   -0.18929180007530170,
   -0.08648117855259871,
    0.00722410154374659,
   -0.00080975594575573,
    0.00010999134432661,
   -0.00001717332998937,
    0.00000298562751447,
   -0.00000056596491457,
    0.00000011526808397,
   -0.00000002495030440,
    0.00000000569232420,
   -0.00000000135995766,
    0.00000000033846628,
   -0.00000000008737853,
    0.00000000002331588,
   -0.00000000000641148,
    0.00000000000181224,
   -0.00000000000052538,
    0.00000000000015592,
   -0.00000000000004729,
    0.00000000000001463,
   -0.00000000000000461,
    0.00000000000000148,
   -0.00000000000000048,
    0.00000000000000016,
   -0.00000000000000005
    ];
const AE14_cs = { length: 25, c: AE14_data, order: 25, a: -1.0, b: 1.0, order_sp: 13 };

// ----------------------------------------------------------------------------

// implementation for E1, allowing for scaling by exp(x)
export function expint_E1_impl(x, scale)
{
    const xmaxt = -GSL_LOG_DBL_MIN;        // XMAXT = -LOG (R1MACH(1))
    const xmax  = xmaxt - Math.log(xmaxt); // XMAX = XMAXT - LOG(XMAXT)

    var s            = 0.0;
    var ln_term      = 0.0;
    var scale_factor = 0.0;
    var result_c     = { val: 0.0, err: 0.0 }; // Result;
    var r            = { val: 0.0, err: 0.0 }; // Result;

    if (x < -xmax && ! scale)
    {
        throw "Overflow"; //SF.OverflowException;
    }
    else if (x <= -10.0)
    {
        s = 1.0 / x;
        if (! scale)
        {
           s = s * Math.exp(-x);
        }
        result_c = cheb_eval_e(AE11_cs, 20.0 / x + 1.0);
        r.val = s * (1.0 + result_c.val);
        r.err = s * result_c.err;
        r.err = r.err + 2.0 * Number.EPSILON * (Math.abs(x) + 1.0) * Math.abs(r.val);
    }
    else if (x <= -4.0)
    {
        s = 1.0 / x;
        if (! scale)
        {
           s = s * Math.exp(-x);
        }
        result_c = cheb_eval_e(AE12_cs, (40.0 / x + 7.0) / 3.0);
        r.val = s * (1.0 + result_c.val);
        r.err = s * result_c.err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else if (x <= -1.0)
    {
        ln_term = -Math.log(Math.abs(x));
        if (scale)
        {
            scale_factor = Math.exp(x);
        }
        else
        {
            scale_factor = 1.0;
        }
        result_c = cheb_eval_e(E11_cs, (2.0 * x + 5.0) / 3.0);
        r.val = scale_factor * (ln_term + result_c.val);
        r.err = scale_factor * (result_c.err + Number.EPSILON * Math.abs(ln_term));
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else if (x == 0.0)
    {
        throw "Domain error"; //SF.DomainException;
    }
    else if (x <= 1.0)
    {
        ln_term = -Math.log(Math.abs(x));
        if (scale)
        {
            scale_factor = Math.exp(x);
        }
        else
        {
            scale_factor = 1.0;
        }
        result_c = cheb_eval_e(E12_cs, x);
        r.val = scale_factor * (ln_term - 0.6875 + x + result_c.val);
        r.err = scale_factor * (result_c.err + Number.EPSILON * Math.abs(ln_term));
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else if (x <= 4.0)
    {
        s = 1.0 / x;
        if (! scale)
        {
           s = s * Math.exp(-x);
        }
        result_c = cheb_eval_e(AE13_cs, (8.0 / x - 5.0) / 3.0);
        r.val = s * (1.0 + result_c.val);
        r.err = s * result_c.err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else if (x <= xmax || scale)
    {
        s = 1.0 / x;
        if (! scale)
        {
           s = s * Math.exp(-x);
        }
        result_c = cheb_eval_e(AE14_cs, 8.0 / x - 1.0);
        r.val = s * (1.0 +  result_c.val);
        r.err = s * (Number.EPSILON + result_c.err);
        r.err = r.err + 2.0 * (x + 1.0) * Number.EPSILON * Math.abs(r.val);
        if (r.val == 0.0)
        {
            throw "SF.UnderflowException";
        }
    }
    else
    {
        throw "SF.UnderflowException";
    }

    return r;

} // expint_E1_impl

// ----------------------------------------------------------------------------

function expint_E2_impl(x, scale)
{
    const xmaxt = -GSL_LOG_DBL_MIN;
    const xmax  = xmaxt - Math.log(xmaxt);
    const c1  = -2.0;
    const c2  =  6.0;
    const c3  = -24.0;
    const c4  =  120.0;
    const c5  = -720.0;
    const c6  =  5040.0;
    const c7  = -40320.0;
    const c8  =  362880.0;
    const c9  = -3628800.0;
    const c10 =  39916800.0;
    const c11 = -479001600.0;
    const c12 =  6227020800.0;
    const c13 = -87178291200.0;
    var ex    = 0.0;
    var s     = 0.0;
    var y     = 0.0;
    var sum6  = 0.0;
    var sum   = 0.0;
    var result_E1 = { val: 0.0, err: 0.0 }; // Result;
    var r         = { val: 0.0, err: 0.0 }; // Result;
  
    if (x < -xmax && !scale)
    {
        throw "Overflow Exception";
    }
    else if (x == 0.0)
    {
        r.val = 1.0; // ??? (scale ? 1.0 : 1.0);
        r.err = 0.0;
    }
    else if (x < 100.0)
    {
        if (scale)
        {
            ex = 1.0;
        }
        else
        {
            ex = Math.exp(-x);
        }
        result_E1 = expint_E1_impl(x, scale);
        r.val = ex - x * result_E1.val;
        r.err = Number.EPSILON * ex + Math.abs(x) * result_E1.err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else if (x < xmax || scale)
    {
        if (scale)
        {
            s = 1.0;
        }
        else
        {
            s = Math.exp(-x);
        }
        y = 1.0 / x;
        sum6 = c6 + y * (c7 + y * (c8 + y * (c9 + y * (c10 + y * (c11 + y * (c12 + y * c13))))));
        sum  = y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * sum6)))));
        r.val = s * (1.0 + sum) / x;
        r.err = 2.0 * (x + 1.0) * Number.EPSILON * r.val;
        if (r.val == 0.0)
        {
            throw "SF.UnderflowException";
        }
    }
    else
    {
        throw "SF.UnderflowException";
    }

    return r;

} // expint_E2_impl

// ----------------------------------------------------------------------------

function expint_En_impl(n, x, scale)
{
    var result_g  = { val: 0.0, err: 0.0 }; // Result;
    var r         = { val: 0.0, err: 0.0 }; // Result;
    var prefactor = 0.0;
    var scale_factor = 0.0;

    if (n < 0)
    {
        throw "SF.DomainException";
    }
    else if (n == 0)
    {
        if (x == 0.0)
        {
            throw "SF.DomainException";
        }
        else
        {
            if (scale)
            {
                r.val = 1.0 / x;
            }
            else
            {
                r.val = Math.exp(-x) / x;
            }
            r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
            //CHECK_UNDERFLOW(result);
        }
    }
    else if (n == 1)
    {
        r = expint_E1_impl(x, scale);
    }
    else if (n == 2)
    {
        r = expint_E2_impl(x, scale);
    }
    else
    {
        if (x < 0.0)
        {
            throw "SF.DomainException";
        }
        if (x == 0.0)
        {
            if (scale)
            {
               r.val = Math.exp(x) * (1.0 / (n - 1.0));
            }
            else
            {
               r.val = 1.0 / (n - 1.0);
            }
            r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
            //CHECK_UNDERFLOW(result);
        }
        else
        {
            prefactor = x ** (n - 1);
            result_g = gsl_sf_gamma_inc_e(1 - n, x);
            if (scale)
            {
                scale_factor = Math.exp(x);
            }
            else
            {
                scale_factor = 1.0;
            }
            r.val = scale_factor * prefactor * result_g.val;
            r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
            r.err = r.err + 2.0 * Math.abs(scale_factor * prefactor * result_g.err);
            //CHECK_UNDERFLOW(result);
        }
    }

    return r;

} // expint_En_impl

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*


export function gsl_sf_expint_E1_e(x)
{ // gsl_sf_expint_E1_e
    return expint_E1_impl(x, false);
} // gsl_sf_expint_E1_e

// ----------------------------------------------------------------------------

export function gsl_sf_expint_E1_scaled_e(x)
{ // gsl_sf_expint_E1_scaled_e
    return expint_E1_impl(x, true);
} // gsl_sf_expint_E1_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_expint_E2_e(x)
{ // gsl_sf_expint_E2_e
    return expint_E2_impl(x, false);
} // gsl_sf_expint_E2_e

// ----------------------------------------------------------------------------

export function gsl_sf_expint_E2_scaled_e(x)
{ // gsl_sf_expint_E2_scaled_e
    return expint_E2_impl(x, true);
} // gsl_sf_expint_E2_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_expint_En_e(n, x)
{ // gsl_sf_expint_En_e
    return expint_En_impl(n, x, false);
} // gsl_sf_expint_En_e

// ----------------------------------------------------------------------------

export function gsl_sf_expint_En_scaled_e(n, x)
{ // gsl_sf_expint_En_scaled_e
    return expint_En_impl(n, x, true);
} // gsl_sf_expint_En_scaled_e

// ----------------------------------------------------------------------------

export function gsl_sf_expint_Ei_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_expint_E1_e(-x);
    r.val = -r.val;
    return r;

} // gsl_sf_expint_Ei_e

// ----------------------------------------------------------------------------

export function gsl_sf_expint_Ei_scaled_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    r = gsl_sf_expint_E1_scaled_e(-x);
    r.val = -r.val;
    return r;

} // gsl_sf_expint_Ei_scaled_e

// ============================================================================

 const expint3_data =
    [
    1.269198414221126014,
   -0.248846446384140982,
    0.80526220717231041e-01,
   -0.25772733251968330e-01,
    0.7599878873073774e-02,
   -0.2030695581940405e-02,
    0.490834586699330e-03,
   -0.107682239142021e-03,
    0.21551726264290e-04,
   -0.3956705137384e-05,
    0.6699240933896e-06,
   -0.105132180807e-06,
    0.15362580199e-07,
   -0.20990960364e-08,
    0.2692109538e-09,
   -0.325195242e-10,
    0.37114816e-11,
   -0.4013652e-12,
    0.412334e-13,
   -0.40338e-14,
    0.3766e-15,
   -0.336e-16,
    0.29e-17,
   -0.2e-18
    ];
const expint3_cs = { length: 23, c: expint3_data, order: 23, a: -1.0, b: 1.0, order_sp: 15 };

const expint3a_data =
    [
    1.9270464955068273729,
   -0.349293565204813805e-01,
    0.14503383718983009e-02,
   -0.8925336718327903e-04,
    0.70542392191184e-05,
   -0.6671727454761e-06,
    0.724267589982e-07,
   -0.87825825606e-08,
    0.11672234428e-08,
   -0.1676631281e-09,
    0.257550158e-10,
   -0.41957888e-11,
    0.7201041e-12,
   -0.1294906e-12,
    0.24287e-13,
   -0.47331e-14,
    0.95531e-15,
   -0.1991e-15,
    0.428e-16,
   -0.94e-17,
    0.21e-17,
   -0.5e-18,
    0.1e-18
    ];
const expint3a_cs = { length: 22, c: expint3a_data, order: 22, a: -1.0, b: 1.0, order_sp: 10 };

// ----------------------------------------------------------------------------

export function gsl_sf_expint_3_e(x)
{
    const val_infinity = 0.892979511569249211;
    var t = 0.0;
    var s = 0.0;
    var c = { val: 0.0, err: 0.0 }; // RESULT;
    var r = { val: 0.0, err: 0.0 }; // RESULT;

    if (x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x < 1.6 * GSL_ROOT3_DBL_EPSILON)
    {
        r.val = x;
        r.err = 0.0;
    }
    else if (x <= 2.0)
    {
        t = x * x * x / 4.0 - 1.0;
        c = cheb_eval_e(expint3_cs, t);
        r.val = x * c.val;
        r.err = x * c.err;
    }
    else if (x < (-GSL_LOG_DBL_EPSILON) ** (1.0/3.0))
    {
        t = 16.0 / (x * x * x) - 1.0;
        s = Math.exp(-x * x * x) / (3.0 * x * x);
        c = cheb_eval_e(expint3a_cs, t);
        r.val = val_infinity - c.val * s;
        r.err = val_infinity * GSL_DBL_EPSILON + s * c.err;
    }
    else
    {
        r.val = val_infinity;
        r.err = val_infinity * GSL_DBL_EPSILON;
    }

    return r;

} // gsl_sf_expint_3_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_expint_E1( x )
{ // gsl_sf_expint_E1
    return EVAL_RESULT_D( gsl_sf_expint_E1_e, x, "gsl_sf_expint_E1" );
} // gsl_sf_expint_E1

export function gsl_sf_expint_E1_scaled( x )
{ // gsl_sf_expint_E1_scaled
    return EVAL_RESULT_D( gsl_sf_expint_E1_scaled_e, x, "gsl_sf_expint_E1_scaled" );
} // gsl_sf_expint_E1_scaled

export function gsl_sf_expint_E2( x )
{ // gsl_sf_expint_E2
    return EVAL_RESULT_D( gsl_sf_expint_E2_e, x, "gsl_sf_expint_E2" );
} // gsl_sf_expint_E2

export function gsl_sf_expint_E2_scaled( x )
{ // gsl_sf_expint_E2_scaled
    return EVAL_RESULT_D( gsl_sf_expint_E2_scaled_e, x, "gsl_sf_expint_E2_scaled" );
} // gsl_sf_expint_E2_scaled

export function gsl_sf_expint_En( n, x )
{ // gsl_sf_expint_En
    return EVAL_RESULT_ID( gsl_sf_expint_En_e, { i: n, x: x }, "gsl_sf_expint_En" );
} // gsl_sf_expint_En;

export function gsl_sf_expint_En_scaled( n, x )
{ // gsl_sf_expint_En_scaled
    return EVAL_RESULT_ID( gsl_sf_expint_En_scaled_e, { i: n, x: x }, "gsl_sf_expint_En_scaled" );
} // gsl_sf_expint_En_scaled;

export function gsl_sf_expint_Ei( x )
{ // gsl_sf_expint_Ei
    return EVAL_RESULT_D( gsl_sf_expint_Ei_e, x, "gsl_sf_expint_Ei" );
} // gsl_sf_expint_Ei

export function gsl_sf_expint_Ei_scaled( x )
{ // gsl_sf_expint_Ei_scaled
    return EVAL_RESULT_D( gsl_sf_expint_Ei_scaled_e, x, "gsl_sf_expint_Ei_scaled" );
} // gsl_sf_expint_Ei_scaled

export function gsl_sf_expint_3( x )
{ // gsl_sf_expint_3
    return EVAL_RESULT_D( gsl_sf_expint_3_e, x, "gsl_sf_expint_3" );
} // gsl_sf_expint_3

// ----------------------------------------------------------------------------
// EOF SF-ExponentialIntegral.mjs
