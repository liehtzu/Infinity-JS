// SF-Erfc.mjs
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

import { cheb_eval_e }           from "./SF-Chebyshev.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { M_SQRTPI }              from "./SF-Math.mjs";
import { M_SQRT2 }               from "./SF-Math.mjs";
import { GSL_ROOT6_DBL_EPSILON } from "./SF-Machine.mjs";
import { gsl_sf_exp_e }          from "./SF-Exponential.mjs";

import { EVAL_RESULT_D }         from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------
// See Hart et al, Computer Approximations, John Wiley and Sons, New York (1968)
// (This applies only to the erfc8 stuff, which is the part
//  of the original code that survives. I have replaced much of
//  the other stuff with Chebyshev fits. These are simpler and
//  more precise than the original approximations. [GJ])
// ----------------------------------------------------------------------------


//     LogRootPi : CONSTANT LONG_FLOAT = 0.57236494292470008706;

// ----------------------------------------------------------------------------

function erfc8_sum(x)
{
    const P =//: CONSTANT ARRAY(0..5) OF LONG_FLOAT =
        [
        2.97886562639399288862,
        7.409740605964741794425,
        6.1602098531096305440906,
        5.019049726784267463450058,
        1.275366644729965952479585264,
        0.5641895835477550741253201704
        ];
    const Q =//: CONSTANT ARRAY(0..6) OF LONG_FLOAT =
        [
        3.3690752069827527677,
        9.608965327192787870698,
        17.08144074746600431571095,
        12.0489519278551290360340491,
        9.396034016235054150430579648,
        2.260528520767326969591866945,
        1.0
        ];
    var num = 0.0;
    var den = 0.0;
    var i   = 0;

    // estimates erfc(x) valid for 8 < x < 100
    // This is based on index 5725 in Hart et al
   
    num = P[5];
    for (i = 4; i >= 0; i--)
    {
        num = x * num + P[i];
    }
    den = Q[6];
    for (i = 5; i >= 0; i--)
    {
        den = x * den + Q[i];
    }

    return num / den;

} // erfc8_sum

// ----------------------------------------------------------------------------

function erfc8(x)
{
    var e = 0.0;

    e = erfc8_sum(x);
    e = e * Math.exp(-x * x);
    return e;

} // erfc8

// ----------------------------------------------------------------------------

function log_erfc8(x)
{
    var e = 0.0;

    e = erfc8_sum(x);
    e = Math.log(e) - x * x;
    return e;

} // log_erfc8

// ----------------------------------------------------------------------------

// Abramowitz+Stegun, 7.1.5
function erfseries(x)
{
    var coef = 0.0;
    var e    = 0.0;
    var del  = 0.0;
    var k    = 0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    coef = x;
    e = coef;
    for (k = 1; k <= 29; k++)
    {
        coef = coef * (-x * x / k);
        del = coef / (2.0 * k + 1.0);
        e = e + del;
    }
    r.val = 2.0 / M_SQRTPI * e;
    r.err = 2.0 / M_SQRTPI * (Math.abs(del) + Number.EPSILON);

    return r;

} // erfseries

// ----------------------------------------------------------------------------

// Chebyshev fit for erfc((t+1)/2), -1 < t < 1

const erfc_xlt1_data =//: CONSTANT Series(0..19) = --[20]
    [
    1.06073416421769980345174155056,
   -0.42582445804381043569204735291,
    0.04955262679620434040357683080,
    0.00449293488768382749558001242,
   -0.00129194104658496953494224761,
   -0.00001836389292149396270416979,
    0.00002211114704099526291538556,
   -5.23337485234257134673693179020e-7,
   -2.78184788833537885382530989578e-7,
    1.41158092748813114560316684249e-8,
    2.72571296330561699984539141865e-9,
   -2.06343904872070629406401492476e-10,
   -2.14273991996785367924201401812e-11,
    2.22990255539358204580285098119e-12,
    1.36250074650698280575807934155e-13,
   -1.95144010922293091898995913038e-14,
   -6.85627169231704599442806370690e-16,
    1.44506492869699938239521607493e-16,
    2.45935306460536488037576200030e-18,
   -9.29599561220523396007359328540e-19
    ];
const erfc_xlt1_cs = { length: 19, c: erfc_xlt1_data, order: 19, a: -1.0, b: 1.0, order_sp: 12 };

// ----------------------------------------------------------------------------

// Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1

const erfc_x15_data =//: CONSTANT Series(0..24) = --[25]
    [
    0.44045832024338111077637466616,
   -0.143958836762168335790826895326,
    0.044786499817939267247056666937,
   -0.013343124200271211203618353102,
    0.003824682739750469767692372556,
   -0.001058699227195126547306482530,
    0.000283859419210073742736310108,
   -0.000073906170662206760483959432,
    0.000018725312521489179015872934,
   -4.62530981164919445131297264430e-6,
    1.11558657244432857487884006422e-6,
   -2.63098662650834130067808832725e-7,
    6.07462122724551777372119408710e-8,
   -1.37460865539865444777251011793e-8,
    3.05157051905475145520096717210e-9,
   -6.65174789720310713757307724790e-10,
    1.42483346273207784489792999706e-10,
   -3.00141127395323902092018744545e-11,
    6.22171792645348091472914001250e-12,
   -1.26994639225668496876152836555e-12,
    2.55385883033257575402681845385e-13,
   -5.06258237507038698392265499770e-14,
    9.89705409478327321641264227110e-15,
   -1.90685978789192181051961024995e-15,
    3.50826648032737849245113757340e-16
    ];
const erfc_x15_cs = { length: 24, c: erfc_x15_data, order: 24, a: -1.0, b: 1.0, oredr_sp: 16 };

// ----------------------------------------------------------------------------

// Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1

const erfc_x510_data =//: CONSTANT Series(0..19) = --[20]
    [
    1.11684990123545698684297865808,
    0.003736240359381998520654927536,
   -0.000916623948045470238763619870,
    0.000199094325044940833965078819,
   -0.000040276384918650072591781859,
    7.76515264697061049477127605790e-6,
   -1.44464794206689070402099225301e-6,
    2.61311930343463958393485241947e-7,
   -4.61833026634844152345304095560e-8,
    8.00253111512943601598732144340e-9,
   -1.36291114862793031395712122089e-9,
    2.28570483090160869607683087722e-10,
   -3.78022521563251805044056974560e-11,
    6.17253683874528285729910462130e-12,
   -9.96019290955316888445830597430e-13,
    1.58953143706980770269506726000e-13,
   -2.51045971047162509999527428316e-14,
    3.92607828989125810013581287560e-15,
   -6.07970619384160374392535453420e-16,
    9.12600607264794717315507477670e-17
    ];
const erfc_x510_cs =  { length: 19, c: erfc_x510_data, order: 19, a: -1.0, b: 1.0, order_sp: 12 };

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_erfc_e(x)
{
    const ax   = Math.abs(x);
    var e_val  = 0.0;
    var e_err  = 0.0;
    var t      = 0.0;
    var ex2    = 0.0;
    var exterm = 0.0;

    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (ax <= 1.0)
    {
        t = 2.0 * ax - 1.0;
        c = cheb_eval_e(erfc_xlt1_cs, t);
        e_val = c.val;
        e_err = c.err;
    }
    else if (ax <= 5.0)
    {
        ex2 = Math.exp(-x * x);
        t = 0.5 * (ax - 3.0);
        c = cheb_eval_e(erfc_x15_cs, t);
        e_val = ex2 * c.val;
        e_err = ex2 * (c.err + 2.0 * Math.abs(x) * Number.EPSILON);
    }
    else if (ax < 10.0)
    {
        exterm = Math.exp(-x * x) / ax;
        t = (2.0 * ax - 15.0) / 5.0;
        c = cheb_eval_e(erfc_x510_cs, t);
        e_val = exterm * c.val;
        e_err = exterm * (c.err + 2.0 * Math.abs(x) * Number.EPSILON + Number.EPSILON);
    }
    else
    {
        e_val = erfc8(ax);
        e_err = (x * x + 1.0) * Number.EPSILON * Math.abs(e_val);
    }
   
    if (x < 0.0)
    {
        r.val = 2.0 - e_val;
        r.err = e_err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else
    {
        r.val = e_val;
        r.err = e_err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }
   
    return r;

} // gsl_sf_erfc_e

// ----------------------------------------------------------------------------

export function gsl_sf_log_erfc_e(x)
{
    var result_erfc = { val: 0.0, err: 0.0 }; // Result;
    var r           = { val: 0.0, err: 0.0 }; // Result;

    if (x * x < 10.0 * GSL_ROOT6_DBL_EPSILON)
    {
        const y = x / M_SQRTPI;
        // series for -1/2 Log[Erfc[Sqrt[Pi] y]] */
        const c3 = (4.0 - M_PI) / 3.0;
        const c4 = 2.0 * (1.0 - M_PI / 3.0);
        const c5 = -0.001829764677455021;  // (96.0 - 40.0*M_PI + 3.0*M_PI*M_PI)/30.0
        const c6 =  0.02629651521057465;   // 2.0*(120.0 - 60.0*M_PI + 7.0*M_PI*M_PI)/45.0
        const c7 = -0.01621575378835404;
        const c8 =  0.00125993961762116;
        const c9 =  0.00556964649138;
        const c10 = -0.0045563339802;
        const c11 =  0.0009461589032;
        const c12 =  0.0013200243174;
        const c13 = -0.00142906;
        const c14 =  0.00048204;
        var series = 0.0;

        series = c8 + y * (c9 + y * (c10 + y * (c11 + y * (c12 + y * (c13 + c14 * y)))));
        series = y * (1.0 + y * (1.0 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * (c7 + y * series)))))));
        r.val = -2.0 * series;
        r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
    //
    //don't like use of log1p(); added above series stuff for small x instead, should be ok [GJ]
    //else if (Math.abs(x) < 1.0) {
    //  gsl_sf_result result_erf;
    //  gsl_sf_erf_e(x, &result_erf);
    //  r.val  = log1p(-result_erf.val);
    //  r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
    //  code = GSL_SUCCESS;
    //}
    //
    }
    else if (x > 8.0)
    {
        r.val = log_erfc8(x);
        r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
    }
    else
    {
        result_erfc = gsl_sf_erfc_e(x);
        r.val = Math.log(result_erfc.val);
        r.err = Math.abs(result_erfc.err / result_erfc.val);
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }

    return r;

} //gsl_sf_log_erfc_e

// ----------------------------------------------------------------------------

export function gsl_sf_erf_e(x)
{
    var result_erfc = { val: 0.0, err: 0.0 }; // Result;
    var r           = { val: 0.0, err: 0.0 }; // Result;

    if (Math.abs(x) < 1.0)
    {
        r = erfseries(x);
    }
    else
    {
        result_erfc = gsl_sf_erfc_e(x);
        r.val = 1.0 - result_erfc.val;
        r.err = result_erfc.err;
        r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_erf_e

// ----------------------------------------------------------------------------

export function gsl_sf_erf_Z_e(x)
{
    var ex2 = 0.0;
    var r   = { val: 0.0, err: 0.0 }; // Result;

    ex2 = Math.exp(-x * x / 2.0);
    r.val = ex2 / (M_SQRT2 * M_SQRTPI);
    r.err = Math.abs(x * r.val) * Number.EPSILON;
    r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);
    //CHECK_UNDERFLOW(result);

    return r;

} // gsl_sf_erf_Z_e

// ----------------------------------------------------------------------------

export function gsl_sf_erf_Q_e(x)
{
    var result_erfc = { val: 0.0, err: 0.0 }; // Result;
    var r           = { val: 0.0, err: 0.0 }; // Result;

    result_erfc = gsl_sf_erfc_e(x / M_SQRT2);
    r.val = 0.5 * result_erfc.val;
    r.err = 0.5 * result_erfc.err;
    r.err = r.err + 2.0 * Number.EPSILON * Math.abs(r.val);

    return r;

} // gsl_sf_erf_Q_e

// ----------------------------------------------------------------------------

export function gsl_sf_hazard_e(x)
{
    var ix2   = 0.0;
    var corrB = 0.0;
    var corrM = 0.0;
    var corrT = 0.0;
    var lnc   = 0.0;
    var arg   = 0.0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < 25.0)
    {
        c = gsl_sf_log_erfc_e(x / M_SQRT2);
        lnc = -0.22579135264472743236; // ln(sqrt(2/pi))
        arg = lnc - 0.5 * x * x - c.val;
        r = gsl_sf_exp_e(arg);
        r.err = r.err + 3.0 * (1.0 + Math.abs(x)) * Number.EPSILON * Math.abs(r.val);
        r.err = r.err + Math.abs(c.err * r.val);
    }
    else
    {
        ix2 = 1.0 / (x * x);
        corrB = 1.0 - 9.0 * ix2 * (1.0 - 11.0 * ix2);
        corrM = 1.0 - 5.0 * ix2 * (1.0 - 7.0 * ix2 * corrB);
        corrT = 1.0 - ix2 * (1.0 - 3.0 * ix2 * corrM);
        r.val = x / corrT;
        r.err = 2.0 * Number.EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_hazard_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_erfc( x )
{ // gsl_sf_erfc
    return EVAL_RESULT_D( gsl_sf_erfc_e, x, "gsl_sf_erfc_e" );
} // gsl_sf_erfc;

export function gsl_sf_log_erfc( x )
{ // gsl_sf_log_erfc
    return EVAL_RESULT_D( gsl_sf_log_erfc_e, x, "gsl_sf_log_erfc_e" );
} // gsl_sf_log_erfc;

export function gsl_sf_erf( x )
{ // gsl_sf_erf
    return EVAL_RESULT_D( gsl_sf_erf_e, x, "gsl_sf_erf_e" );
} // gsl_sf_erf;

export function gsl_sf_erf_Z( x )
{ // gsl_sf_erf_Z
    return EVAL_RESULT_D( gsl_sf_erf_Z_e, x, "gsl_sf_erf_Z_e" );
} // gsl_sf_erf_Z;

export function gsl_sf_erf_Q( x )
{ // gsl_sf_erf_Q
    return EVAL_RESULT_D( gsl_sf_erf_Q_e, x, "gsl_sf_erf_Q_e" );
} // gsl_sf_erf_Q;

export function gsl_sf_hazard( x )
{ // gsl_sf_hazard
    return EVAL_RESULT_D( gsl_sf_hazard_e, x, "gsl_sf_hazard_e" );
} // gsl_sf_hazard;

// ----------------------------------------------------------------------------
// EOF SF-Erfc.mjs
