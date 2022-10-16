// SF-Logarithmic.mjs
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

import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_ROOT5_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_ROOT6_DBL_EPSILON } from "./SF-Machine.mjs";
import { cheb_eval_e }           from "./SF-Chebyshev.mjs";

import { EVAL_RESULT_D }         from "./SF-Evaluate.mjs";

// Chebyshev expansion for log(1 + x(t))/x(t)
//
// x(t) = (4t-1)/(2(4-t))
// t(x) = (8x+1)/(2(x+2))
// -1/2 < x < 1/2
// -1 < t < 1

const lopx_data =//: CONSTANT Series(0..20) = --[21]
    [
    2.16647910664395270521272590407,
   -0.28565398551049742084877469679,
    0.01517767255690553732382488171,
   -0.00200215904941415466274422081,
    0.00019211375164056698287947962,
   -0.00002553258886105542567601400,
    2.9004512660400621301999384544e-06,
   -3.8873813517057343800270917900e-07,
    4.7743678729400456026672697926e-08,
   -6.4501969776090319441714445454e-09,
    8.2751976628812389601561347296e-10,
   -1.1260499376492049411710290413e-10,
    1.4844576692270934446023686322e-11,
   -2.0328515972462118942821556033e-12,
    2.7291231220549214896095654769e-13,
   -3.7581977830387938294437434651e-14,
    5.1107345870861673561462339876e-15,
   -7.0722150011433276578323272272e-16,
    9.7089758328248469219003866867e-17,
   -1.3492637457521938883731579510e-17,
    1.8657327910677296608121390705e-18
    ];
const lopx_cs = { length: 20, c: lopx_data, order: 20, a: -1.0, b: 1.0, order_sp: 10 };

// Chebyshev expansion for (log(1 + x(t)) - x(t))/x(t)^2
//
// x(t) = (4t-1)/(2(4-t))
// t(x) = (8x+1)/(2(x+2))
// -1/2 < x < 1/2
// -1 < t < 1

const lopxmx_data =//: CONSTANT Series(0..19) = --[20]
    [
   -1.12100231323744103373737274541,
    0.19553462773379386241549597019,
   -0.01467470453808083971825344956,
    0.00166678250474365477643629067,
   -0.00018543356147700369785746902,
    0.00002280154021771635036301071,
   -2.8031253116633521699214134172e-06,
    3.5936568872522162983669541401e-07,
   -4.6241857041062060284381167925e-08,
    6.0822637459403991012451054971e-09,
   -8.0339824424815790302621320732e-10,
    1.0751718277499375044851551587e-10,
   -1.4445310914224613448759230882e-11,
    1.9573912180610336168921438426e-12,
   -2.6614436796793061741564104510e-13,
    3.6402634315269586532158344584e-14,
   -4.9937495922755006545809120531e-15,
    6.8802890218846809524646902703e-16,
   -9.5034129794804273611403251480e-17,
    1.3170135013050997157326965813e-17
    ];
const lopxmx_cs = { length: 19, c: lopxmx_data, order: 19, a: -1.0, b: 1.0, order_sp: 9 };

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_log_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= 0.0)
    {
        throw "SF.DomainException";
    }
    else
    {
        r.val = Math.log(x);
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_log_e

// ----------------------------------------------------------------------------

export function gsl_sf_log_abs_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x == 0.0)
    {
        throw "SF.DomainException";
    }
    else
    {
        r.val = Math.log(Math.abs(x));
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_log_abs_e

// ----------------------------------------------------------------------------

export function gsl_sf_log_1plusx_e(x)
{
    const c1 = -0.5;
    const c2 =  1.0/3.0;
    const c3 = -1.0/4.0;
    const c4 =  1.0/5.0;
    const c5 = -1.0/6.0;
    const c6 =  1.0/7.0;
    const c7 = -1.0/8.0;
    const c8 =  1.0/9.0;
    const c9 = -1.0/10.0;
    var t = 0.0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= -1.0)
    {
        throw "SF.DomainException";
    }
    else if (Math.abs(x) < GSL_ROOT6_DBL_EPSILON)
    {
        t =  c5 + x * (c6 + x * (c7 + x * (c8 + x * c9)));
        r.val = x * (1.0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * t)))));
        r.err = GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (Math.abs(x) < 0.5)
    {
        t = 0.5 * (8.0 * x + 1.0) / (x + 2.0);
        c = cheb_eval_e(lopx_cs, t);
        r.val = x * c.val;
        r.err = Math.abs(x * c.err);
    }
    else
    {
        r.val = Math.log(1.0 + x);
        r.err = GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // gsl_sf_log_1plusx_e

// ----------------------------------------------------------------------------

export function gsl_sf_log_1plusx_mx_e(x)
{
    const c1 = -0.5;
    const c2 =  1.0 / 3.0;
    const c3 = -1.0 / 4.0;
    const c4 =  1.0 / 5.0;
    const c5 = -1.0 / 6.0;
    const c6 =  1.0 / 7.0;
    const c7 = -1.0 / 8.0;
    const c8 =  1.0 / 9.0;
    const c9 = -1.0 / 10.0;
    var t = 0.0;
    var lterm = 0.0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x <= -1.0)
    {
        throw "SF.DomainException";
    }
    else if (Math.abs(x) < GSL_ROOT5_DBL_EPSILON)
    {
        t =  c5 + x * (c6 + x * (c7 + x * (c8 + x * c9)));
        r.val = x * x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * t))));
        r.err = GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (Math.abs(x) < 0.5)
    {
        t = 0.5 * (8.0 * x + 1.0) / (x + 2.0);
        c = cheb_eval_e(lopxmx_cs, t);
        r.val = x * x * c.val;
        r.err = x * x * c.err;
    }
    else
    {
        lterm = Math.log(1.0 + x);
        r.val = lterm - x;
        r.err = GSL_DBL_EPSILON * (Math.abs(lterm) + Math.abs(x));
    }

    return r;

} // gsl_sf_log_1plusx_mx_e

// ----------------------------------------------------------------------------

export function gsl_sf_complex_log_e( zr, zi, lnr, theta )
{
    var ax  = 0.0;
    var ay  = 0.0;
    var min = 0.0;
    var max = 0.0;

    if (zr != 0.0 || zi != 0.0)
    {
        ax = Math.abs( zr );
        ay = Math.abs( zi );
        min = Math.min( ax, ay );
        max = Math.max( ax, ay );
        lnr.val = Math.log( max ) + 0.5 * Math.log( 1.0 + (min / max) * (min / max) );
        lnr.err = 2.0 * GSL_DBL_EPSILON * Math.abs( lnr.val );
        theta.val = Math.atan2( zi, zr );
        theta.err = GSL_DBL_EPSILON * Math.abs( lnr.val );
    }
    else
    {
        throw "SF.DomainException";
    }

} // gsl_sf_complex_log_e

//-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_log( x )
{ // gsl_sf_log
    return EVAL_RESULT_D( gsl_sf_log_e, x, "gsl_sf_log" );
} // gsl_sf_log

export function gsl_sf_log_abs( x )
{ // gsl_sf_log_abs
    return EVAL_RESULT_D( gsl_sf_log_abs_e, x, "gsl_sf_log_abs" );
} // gsl_sf_log_abs

export function gsl_sf_log_1plusx( x )
{ // gsl_sf_log_1plusx
    return EVAL_RESULT_D( gsl_sf_log_1plusx_e, x, "gsl_sf_log_1plusx" );
} // gsl_sf_log_1plusx

export function gsl_sf_log_1plusx_mx( x )
{ // gsl_sf_log_1plusx_mx
    return EVAL_RESULT_D( gsl_sf_log_1plusx_mx_e, x, "gsl_sf_log_1plusx_mx" );
} // gsl_sf_log_1plusx_mx

// ----------------------------------------------------------------------------
// EOF SF-Logarithmic.mjs
