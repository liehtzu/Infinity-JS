// SF-Synchrotron.mjs
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

import { GSL_DBL_EPSILON }        from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }   from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MIN }        from "./SF-Machine.mjs";
import { M_PI }                   from "./SF-Math.mjs";
import { M_SQRT2 }                from "./SF-Math.mjs";
import { M_SQRT3 }                from "./SF-Math.mjs";
import { cheb_eval_e }            from "./SF-Chebyshev.mjs";

import { EVAL_RESULT_D }          from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const synchrotron1_data =
    [
   30.364682982501076273,
   17.079395277408394574,
    4.560132133545072889,
    0.549281246730419979,
    0.372976075069301172e-01,
    0.161362430201041242e-02,
    0.481916772120371e-04,
    0.10512425288938e-05,
    0.174638504670e-07,
    0.22815486544e-09,
    0.240443082e-11,
    0.2086588e-13,
    0.15167e-15
    ];
const synchrotron1_cs = { length: 12, c: synchrotron1_data, order: 12, a: -1.0, b: 1.0, order_sp: 9 };

const synchrotron2_data =
    [
    0.4490721623532660844,
    0.898353677994187218e-01,
    0.81044573772151290e-02,
    0.4261716991089162e-03,
    0.147609631270746e-04,
    0.3628633615300e-06,
    0.66634807498e-08,
    0.949077166e-10,
    0.1079125e-11,
    0.10022e-13,
    0.77e-16,
    0.5e-18
    ];
const synchrotron2_cs = { length: 11, c: synchrotron2_data, order: 11, a: -1.0, b: 1.0, order_sp: 7 };

const synchrotron1a_data =
    [
    2.1329305161355000985,
    0.741352864954200240e-01,
    0.86968099909964198e-02,
    0.11703826248775692e-02,
    0.1645105798619192e-03,
    0.240201021420640e-04,
    0.35827756389389e-05,
    0.5447747626984e-06,
    0.838802856196e-07,
    0.13069882684e-07,
    0.2053099071e-08,
    0.325187537e-09,
    0.517914041e-10,
    0.83002988e-11,
    0.13352728e-11,
    0.2159150e-12,
    0.349967e-13,
    0.56994e-14,
    0.9291e-15,
    0.152e-15,
    0.249e-16,
    0.41e-17,
    0.7e-18
    ];
const synchrotron1a_cs = { length: 22, c: synchrotron1a_data, order: 22, a: -1.0, b: 1.0, order_sp: 11 };

const synchrotron21_data =
    [
    38.617839923843085480,
    23.037715594963734597,
    5.3802499868335705968,
    0.6156793806995710776,
    0.406688004668895584e-01,
    0.17296274552648414e-02,
    0.51061258836577e-04,
    0.110459595022e-05,
    0.18235530206e-07,
    0.2370769803e-09,
    0.24887296e-11,
    0.21529e-13,
    0.156e-15
    ];
const synchrotron21_cs = { length: 12, c: synchrotron21_data, order: 12, a: -1.0, b: 1.0, order_sp: 9 };

const synchrotron22_data =
    [
    7.9063148270660804288,
    3.1353463612853425684,
    0.4854879477453714538,
    0.394816675827237234e-01,
    0.19661622334808802e-02,
    0.659078932293042e-04,
    0.15857561349856e-05,
    0.286865301123e-07,
    0.4041202360e-09,
    0.45568444e-11,
    0.420459e-13,
    0.3232e-15,
    0.21e-17
    ];
const synchrotron22_cs = { length: 12, c: synchrotron22_data, order: 12, a: -1.0, b:1.0, order_sp: 8 };

const synchrotron2a_data =
    [
    2.020337094170713600,
    0.10956237121807404e-01,
    0.8542384730114676e-03,
    0.723430242132822e-04,
    0.63124427962699e-05,
    0.5648193141174e-06,
    0.512832480138e-07,
    0.47196532914e-08,
    0.4380744214e-09,
    0.410268149e-10,
    0.38623072e-11,
    0.3661323e-12,
    0.348023e-13,
    0.33301e-14,
    0.319e-15,
    0.307e-16,
    0.3e-17
    ];
const synchrotron2a_cs = { length: 16, c: synchrotron2a_data, order: 16, a: -1.0, b: 1.0, order_sp: 8 };


//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_synchrotron_1_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0) 
    {
        throw "SF.DomainException";
    }
    else if (x < 2.0 * M_SQRT2 * GSL_SQRT_DBL_EPSILON) 
    {
        // BJG: added first order correction term.  The taylor series
        // is  S1(x) = ((4pi)/(sqrt(3)gamma(1/3))) * (x/2)^(1/3) 
        // * (1 - (gamma(1/3)/2)*(x/2)^2/3 + (3/4) * (x/2)^2 ....)
        const z  = Math.pow(x, (1.0 / 3.0));
        const cf = 1.0 - 8.43812762813205e-01 * z * z;

        r.val = 2.14952824153447863671 * z * cf;
        r.err = GSL_DBL_EPSILON * r.val;
        return r;
    }
    else if (x <= 4.0)
    {
        const c0   = M_PI / M_SQRT3;
        const px   = Math.pow(x, (1.0 / 3.0));
        const px11 = Math.pow(px, 11); //gsl_sf_pow_int(px,11);
        const t    = x * x / 8.0 - 1.0;
        let result_c1 = { val: 0.0, err: 0.0 }; // Result;
        let result_c2 = { val: 0.0, err: 0.0 }; // Result;

        result_c1 = cheb_eval_e(synchrotron1_cs, t);
        result_c2 = cheb_eval_e(synchrotron2_cs, t);
        r.val = px * result_c1.val - px11 * result_c2.val - c0 * x;
        r.err = px * result_c1.err + px11 * result_c2.err + c0 * x * GSL_DBL_EPSILON;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }
    else if (x < -8.0 * GSL_LOG_DBL_MIN / 7.0)
    {
        const c0 = 0.2257913526447274323630976; // log(sqrt(pi/2))
        const t  = (12.0 - x) / (x + 4.0);
        let result_c1 = { val: 0.0, err: 0.0 }; // Result;

        result_c1 = cheb_eval_e(synchrotron1a_cs, t);
        r.val = Math.sqrt(x) * result_c1.val * Math.exp(c0 - x);
        r.err = 2.0 * GSL_DBL_EPSILON * r.val * (Math.abs(c0 - x) + 1.0);
        return r;
    }
    else
    {
        throw "SF.UnderflowException";
    }

} // gsl_sf_synchrotron_1_e

// ----------------------------------------------------------------------------

export function gsl_sf_synchrotron_2_e(x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0) 
    {
        throw "SF.DomainException";
    }
    else if (x < 2.0 * M_SQRT2 * GSL_SQRT_DBL_EPSILON) 
    {
        // BJG: added first order correction term.  The taylor series
        // is  S2(x) = ((2pi)/(sqrt(3)*gamma(1/3))) * (x/2)^(1/3) 
        // * (1 - (gamma(1/3)/gamma(4/3))*(x/2)^(4/3) + (gamma(1/3)/gamma(4/3))*(x/2)^2...)
        const z  = Math.pow(x, (1.0 / 3.0));
        const cf = 1.0 - 1.17767156510235e+00 * z * x;

        r.val = 1.07476412076723931836 * z * cf ;
        r.err = 2.0 * GSL_DBL_EPSILON * r.val;
        return r;
    }
    else if (x <= 4.0) 
    {
        const px  = Math.pow(x, (1.0 / 3.0));
        const px5 = Math.pow(px, 5); //gsl_sf_pow_int(px,5);
        const t   = x * x / 8.0 - 1.0;
        let cheb1 = { val: 0.0, err: 0.0 }; // Result;
        let cheb2 = { val: 0.0, err: 0.0 }; // Result;

        cheb1 = cheb_eval_e(synchrotron21_cs, t);
        cheb2 = cheb_eval_e(synchrotron22_cs, t);
        r.val = px * cheb1.val - px5 * cheb2.val;
        r.err = px * cheb1.err + px5 * cheb2.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        return r;
    }
    else if (x < -8.0*GSL_LOG_DBL_MIN / 7.0)
    { 
        const c0 = 0.22579135264472743236;   // log(sqrt(pi/2))
        const t  = (10.0 - x) / (x + 2.0);
        let cheb1 = { val: 0.0, err: 0.0 }; // Result;

        cheb1 = cheb_eval_e(synchrotron2a_cs, t);
        r.val = Math.sqrt(x) * Math.exp(c0 - x) * cheb1.val;
        r.err = GSL_DBL_EPSILON * r.val * (Math.abs(c0 - x) + 1.0);
        return r;
    }
    else
    {
        throw "SF.UnderflowException";
    }

} // gsl_sf_synchrotron_2_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_synchrotron_1( x )
{ // gsl_sf_synchrotron_1
    return EVAL_RESULT_D( gsl_sf_synchrotron_1_e, x, "gsl_sf_synchrotron_1" );
} // gsl_sf_synchrotron_1

export function gsl_sf_synchrotron_2( x )
{ // gsl_sf_synchrotron_2
    return EVAL_RESULT_D( gsl_sf_synchrotron_2_e, x, "gsl_sf_synchrotron_2" );
} // gsl_sf_synchrotron_2

// ----------------------------------------------------------------------------
// EOF SF-Synchrotron.mjs
