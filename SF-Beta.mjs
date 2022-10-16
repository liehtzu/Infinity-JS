// SF-Beta.mjs
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

import { M_SQRT2 }              from "./SF-Math.mjs";
import { M_SQRTPI }             from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";
import { gsl_sf_lngamma_sgn_e } from "./SF-Gamma.mjs";
import { gsl_sf_gammastar_e }   from "./SF-Gamma.mjs";
import { gsl_sf_gamma_e }       from "./SF-Gamma.mjs";
import { gsl_sf_log_1plusx_e }  from "./SF-Logarithmic.mjs";
import { gsl_sf_exp_err_e }     from "./SF-Exponential.mjs";

import { EVAL_RESULT_DD }       from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

function isnegint(x)
{
    return (x < 0.0) && (x == Math.floor(x));
} // isnegint

// ----------------------------------------------------------------------------

export function gsl_sf_lnbeta_e(x, y)
{
    var sgn = 0.0;

    var rs = { val: 0.0, err: 0.0, sign: 0.0 };
    var r  = { val: 0.0, err: 0.0 };

    rs = gsl_sf_lnbeta_sgn_e(x, y); //, r, sgn);
    if (rs.sign == -1.0)
    {
        throw "SF.DomainException";
    }
    r.val = rs.val;
    r.err = rs.err;
    return r;

} // gsl_sf_lnbeta_e

// ----------------------------------------------------------------------------

export function gsl_sf_lnbeta_sgn_e(x, y) //, r, sgn)
{
    var max = 0.0;
    var min = 0.0;
    var rat = 0.0;

    var r = { val: 0.0, err: 0.0, sign: 0.0 };

    if (x == 0.0 || y == 0.0)
    {
        //gn = 0.0;
        throw "SF.DomainException";
    }
    else if (isnegint(x) || isnegint(y))
    {
        //sgn = 0.0;
        throw "SF.DomainException"; // not defined for negative integers
    }
   
    // See if we can handle the postive case with min/max < 0.2
   
    if (x > 0.0 && y > 0.0)
    {
        max = Math.max(x, y);
        min = Math.min(x, y);
        rat = min / max;
        
        if (rat < 0.2)
        {
            // min << max, so be careful
            // with the subtraction
            //
            var lnpre_val = 0.0;
            var lnpre_err = 0.0;
            var lnpow_val = 0.0;
            var lnpow_err = 0.0;
            var t1 = 0.0;
            var t2 = 0.0;
            var t3 = 0.0;
            var lnopr = { val: 0.0, err: 0.0 }; // Result;
            var gsx   = { val: 0.0, err: 0.0 }; // Result;
            var gsy   = { val: 0.0, err: 0.0 }; // Result;
            var gsxy  = { val: 0.0, err: 0.0 }; // Result;

            gsx = gsl_sf_gammastar_e(x);
            gsy = gsl_sf_gammastar_e(y);
            gsxy = gsl_sf_gammastar_e(x + y);
            lnopr = gsl_sf_log_1plusx_e(rat);
            lnpre_val = Math.log(gsx.val * gsy.val / gsxy.val * M_SQRT2 * M_SQRTPI);
            lnpre_err = gsx.err / gsx.val + gsy.err / gsy.val + gsxy.err / gsxy.val;
            t1 = min * Math.log(rat);
            t2 = 0.5 * Math.log(min);
            t3 = (x + y - 0.5) * lnopr.val;
            lnpow_val = t1 - t2 - t3;
            lnpow_err = GSL_DBL_EPSILON * (Math.abs(t1) + Math.abs(t2) + Math.abs(t3));
            lnpow_err = lnpow_err + Math.abs(x + y - 0.5) * lnopr.err;
            r.val = lnpre_val + lnpow_val;
            r.err = lnpre_err + lnpow_err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
            r.sign = 1.0;
            return r;
        }
    }

    // General case - Fallback
    var lgx  = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var lgy  = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    var lgxy = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
    //var sgx  = 0.0;
    //var sgy  = 0.0;
    //var sgxy = 0.0;
    var xy   = 0.0;

    xy = x + y;
    lgx = gsl_sf_lngamma_sgn_e( x); //, lgx,  sgx);
    lgy = gsl_sf_lngamma_sgn_e( y); //, lgy,  sgy);
    lgxy = gsl_sf_lngamma_sgn_e(xy); //, lgxy, sgxy);
    r.sign = lgx.sign * lgy.sign * lgxy.sign;
    r.val = lgx.val + lgy.val - lgxy.val;
    r.err = lgx.err + lgy.err + lgxy.err;
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * (Math.abs(lgx.val) + Math.abs(lgy.val) + Math.abs(lgxy.val));
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    return r;

} // gsl_sf_lnbeta_sgn_e

// ----------------------------------------------------------------------------

export function gsl_sf_beta_e(x, y)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ((x > 0.0 && y > 0.0) && x < 50.0 && y < 50.0)
    {
        // Handle the easy case
        var gx  = { val: 0.0, err: 0.0 }; // Result;
        var gy  = { val: 0.0, err: 0.0 }; // Result;
        var gxy = { val: 0.0, err: 0.0 }; // Result;

        gx  = gsl_sf_gamma_e(x);
        gy  = gsl_sf_gamma_e(y);
        gxy = gsl_sf_gamma_e(x + y);
        r.val = (gx.val * gy.val) / gxy.val;
        r.err = gx.err * Math.abs(gy.val / gxy.val);
        r.err = r.err + gy.err * Math.abs(gx.val / gxy.val);
        r.err = r.err + Math.abs((gx.val * gy.val) / (gxy.val * gxy.val)) * gxy.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);

    }
    else if (isnegint(x) || isnegint(y))
    {
        throw "SF.DomainException";
    }
    else if (isnegint(x + y))  // infinity in the denominator
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else
    {
        var lb  = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;
        //var sgn = 0.0;

        lb = gsl_sf_lnbeta_sgn_e(x, y); //, lb, sgn);
        r = gsl_sf_exp_err_e(lb.val, lb.err);
        r.val = r.val * lb.sign;
    }

    return r;

} // gsl_sf_beta_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_lnbeta( x, y )
{ // gsl_sf_lnbeta
    return EVAL_RESULT_DD( gsl_sf_lnbeta_e, { x: x, y: y }, "gsl_sf_lnbeta" );
} // gsl_sf_lnbeta

export function gsl_sf_beta( x, y )
{ // gsl_sf_beta
    return EVAL_RESULT_DD( gsl_sf_beta_e, { x: x, y: y }, "gsl_sf_beta" );
} // gsl_sf_beta

// ----------------------------------------------------------------------------
// EOF SF-Beta.mjs
