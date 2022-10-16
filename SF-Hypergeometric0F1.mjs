// SF-Hypergeometric0F1.mjs
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

import { M_PI }                   from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }        from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }        from "./SF-Machine.mjs";
import { gsl_sf_lngamma_sgn_e }   from "./SF-Gamma.mjs";
import { gsl_sf_exp_mult_err_e }  from "./SF-Exponential.mjs";
import { gsl_sf_bessel_Inu_scaled_e } from "./SF-BesselInu.mjs";
import { gsl_sf_bessel_Knu_scaled_e } from "./SF-BesselKnu.mjs";
import { gsl_sf_bessel_Jnu_e }   from "./SF-BesselJnu.mjs";
import { gsl_sf_bessel_Ynu_e }   from "./SF-BesselYnu.mjs";

import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

const locEPS = 1000.0 * GSL_DBL_EPSILON;


// Evaluate bessel_I(nu, x), allowing nu < 0.
// This is fine here because we do not not allow
// nu to be a negative integer.
// x > 0.
//
function hyperg_0F1_bessel_I(nu, x)
{
    var anu = 0.0;
    var s   = 0.0;
    var ex  = 0.0;
    var I   = { val: 0.0, err: 0.0 }; // Result;
    var K   = { val: 0.0, err: 0.0 }; // Result;
    var r   = { val: 0.0, err: 0.0 }; // Result;

    if (x > GSL_LOG_DBL_MAX)
    {
        throw "SF.OverflowException";
    }
   
    if (nu < 0.0)
    {
        anu = -nu;
        s   = 2.0 / M_PI * Math.sin(anu * M_PI);
        ex  = Math.exp(x);
        I = gsl_sf_bessel_Inu_scaled_e(anu, x);
        K = gsl_sf_bessel_Knu_scaled_e(anu, x);
        r.val = ex * I.val + s * (K.val / ex);
        r.err = ex * I.err + Math.abs(s * K.err / ex);
        r.err = r.err + Math.abs(s * (K.val / ex)) * GSL_DBL_EPSILON * anu * M_PI;
    }
    else
    {
        ex = Math.exp(x);
        I = gsl_sf_bessel_Inu_scaled_e(nu, x);
        r.val = ex * I.val;
        r.err = ex * I.err + GSL_DBL_EPSILON * Math.abs(r.val);
    }

    return r;

} // hyperg_0F1_bessel_I

// ----------------------------------------------------------------------------

// Evaluate bessel_J(nu, x), allowing nu < 0.
// This is fine here because we do not not allow
// nu to be a negative integer.
// x > 0.
//
function hyperg_0F1_bessel_J(nu, x)
{
    var anu = 0.0;
    var s   = 0.0;
    var c   = 0.0;
    var J   = { val: 0.0, err: 0.0 }; // Result;
    var Y   = { val: 0.0, err: 0.0 }; // Result;
    var r   = { val: 0.0, err: 0.0 }; // Result;

  if (nu < 0.0)
  {
      anu = -nu;
      s   = Math.sin(anu * M_PI);
      c   = Math.cos(anu * M_PI);
      J = gsl_sf_bessel_Jnu_e(anu, x);
      Y = gsl_sf_bessel_Ynu_e(anu, x);
      r.val = c * J.val - s * Y.val;
      r.err = Math.abs(c * J.err) + Math.abs(s * Y.err);
      r.err = r.err + Math.abs(anu * M_PI) * GSL_DBL_EPSILON * Math.abs(J.val + Y.val);
  }
  else
  {
      r = gsl_sf_bessel_Jnu_e(nu, x);
  }

  return r;

} // hyperg_0F1_bessel_J

//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hyperg_0F1_e(c, x)
{
    var rintc  = 0.0;
    var sgn    = 0.0;
    var tl     = 0.0;
    var ln_pre_val = 0.0;
    var ln_pre_err = 0.0;
    var c_neg_integer = false;
    var Jcm1   = { val: 0.0, err: 0.0 }; // Result;
    var Icm1   = { val: 0.0, err: 0.0 }; // Result;
    var r      = { val: 0.0, err: 0.0 }; // Result;
    var lg_c   = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;

    rintc = Math.floor(c + 0.5);
    c_neg_integer = (c < 0.0 && Math.abs(c - rintc) < locEPS);
   
    if (c == 0.0 || c_neg_integer)
    {
        throw "SF.DomainException";
    }
    else if (x < 0.0)
    {
        lg_c = gsl_sf_lngamma_sgn_e(c); //, lg_c, sgn);
        sgn = lg_c.sign;
        Jcm1 = hyperg_0F1_bessel_J(c - 1.0, 2.0 * Math.sqrt(-x));
        tl = Math.log(-x) * 0.5 * (1.0 - c);
        ln_pre_val = lg_c.val + tl;
        ln_pre_err = lg_c.err + 2.0 * GSL_DBL_EPSILON * Math.abs(tl);
        r = gsl_sf_exp_mult_err_e(ln_pre_val, ln_pre_err, sgn * Jcm1.val, Jcm1.err);
    }
    else if (x == 0.0)
    {
        r.val = 1.0;
        r.err = 1.0;
    }
    else
    {
        lg_c = gsl_sf_lngamma_sgn_e(c); //, lg_c, sgn);
        sgn = lg_c.sign;
        Icm1 = hyperg_0F1_bessel_I(c - 1.0, 2.0 * Math.sqrt(x));
        tl = Math.log(x) * 0.5 * (1.0 - c);
        ln_pre_val = lg_c.val + tl;
        ln_pre_err = lg_c.err + 2.0 * GSL_DBL_EPSILON * Math.abs(tl);
        r = gsl_sf_exp_mult_err_e(ln_pre_val, ln_pre_err, sgn * Icm1.val, Icm1.err);
    }

    return r;

} // gsl_sf_hyperg_0F1_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hyperg_0F1( c, x )
{ // gsl_sf_hyperg_0F1
    return EVAL_RESULT_DD( gsl_sf_hyperg_0F1_e, { x: c, y: x }, "gsl_sf_hyperg_0F1" );
} // gsl_sf_hyperg_0F1

// ----------------------------------------------------------------------------
// EOF SF-Hypergeometric0F1.mjs
