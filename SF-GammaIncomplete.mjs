// SF-GammaIncomplete.mjs
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

import { M_PI } from "./SF-Math.mjs";
import { M_SQRT2 } from "./SF-Math.mjs";
import { M_SQRTPI } from "./SF-Math.mjs";
import { M_EULER } from "./SF-Math.mjs";
import { GSL_IS_ODD } from "./SF-Math.mjs";
import { GSL_SIGN } from "./SF-Math.mjs";
import { cheb_eval_e }      from "./SF-Chebyshev.mjs";
import { GSL_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_ROOT5_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MIN }  from "./SF-Machine.mjs";
import { gsl_sf_expint_E1_e }from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_gamma_e } from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_e } from "./SF-Gamma.mjs";
import { gsl_sf_gammastar_e } from "./SF-Gamma.mjs";
import { gsl_sf_log_1plusx_mx_e } from "./SF-Logarithmic.mjs";
import { gsl_sf_exprel_n_CF_e } from "./SF-Exponential.mjs";
import { gsl_sf_exp_err_e } from "./SF-Exponential.mjs";
import { gsl_sf_erfc_e } from "./SF-Erfc.mjs";

// ************************************************************************* --
// The dominant part,
// D(a,x) = x^a e^(-x) / Gamma(a+1)
//
//static
//int
//gamma_inc_D(const double a, const double x, gsl_sf_result * result)
//{
//  if(a < 10.0) {
//    double lnr;
//    gsl_sf_result lg;
//    gsl_sf_lngamma_e(a+1.0, &lg);
//    lnr = a * log(x) - x - lg.val;
//    result->val = exp(lnr);
//    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lnr) + 1.0) * fabs(result->val);
//    return GSL_SUCCESS;
//  }
//  else {
//    gsl_sf_result gstar;
//    gsl_sf_result ln_term;
//    double term1;
//    if (x < 0.5*a) {
//      double u = x/a;   
//      double ln_u = log(u);
//      ln_term.val = ln_u - u + 1.0;
//      ln_term.err = (fabs(ln_u) + fabs(u) + 1.0) * GSL_DBL_EPSILON;
//    } else {
//      double mu = (x-a)/a;
//      gsl_sf_log_1plusx_mx_e(mu, &ln_term);  /* log(1+mu) - mu */
//      /* Propagate cancellation error from x-a, since the absolute
//         error of mu=x-a is DBL_EPSILON */
//      ln_term.err += GSL_DBL_EPSILON * fabs(mu);
//    };
//    gsl_sf_gammastar_e(a, &gstar);
//    term1 = exp(a*ln_term.val)/sqrt(2.0*M_PI*a);
//    result->val  = term1/gstar.val;
//    result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(a*ln_term.val) + 1.0) * fabs(result->val);
//    /* Include propagated error from log term */
//    result->err += fabs(a) * ln_term.err * fabs(result->val);
//    result->err += gstar.err/fabs(gstar.val) * fabs(result->val);
//    return GSL_SUCCESS;
//  }
//
//}

// ----------------------------------------------------------------------------

// Continued fraction which occurs in evaluation
// of Q(a,x) or Gamma(a,x).
//
//              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
//   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
//             1 +   1 +     1 +   1 +      1 +   1 +
//
// Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
//
// Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
// See gamma_inc_Q_CF() below.
//
//
//static int
//gamma_inc_F_CF(const double a, const double x, gsl_sf_result * result)
//{
//  const int    nmax  =  5000;
//  const double small =  gsl_pow_3 (GSL_DBL_EPSILON);
//
//  double hn = 1.0;           /* convergent */
//  double Cn = 1.0 / small;
//  double Dn = 1.0;
//  int n;
//
//  /* n == 1 has a_1, b_1, b_0 independent of a,x,
//     so that has been done by hand                */
//  for ( n = 2 ; n < nmax ; n++ )
//  {
//    double an;
//    double delta;
//
//    if(GSL_IS_ODD(n))
//      an = 0.5*(n-1)/x;
//    else
//      an = (0.5*n-a)/x;
//
//    Dn = 1.0 + an * Dn;
//    if ( fabs(Dn) < small )
//      Dn = small;
//    Cn = 1.0 + an/Cn;
//    if ( fabs(Cn) < small )
//      Cn = small;
//    Dn = 1.0 / Dn;
//    delta = Cn * Dn;
//    hn *= delta;
//    if(fabs(delta-1.0) < GSL_DBL_EPSILON) break;
//  }
//
//  result->val = hn;
//  result->err = 2.0*GSL_DBL_EPSILON * fabs(hn);
//  result->err += GSL_DBL_EPSILON * (2.0 + 0.5*n) * fabs(result->val);
//
//  if(n == nmax)
//    GSL_ERROR ("error in CF for F(a,x)", GSL_EMAXITER);
//  else
//    return GSL_SUCCESS;
//}

// ----------------------------------------------------------------------------

//static int
//gamma_inc_CF(double a, double x, gsl_sf_result * result)
//{
//  gsl_sf_result F;
//  gsl_sf_result pre;
//  const double am1lgx = (a-1.0)*log(x);
//  const int stat_F = gamma_inc_F_CF(a, x, &F);
//  const int stat_E = gsl_sf_exp_err_e(am1lgx - x, GSL_DBL_EPSILON*fabs(am1lgx), &pre);
//
//  result->val = F.val * pre.val;
//  result->err = fabs(F.err * pre.val) + fabs(F.val * pre.err);
//  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
//
//  return GSL_ERROR_SELECT_2(stat_F, stat_E);
//}

// ************************************************************************* --

// The dominant part,
// D(a,x) := x^a e^(-x) / Gamma(a+1)
//
function gamma_inc_D(a, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (a < 10.0)
    {
        var lnr = 0.0;
        var lg  = { val: 0.0, err: 0.0 }; // Result;

        lg = gsl_sf_lngamma_e(a + 1.0);
        lnr = a * Math.log(x) - x - lg.val;
        r.val = Math.exp(lnr);
        r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs(lnr) + 1.0) * Math.abs(r.val);
    }
    else
    {
        var gstar   = { val: 0.0, err: 0.0 }; // Result;
        var ln_term = { val: 0.0, err: 0.0 }; // Result;
        var term1   = 0.0;
        var u       = 0.0;
        var ln_u    = 0.0;
        var mu      = 0.0;

        if (x < 0.5 * a)
        {
            u = x / a;   
            ln_u = Math.log(u);
            ln_term.val = ln_u - u + 1.0;
            ln_term.err = (Math.abs(ln_u) + Math.abs(u) + 1.0) * GSL_DBL_EPSILON;
        }
        else
        {
            mu = (x - a) / a;
            ln_term = gsl_sf_log_1plusx_mx_e(mu);  // log(1+mu) - mu
            // Propagate cancellation error from x-a, since the absolute
            // error of mu=x-a is DBL_EPSILON
            ln_term.err = ln_term.err + GSL_DBL_EPSILON * Math.abs(mu);
        }
        gstar = gsl_sf_gammastar_e(a);
        term1 = Math.exp(a * ln_term.val) / Math.sqrt(2.0 * M_PI * a);
        r.val = term1 / gstar.val;
        r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs(a*ln_term.val) + 1.0) * Math.abs(r.val);
        // Include propagated error from log term
        r.err = r.err + Math.abs(a) * ln_term.err * Math.abs(r.val);
        r.err = r.err + gstar.err / Math.abs(gstar.val) * Math.abs(r.val);
    }

    return r;

} // gamma_inc_D

// ----------------------------------------------------------------------------

// P series representation.
//
function gamma_inc_P_series(a, x)
{
    const nmax = 10000;
    var n    = 0;
    var nlow = 0;
    var sum       = 0.0;
    var term      = 0.0;
    var remainder = 0.0;
    var tnp1      = 0.0;
    var c = { val: 0.0, err: 0.0 }; // Result;
    var D = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    D = gamma_inc_D(a, x);

    // Approximating the terms of the series using Stirling's
    // approximation gives t_n = (x/a)^n * exp(-n(n+1)/(2a)), so the
    // convergence condition is n^2 / (2a) + (1-(x/a) + (1/2a)) n >>
    // -log(GSL_DBL_EPS) if we want t_n < O(1e-16) t_0. The condition
    // below detects cases where the minimum value of n is > 5000 */

    if (x > 0.995 * a && a > 1.0e5) // Difficult case: try continued fraction
    {
        c = gsl_sf_exprel_n_CF_e(a, x);
        r.val = D.val * c.val;
        r.err = Math.abs(D.val * c.err) + Math.abs(D.err * c.val);
        return r;
    }

    // Series would require excessive number of terms

    if (x > (a + (nmax)))
    {
        throw "SF.MaxIterationsException"; // WITH "gamma_inc_P_series x>>a exceeds range";
    }

    // Normal case: sum the series

    sum  = 1.0;
    term = 1.0;

    // Handle lower part of the series where t_n is increasing, |x| > a+n

    if (x > a)
    {
        nlow = Math.trunc(x - a);
    }
    else
    {
        nlow = 0;
    }

    n = 1;
    while (n < nlow)
    {
        term = term * x / (a + (n));
        sum = sum + term;
        n = n + 1;
    }

    // Handle upper part of the series where t_n is decreasing, |x| < a+n

    // n = previous n
    while (n < nmax)
    {
        term = term * x / (a + (n));
        sum = sum + term;
        if (Math.abs(term / sum) < GSL_DBL_EPSILON)
        {
            break;
        }
        n = n + 1;
    }

    //  Estimate remainder of series ~ t_(n+1)/(1-x/(a+n+1))
    tnp1 = (x / (a + (n))) * term;
    remainder =  tnp1 / (1.0 - x / (a + (n) + 1.0));

    r.val = D.val * sum;
    r.err = D.err * Math.abs(sum) + Math.abs(D.val * remainder);
    r.err = r.err + (1 + n) * GSL_DBL_EPSILON * Math.abs(r.val);

    if (n >= nmax && Math.abs(remainder / sum) > GSL_SQRT_DBL_EPSILON)
    {
        throw "SF.MaxIterationsException"; //WITH "gamma_inc_P_series failed to converge";
    }

    return r;

} // gamma_inc_P_series

// ----------------------------------------------------------------------------

// Q large x asymptotic
//
// function gamma_inc_Q_large_x(a: LONG_FLOAT; x: LONG_FLOAT) return Result IS

//     nmax : CONSTANT INTEGER = 5000;
//     n    : INTEGER = 0;

//     sum  : LONG_FLOAT = 0.0;
//     term : LONG_FLOAT = 0.0;
//     last : LONG_FLOAT = 0.0;

//     D : Result;
//     r : Result;

// BEGIN -- gamma_inc_Q_large_x

//     D = gamma_inc_D(a, x);

//     sum  = 1.0;
//     term = 1.0;
//     last = 1.0;
//     n = 1;
//     WHILE (n < nmax) LOOP
//         term = term * (a - LONG_FLOAT(n)) / x;
//         EXIT WHEN (Math.abs(term / last) > 1.0);
//         EXIT WHEN (Math.abs(term / sum) < GSL_DBL_EPSILON);
//         sum  = sum + term;
//         last = term;
//         n = n + 1;
//     END LOOP;

//     IF (n >= nmax)
//         RAISE SF.MaxIterationsException WITH "error in large x asymptotic";
//     END IF;

//     r.val = D.val * (a / x) * sum;
//     r.err = D.err * Math.abs((a / x) * sum);
//     r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);

//     return r;

// END gamma_inc_Q_large_x;

// ----------------------------------------------------------------------------

// Uniform asymptotic for x near a, a and x large.
// See [Temme, p. 285]
//
function gamma_inc_Q_asymp_unif(a, x)
{
    var rta = 0.0;
    var eps = 0.0;
    var eta = 0.0;
    var c0  = 0.0;
    var c1  = 0.0;
    var R1  = 0.0;
    var lam = 0.0;
    var rt  = 0.0;

    var c = { val: 0.0, err: 0.0 }; // Result;
    var t = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    rta = Math.sqrt(a);
    eps = (x - a) / a;

    t = gsl_sf_log_1plusx_mx_e(eps);  // log(1+eps) - eps
    eta = GSL_SIGN(eps) * Math.sqrt(-2.0 * t.val);

    // This used to say erfc(eta*M_SQRT2*rta), which is wrong.
    // The sqrt(2) is in the denominator. Oops.
    // Fixed: [GJ] Mon Nov 15 13:25:32 MST 2004

    c = gsl_sf_erfc_e(eta * rta / M_SQRT2);

    if (Math.abs(eps) < GSL_ROOT5_DBL_EPSILON)
    {
        c0 = -1.0 / 3.0 + eps * (1.0 / 12.0 - eps * (23.0 / 540.0 - eps * (353.0 / 12960.0 - eps * 589.0 / 30240.0)));
        c1 = -1.0 / 540.0 - eps / 288.0;
    }
    else
    {
        rt = Math.sqrt(-2.0 * t.val / (eps * eps));
        lam = x / a;
        c0 = (1.0 - 1.0 / rt) / eps;
        c1 = -(eta * eta * eta * (lam * lam + 10.0 * lam + 1.0) - 12.0 * eps * eps * eps) / (12.0 * eta * eta * eta * eps * eps * eps);
    }

    R1 = Math.exp(-0.5 * a * eta * eta) / (M_SQRT2 * M_SQRTPI * rta) * (c0 + c1 / a);

    r.val = 0.5 * c.val + R1;
    r.err = GSL_DBL_EPSILON * Math.abs(R1 * 0.5 * a * eta * eta) + 0.5 * c.err;
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);

    return r;

} // gamma_inc_Q_asymp_unif

// ----------------------------------------------------------------------------

// Continued fraction which occurs in evaluation
// of Q(a,x) or Gamma(a,x).
//
//              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
//   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
//             1 +   1 +     1 +   1 +      1 +   1 +
//
// Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
//
// Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
// See gamma_inc_Q_CF() below.
//
//
function gamma_inc_F_CF(a, x)
{
    const nmax = 5000;
    var n = 0;
    var m = 0;

    var small = 0.0;
    var hn = 0.0;
    var Cn = 0.0;
    var Dn = 0.0;
    var an = 0.0;
    var delta1 = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    small = Math.pow(GSL_DBL_EPSILON, 3); //gsl_pow_3 (GSL_DBL_EPSILON);
    hn = 1.0;           // convergent
    Cn = 1.0 / small;
    Dn = 1.0;
   
    // n = 1 has a_1, b_1, b_0 independent of a,x,
    // so that has been done by hand
    m = 2;
    for (n = 2; n <= nmax - 1; n++)
    {
        if (GSL_IS_ODD(n))
        {
            an = 0.5 * (n - 1) / x;
        }
        else
        {
            an = (0.5 * (n) - a) / x;
        }
       
        Dn = 1.0 + an * Dn;
        if (Math.abs(Dn) < small)
        {
            Dn = small;
        }
        Cn = 1.0 + an / Cn;
        if (Math.abs(Cn) < small)
        {
            Cn = small;
        }
        Dn = 1.0 / Dn;
        delta1 = Cn * Dn;
        hn = hn * delta1;
        m  = m + 1;
        if (Math.abs(delta1 - 1.0) < GSL_DBL_EPSILON)
        {
            break;
        }
    }
   
    r.val = hn;
    r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(hn);
    r.err = r.err + GSL_DBL_EPSILON * (2.0 + 0.5 * (n)) * Math.abs(r.val);
   
    if (m >= nmax)
    {
        throw "SF.MaxIterationsException"; // WITH "error in CF for F(a,x)";
    }

    return r;

} // gamma_inc_F_CF

// ----------------------------------------------------------------------------

// Continued fraction for Q.
//
// Q(a,x) = D(a,x) a/x F(a,x)
//
// Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no):
//
// Since the Gautschi equivalent series method for CF evaluation may lead
// to singularities, I have replaced it with the modified Lentz algorithm
// given in
//
// I J Thompson and A R Barnett
// Coulomb and Bessel Functions of Complex Arguments and Order
// J Computational Physics 64:490-509 (1986)
//
// In consequence, gamma_inc_Q_CF_protected() is now obsolete and has been
// removed.
//
// Identification of terms between the above equation for F(a, x) and
// the first equation in the appendix of Thompson&Barnett is as follows:
//
//    b_0 = 0, b_n = 1 for all n > 0
//
//    a_1 = 1
//    a_n = (n/2-a)/x    for n even
//    a_n = (n-1)/(2x)   for n odd
//
//
function gamma_inc_Q_CF(a, x)
{
    var D = { val: 0.0, err: 0.0 }; // Result;
    var F = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    D = gamma_inc_D(a, x);
    F = gamma_inc_F_CF(a, x);

    r.val = D.val * (a / x) * F.val;
    r.err = D.err * Math.abs((a / x) * F.val) + Math.abs(D.val * a / x * F.err);

    return r;

} // gamma_inc_Q_CF

// ----------------------------------------------------------------------------

// Useful for small a and x. Handles the subtraction analytically.
//
function gamma_inc_Q_series(a, x)
{
    const nmax = 5000;
    var n = 0;

    var sum   = 0.0; // 1 + (a+1)/(a+2)(-x)/2! + (a+1)/(a+3)(-x)^2/3! + ...
    var term1 = 0.0; // 1 - x^a/Gamma(a+1)
    var term2 = 0.0; // a temporary variable used at the end
    var t     = 0.0;

    var r = { val: 0.0, err: 0.0 }; // Result;

    // Evaluate series for 1 - x^a/Gamma(a+1), small a

    var pg21 = -2.404113806319188570799476;  // PolyGamma[2,1]
    var lnx  = Math.log(x);
    var el   = M_EULER + lnx;
    var c1   = -el;
    var c2   = M_PI * M_PI / 12.0 - 0.5 * el * el;
    var c3   = el * (M_PI * M_PI / 12.0 - el * el / 6.0) + pg21 / 6.0;
    var c4   = -0.04166666666666666667
            * (-1.758243446661483480 + lnx)
            * (-0.764428657272716373 + lnx)
            * ( 0.723980571623507657 + lnx)
            * ( 4.107554191916823640 + lnx);
    var c5   = -0.0083333333333333333
            * (-2.06563396085715900 + lnx)
            * (-1.28459889470864700 + lnx)
            * (-0.27583535756454143 + lnx)
            * ( 1.33677371336239618 + lnx)
            * ( 5.17537282427561550 + lnx);
    var c6   = -0.0013888888888888889
            * (-2.30814336454783200 + lnx)
            * (-1.65846557706987300 + lnx)
            * (-0.88768082560020400 + lnx)
            * ( 0.17043847751371778 + lnx)
            * ( 1.92135970115863890 + lnx)
            * ( 6.22578557795474900 + lnx);
    var c7   = -0.00019841269841269841
            * (-2.5078657901291800 + lnx)
            * (-1.9478900888958200 + lnx)
            * (-1.3194837322612730 + lnx)
            * (-0.5281322700249279 + lnx)
            * ( 0.5913834939078759 + lnx)
            * ( 2.4876819633378140 + lnx)
            * ( 7.2648160783762400 + lnx);
    var c8   = -0.00002480158730158730
            * (-2.677341544966400 + lnx)
            * (-2.182810448271700 + lnx)
            * (-1.649350342277400 + lnx)
            * (-1.014099048290790 + lnx)
            * (-0.191366955370652 + lnx)
            * ( 0.995403817918724 + lnx)
            * ( 3.041323283529310 + lnx)
            * ( 8.295966556941250 + lnx);
    var c9   = -2.75573192239859e-6
            * (-2.8243487670469080 + lnx)
            * (-2.3798494322701120 + lnx)
            * (-1.9143674728689960 + lnx)
            * (-1.3814529102920370 + lnx)
            * (-0.7294312810261694 + lnx)
            * ( 0.1299079285269565 + lnx)
            * ( 1.3873333251885240 + lnx)
            * ( 3.5857258865210760 + lnx)
            * ( 9.3214237073814600 + lnx);
    var c10  = -2.75573192239859e-7
            * (-2.9540329644556910 + lnx)
            * (-2.5491366926991850 + lnx)
            * (-2.1348279229279880 + lnx)
            * (-1.6741881076349450 + lnx)
            * (-1.1325949616098420 + lnx)
            * (-0.4590034650618494 + lnx)
            * ( 0.4399352987435699 + lnx)
            * ( 1.7702236517651670 + lnx)
            * ( 4.1231539047474080 + lnx)
            * ( 10.342627908148680 + lnx);

    term1 = a * (c1 + a * (c2 + a * (c3 + a * (c4 + a * (c5 + a * (c6 + a * (c7 + a * (c8 + a * (c9 + a * c10)))))))));

    // Evaluate the sum.

    t = 1.0;
    sum = 1.0;

    n = 1;
    while (n < nmax)
    {
        t = -t * x / (n + 1);
        sum = sum + (a + 1.0) / (a + (n) + 1.0) * t;
        if (Math.abs(t / sum) < GSL_DBL_EPSILON)
        {
            break;
        }
        n = n + 1;
    }

    if (n >= nmax)
    {
       throw "SF.MaxIterationsException";
    }

    term2 = (1.0 - term1) * a / (a + 1.0) * x * sum;
    r.val = term1 + term2;
    r.err = GSL_DBL_EPSILON * (Math.abs(term1) + 2.0 * Math.abs(term2));
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);

    return r;

} // gamma_inc_Q_series

// ----------------------------------------------------------------------------

// series for small a and x, but not defined for a = 0
function gamma_inc_series(a, x)
{
    var Q = { val: 0.0, err: 0.0 }; // Result;
    var G = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    Q = gamma_inc_Q_series(a, x);
    G = gsl_sf_gamma_e(a);
    r.val = Q.val * G.val;
    r.err = Math.abs(Q.val * G.err) + Math.abs(Q.err * G.val);
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);

    return r;

} // gamma_inc_series

// ----------------------------------------------------------------------------

function gamma_inc_a_gt_0(a, x)
{
    var Q = { val: 0.0, err: 0.0 }; // Result;
    var G = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    // x > 0 and a > 0; use result for Q
    Q = gsl_sf_gamma_inc_Q_e(a, x);
    G = gsl_sf_gamma_e(a);

    r.val = G.val * Q.val;
    r.err = Math.abs(G.val * Q.err) + Math.abs(G.err * Q.val);
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);

    return r;

} // gamma_inc_a_gt_0

// ----------------------------------------------------------------------------

function gamma_inc_CF(a, x)
{
    var am1lgx = 0.0;

    var r   = { val: 0.0, err: 0.0 }; // Result;
    var F   = { val: 0.0, err: 0.0 }; // Result;
    var pre = { val: 0.0, err: 0.0 }; // Result;

    am1lgx = (a - 1.0) * Math.log(x);
    F = gamma_inc_F_CF(a, x);
    pre = gsl_sf_exp_err_e(am1lgx - x, GSL_DBL_EPSILON * Math.abs(am1lgx));
   
    r.val = F.val * pre.val;
    r.err = Math.abs(F.err * pre.val) + Math.abs(F.val * pre.err);
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
   
    return r;

} // gamma_inc_CF

// ----------------------------------------------------------------------------

// evaluate Gamma(0,x), x > 0
function GAMMA_INC_A_0(x)
{
    return gsl_sf_expint_E1_e(x);
} // GAMMA_INC_A_0

// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_gamma_inc_Q_e(a, x)
{
    var P = { val: 0.0, err: 0.0 }; // Result;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (a < 0.0 || x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        r.val = 1.0;
        r.err = 0.0;
    }
    else if (a == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (x <= 0.5 * a)
    {
        // If the series is quick, do that. It is
        // robust and simple.
        P = gamma_inc_P_series(a, x);
        r.val = 1.0 - P.val;
        r.err = P.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (a >= 1.0e+06 && (x - a) * (x - a) < a)
    {
        // Then try the difficult asymptotic regime.
        // This is the only way to do this region.
        r = gamma_inc_Q_asymp_unif(a, x);
    }
    else if (a < 0.2 && x < 5.0)
    {
        // Cancellations at small a must be handled
        // analytically; x should not be too big
        // either since the series terms grow
        // with x and log(x).
        r = gamma_inc_Q_series(a, x);
    }
    else if (a <= x)
    {
        if (x <= 1.0e+06)
        {
            // Continued fraction is excellent for x >~ a.
            // We do not let x be too large when x > a since
            // it is somewhat pointless to try this there;
            // the function is rapidly decreasing for
            // x large and x > a, and it will just
            // underflow in that region anyway. We
            // catch that case in the standard
            // large-x method.
            r = gamma_inc_Q_CF(a, x);
        }
        else
        {
            r = gamma_inc_Q_large_x(a, x);
        }
    }
    else
    {
        if (x > a - Math.sqrt(a))
        {
            // Continued fraction again. The convergence
            // is a little slower here, but that is fine.
            // We have to trade that off against the slow
            // convergence of the series, which is the
            // only other option.
            r = gamma_inc_Q_CF(a, x);
        }
        else
        {
            P = gamma_inc_P_series(a, x);
            r.val = 1.0 - P.val;
            r.err = P.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
    }

    return r;

} // gsl_sf_gamma_inc_Q_e

// ----------------------------------------------------------------------------

export function gsl_sf_gamma_inc_P_e(a, x)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    var Q = { val: 0.0, err: 0.0 }; // Result;

    if (a <= 0.0 || x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (x < 20.0 || x < 0.5 * a)
    {
        // Do the easy series cases. Robust and quick.
        //
        r = gamma_inc_P_series(a, x);
    }
    else if (a > 1.0e+06 && (x - a) * (x - a) < a)
    {
        // Crossover region. Note that Q and P are
        // roughly the same order of magnitude here,
        // so the subtraction is stable.
        //
        Q = gamma_inc_Q_asymp_unif(a, x);
        r.val = 1.0 - Q.val;
        r.err = Q.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (a <= x)
    {
        // Q <~ P in this area, so the
        // subtractions are stable.
        //
        if (a > 0.2 * x)
        {
            Q = gamma_inc_Q_CF(a, x);
        }
        else
        {
            Q = gamma_inc_Q_large_x(a, x);
        }
        r.val = 1.0 - Q.val;
        r.err = Q.err;
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        if ((x - a) * (x - a) < a)
        {
            // This condition is meant to insure
            // that Q is not very close to 1,
            // so the subtraction is stable.
            //
            Q = gamma_inc_Q_CF(a, x);
            r.val = 1.0 - Q.val;
            r.err = Q.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
        else
        {
            r = gamma_inc_P_series(a, x);
        }
    }

    return r;

} // gsl_sf_gamma_inc_P_e

// ----------------------------------------------------------------------------

export function gsl_sf_gamma_inc_e(a, x)
{
    var fa    = 0.0;
    var da    = 0.0;
    var alpha = 0.0;
    var gax   = 0.0;
    var shift = 0.0;

    var r    = { val: 0.0, err: 0.0 }; // Result;
    var g_da = { val: 0.0, err: 0.0 }; // Result;

    if (x < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x == 0.0)
    {
        r = gsl_sf_gamma_e(a);
    }
    else if (a == 0.0)
    {
        r = GAMMA_INC_A_0(x);
    }
    else if (a > 0.0)
    {
        r = gamma_inc_a_gt_0(a, x);
    }
    else if (x > 0.25)
    {
        // continued fraction seems to fail for x too small; otherwise
        // it is ok, independent of the value of |x/a|, because of the
        // non-oscillation in the expansion, i.e. the CF is
        // un-conditionally convergent for a < 0 and x > 0
        r = gamma_inc_CF(a, x);
    }
    else if (Math.abs(a) < 0.5)
    {
        r = gamma_inc_series(a, x);
    }
    else
    {
        // a = fa + da; da >= 0
        fa = Math.floor(a);
        da = a - fa;

        if (da > 0.0)
        {
            g_da = gamma_inc_a_gt_0(da, x);
        }
        else
        {
            g_da = GAMMA_INC_A_0(x);
        }

        alpha = da;
        gax = g_da.val;

        // Gamma(alpha-1,x) = 1/(alpha-1) (Gamma(a,x) - x^(alpha-1) e^-x)
        while (true)
        {
            shift = Math.exp(-x + (alpha - 1.0) * Math.log(x));
            gax = (gax - shift) / (alpha - 1.0);
            alpha = alpha - 1.0;
            //EXIT WHEN NOT (alpha > a);
            if (alpha <= a)
            {
                break;
            }
        }

        r.val = gax;
        r.err = 2.0 * (1.0 + Math.abs(a)) * GSL_DBL_EPSILON * Math.abs(gax);
    }

    return r;

 } // gsl_sf_gamma_inc_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

// function gsl_sf_gamma_inc_P(a: LONG_FLOAT; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_gamma_inc_P
//     return EVAL_RESULT(gsl_sf_gamma_inc_P_e'Access, (a, x), "gsl_sf_gamma_inc_P");
// END gsl_sf_gamma_inc_P;

// function gsl_sf_gamma_inc_Q(a: LONG_FLOAT; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_gamma_inc_Q
//     return EVAL_RESULT(gsl_sf_gamma_inc_Q_e'Access, (a, x), "gsl_sf_gamma_inc_Q");
// END gsl_sf_gamma_inc_Q;

// function gsl_sf_gamma_inc(a: LONG_FLOAT; x: LONG_FLOAT) return LONG_FLOAT IS
// BEGIN -- gsl_sf_gamma_inc
//      return EVAL_RESULT(gsl_sf_gamma_inc_e'Access, (a, x), "gsl_sf_gamma_inc");
// END gsl_sf_gamma_inc;

// ----------------------------------------------------------------------------
// EOF SF-GammaIncomplete.mjs
