// SF-BesselTemme.mjs
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

import { M_PI }                  from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { cheb_eval_e }           from "./SF-Chebyshev.mjs";

// ----------------------------------------------------------------------------

// nu = (x+1)/4, -1<x<1, 1/(2nu)(1/Gamma[1-nu]-1/Gamma[1+nu])
const g1_dat =
    [
   -1.14516408366268311786898152867,
    0.00636085311347084238122955495,
    0.00186245193007206848934643657,
    0.000152833085873453507081227824,
    0.000017017464011802038795324732,
   -6.4597502923347254354668326451e-07,
   -5.1819848432519380894104312968e-08,
    4.5189092894858183051123180797e-10,
    3.2433227371020873043666259180e-11,
    6.8309434024947522875432400828e-13,
    2.8353502755172101513119628130e-14,
   -7.9883905769323592875638087541e-16,
   -3.3726677300771949833341213457e-17,
   -3.6586334809210520744054437104e-20
    ];
const g1_cs = { length: 13, c: g1_dat, order: 13, a: -1.0, b: 1.0, order_sp: 7 };

// nu = (x+1)/4, -1<x<1,  1/2 (1/Gamma[1-nu]+1/Gamma[1+nu])
const g2_dat =
    [
    1.882645524949671835019616975350,
   -0.077490658396167518329547945212,  
   -0.018256714847324929419579340950,
    0.0006338030209074895795923971731,
    0.0000762290543508729021194461175,
   -9.5501647561720443519853993526e-07,
   -8.8927268107886351912431512955e-08,
   -1.9521334772319613740511880132e-09,
   -9.4003052735885162111769579771e-11,
    4.6875133849532393179290879101e-12,
    2.2658535746925759582447545145e-13,
   -1.1725509698488015111878735251e-15,
   -7.0441338200245222530843155877e-17,
   -2.4377878310107693650659740228e-18,
   -7.5225243218253901727164675011e-20
    ];
const g2_cs = { length: 14, c: g2_dat, order: 14, a: -1.0, b: 1.0, order_sp: 8 };

// ----------------------------------------------------------------------------

function gsl_sf_temme_gamma( nu )//: LONG_FLOAT; g_1pnu: IN OUT LONG_FLOAT; g_1mnu: IN OUT LONG_FLOAT; g1: IN OUT LONG_FLOAT; g2: IN OUT LONG_FLOAT) IS
{
    var anu  = Math.abs( nu ); // functions are even
    var x    = 4.0 * anu - 1.0;
    var r_g1 = { val: 0.0, err: 0.0 }; // Result;
    var r_g2 = { val: 0.0, err: 0.0 }; // Result;
    var r = { g1: 0.0, g2: 0.0, g_1mnu: 0.0, g_1pnu: 0.0 };

    r_g1 = cheb_eval_e( g1_cs, x );
    r_g2 = cheb_eval_e( g2_cs, x );
    r.g1 = r_g1.val;
    r.g2 = r_g2.val;
    r.g_1mnu = 1.0 / (r_g2.val + nu * r_g1.val);
    r.g_1pnu = 1.0 / (r_g2.val - nu * r_g1.val);

    return r;

} // gsl_sf_temme_gamma

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_Y_temme( nu, x, Ynu, Ynup1 )
{
    const max_iter = 15000;
    var k = 0;
    
    var half_x    = 0.0;
    var ln_half_x = 0.0;
    var half_x_nu = 0.0;
    var pi_nu     = 0.0;
    var alpha     = 0.0;
    var sigma     = 0.0;
    var sinrat    = 0.0;
    var sinhrat   = 0.0;
    var sinhalf   = 0.0;
    var sin_sqr   = 0.0;
    var sum0      = 0.0;
    var sum1      = 0.0;
    var fk        = 0.0;
    var pk        = 0.0;
    var qk        = 0.0;
    var hk        = 0.0;
    var ck        = 0.0;
    var g_1pnu    = 0.0;
    var g_1mnu    = 0.0;
    var g1        = 0.0;
    var g2        = 0.0;
    var del0      = 0.0;
    var del1      = 0.0;
    var gk        = 0.0;

    var rg = { g1: 0.0, g2: 0.0, g_1mnu: 0.0, g_1pnu: 0.0 };

    half_x = 0.5 * x;
    ln_half_x = Math.log( half_x );
    half_x_nu = Math.exp( nu * ln_half_x );
    pi_nu   = M_PI * nu;
    alpha   = pi_nu / 2.0;
    sigma   = -nu * ln_half_x;
    if ( Math.abs( pi_nu ) < GSL_DBL_EPSILON )
    {
        sinrat = 1.0;
    }
    else
    {
        sinrat = pi_nu / Math.sin( pi_nu );
    }
    if ( Math.abs( sigma ) < GSL_DBL_EPSILON )
    {
        sinhrat = 1.0;
    }
    else
    {
        sinhrat = Math.sinh( sigma ) / sigma;
    }
    if ( Math.abs( alpha ) < GSL_DBL_EPSILON )
    {
        sinhalf = 1.0;
    }
    else
    {
        sinhalf = Math.sin( alpha ) / alpha;
    }
    sin_sqr = nu * M_PI * M_PI * 0.5 * sinhalf * sinhalf;
    
    rg = gsl_sf_temme_gamma( nu );//, g_1pnu, g_1mnu, g1, g2 );
    g_1pnu = rg.g_1pnu;
    g_1mnu = rg.g_1mnu;
    g1 = rg.g1;
    g2 = rg.g2;
   
    fk = 2.0 / M_PI * sinrat * (Math.cosh(sigma) * g1 - sinhrat * ln_half_x * g2);
    pk = 1.0 / M_PI / half_x_nu * g_1pnu;
    qk = 1.0 / M_PI * half_x_nu * g_1mnu;
    hk = pk;
    ck = 1.0;
   
    sum0 = fk + sin_sqr * qk;
    sum1 = pk;
   
    k = 0;
    while ( k < max_iter )
    {
        k = k + 1;
        fk  = ((k) * fk + pk + qk) / ((k) * (k) - nu * nu);
        ck  = ck * (-half_x * half_x / (k));
        pk  = pk / ((k) - nu);
        qk  = qk / ((k) + nu);
        gk  = fk + sin_sqr * qk;
        hk  = -(k) * gk + pk; 
        del0 = ck * gk;
        del1 = ck * hk;
        sum0 = sum0 + del0;
        sum1 = sum1 + del1;
        if ( Math.abs( del0 ) < 0.5 * (1.0 + Math.abs( sum0 )) * GSL_DBL_EPSILON ) break;
    }
   
    Ynu.val   = -sum0;
    Ynu.err   = (2.0 + 0.5 * (k)) * GSL_DBL_EPSILON * Math.abs ( Ynu.val );
    Ynup1.val = -sum1 * 2.0 / x;
    Ynup1.err = (2.0 + 0.5 * (k)) * GSL_DBL_EPSILON * Math.abs ( Ynup1.val );
   
    if ( k >= max_iter )
    {
        throw "SF.MaxIterationsException";
    }

} // gsl_sf_bessel_Y_temme

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_K_scaled_temme( nu, x )//, K_nu: IN OUT LONG_FLOAT; K_nup1: IN OUT LONG_FLOAT; Kp_nu: IN OUT LONG_FLOAT) IS
{
    var rk = { K_nu: 0.0, K_nup1: 0.0, Kp_nu: 0.0 };

    const max_iter = 15000;
  
    var half_x    = 0.0;
    var ln_half_x = 0.0;
    var half_x_nu = 0.0;
    var pi_nu     = 0.0;
    var sigma     = 0.0;
    var sinrat    = 0.0;
    var sinhrat   = 0.0;
    var ex        = 0.0;
    var sum0      = 0.0;
    var sum1      = 0.0;
    var fk        = 0.0;
    var pk        = 0.0;
    var qk        = 0.0;
    var hk        = 0.0;
    var ck        = 0.0;
    var g_1pnu    = 0.0;
    var g_1mnu    = 0.0;
    var g1        = 0.0;
    var g2        = 0.0;
    var del0      = 0.0;
    var del1      = 0.0;
    var k = 0;

    var rg = { g1: 0.0, g2: 0.0, g_1mnu: 0.0, g_1pnu: 0.0 };
  
    half_x    = 0.5 * x;
    ln_half_x = Math.log( half_x );
    half_x_nu = Math.exp( nu * ln_half_x );
    pi_nu   = M_PI * nu;
    sigma   = -nu * ln_half_x;
    if ( Math.abs( pi_nu ) < GSL_DBL_EPSILON )
    {
        sinrat = 1.0;
    }
    else
    {
        sinrat = pi_nu / Math.sin( pi_nu );
    }
    if ( Math.abs( sigma ) < GSL_DBL_EPSILON )
    {
        sinhrat = 1.0;
    }
    else
    {
        sinhrat = Math.sinh( sigma ) / sigma;
    }
    ex = Math.exp( x );
  
    k = 0;
  
    rg = gsl_sf_temme_gamma( nu ); //, g_1pnu, g_1mnu, g1, g2);
    g_1pnu = rg.g_1pnu;
    g_1mnu = rg.g_1mnu;
    g1 = rg.g1;
    g2 = rg.g2;
  
    fk = sinrat * (Math.cosh( sigma ) * g1 - sinhrat * ln_half_x * g2);
    pk = 0.5 / half_x_nu * g_1pnu;
    qk = 0.5 * half_x_nu * g_1mnu;
    hk = pk;
    ck = 1.0;
    sum0 = fk;
    sum1 = hk;
    while ( k < max_iter )
    {
        k = k + 1;
        fk = ((k) * fk + pk + qk) / ((k) * (k) - nu * nu);
        ck = ck * half_x * half_x / (k);
        pk = pk / ((k) - nu);
        qk = qk / ((k) + nu);
        hk  = -(k) * fk + pk;
        del0 = ck * fk;
        del1 = ck * hk;
        sum0 = sum0 + del0;
        sum1 = sum1 + del1;
        if ( Math.abs( del0 ) < 0.5 * Math.abs( sum0 ) * GSL_DBL_EPSILON ) break;
    }
    
    rk.K_nu   = sum0 * ex;
    rk.K_nup1 = sum1 * 2.0 / x * ex;
    rk.Kp_nu  = -rk.K_nup1 + nu / x * rk.K_nu;
  
    if ( k >= max_iter )
    {
        throw "SF.MaxIterationsException";
    }

    return rk;

} // gsl_sf_bessel_K_scaled_temme

// ----------------------------------------------------------------------------
// EOF SF-BesselTemme.mjs
