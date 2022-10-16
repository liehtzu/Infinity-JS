// SF-Coulomb.mjs
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
// Translation to Ada: Igor Izvarin

import { M_PI }                 from "./SF-Math.mjs";
import { M_EULER }              from "./SF-Math.mjs";
import { M_LN2 }                from "./SF-Math.mjs";
import { GSL_SIGN }             from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MIN }      from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_DBL_MAX }          from "./SF-Machine.mjs";
import { gsl_sf_lngamma_e }     from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_complex_e } from "./SF-Gamma.mjs";
import { gsl_sf_exp_err_e }     from "./SF-Exponential.mjs";
import { gsl_sf_expm1_e }       from "./SF-Exponential.mjs";
import { gsl_sf_psi_1piy_e }    from "./SF-Psi.mjs";
import { GSL_MODE_DEFAULT }     from "./SF-Mode.mjs";
import { gsl_sf_hypot }         from "./SF-Trigonometric.mjs";
import { gsl_sf_airy_Ai_scaled_e } from "./SF-Airy.mjs";
import { gsl_sf_airy_Bi_scaled_e } from "./SF-Airy.mjs";

// ----------------------------------------------------------------------------

// Evaluation of Coulomb wave functions F_L(eta, x), G_L(eta, x),
// and their derivatives. A combination of Steed's method, asymptotic
// results, and power series.
//
// Steed's method:
//  [Barnett, CPC 21, 297 (1981)]
// Power series and other methods:
//  [Biedenharn et al., PR 97, 542 (1954)]
//  [Bardin et al., CPC 3, 73 (1972)]
//  [Abad+Sesma, CPC 71, 110 (1992)]
//

// the L=0 normalization constant
// [Abramowitz+Stegun 14.1.8]
//
function C0sq( eta )
{
    var twopieta = 2.0 * M_PI * eta;

    if ( Math.abs( eta ) < GSL_DBL_EPSILON )
    {
        return 1.0;
    }
    else if ( twopieta > GSL_LOG_DBL_MAX )
    {
        return 0.0;
    }
    else
    {
        let scale= { val: 0.0, err: 0.0 };
        scale = gsl_sf_expm1_e( twopieta);
        return twopieta / scale.val;
    }
}

// ----------------------------------------------------------------------------

// the full definition of C_L(eta) for any valid L and eta
// [Abramowitz and Stegun 14.1.7]
// This depends on the complex gamma function. For large
// arguments the phase of the complex gamma function is not
// very accurately determined. However the modulus is, and that
// is all that we need to calculate C_L.
//
// This is not valid for L <= -3/2  or  L = -1.
//
function CLeta( /* IN Double */ L, /* IN Double */ eta )
{

    var sgn     = 1.0;
    var arg_val = 0.0;
    var arg_err = 0.0;
    var ln1     = { val: 0.0, err: 0.0 }; // Result; // log of numerator Gamma function
    var ln2     = { val: 0.0, err: 0.0 }; // Result; // log of denominator Gamma function
    var p1      = { val: 0.0, err: 0.0 }; // Result; // phase of numerator Gamma    not used
  
    if ( Math.abs( eta / (L + 1.0) ) < GSL_DBL_EPSILON )
    {
        ln1 = gsl_sf_lngamma_e( L + 1.0 );
    }
    else
    {
        gsl_sf_lngamma_complex_e( L + 1.0, eta, ln1, p1 ); // should be ok
    }
  
    ln2 = gsl_sf_lngamma_e( 2.0 * (L + 1.0) );
    if ( L < -1.0 )
    {
        sgn = -sgn;
    }
  
    arg_val = L * M_LN2 - 0.5 * eta * M_PI + ln1.val - ln2.val;
    arg_err = ln1.err + ln2.err;
    arg_err = arg_err + GSL_DBL_EPSILON * (Math.abs( L * M_LN2 ) + Math.abs( 0.5 * eta * M_PI ));
    return gsl_sf_exp_err_e( arg_val, arg_err );

} // CLeta

// ----------------------------------------------------------------------------

function gsl_sf_coulomb_CL_e( lam, eta )
{
    var result = { val: 0.0, err: 0.0 };

    if ( lam <= -1.0 )
    {
        throw "SF.DomainException";
    }
    else if ( Math.abs(lam) < GSL_DBL_EPSILON )
    {
        /* saves a calculation of complex_lngamma(), otherwise not necessary */
        result.val = Math.sqrt( C0sq( eta ) );
        result.err = 2.0 * GSL_DBL_EPSILON * result.val;
        return result;
    }
    else
    {
        return CLeta( lam, eta, result );
    }
}

// ----------------------------------------------------------------------------

// cl[0] .. cl[kmax] = C_{lam_min}(eta) .. C_{lam_min+kmax}(eta)
//
export function gsl_sf_coulomb_CL_array( lam_min, kmax, eta, cl )
{
    var cl_0 = { val: 0.0, err: 0.0 };
    cl_0 = gsl_sf_coulomb_CL_e( lam_min, eta );
    cl[0] = cl_0.val;

    for ( let k = 1; k <= kmax; k++ )
    {
        let L = lam_min + k;
        cl[k] = cl[k-1] * gsl_sf_hypot( L, eta ) / (L * (2.0 * L + 1.0));
    }
}

// ----------------------------------------------------------------------------

// Determine the connection phase, phi_lambda.
// See coulomb_FG_series() below. We have
// to be careful about sin(phi).0. Note that
// there is an underflow condition for large 
// positive eta in any case.
//
function coulomb_connection
    (
    lam,     // IN     Double
    eta,     // IN     Double
    cos_phi, // IN OUT Double
    sin_phi  // IN OUT Double
    )
{
    var eps = 0.0;
    var tpl = 0.0;
    var dth = 0.0;
    var phi = 0.0;
    var X   = 0.0;

    if ( eta > -GSL_LOG_DBL_MIN / 2.0 * M_PI - 1.0 )
    {
        //cos_phi = 1.0;
        //sin_phi = 0.0;
        throw "SF.UnderflowException";
    }
    else if ( eta > -GSL_LOG_DBL_EPSILON / (4.0 * M_PI) )
    {
        eps = 2.0 * Math.exp( -2.0 * M_PI * eta );
        tpl = Math.tan( M_PI * lam );
        dth = eps * tpl / (tpl * tpl + 1.0);
        cos_phi.Double = -1.0 + 0.5 * dth * dth;
        sin_phi.Double = -dth;
        return;
    }
    else
    {
        X   = Math.tanh( M_PI * eta ) / Math.tan( M_PI * lam );
        phi = -Math.atan( X ) - (lam + 0.5) * M_PI;
        cos_phi.Double = Math.cos( phi );
        sin_phi.Double = Math.sin( phi );
        return;
    }

} // coulomb_connection

// ----------------------------------------------------------------------------

// Evaluate the Frobenius series for F_lam(eta,x) and G_lam(eta,x).
// Homegrown algebra. Evaluates the series for F_{lam} and
// F_{-lam-1}, then uses
//    G_{lam} = (F_{lam} cos(phi) - F_{-lam-1}) / sin(phi)
// where
//    phi = Arg[Gamma[1+lam+I eta]] - Arg[Gamma[-lam + I eta]] - (lam+1/2)Pi
//        = Arg[Sin[Pi(-lam+I eta)] - (lam+1/2)Pi
//        = atan2(-cos(lam Pi)sinh(eta Pi), -sin(lam Pi)cosh(eta Pi)) - (lam+1/2)Pi
//
//        = -atan(X) - (lam+1/2) Pi,  X = tanh(eta Pi)/tan(lam Pi)
//
// Not appropriate for lam <= -1/2, lam = 0, or lam >= 1/2.
//
function coulomb_FG_series
    (
    lam, // IN    Double
    eta, // IN    Double
    x,   // IN    Double
    F,   // IN OUT Result
    G    // IN OUT Result
    )
{
    const max_iter = 800;
    const tlp1     = 2.0 * lam + 1.0;
    const pow_x    = Math.pow( x, lam );
    var m         = 0;
    var cos_phi_lam = { Double: 0.0 };
    var sin_phi_lam = { Double: 0.0 };
    var uA_mm2      = 0.0;
    var uA_mm1      = 0.0;
    var uA_m        = 0.0;
    var uB_mm2      = 0.0;
    var uB_mm1      = 0.0;
    var uB_m        = 0.0;
    var A_sum       = 0.0;
    var B_sum       = 0.0;
    var abs_dA      = 0.0;
    var abs_dB      = 0.0;
    var max_abs_dA  = 0.0;
    var max_abs_dB  = 0.0;
    var abs_A       = 0.0;
    var abs_B       = 0.0;
    var A_abs_del_prev = 0.0;
    var B_abs_del_prev = 0.0;
    var FA    = { val: 0.0, err: 0.0 }; // Result;
    var FB    = { val: 0.0, err: 0.0 }; // Result;
    var ClamA = { val: 0.0, err: 0.0 }; // Result;
    var ClamB = { val: 0.0, err: 0.0 }; // Result;

    ClamA = CLeta( lam, eta );
    ClamB = CLeta( -lam - 1.0, eta );
  
    uA_mm2 = 1.0;                  // uA sum is for F_{lam}
    uA_mm1 = x * eta / (lam + 1.0);
    uB_mm2 = 1.0;                  // uB sum is for F_{-lam-1}
    uB_mm1 = -x * eta / lam;
    A_sum = uA_mm2 + uA_mm1;
    B_sum = uB_mm2 + uB_mm1;
    A_abs_del_prev = Math.abs( A_sum );
    B_abs_del_prev = Math.abs( B_sum );
    m = 2;
  
    coulomb_connection( lam, eta, cos_phi_lam, sin_phi_lam );
  
    //if (stat_conn = GSL_EUNDRFLW)
        //F.val = 0.0;  // FIXME: should this be set to Inf too like G?
        //F.err = 0.0;
    //    throw SF.OverflowException;
    //END if;
  
    while ( m < max_iter )
    {
        uA_m = x * (2.0 * eta * uA_mm1 - x * uA_mm2) / ((m) * ((m) + tlp1));
        uB_m = x * (2.0 * eta * uB_mm1 - x * uB_mm2) / ((m) * ((m) - tlp1));
        A_sum = A_sum + uA_m;
        B_sum = B_sum + uB_m;
        abs_dA = Math.abs( uA_m );
        abs_dB = Math.abs( uB_m );
        if ( m > 15 )
        {
            // Don't bother checking until we have gone out a little ways;
            // a minor optimization. Also make sure to check both the
            // current and the previous increment because the odd and even
            // terms of the sum can have very different behaviour, depending
            // on the value of eta.
            //
            max_abs_dA = Math.max( abs_dA, A_abs_del_prev );
            max_abs_dB = Math.max( abs_dB, B_abs_del_prev );
            abs_A = Math.abs( A_sum );
            abs_B = Math.abs( B_sum );
            if ( max_abs_dA / (max_abs_dA + abs_A) < 4.0 * GSL_DBL_EPSILON
                && max_abs_dB / (max_abs_dB + abs_B) < 4.0 * GSL_DBL_EPSILON ) break;
        }
        A_abs_del_prev = abs_dA;
        B_abs_del_prev = abs_dB;
        uA_mm2 = uA_mm1;
        uA_mm1 = uA_m;
        uB_mm2 = uB_mm1;
        uB_mm1 = uB_m;
        m = m + 1;
    }
  
    FA.val = A_sum * ClamA.val * pow_x * x;
    FA.err = Math.abs( A_sum ) * ClamA.err * pow_x * x + 2.0 * GSL_DBL_EPSILON * Math.abs( FA.val );
    FB.val = B_sum * ClamB.val / pow_x;
    FB.err = Math.abs( B_sum ) * ClamB.err / pow_x + 2.0 * GSL_DBL_EPSILON * Math.abs( FB.val );
  
    F.val = FA.val;
    F.err = FA.err;
  
    G.val = (FA.val * cos_phi_lam.Double - FB.val) / sin_phi_lam.Double;
    G.err = (FA.err * Math.abs( cos_phi_lam.Double ) + FB.err) / Math.abs( sin_phi_lam.Double );
  
    if ( m >= max_iter )
    {
        throw "SF.MaxIterationsException";
    }

} // coulomb_FG_series

// ----------------------------------------------------------------------------

// Evaluate the Frobenius series for F_0(eta,x) and G_0(eta,x).
// See [Bardin et al., CPC 3, 73 (1972), (14)-(17)];
// note the misprint in (17): nu_0=1 is correct, not nu_0=0.
//
function coulomb_FG0_series
    (
    eta, // IN    Double
    x,   // IN    Double
    F,   // IN OUT Result
    G    // IN OUT Result
    )
{
    const max_iter   = 800;
    const x2         = x * x;
    const tex        = 2.0 * eta * x;
    var m          = 0;
    //var stat_CL    = 0;
    //var psi_stat   = 0;
    var u_mm2      = 0.0;
    var u_mm1      = 0.0;
    var u_m        = 0.0;
    var v_mm2      = 0.0;
    var v_mm1      = 0.0;
    var v_m        = 0.0;
    var u_sum      = 0.0;
    var v_sum      = 0.0;
    var u_abs_del_prev = 0.0;
    var v_abs_del_prev = 0.0;
    var u_sum_err  = 0.0;
    var v_sum_err  = 0.0;
    var ln2x       = 0.0;
    var abs_du     = 0.0;
    var abs_dv     = 0.0;
    var m_mm1      = 0.0;
    var max_abs_du = 0.0;
    var max_abs_dv = 0.0;
    var abs_u      = 0.0;
    var abs_v      = 0.0;
    var C0         = { val: 0.0, err: 0.0 }; // Result;
    var r1pie      = { val: 0.0, err: 0.0 }; // Result;

    C0 = CLeta( 0.0, eta );
    r1pie = gsl_sf_psi_1piy_e( eta );
    u_mm2 = 0.0;  // u_0
    u_mm1 = x;    // u_1
    v_mm2 = 1.0;                                       // nu_0
    v_mm1 = tex * (2.0 * M_EULER - 1.0 + r1pie.val);   // nu_1
    u_sum = u_mm2 + u_mm1;
    v_sum = v_mm2 + v_mm1;
    u_abs_del_prev = Math.abs( u_sum );
    v_abs_del_prev = Math.abs( v_sum );
    m = 2;
    u_sum_err = 2.0 * GSL_DBL_EPSILON * Math.abs( u_sum );
    v_sum_err = 2.0 * GSL_DBL_EPSILON * Math.abs( v_sum );
    ln2x = Math.log( 2.0 * x );
  
    while ( m < max_iter )
    {
        m_mm1 = (m) * (m - 1);
        u_m = (tex * u_mm1 - x2 * u_mm2) / m_mm1;
        v_m = (tex * v_mm1 - x2 * v_mm2 - 2.0 * eta * (2 * m - 1) * u_m) / m_mm1;
        u_sum = u_sum + u_m;
        v_sum = v_sum + v_m;
        abs_du = Math.abs( u_m );
        abs_dv = Math.abs( v_m );
        u_sum_err = u_sum_err + 2.0 * GSL_DBL_EPSILON * abs_du;
        v_sum_err = v_sum_err + 2.0 * GSL_DBL_EPSILON * abs_dv;
        if ( m > 15 )
        {
            // Don't bother checking until we have gone out a little ways;
            // a minor optimization. Also make sure to check both the
            // current and the previous increment because the odd and even
            // terms of the sum can have very different behaviour, depending
            // on the value of eta.
            //
            max_abs_du = Math.max( abs_du, u_abs_del_prev );
            max_abs_dv = Math.max( abs_dv, v_abs_del_prev );
            abs_u = Math.abs( u_sum );
            abs_v = Math.abs( v_sum );
            if ( max_abs_du / (max_abs_du + abs_u) < 40.0 * GSL_DBL_EPSILON
                && max_abs_dv / (max_abs_dv + abs_v) < 40.0 * GSL_DBL_EPSILON ) break;
        }
        u_abs_del_prev = abs_du;
        v_abs_del_prev = abs_dv;
        u_mm2 = u_mm1;
        u_mm1 = u_m;
        v_mm2 = v_mm1;
        v_mm1 = v_m;
        m = m + 1;
    }
  
    F.val = C0.val * u_sum;
    F.err = C0.err * Math.abs( u_sum );
    F.err = F.err + Math.abs( C0.val ) * u_sum_err;
    F.err = F.err + 2.0 * GSL_DBL_EPSILON * Math.abs( F.val );
  
    G.val = (v_sum + 2.0 * eta * u_sum * ln2x) / C0.val;
    G.err = (Math.abs( v_sum ) + Math.abs( 2.0 * eta * u_sum * ln2x )) / Math.abs( C0.val ) * Math.abs( C0.err / C0.val );
    G.err = G.err + (v_sum_err + Math.abs( 2.0 * eta * u_sum_err * ln2x )) / Math.abs( C0.val );
    G.err = G.err + 2.0 * GSL_DBL_EPSILON * Math.abs( G.val );
  
    if ( m >= max_iter )
    {
        throw "SF.MaxIterationsException";
    }

} // coulomb_FG0_series

// ----------------------------------------------------------------------------

// Evaluate the Frobenius series for F_{-1/2}(eta,x) and G_{-1/2}(eta,x).
// Homegrown algebra.
//
function coulomb_FGmhalf_series
    (
    eta, // IN    Double
    x,   // IN    Double
    F,   // IN OUT Result
    G    // IN OUT Result
    )
{
    const max_iter = 800;
    var rx    = Math.sqrt(x);
    var x2    = x * x;
    var tex   = 2.0 * eta * x;
    var m     = 0;
    var u_mm2 = 0.0;
    var u_mm1 = 0.0;
    var u_m   = 0.0;
    var v_mm2 = 0.0;
    var v_mm1 = 0.0;
    var v_m   = 0.0;
    var f_sum = 0.0;
    var g_sum = 0.0;
    var tmp1  = 0.0;
    var m2    = 0.0;
    var Cmhalf    = { val: 0.0, err: 0.0 }; // Result;
    var rpsi_1pe  = { val: 0.0, err: 0.0 }; // Result;
    var rpsi_1p2e = { val: 0.0, err: 0.0 }; // Result;

    Cmhalf = CLeta( -0.5, eta );
    u_mm2 = 1.0;                      // u_0
    u_mm1 = tex * u_mm2;              // u_1
    m = 2;
  
    rpsi_1pe  = gsl_sf_psi_1piy_e( eta );
    rpsi_1p2e = gsl_sf_psi_1piy_e( 2.0 * eta );
  
    v_mm2 = 2.0 * M_EULER - M_LN2 - rpsi_1pe.val + 2.0 * rpsi_1p2e.val;
    v_mm1 = tex * (v_mm2 - 2.0 * u_mm2);
  
    f_sum = u_mm2 + u_mm1;
    g_sum = v_mm2 + v_mm1;
  
    while ( m < max_iter )
    {
        m2 = (m) * (m);
        u_m = (tex * u_mm1 - x2 * u_mm2) / m2;
        v_m = (tex * v_mm1 - x2 * v_mm2 - 2.0 * (m) * u_m) / m2;
        f_sum = f_sum + u_m;
        g_sum = g_sum + v_m;
        if ( f_sum != 0.0 && g_sum != 0.0 && (Math.abs( u_m / f_sum ) + Math.abs( v_m / g_sum ) < 10.0 * GSL_DBL_EPSILON) ) break;
        u_mm2 = u_mm1;
        u_mm1 = u_m;
        v_mm2 = v_mm1;
        v_mm1 = v_m;
        m = m + 1;
    }
    
    F.val = Cmhalf.val * rx * f_sum;
    F.err = Cmhalf.err * Math.abs( rx * f_sum ) + 2.0 * GSL_DBL_EPSILON * Math.abs( F.val );
  
    tmp1  = f_sum * Math.log( x );
    G.val = -rx * (tmp1 + g_sum) / Cmhalf.val;
    G.err = Math.abs( rx ) * (Math.abs( tmp1 ) + Math.abs( g_sum )) / Math.abs( Cmhalf.val ) * Math.abs( Cmhalf.err / Cmhalf.val );
  
    if ( m >= max_iter )
    {
        throw "SF.MaxIterationsException";
    }

} // coulomb_FGmhalf_series

// ----------------------------------------------------------------------------

// Evolve the backwards recurrence for F,F'.
//
//    F_{lam-1}  = (S_lam F_lam + F_lam') / R_lam
//    F_{lam-1}' = (S_lam F_{lam-1} - R_lam F_lam)
// where
//    R_lam = sqrt(1 + (eta/lam)^2)
//    S_lam = lam/x + eta/lam
//
//
function coulomb_F_recur
    (
    lam_min,    // IN     Double
    kmax,       // IN     Integer
    eta,        // IN     Double
    x,          // IN     Double
    F_lam_max,  // IN     Double
    Fp_lam_max, // IN     Double
    F_lam_min,  // IN OUT Double
    Fp_lam_min  // IN OUT Double
    )
{
    var x_inv   = 0.0;
    var fcl     = 0.0;
    var fpl     = 0.0;
    var lam_max = 0.0;
    var lam     = 0.0;
    var el      = 0.0;
    var rl      = 0.0;
    var sl      = 0.0;
    var fc_lm1  = 0.0;

    x_inv   = 1.0 / x;
    fcl     = F_lam_max;
    fpl     = Fp_lam_max;
    lam_max = lam_min + (kmax);
    lam     = lam_max;
  
    for ( let k = kmax - 1; k >= 0; k-- )
    {
        el     = eta / lam;
        rl     = gsl_sf_hypot( 1.0, el );
        sl     = el  + lam * x_inv;
        fc_lm1 = (fcl * sl + fpl) / rl;
        fpl    =  fc_lm1 * sl - fcl * rl;
        fcl    =  fc_lm1;
        lam    = lam - 1.0;
    }
  
    F_lam_min.Double  = fcl;
    Fp_lam_min.Double = fpl;  

} // coulomb_F_recur

// ----------------------------------------------------------------------------

// Evolve the forward recurrence for G,G'.
//
//   G_{lam+1}  = (S_lam G_lam - G_lam')/R_lam
//   G_{lam+1}' = R_{lam+1} G_lam - S_lam G_{lam+1}
//
// where S_lam and R_lam are as above in the F recursion.
//
function coulomb_G_recur
    (
    lam_min,    // IN     Double
    kmax,       // IN     Integer
    eta,        // IN     Double
    x,          // IN     Double
    G_lam_min,  // IN     Double
    Gp_lam_min, // IN     Double
    G_lam_max,  // IN OUT Double
    Gp_lam_max  // IN OUT Double
    )
{
    var x_inv = 0.0;
    var gcl   = 0.0;
    var gpl   = 0.0;
    var lam   = 0.0;
    var el    = 0.0;
    var rl    = 0.0;
    var sl    = 0.0;
    var gcl1  = 0.0;

    x_inv = 1.0 / x;
    gcl   = G_lam_min;
    gpl   = Gp_lam_min;
    lam   = lam_min + 1.0;
  
    for ( let k = 1; k <= kmax; k++ )
    {
        el   = eta / lam;
        rl   = gsl_sf_hypot( 1.0, el );
        sl   = el + lam * x_inv;
        gcl1 = (sl * gcl - gpl) / rl;
        gpl  = rl * gcl - sl * gcl1;
        gcl  = gcl1;
        lam  = lam + 1.0;
    }
    
    G_lam_max.Double  = gcl;
    Gp_lam_max.Double = gpl;

} // coulomb_G_recur

// ----------------------------------------------------------------------------

// Evaluate the first continued fraction, giving
// the ratio F'/F at the upper lambda value.
// We also determine the sign of F at that point,
// since it is the sign of the last denominator
// in the continued fraction.
//
function coulomb_CF1
        (
        lambda,   // IN     Double
        eta,      // IN     Double
        x,        // IN     Double
        fcl_sign, // IN OUT Double
        r,        // IN OUT Double
        count     // IN OUT Integer
        )
{
    var CF1_small = 1.0e-30;
    var CF1_abort = 1.0e+05;
    var CF1_acc   = 2.0 * GSL_DBL_EPSILON;
    var x_inv     = 1.0 / x;
    var px        = lambda + 1.0 + CF1_abort;
  
    var pk = lambda + 1.0;
    var F  = eta / pk + pk * x_inv;
    var D  = 0.0;
    var C  = 0.0;
    var df = 0.0;
 
    var pk1 = 0.0;
    var ek  = 0.0;
    var rk2 = 0.0;
    var tk  = 0.0;
  
    fcl_sign.Double = 1.0;
    count.Integer = 0;
  
    if ( Math.abs( F ) < CF1_small )
    {
        F = CF1_small;
    }
    D = 0.0;
    C = F;
  
    while ( true )
    {
        pk1 = pk + 1.0;
        ek  = eta / pk;
        rk2 = 1.0 + ek * ek;
        tk  = (pk + pk1) * (x_inv + ek / pk1);
        D   =  tk - rk2 * D;
        C   =  tk - rk2 / C;
        if ( Math.abs( C ) < CF1_small )
        {
            C = CF1_small;
        }
        if ( Math.abs( D ) < CF1_small )
        {
            D = CF1_small;
        }
        D = 1.0 / D;
        df = D * C;
        F  = F * df;
        if ( D < 0.0 )
        {
            // sign of result depends on sign of denominator
            fcl_sign.Double = -fcl_sign.Double;
        }
        pk = pk1;
        if ( pk > px )
        {
            //r = F;
            //GSL_ERROR("error", GSL_ERUNAWAY);
            throw "SF.IterationException";
        }
        count.Integer = count.Integer + 1;
        if ( ! (Math.abs( df - 1.0 ) > CF1_acc) ) break;
    }
    
    r.Double = F;

} // coulomb_CF1


// ----------------------------------------------------------------------------

// Evaluate the second continued fraction to 
// obtain the ratio
//    (G' + i F')/(G + i F) = P + i Q
// at the specified lambda value.
//
function coulomb_CF2
    (
    lambda, // in Double
    eta,    // in Double
    x,      // in Double
    result_P, // IN OUT Double
    result_Q, // IN OUT Double
    count     // IN OUT INTEGER
    )
{
    const CF2_acc   = 4.0 * GSL_DBL_EPSILON;
    const CF2_abort = 2.0e+05;
    const wi    = 2.0 * eta;
    const x_inv = 1.0 / x;
    const e2mm1 = eta * eta + lambda * (lambda + 1.0);
    var ar = 0.0;
    var ai = 0.0;
    var br = 0.0;
    var bi = 0.0;
    var dr = 0.0;
    var di = 0.0;
    var dp = 0.0;
    var dq = 0.0;
    var pk = 0.0;
    var P  = 0.0;
    var Q  = 0.0;
    var A  = 0.0;
    var B  = 0.0;
    var C  = 0.0;
    var D  = 0.0;

    ar = -e2mm1;
    ai =  eta;
  
    br =  2.0 * (x - eta);
    bi =  2.0;
  
    dr =  br / (br * br + bi * bi);
    di = -bi / (br * br + bi * bi);
  
    dp = -x_inv * (ar * di + ai * dr);
    dq =  x_inv * (ar * dr - ai * di);
  
    pk =  0.0;
    P  =  0.0;
    Q  =  1.0 - eta * x_inv;
  
    count.Integer = 0;
   
    while ( true )
    {
        P = P + dp;
        Q = Q + dq;
        pk = pk + 2.0;
        ar = ar + pk;
        ai = ai + wi;
        bi = bi + 2.0;
        D  = ar * dr - ai * di + br;
        di = ai * dr + ar * di + bi;
        C  = 1.0 / (D * D + di * di);
        dr =  C * D;
        di = -C * di;
        A  = br * dr - bi * di - 1.0;
        B  = bi * dr + br * di;
        C  = dp * A  - dq * B;
        dq = dp * B  + dq * A;
        dp = C;
        if ( pk > CF2_abort )
        {
            throw "SF.IterationException";
            //EXIT;
        }
        count.Integer = count.Integer + 1;
        if ( ! (Math.abs( dp ) + Math.abs( dq ) > (Math.abs( P ) + Math.abs( Q )) * CF2_acc) ) break;
    }
  
    if ( Q < CF2_abort * GSL_DBL_EPSILON * Math.abs( P ) )
    {
        throw "SF.AccuracyLossException";
    }
  
    result_P.Double = P;
    result_Q.Double = Q;

} // coulomb_CF2

// ----------------------------------------------------------------------------

// WKB evaluation of F, G. Assumes  0 < x < turning point.
// Overflows are trapped, GSL_EOVRFLW is signalled,
// and an exponent is returned such that:
//
//   result_F = fjwkb * exp(-exponent)
//   result_G = gjwkb * exp( exponent)
//
// See [Biedenharn et al. Phys. Rev. 97, 542-554 (1955), Section IV]
//
// Unfortunately, this is not very accurate in general. The
// test cases typically have 3-4 digits of precision. One could
// argue that this is ok for general use because, for instance,
// F is exponentially small in this region and so the absolute
// accuracy is still roughly acceptable. But it would be better
// to have a systematic method for improving the precision. See
// the Abad+Sesma method discussion below.
//
function coulomb_jwkb(
    lam, // in Double
    eta, // in Double
    x,   // in Double
    fjwkb, // IN OUT Result
    gjwkb, // IN OUT Result
    exponent // IN OUT Double
    )
{
    const llp1      = lam * (lam + 1.0) + 6.0 / 35.0;
    const llp1_eff  = Math.max( llp1, 0.0 );
    const rho_ghalf = Math.sqrt( x * (2.0 * eta - x) + llp1_eff );
    const sinh_arg  = Math.sqrt( llp1_eff / (eta * eta + llp1_eff) ) * rho_ghalf / x;
    const sinh_inv  = Math.log( sinh_arg + gsl_sf_hypot( 1.0, sinh_arg ) );
    const phi       = Math.abs(rho_ghalf - eta * Math.atan2( rho_ghalf, x - eta ) - Math.sqrt( llp1_eff ) * sinh_inv);
    const zeta_half = (3.0 * phi / 2.0) ** (1.0 / 3.0);
    const prefactor = Math.sqrt( M_PI * phi * x / (6.0 * rho_ghalf) );
    var F         = 0.0;
    var G         = 0.0;
    var F_exp     = 0.0;
    var G_exp     = 0.0;
    var airy_scale_exp = 0.0;
    var ai = { val: 0.0, err: 0.0 }; // Result;
    var bi = { val: 0.0, err: 0.0 }; // Result;

    F = prefactor * 3.0 / zeta_half;
    G = prefactor * 3.0 / zeta_half; // Note the sqrt(3) from Bi normalization
    
    airy_scale_exp = phi;
    ai = gsl_sf_airy_Ai_scaled_e( zeta_half * zeta_half, GSL_MODE_DEFAULT );
    bi = gsl_sf_airy_Bi_scaled_e( zeta_half * zeta_half, GSL_MODE_DEFAULT );
    F = F * ai.val;
    G = G * bi.val;
    F_exp = Math.log( F ) - airy_scale_exp;
    G_exp = Math.log( G ) + airy_scale_exp;
  
    if ( G_exp >= GSL_LOG_DBL_MAX )
    {
        //fjwkb.val = F;
        //gjwkb.val = G;
        //fjwkb.err = 1.0e-3 * Math.abs(F); // FIXME: real error here ... could be smaller
        //gjwkb.err = 1.0e-3 * Math.abs(G);
        //exponent  = airy_scale_exp;
        throw "SF.OverflowException";
    }
    else
    {
        fjwkb.val = Math.exp( F_exp );
        gjwkb.val = Math.exp( G_exp );
        fjwkb.err = 1.0e-3 * Math.abs( fjwkb.val );
        gjwkb.err = 1.0e-3 * Math.abs( gjwkb.val );
        exponent.Double  = 0.0;
    }

} // coulomb_jwkb


//*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_coulomb_wave_FG_e(
          eta,
          x,
          lam_F,
          k_lam_G, // lam_G = lam_F - k_lam_G
          F,       // IN OUT Result
          Fp,      // IN OUT Result
          G,       // IN OUT Result
          Gp,      // IN OUT Result
          exp_F,   // IN OUT Double
          exp_G    // IN OUT Double
          )
{

    var lam_G = 0.0;

    lam_G = lam_F - (k_lam_G);

    if ( x < 0.0 || lam_F <= -0.5 || lam_G <= -0.5 )
    {
        throw "SF.DomainException";
    }
    else if ( x == 0.0 )
    {
        throw "SF.DomainException";
        // After all, since we are asking for G, this is a domain error...
    }
    else if ( x < 1.2 && 2.0 * M_PI * eta < 0.9 * (-GSL_LOG_DBL_MIN) && Math.abs( eta * x ) < 10.0 )
    {
        // Reduce to a small lambda value and use the series
        // representations for F and G. We cannot allow eta to
        // be large and positive because the connection formula
        // for G_lam is badly behaved due to an underflow in sin(phi_lam) 
        // [see coulomb_FG_series() and coulomb_connection() above].
        // Note that large negative eta is ok however.
        //
        const SMALL = GSL_SQRT_DBL_EPSILON;
        const N     = Math.trunc( lam_F + 0.5 );
        const span  = Math.max( k_lam_G, N );
        let lam_min           = lam_F - (N);    // -1/2 <= lam_min < 1/2
        let F_lam_F           = 0.0;
        let Fp_lam_F          = 0.0;
        let G_lam_G           = { Double: 0.0 };
        let Gp_lam_G          = { Double: 0.0 };
        let F_lam_F_err       = 0.0;
        let Fp_lam_F_err      = 0.0;
        let Fp_over_F_lam_F   = { Double: 0.0 };
        let F_sign_lam_F      = { Double: 0.0 };
        let F_lam_min_unnorm  = { Double: 0.0 };
        let Fp_lam_min_unnorm = { Double: 0.0 };
        let Fp_over_F_lam_min = 0.0;
        let F_scale           = 0.0;
        let Gerr_frac         = 0.0;
        let F_scale_frac_err  = 0.0;
        let F_unnorm_frac_err = 0.0;
        let F_lam_min  = { val: 0.0, err: 0.0 }; // Result;
        let G_lam_min  = { val: 0.0, err: 0.0 }; // Result;
        let Gp_lam_min = { val: 0.0, err: 0.0 }; // Result;
           
        // Determine F'/F at lam_F.
        let CF1_count = { Integer: 0 };

        coulomb_CF1( lam_F, eta, x, F_sign_lam_F, Fp_over_F_lam_F, CF1_count ); //*****!!!
        // Recurse down with unnormalized F,F' values.
        F_lam_F  = SMALL;
        Fp_lam_F = Fp_over_F_lam_F.Double * F_lam_F;
        if ( span != 0 )
        {
            coulomb_F_recur( lam_min, span, eta, x, F_lam_F, Fp_lam_F, F_lam_min_unnorm, Fp_lam_min_unnorm );
        }
        else
        {
            F_lam_min_unnorm.Double  =  F_lam_F;
            Fp_lam_min_unnorm.Double = Fp_lam_F;
        }
        
        // Determine F and G at lam_min.
        if ( lam_min == -0.5 )
        {
            coulomb_FGmhalf_series( eta, x, F_lam_min, G_lam_min );
        }
        else if ( lam_min == 0.0 )
        {
            coulomb_FG0_series( eta, x, F_lam_min, G_lam_min );
        }
        else if ( lam_min == 0.5 )
        {
            // This cannot happen.
            throw "SF.SanityException";
        }
        else
        {
            coulomb_FG_series( lam_min, eta, x, F_lam_min, G_lam_min );
        }
        
        // Determine remaining quantities.
        Fp_over_F_lam_min = Fp_lam_min_unnorm.Double / F_lam_min_unnorm.Double;
        Gp_lam_min.val  = Fp_over_F_lam_min * G_lam_min.val - 1.0 / F_lam_min.val;
        Gp_lam_min.err  = Math.abs( Fp_over_F_lam_min ) * G_lam_min.err;
        Gp_lam_min.err = Gp_lam_min.err + Math.abs( 1.0 / F_lam_min.val ) * Math.abs( F_lam_min.err / F_lam_min.val );
        F_scale     = F_lam_min.val / F_lam_min_unnorm.Double;
        
        // Apply scale to the original F,F' values.
        F_scale_frac_err  = Math.abs( F_lam_min.err / F_lam_min.val );
        F_unnorm_frac_err = 2.0 * GSL_DBL_EPSILON * (CF1_count.Integer + span + 1);
        F_lam_F      = F_lam_F * F_scale;
        F_lam_F_err  = Math.abs( F_lam_F ) * (F_unnorm_frac_err + F_scale_frac_err);
        Fp_lam_F     = Fp_lam_F * F_scale;
        Fp_lam_F_err = Math.abs( Fp_lam_F ) * (F_unnorm_frac_err + F_scale_frac_err);
        
        // Recurse up to get the required G,G' values.
        coulomb_G_recur( lam_min, Math.max( N - k_lam_G, 0 ), eta, x,
                                    G_lam_min.val, Gp_lam_min.val,
                                    G_lam_G, Gp_lam_G );
        
        F.val = F_lam_F;
        F.err = F_lam_F_err;
        F.err = F.err + 2.0 * GSL_DBL_EPSILON * Math.abs( F_lam_F );
        
        Fp.val = Fp_lam_F;
        Fp.err = Fp_lam_F_err;
        Fp.err = Fp.err + 2.0 * GSL_DBL_EPSILON * Math.abs( Fp_lam_F );
        
        Gerr_frac = Math.abs( G_lam_min.err / G_lam_min.val ) + Math.abs( Gp_lam_min.err / Gp_lam_min.val );
        
        G.val = G_lam_G.Double;
        G.err = Gerr_frac * Math.abs( G_lam_G.Double );
        G.err = G.err + 2.0 * (CF1_count.Integer + 1) * GSL_DBL_EPSILON * Math.abs( G.val );
        
        Gp.val = Gp_lam_G.Double;
        Gp.err = Gerr_frac * Math.abs( Gp.val );
        Gp.err = Gp.err + 2.0 * (CF1_count.Integer + 1) * GSL_DBL_EPSILON * Math.abs( Gp.val );
        
        exp_F.Double = 0.0;
        exp_G.Double = 0.0;
        
        return;
    }
    else if ( x < 2.0 * eta )
    {
        // Use WKB approximation to obtain F and G at the two
        // lambda values, and use the Wronskian and the
        // continued fractions for F'/F to obtain F' and G'.
        //
        let exp_lam_F       = { Double: 0.0 };
        let exp_lam_G       = { Double: 0.0 };
        let Fp_over_F_lam_F = { Double: 0.0 };
        let Fp_over_F_lam_G = { Double: 0.0 };
        let F_sign_lam_F    = { Double: 0.0 };
        let F_sign_lam_G    = { Double: 0.0 };
        let F_lam_F = { val: 0.0, err: 0.0 }; // Result;
        let G_lam_F = { val: 0.0, err: 0.0 }; // Result;
        let F_lam_G = { val: 0.0, err: 0.0 }; // Result;
        let G_lam_G = { val: 0.0, err: 0.0 }; // Result;

        let CF1_count = { Integer: 0 };

        coulomb_jwkb( lam_F, eta, x, F_lam_F, G_lam_F, exp_lam_F );
        if ( k_lam_G == 0 )
        {
            F_lam_G = F_lam_F;
            G_lam_G = G_lam_F;
            exp_lam_G = exp_lam_F;
        }
        else
        {
            coulomb_jwkb( lam_G, eta, x, F_lam_G, G_lam_G, exp_lam_G );
        }
        
        coulomb_CF1( lam_F, eta, x, F_sign_lam_F, Fp_over_F_lam_F, CF1_count );
        if ( k_lam_G == 0 )
        {
            F_sign_lam_G    = F_sign_lam_F;
            Fp_over_F_lam_G = Fp_over_F_lam_F;
        }
        else
        {
            coulomb_CF1( lam_G, eta, x, F_sign_lam_G, Fp_over_F_lam_G, CF1_count );
        }
        
        F.val = F_lam_F.val;
        F.err = F_lam_F.err;
        
        G.val = G_lam_G.val;
        G.err = G_lam_G.err;
        
        Fp.val = Fp_over_F_lam_F.Double * F_lam_F.val;
        Fp.err = Math.abs( Fp_over_F_lam_F.Double ) * F_lam_F.err;
        Fp.err = Fp.err + 2.0 * GSL_DBL_EPSILON * Math.abs( Fp.val );
        
        Gp.val = Fp_over_F_lam_G.Double * G_lam_G.val - 1.0 / F_lam_G.val;
        Gp.err = Math.abs( Fp_over_F_lam_G.Double ) * G_lam_G.err;
        Gp.err = Gp.err + Math.abs( 1.0 / F_lam_G.val ) * Math.abs( F_lam_G.err / F_lam_G.val );
        
        exp_F.Double = exp_lam_F.Double;
        exp_G.Double = exp_lam_G.Double;
        
        return;
    }
    else
    {
        // x > 2 eta, so we know that we can find a lambda value such
        // that x is above the turning point. We do this, evaluate
        // using Steed's method at that oscillatory point, then
        // use recursion on F and G to obtain the required values.
        //
        // lam_0   = a value of lambda such that x is below the turning point
        // lam_min = minimum of lam_0 and the requested lam_G, since
        //           we must go at least as low as lam_G
        //
        const SMALL   = GSL_SQRT_DBL_EPSILON;
        const C       = Math.sqrt( 1.0 + 4.0 * x * (x - 2.0 * eta) );
        const N       = Math.trunc( Math.ceil( lam_F - C + 0.5 ) );
        const lam_0   = lam_F - (Math.max( N, 0 ));
        const lam_min = Math.min( lam_0, lam_G );
        let F_lam_F           = 0.0;
        let Fp_lam_F          = 0.0;
        let G_lam_G           = { Double: 0.0 };
        let Gp_lam_G          = { Double: 0.0 };
        let F_lam_min_unnorm  = { Double: 0.0 };
        let Fp_lam_min_unnorm = { Double: 0.0 };
        let F_lam_min         = 0.0;
        let Fp_lam_min        = 0.0;
        let G_lam_min         = 0.0;
        let Gp_lam_min        = 0.0;
        let Fp_over_F_lam_F   = { Double: 0.0 };
        let Fp_over_F_lam_min = 0.0;
        let F_sign_lam_F      = { Double: 0.0 };
        let F_sign_lam_min    = 0.0;
        let P_lam_min         = { Double: 0.0 };
        let Q_lam_min         = { Double: 0.0 };
        let alpha             = 0.0;
        let gamma             = 0.0;
        let F_scale           = 0.0;
        
        let CF1_count = { Integer: 0 };
        let CF2_count = { Integer: 0 };
        
        let F_recur_count = 0;
        let G_recur_count = 0;
        
        let err_amplify = 0.0;

        coulomb_CF1( lam_F, eta, x, F_sign_lam_F, Fp_over_F_lam_F, CF1_count );
        F_lam_F  = F_sign_lam_F.Double * SMALL;  // unnormalized
        Fp_lam_F = Fp_over_F_lam_F.Double * F_lam_F;
        
        // Backward recurrence to get F,Fp at lam_min
        F_recur_count = Math.max( k_lam_G, N );
        coulomb_F_recur( lam_min, F_recur_count, eta, x,
                                    F_lam_F, Fp_lam_F,
                                    F_lam_min_unnorm, Fp_lam_min_unnorm
                                    );
        Fp_over_F_lam_min = Fp_lam_min_unnorm.Double / F_lam_min_unnorm.Double;
        
        // Steed evaluation to complete evaluation of F,Fp,G,Gp at lam_min
        coulomb_CF2( lam_min, eta, x, P_lam_min, Q_lam_min, CF2_count );
        alpha = Fp_over_F_lam_min - P_lam_min.Double;
        gamma = alpha / Q_lam_min.Double;
        
        F_sign_lam_min = GSL_SIGN( F_lam_min_unnorm.Double );
        
        F_lam_min  = F_sign_lam_min / Math.sqrt( alpha * alpha / Q_lam_min.Double + Q_lam_min.Double );
        Fp_lam_min = Fp_over_F_lam_min * F_lam_min;
        G_lam_min  = gamma * F_lam_min;
        Gp_lam_min = (P_lam_min.Double * gamma - Q_lam_min.Double) * F_lam_min;
        
        // Apply scale to values of F,Fp at lam_F (the top).
        F_scale  = F_lam_min / F_lam_min_unnorm.Double;    
        F_lam_F  = F_lam_F * F_scale;
        Fp_lam_F = Fp_lam_F * F_scale;
        
        // Forward recurrence to get G,Gp at lam_G (the top).
        G_recur_count = Math.max( N - k_lam_G, 0 );
        coulomb_G_recur( lam_min, G_recur_count, eta, x,
                                    G_lam_min, Gp_lam_min,
                                    G_lam_G, Gp_lam_G
                                    );
        
        err_amplify = (CF1_count.Integer + CF2_count.Integer + F_recur_count + G_recur_count + 1);
        
        F.val  = F_lam_F;
        F.err  = 8.0 * err_amplify * GSL_DBL_EPSILON * Math.abs( F.val );
        
        Fp.val = Fp_lam_F;
        Fp.err = 8.0 * err_amplify * GSL_DBL_EPSILON * Math.abs( Fp.val );
        
        G.val  = G_lam_G.Double;
        G.err  = 8.0 * err_amplify * GSL_DBL_EPSILON * Math.abs( G.val );
        
        Gp.val = Gp_lam_G.Double;
        Gp.err = 8.0 * err_amplify * GSL_DBL_EPSILON * Math.abs( Gp.val );
        
        exp_F.Double = 0.0;
        exp_G.Double = 0.0;
        
        return;
    }

} // gsl_sf_coulomb_wave_FG_e

// ----------------------------------------------------------------------------

export function gsl_sf_coulomb_wave_F_array
    (
    lam_min,
    kmax,
    eta,
    x, 
    fc_array, //OUT array
    F_exp     // OUT Double
    )
{

    if ( x == 0.0 )
    {
        F_exp.Double = 0.0;
        for ( let k = 0; k <= kmax; k++ )
        {
            fc_array[k] = 0.0;
        }
        if ( lam_min == 0.0 )
        {
            let f_0 = { val: 0.0, err: 0.0 };
            f_0 = CLeta( 0.0, eta);
            fc_array[0] = f_0.val;
        }
        return;
    }
    else
    {
        const x_inv = 1.0 / x;
        const lam_max = lam_min + kmax;
        let F  = { val: 0.0, err: 0.0 };
        let Fp = { val: 0.0, err: 0.0 };
        let G  = { val: 0.0, err: 0.0 };
        let Gp = { val: 0.0, err: 0.0 };
        let G_exp = { Double: 0.0 };

        gsl_sf_coulomb_wave_FG_e( eta, x, lam_max, 0, F, Fp, G, Gp, F_exp, G_exp );

        let fcl  = F.val;
        let fpl = Fp.val;
        let lam = lam_max;

        fc_array[kmax] = F.val;

        for ( let k = kmax - 1; k >= 0; k-- )
        {
            let el = eta / lam;
            let rl = gsl_sf_hypot( 1.0, el );
            let sl = el  + lam * x_inv;
            let fc_lm1 = (fcl * sl + fpl) / rl;
            fc_array[k] = fc_lm1;
            fpl         = fc_lm1 * sl - fcl * rl;
            fcl         = fc_lm1;
            lam -= 1.0;
        }

        return;
    }

}

// ----------------------------------------------------------------------------

export function gsl_sf_coulomb_wave_FG_array
    (
    lam_min,
    kmax,
    eta,
    x,
    fc_array, // OUT array
    gc_array, // OUT array
    F_exp,    // OUT Double
    G_exp     // OUT Double
    )
{

    const x_inv = 1.0 / x;
    const lam_max = lam_min + kmax;
    var F  = { val: 0.0, err: 0.0 };
    var Fp = { val: 0.0, err: 0.0 };
    var G  = { val: 0.0, err: 0.0 };
    var Gp = { val: 0.0, err: 0.0 };

    gsl_sf_coulomb_wave_FG_e( eta, x, lam_max, kmax, F, Fp, G, Gp, F_exp, G_exp );

    var fcl = F.val;
    var fpl = Fp.val;
    var lam = lam_max;

    var gcl = 0.0;
    var gpl = 0.0;

    fc_array[kmax] = F.val;

    for ( let k = kmax - 1; k >= 0; k-- )
    {
        let el = eta / lam;
        let rl = gsl_sf_hypot( 1.0, el );
        let sl = el + lam * x_inv;
        let fc_lm1;
        fc_lm1 = (fcl * sl + fpl) / rl;
        fc_array[k] = fc_lm1;
        fpl         = fc_lm1 * sl - fcl * rl;
        fcl         = fc_lm1;
        lam -= 1.0;
    }

    gcl = G.val;
    gpl = Gp.val;
    lam = lam_min + 1.0;

    gc_array[0] = G.val;

    for ( let k = 1; k <= kmax; k++ )
    {
        let el = eta / lam;
        let rl = gsl_sf_hypot( 1.0, el );
        let sl = el + lam * x_inv;
        let gcl1 = (sl * gcl - gpl) / rl;
        gc_array[k] = gcl1;
        gpl         = rl * gcl - sl * gcl1;
        gcl         = gcl1;
        lam += 1.0;
    }

}

// ----------------------------------------------------------------------------

export function gsl_sf_coulomb_wave_FGp_array
    (
    lam_min,
    kmax,
    eta,
    x,
    fc_array,  // OUT array
    fcp_array, // OUT array
    gc_array,  // OUT array
    gcp_array, // OUT array
    F_exp,     // OUT Double
    G_exp      // OUT Double
    )
{

    const x_inv = 1.0 / x;
    const lam_max = lam_min + kmax;
    var F  = { val: 0.0, err: 0.0 };
    var Fp = { val: 0.0, err: 0.0 };
    var G  = { val: 0.0, err: 0.0 };
    var Gp = { val: 0.0, err: 0.0 };

    gsl_sf_coulomb_wave_FG_e( eta, x, lam_max, kmax, F, Fp, G, Gp, F_exp, G_exp );

    var fcl = F.val;
    var fpl = Fp.val;
    var lam = lam_max;

    var gcl = 0.0;
    var gpl = 0.0;

    fc_array[kmax]  = F.val;
    fcp_array[kmax] = Fp.val;

    for ( let k = kmax - 1; k >= 0; k-- )
    {
        let el = eta / lam;
        let rl = gsl_sf_hypot( 1.0, el );
        let sl = el  + lam * x_inv;
        let fc_lm1 = 0.0;
        fc_lm1 = (fcl * sl + fpl) / rl;
        fc_array[k]  = fc_lm1;
        fpl          = fc_lm1 * sl - fcl * rl;
        fcp_array[k] = fpl;
        fcl          = fc_lm1;
        lam -= 1.0;
    }

    gcl = G.val;
    gpl = Gp.val;
    lam = lam_min + 1.0;

    gc_array[0]  = G.val;
    gcp_array[0] = Gp.val;

    for ( let k = 1; k <= kmax; k++ )
    {
        let el = eta / lam;
        let rl = gsl_sf_hypot( 1.0, el );
        let sl = el + lam * x_inv;
        let gcl1 = (sl * gcl - gpl) / rl;
        gc_array[k]  = gcl1;
        gpl          = rl * gcl - sl * gcl1;
        gcp_array[k] = gpl;
        gcl          = gcl1;
        lam += 1.0;
    }

}

// ----------------------------------------------------------------------------

export function gsl_sf_coulomb_wave_sphF_array
    (
    lam_min,
    kmax,
    eta,
    x,
    fc_array, // OUT array
    F_exp     // OUT Double
    )
{

    if ( x < 0.0 || lam_min < -0.5 )
    {
        throw "SF.DomainException";
    }
    else if ( x < 10.0 / GSL_DBL_MAX )
    {
        for ( let k = 0; k <= kmax; k++ )
        {
            fc_array[k] = 0.0;
        }
        if ( lam_min == 0.0 )
        {
            fc_array[0] = Math.sqrt( C0sq( eta ) );
        }
        F_exp.Double = 0.0;
        if ( x == 0.0 )
        {
            return;
        }
        else
        {
            throw "SF.UnderflowException"; // ***!!!
        }
    }
    else
    {
        gsl_sf_coulomb_wave_F_array( lam_min, kmax, eta, x, fc_array, F_exp );
        for ( let k = 0; k <= kmax; k++ )
        {
            fc_array[k] = fc_array[k] / x;
        }
        return;
    }

}

// ----------------------------------------------------------------------------
// EOF SF-Coulomb.mjs
