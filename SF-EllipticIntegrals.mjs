// SF-Trigonometric.mjs
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

import { GSL_DBL_MIN }          from "./SF-Machine.mjs";
import { GSL_DBL_MAX }          from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { M_PI }                 from "./SF-Math.mjs";
import { GSL_SIGN }             from "./SF-Math.mjs";
import { GSL_MODE_PREC }        from "./SF-Mode.mjs";
import { GSL_PREC_DOUBLE }      from "./SF-Mode.mjs";
import { gsl_prec_eps }         from "./SF-Mode.mjs";

import { EVAL_RESULT_DM }       from "./SF-Evaluate.mjs";
import { EVAL_RESULT_DDM }      from "./SF-Evaluate.mjs";
import { EVAL_RESULT_3DM }      from "./SF-Evaluate.mjs";
import { EVAL_RESULT_4DM }      from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// *-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*

function locMAX3(x, y, z)
{
    return Math.max(Math.max(x, y), z);
} // locMAX3

function locMAX4(x, y, z, w)
{
    return Math.max(Math.max(Math.max(x, y), z), w);
} // locMAX4


// *-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*


// based on Carlson's algorithms:
// [B. C. Carlson Numer. Math. 33, 1 (1979)]
// 
// see also:
// [B.C. Carlson, Special Functions of Applied Mathematics (1977)]
//

// According to Carlson's algorithm, the errtol parameter
// typically effects the relative error in the following way:
//
// relative error < 16 errtol^6 / (1 - 2 errtol)
//
//   errtol     precision
//   ------     ----------
//   0.001       1.0e-17
//   0.003       2.0e-14 
//   0.01        2.0e-11
//   0.03        2.0e-8
//   0.1         2.0e-5
//


export function gsl_sf_ellint_RC_e(x, y, mode)
{
    const lolim  = 5.0 * GSL_DBL_MIN;
    const uplim  = 0.2 * GSL_DBL_MAX;
    var goal   = GSL_MODE_PREC(mode);
    var errtol = 0.0;
    var prec   = 0.0;
    var r      = { val: 0.0, err: 0.0 }; // Result;

    if (goal == GSL_MODE_PREC(GSL_PREC_DOUBLE))
    {
        errtol = 0.001;
    }
    else
    {
        errtol = 0.03;
    }
    prec = gsl_prec_eps[goal];

    if (x < 0.0 || y < 0.0 || x + y < lolim)
    {
        throw "SF.DomainException";
    }
    else if (Math.max(x, y) < uplim) 
    {
        const c1 = 1.0 / 7.0;
        const c2 = 9.0 / 22.0;
        var xn = x;
        var yn = y;
        var mu = 0.0;
        var sn = 0.0;
        var lamda = 0.0;
        var s  = 0.0;

        while (true)
        {
            mu = (xn + yn + yn) / 3.0;
            sn = (yn + mu) / mu - 2.0;
            if (Math.abs(sn) < errtol) break;
            lamda = 2.0 * Math.sqrt(xn) * Math.sqrt(yn) + yn;
            xn = (xn + lamda) * 0.25;
            yn = (yn + lamda) * 0.25;
        }
        s = sn * sn * (0.3 + sn * (c1 + sn * (0.375 + sn * c2)));
        r.val = (1.0 + s) / Math.sqrt(mu);
        r.err = prec * Math.abs(r.val);
    }
    else
    {
        throw "SF.DomainException";
    }

    return r;

} // gsl_sf_ellint_RC_e

// ----------------------------------------------------------------------------

export function gsl_sf_ellint_RD_e(x, y, z, mode)
{

    const lolim  = 2.0 / (GSL_DBL_MAX ** (2.0 / 3.0));
    var uplim  = 0.0;
    var goal   = GSL_MODE_PREC(mode);
    var errtol = 0.0;
    var prec   = 0.0;
    var r      = { val: 0.0, err: 0.0 }; // Result;

    if (goal == GSL_MODE_PREC(GSL_PREC_DOUBLE))
    {
        errtol = 0.001;
    }
    else
    {
        errtol = 0.03;
    }
    prec = gsl_prec_eps[goal];
    uplim = (0.1 * errtol / GSL_DBL_MIN) ** (2.0 / 3.0);

    if (Math.min(x, y) < 0.0 || Math.min(x + y, z) < lolim)
    {
        throw "SF.DomainException";
    }
    else if (locMAX3(x,y,z) < uplim)
    {
        const c1 = 3.0 / 14.0;
        const c2 = 1.0 /  6.0;
        const c3 = 9.0 / 22.0;
        const c4 = 3.0 / 26.0;
        var xn = x;
        var yn = y;
        var zn = z;
        var sigma  = 0.0;
        var power4 = 1.0;
        var ea = 0.0;
        var eb = 0.0;
        var ec = 0.0;
        var ed = 0.0;
        var ef = 0.0;
        var s1 = 0.0;
        var s2 = 0.0;
        var mu = 0.0;
        var xndev = 0.0;
        var yndev = 0.0;
        var zndev = 0.0;
        var xnroot = 0.0;
        var ynroot = 0.0;
        var znroot = 0.0;
        var lamda = 0.0;
        var epslon = 0.0;

        while (true)
        {
            mu = (xn + yn + 3.0 * zn) * 0.2;
            xndev = (mu - xn) / mu;
            yndev = (mu - yn) / mu;
            zndev = (mu - zn) / mu;
            epslon = locMAX3(Math.abs(xndev), Math.abs(yndev), Math.abs(zndev));
            if (epslon < errtol) break;
            xnroot = Math.sqrt(xn);
            ynroot = Math.sqrt(yn);
            znroot = Math.sqrt(zn);
            lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
            sigma = sigma + power4 / (znroot * (zn + lamda));
            power4 = power4 * 0.25;
            xn = (xn + lamda) * 0.25;
            yn = (yn + lamda) * 0.25;
            zn = (zn + lamda) * 0.25;
        }
        ea = xndev * yndev;
        eb = zndev * zndev;
        ec = ea - eb;
        ed = ea - 6.0 * eb;
        ef = ed + ec + ec;
        s1 = ed * (- c1 + 0.25 * c3 * ed - 1.5 * c4 * zndev * ef);
        s2 = zndev * (c2 * ef + zndev * (- c3 * ec + zndev * c4 * ea));
        r.val = 3.0 * sigma + power4 * (1.0 + s1 + s2) / (mu * Math.sqrt(mu));
        r.err = prec * Math.abs(r.val);
    }
    else
    {
        throw "SF.DomainException";
    }

    return r;

} // gsl_sf_ellint_RD_e

// ----------------------------------------------------------------------------

export function gsl_sf_ellint_RF_e(x, y, z, mode)
{
    const lolim  = 5.0 * GSL_DBL_MIN;
    const uplim  = 0.2 * GSL_DBL_MAX;
    var goal   = GSL_MODE_PREC(mode);
    var errtol = 0.0;
    var prec   = 0.0;
    var r      = { val: 0.0, err: 0.0 }; // Result;

    if (goal == GSL_MODE_PREC(GSL_PREC_DOUBLE))
    {
        errtol = 0.001;
    }
    else
    {
        errtol = 0.03;
    }
    prec = gsl_prec_eps[goal];

    if (x < 0.0 || y < 0.0 || z < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x + y < lolim || x + z < lolim || y + z < lolim)
    {
        throw "SF.DomainException";
    }
    else if (locMAX3(x, y, z) < uplim)
    {
        const c1 = 1.0 / 24.0;
        const c2 = 3.0 / 44.0;
        const c3 = 1.0 / 14.0;
        var xn = x;
        var yn = y;
        var zn = z;
        var mu = 0.0;
        var xndev = 0.0;
        var yndev = 0.0;
        var zndev = 0.0;
        var e2 = 0.0;
        var e3 = 0.0;
        var s = 0.0;
        var xnroot = 0.0;
        var ynroot = 0.0;
        var znroot = 0.0;
        var lamda = 0.0;
        var epslon = 0.0;

        while (true)
        {
            mu = (xn + yn + zn) / 3.0;
            xndev = 2.0 - (mu + xn) / mu;
            yndev = 2.0 - (mu + yn) / mu;
            zndev = 2.0 - (mu + zn) / mu;
            epslon = locMAX3(Math.abs(xndev), Math.abs(yndev), Math.abs(zndev));
            if (epslon < errtol) break;
            xnroot = Math.sqrt(xn);
            ynroot = Math.sqrt(yn);
            znroot = Math.sqrt(zn);
            lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
            xn = (xn + lamda) * 0.25;
            yn = (yn + lamda) * 0.25;
            zn = (zn + lamda) * 0.25;
        }
        e2 = xndev * yndev - zndev * zndev;
        e3 = xndev * yndev * zndev;
        s = 1.0 + (c1 * e2 - 0.1 - c2 * e3) * e2 + c3 * e3;
        r.val = s / Math.sqrt(mu);
        r.err = prec * Math.abs(r.val);
    }
    else
    {
        throw "SF.DomainException";
    }

    return r;

} // gsl_sf_ellint_RF_e

// ----------------------------------------------------------------------------

export function gsl_sf_ellint_RJ_e(x, y, z, p, mode)
{
    const lolim  = (5.0 * GSL_DBL_MIN) ** (1.0 / 3.0);
    const uplim  = 0.3 * ((0.2 * GSL_DBL_MAX) ** (1.0 / 3.0));
    var goal   = GSL_MODE_PREC(mode);
    var errtol = 0.0;
    var prec   = 0.0;
    var r      = { val: 0.0, err: 0.0 }; // Result;

    if (goal == GSL_MODE_PREC(GSL_PREC_DOUBLE))
    {
        errtol = 0.001;
    }
    else
    {
        errtol = 0.03;
    }
    prec = gsl_prec_eps[goal];

    if (x < 0.0 || y < 0.0 || z < 0.0)
    {
        throw "SF.DomainException";
    }
    else if (x + y < lolim || x + z < lolim || y + z < lolim || p < lolim)
    {
        throw "SF.DomainException";
    }
    else if (locMAX4(x, y, z, p) < uplim)
    {
        const c1 = 3.0 / 14.0;
        const c2 = 1.0 /  3.0;
        const c3 = 3.0 / 22.0;
        const c4 = 3.0 / 26.0;
        var xn = x;
        var yn = y;
        var zn = z;
        var pn = p;
        var sigma = 0.0;
        var power4 = 1.0;
        var mu = 0.0;
        var xndev = 0.0;
        var yndev = 0.0;
        var zndev = 0.0;
        var pndev = 0.0;
        var ea = 0.0;
        var eb = 0.0;
        var ec = 0.0;
        var e2 = 0.0;
        var e3 = 0.0;
        var s1 = 0.0;
        var s2 = 0.0;
        var s3 = 0.0;
        var xnroot = 0.0;
        var ynroot = 0.0;
        var znroot = 0.0;
        var lamda = 0.0;
        var alfa = 0.0;
        var beta = 0.0;
        var epslon = 0.0;
        var rcresult = { val: 0.0, err: 0.0 }; // Result;

        while (true)
        {
            mu = (xn + yn + zn + pn + pn) * 0.2;
            xndev = (mu - xn) / mu;
            yndev = (mu - yn) / mu;
            zndev = (mu - zn) / mu;
            pndev = (mu - pn) / mu;
            epslon = locMAX4(Math.abs(xndev), Math.abs(yndev), Math.abs(zndev), Math.abs(pndev));
            if (epslon < errtol) break;
            xnroot = Math.sqrt(xn);
            ynroot = Math.sqrt(yn);
            znroot = Math.sqrt(zn);
            lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
            alfa = pn * (xnroot + ynroot + znroot) + xnroot * ynroot * znroot;
            alfa = alfa * alfa;
            beta = pn * (pn + lamda) * (pn + lamda);
            rcresult = gsl_sf_ellint_RC_e(alfa, beta, mode);
            sigma  = sigma + power4 * rcresult.val;
            power4 = power4 * 0.25;
            xn = (xn + lamda) * 0.25;
            yn = (yn + lamda) * 0.25;
            zn = (zn + lamda) * 0.25;
            pn = (pn + lamda) * 0.25;
        }
            
        ea = xndev * (yndev + zndev) + yndev * zndev;
        eb = xndev * yndev * zndev;
        ec = pndev * pndev;
        e2 = ea - 3.0 * ec;
        e3 = eb + 2.0 * pndev * (ea - ec);
        s1 = 1.0 + e2 * (- c1 + 0.75 * c3 * e2 - 1.5 * c4 * e3);
        s2 = eb * (0.5 * c2 + pndev * (- c3 - c3 + pndev * c4));
        s3 = pndev * ea * (c2 - pndev * c3) - c2 * pndev * ec;
        r.val = 3.0 * sigma + power4 * (s1 + s2 + s3) / (mu * Math.sqrt(mu));
        r.err = prec * Math.abs(r.val);
    }
    else
    {
        throw "SF.DomainException";
    }

    return r;

} // gsl_sf_ellint_RJ_e

// ----------------------------------------------------------------------------

// [Carlson, Numer. Math. 33 (1979) 1, (4.1)]
export function gsl_sf_ellint_F_e(phi0, k, mode)
{
    var nc       = 0.0;
    var phi_red  = 0.0;
    var sin_phi  = 0.0;
    var sin2_phi = 0.0;
    var x        = 0.0;
    var y        = 0.0;
    var status   = 0;
    var rkstatus = 0;
    var rf  = { val: 0.0, err: 0.0 }; // Result;
    var rk  = { val: 0.0, err: 0.0 }; // Result;
    var r   = { val: 0.0, err: 0.0 }; // Result;
    var phi = phi0;

    // Angular reduction to -pi/2 < phi < pi/2 (we should really use an
    // exact reduction but this will have to do for now) BJG

    nc = Math.floor(phi / M_PI + 0.5);
    phi_red = phi - nc * M_PI;
    phi = phi_red;
  
    sin_phi  = Math.sin(phi);
    sin2_phi = sin_phi * sin_phi;
    x = 1.0 - sin2_phi;
    y = 1.0 - k * k * sin2_phi;
    rf = gsl_sf_ellint_RF_e(x, y, 1.0, mode);
    r.val = sin_phi * rf.val;
    r.err = GSL_DBL_EPSILON * Math.abs(r.val) + Math.abs(sin_phi * rf.err);
    if (nc != 0.0)
    {
        rk = gsl_sf_ellint_Kcomp_e(k, mode);  
        r.val = r.val + 2.0 * nc * rk.val;
        r.err = r.err + 2.0 * Math.abs(nc) * rk.err;
    }

    return r;

} // gsl_sf_ellint_F_e

// ----------------------------------------------------------------------------

// [Carlson, Numer. Math. 33 (1979) 1, (4.2)]
export function gsl_sf_ellint_E_e(phi0, k, mode)
{
    var nc       = 0.0;
    var phi_red  = 0.0;
    var sin_phi  = 0.0;
    var sin2_phi = 0.0;
    var sin3_phi = 0.0;
    var x        = 0.0;
    var y        = 0.0;
    var re = { val: 0.0, err: 0.0 }; // Result;
    var rf = { val: 0.0, err: 0.0 }; // Result;
    var rd = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;
    var phi = phi0;

    // Angular reduction to -pi/2 < phi < pi/2 (we should really use an
    // exact reduction but this will have to do for now) BJG

    nc = Math.floor(phi / M_PI + 0.5);
    phi_red = phi - nc * M_PI;
    phi = phi_red;

    sin_phi  = Math.sin(phi);
    sin2_phi = sin_phi  * sin_phi;
    x = 1.0 - sin2_phi;
    y = 1.0 - k * k * sin2_phi;

    if (x < GSL_DBL_EPSILON)
    {
        re = gsl_sf_ellint_Ecomp_e(k, mode);  
        // could use A&S 17.4.14 to improve the value below
        r.val = 2.0 * nc * re.val + GSL_SIGN(sin_phi) * re.val;
        r.err = 2.0 * Math.abs(nc) * re.err + re.err;
    }
    else
    {
        sin3_phi = sin2_phi * sin_phi;
        rf = gsl_sf_ellint_RF_e(x, y, 1.0, mode);
        rd = gsl_sf_ellint_RD_e(x, y, 1.0, mode);
        r.val = sin_phi * rf.val - k * k/3.0 * sin3_phi * rd.val;
        r.err = GSL_DBL_EPSILON * Math.abs(sin_phi * rf.val);
        r.err = r.err + Math.abs(sin_phi * rf.err);
        r.err = r.err + k * k / 3.0 * GSL_DBL_EPSILON * Math.abs(sin3_phi * rd.val);
        r.err = r.err + k * k / 3.0 * Math.abs(sin3_phi * rd.err);
        if (nc != 0.0)
        {
            re = gsl_sf_ellint_Ecomp_e(k, mode);  
            r.val = r.val + 2.0 * nc * re.val;
            r.err = r.err + 2.0 * Math.abs(nc) * re.err;
        }
    }

    return r;

} // gsl_sf_ellint_E_e

// ----------------------------------------------------------------------------

// [Carlson, Numer. Math. 33 (1979) 1, (4.3)]
export function gsl_sf_ellint_P_e(phi0, k, n, mode)
{
    var nc       = 0.0;
    var phi_red  = 0.0;
    var sin_phi  = 0.0;
    var sin2_phi = 0.0;
    var sin3_phi = 0.0;
    var x        = 0.0;
    var y        = 0.0;
    var rf = { val: 0.0, err: 0.0 }; // Result;
    var rj = { val: 0.0, err: 0.0 }; // Result;
    var rp = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;
    var phi = phi0;

    // Angular reduction to -pi/2 < phi < pi/2 (we should really use an
    // exact reduction but this will have to do for now) BJG

    nc = Math.floor(phi / M_PI + 0.5);
    phi_red = phi - nc * M_PI;
    phi = phi_red;

    // FIXME: need to handle the case of small x, as for E,F

    sin_phi  = Math.sin(phi);
    sin2_phi = sin_phi  * sin_phi;
    sin3_phi = sin2_phi * sin_phi;
    x = 1.0 - sin2_phi;
    y = 1.0 - k * k * sin2_phi;
    rf = gsl_sf_ellint_RF_e(x, y, 1.0, mode);
    rj = gsl_sf_ellint_RJ_e(x, y, 1.0, 1.0 + n * sin2_phi, mode);
    r.val = sin_phi * rf.val - n/3.0*sin3_phi * rj.val;
    r.err = GSL_DBL_EPSILON * Math.abs(sin_phi * rf.val);
    r.err = r.err + Math.abs(sin_phi * rf.err);
    r.err = r.err + n / 3.0 * GSL_DBL_EPSILON * Math.abs(sin3_phi * rj.val);
    r.err = r.err + n / 3.0 * Math.abs(sin3_phi * rj.err);
    if (nc != 0.0)
    {
        rp = gsl_sf_ellint_Pcomp_e(k, n, mode);  
        r.val = r.val + 2.0 * nc * rp.val;
        r.err = r.err + 2.0 * Math.abs(nc) * rp.err;
    }

    return r;

} // gsl_sf_ellint_P_e

// ----------------------------------------------------------------------------

// [Carlson, Numer. Math. 33 (1979) 1, (4.4)]
export function gsl_sf_ellint_D_e(phi0, k, n, mode)
{
    var nc       = 0.0;
    var phi_red  = 0.0;
    var sin_phi  = 0.0;
    var sin2_phi = 0.0;
    var sin3_phi = 0.0;
    var x        = 0.0;
    var y        = 0.0;
    var rd = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;
    var phi = phi0;

    // Angular reduction to -pi/2 < phi < pi/2 (we should really use an
    // exact reduction but this will have to do for now) BJG

    nc = Math.floor(phi / M_PI + 0.5);
    phi_red = phi - nc * M_PI;
    phi = phi_red;
   
    // FIXME: need to handle the case of small x, as for E,F
    sin_phi  = Math.sin(phi);
    sin2_phi = sin_phi  * sin_phi;
    sin3_phi = sin2_phi * sin_phi;
    x = 1.0 - sin2_phi;
    y = 1.0 - k * k * sin2_phi;
    rd = gsl_sf_ellint_RD_e(x, y, 1.0, mode);
    r.val = sin3_phi / 3.0 * rd.val;
    r.err = GSL_DBL_EPSILON * Math.abs(r.val) + Math.abs(sin3_phi / 3.0 * rd.err);
    if (nc != 0.0)
    {
        rd = gsl_sf_ellint_Dcomp_e(k, mode);  
        r.val = r.val + 2.0 * nc * rd.val;
        r.err = r.err + 2.0 * Math.abs(nc) * rd.err;
    }

    return r;

} // gsl_sf_ellint_D_e

// ----------------------------------------------------------------------------

export function gsl_sf_ellint_Dcomp_e(k, mode)
{
    const y  = 1.0 - k * k; // FIXME: still need to handle k~=~1
    var rd = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (k * k >= 1.0)
    {
        throw "SF.DomainException";
    }
    else
    {
        rd = gsl_sf_ellint_RD_e(0.0, y, 1.0, mode);
        r.val = (1.0 / 3.0) * rd.val;
        r.err = GSL_DBL_EPSILON * Math.abs(r.val) + Math.abs(1.0 / 3.0 * rd.err);
    }

    return r;

} // gsl_sf_ellint_Dcomp_e

// ----------------------------------------------------------------------------

// [Carlson, Numer. Math. 33 (1979) 1, (4.5)]
export function gsl_sf_ellint_Kcomp_e(k, mode)
{
    var y = 1.0 - k * k;
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (k * k >= 1.0)
    {
        throw "SF.DomainException";
    }
    else if (k * k >= 1.0 - GSL_SQRT_DBL_EPSILON)
    {
        // [Abramowitz+Stegun, 17.3.34]
        const a  = [ 1.38629436112, 0.09666344259, 0.03590092383 ];
        const b  = [ 0.5, 0.12498593597, 0.06880248576 ];
        const  ta = a[0] + y * (a[1] + y * a[2]);
        const  tb = -Math.log(y) * (b[0] + y * (b[1] + y * b[2]));

        r.val = ta + tb;
        r.err = 2.0 * GSL_DBL_EPSILON * (Math.abs(r.val) + Math.abs(k / y));
    }
    else
    {
        // This was previously computed as,
        //
        //   code = gsl_sf_ellint_RF_e(0.0, 1.0 - k*k, 1.0, mode, result);
        //
        // but this underestimated the total error for small k, since the 
        // argument y=1-k^2 is not exact (there is an absolute error of
        // GSL_DBL_EPSILON near y=0 due to cancellation in the subtraction).
        // Taking the singular behavior of -log(y) above gives an error
        // of 0.5*epsilon/y near y=0. (BJG)
        r = gsl_sf_ellint_RF_e(0.0, y, 1.0, mode);
        r.err = r.err + 0.5 * GSL_DBL_EPSILON / y;
    }

    return r;

} // gsl_sf_ellint_Kcomp_e

// ----------------------------------------------------------------------------

// [Carlson, Numer. Math. 33 (1979) 1, (4.6)]
export function gsl_sf_ellint_Ecomp_e(k, mode)
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if (k * k >= 1.0)
    {
        throw "SF.DomainException";
    }
    else if (k * k >= 1.0 - GSL_SQRT_DBL_EPSILON)
    {
        // [Abramowitz+Stegun, 17.3.36]
        const a  = [ 0.44325141463, 0.06260601220, 0.04757383546 ];
        const b  = [ 0.24998368310, 0.09200180037, 0.04069697526 ];
        const y  = 1.0 - k * k;
        const ta = 1.0 + y * (a[0] + y * (a[1] + a[2] * y));
        const tb = -y * Math.log(y) * (b[0] + y * (b[1] + b[2] * y));

        r.val = ta + tb;
        r.err = 2.0 * GSL_DBL_EPSILON * r.val;
    }
    else
    {
        var rf = { val: 0.0, err: 0.0 }; // Result;
        var rd = { val: 0.0, err: 0.0 }; // Result;
        const y = 1.0 - k * k;
        var rfstatus = 0;
        var rdstatus = 0;

        rf = gsl_sf_ellint_RF_e(0.0, y, 1.0, mode);
        rd = gsl_sf_ellint_RD_e(0.0, y, 1.0, mode);
        r.val = rf.val - k * k / 3.0 * rd.val;
        r.err = rf.err + k * k / 3.0 * rd.err;
    }

    return r;

} // gsl_sf_ellint_Ecomp_e

// ----------------------------------------------------------------------------

// [Carlson, Numer. Math. 33 (1979) 1, (4.3) phi=pi/2]
export function gsl_sf_ellint_Pcomp_e(k, n, mode)
{
    const y  = 1.0 - k * k;
    var rf = { val: 0.0, err: 0.0 }; // Result;
    var rj = { val: 0.0, err: 0.0 }; // Result;
    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (k * k >= 1.0)
    {
        throw "SF.DomainException";
    }
    // FIXME: need to handle k ~=~ 1  cancellations
    else
    {
        rf = gsl_sf_ellint_RF_e(0.0, y, 1.0, mode);
        rj = gsl_sf_ellint_RJ_e(0.0, y, 1.0, 1.0 + n, mode);
        r.val = rf.val - (n / 3.0) * rj.val;
        r.err = rf.err + Math.abs(n / 3.0) * rj.err;
    }

    return r;

} // gsl_sf_ellint_Pcomp_e

// *-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_ellint_Kcomp( k, mode )
{ // gsl_sf_ellint_Kcomp
    return EVAL_RESULT_DM( gsl_sf_ellint_Kcomp_e, { x: k, m: mode }, "gsl_sf_ellint_Kcomp" );
} // gsl_sf_ellint_Kcomp;

export function gsl_sf_ellint_Ecomp( k, mode )
{ // gsl_sf_ellint_Ecomp
    return EVAL_RESULT_DM( gsl_sf_ellint_Ecomp_e, { x: k, m: mode }, "gsl_sf_ellint_Ecomp" );
} // gsl_sf_ellint_Ecomp;

export function gsl_sf_ellint_Pcomp( k, n, mode )
{ // gsl_sf_ellint_Pcomp
    return EVAL_RESULT_DDM( gsl_sf_ellint_Pcomp_e, { x: k, y: n, m: mode }, "gsl_sf_ellint_Pcomp" );
} // gsl_sf_ellint_Pcomp;

export function gsl_sf_ellint_Dcomp( k, mode )
{ // gsl_sf_ellint_Dcomp
    return EVAL_RESULT_DM( gsl_sf_ellint_Dcomp_e, { x: k, m: mode }, "gsl_sf_ellint_Dcomp" );
} // gsl_sf_ellint_Dcomp;

export function gsl_sf_ellint_F( phi, k, mode )
{ // gsl_sf_ellint_F
    return EVAL_RESULT_DDM( gsl_sf_ellint_F_e, { x: phi, y: k, m: mode }, "gsl_sf_ellint_F" );
} // gsl_sf_ellint_F;

export function gsl_sf_ellint_E( phi, k, mode )
{ // gsl_sf_ellint_E
    return EVAL_RESULT_DDM( gsl_sf_ellint_E_e, { x: phi, y: k, m: mode }, "gsl_sf_ellint_E" );
} // gsl_sf_ellint_E;

export function gsl_sf_ellint_P( phi, k, n, mode )
{ // gsl_sf_ellint_P
    return EVAL_RESULT_3DM( gsl_sf_ellint_P_e, { x: phi, y: k, z: n, m: mode }, "gsl_sf_ellint_P" );
} // gsl_sf_ellint_P;

export function gsl_sf_ellint_D( phi, k, n, mode )
{ // gsl_sf_ellint_D
    return EVAL_RESULT_3DM( gsl_sf_ellint_D_e, { x: phi, y: k, z: n, m: mode }, "gsl_sf_ellint_D" );
} // gsl_sf_ellint_D;

export function gsl_sf_ellint_RC( x, y, mode )
{ // gsl_sf_ellint_RC
    return EVAL_RESULT_DDM( gsl_sf_ellint_RC_e, { x: x, y: y, m: mode }, "gsl_sf_ellint_RC" );
} // gsl_sf_ellint_RC;

export function gsl_sf_ellint_RD( x, y, z, mode )
{ // gsl_sf_ellint_RD
    return EVAL_RESULT_3DM( gsl_sf_ellint_RD_e, { x: x, y: y, z: z, m: mode }, "gsl_sf_ellint_RD" );
} // gsl_sf_ellint_RD;

export function gsl_sf_ellint_RF( x, y, z, mode )
{ // gsl_sf_ellint_RF
    return EVAL_RESULT_3DM( gsl_sf_ellint_RF_e, { x: x, y: y, z: z, m: mode }, "gsl_sf_ellint_RF" );
} // gsl_sf_ellint_RF;

export function gsl_sf_ellint_RJ( x, y, z, p, mode )
{ // gsl_sf_ellint_RJ
    return EVAL_RESULT_4DM( gsl_sf_ellint_RJ_e, { x: x, y: y, z: z, p: p, m: mode }, "gsl_sf_ellint_RJ" );
} // gsl_sf_ellint_RJ;

// ----------------------------------------------------------------------------
// EOF SF-EllipticIntegrals.mjs
