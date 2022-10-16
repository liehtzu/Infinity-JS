// SF-EllipticJacobi.mjs
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

import { GSL_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_SIGN }         from "./SF-Math.mjs";
import { gsl_sf_hypot }     from "./SF-Trigonometric.mjs";

// ----------------------------------------------------------------------------
// See Hart et al, Computer Approximations, John Wiley and Sons, New York (1968)
// (This applies only to the erfc8 stuff, which is the part
//  of the original code that survives. I have replaced much of
//  the other stuff with Chebyshev fits. These are simpler and
//  more precise than the original approximations. [GJ])
// ----------------------------------------------------------------------------

// GJ: See [Thompson, Atlas for Computing Mathematical Functions]

// BJG 2005-07: New algorithm based on Algorithm 5 from Numerische
// Mathematik 7, 78-90 (1965) "Numerical Calculation of Elliptic
// Integrals and Elliptic Functions" R. Bulirsch.
//
// Minor tweak is to avoid division by zero when sin(x u_l) = 0 by
// computing reflected values sn(K-u) cn(K-u) dn(K-u) and using
// transformation from Abramowitz & Stegun table 16.8 column "K-u"*/

export function gsl_sf_elljac_e(u, m, sn, cn, dn)
{

    if (Math.abs(m) > 1.0)
    {
        //sn = 0.0;
        //cn = 0.0;
        //dn = 0.0;
        throw "SF.DomainException"; // WITH "|m| > 1.0";
    }
    else if (Math.abs(m) < 2.0 * GSL_DBL_EPSILON)
    {
        sn.Double = Math.sin(u);
        cn.Double = Math.cos(u);
        dn.Double = 1.0;
    }
    else if (Math.abs(m - 1.0) < 2.0 * GSL_DBL_EPSILON)
    {
        sn.Double = Math.tanh(u);
        cn.Double = 1.0 / Math.cosh(u);
        dn.Double = cn.Double;
    }
    else
    {
        const NMAX = 16;
        var mu = [];
        var nu = [];
        var c = [];
        var d = [];
        var sin_umu = 0.0;
        var cos_umu = 0.0;
        var t  = 0.0;
        var r  = 0.0;
        var n = 0;

        mu[0] = 1.0;
        nu[0] = Math.sqrt(1.0 - m);
        
        while (Math.abs(mu[n] - nu[n]) > 4.0 * GSL_DBL_EPSILON * Math.abs(mu[n] + nu[n]))
        {
            mu[n+1] = 0.5 * (mu[n] + nu[n]);
            nu[n+1] = Math.sqrt(mu[n] * nu[n]);
            n = n + 1;
            if (n >= NMAX - 1)
            {
                //break;
                throw "SF.MaxIterationsException";
            }
        }
        
        sin_umu = Math.sin(u * mu[n]);
        cos_umu = Math.cos(u * mu[n]);
        
        // Since sin(u*mu(n)) can be zero we switch to computing sn(K-u),
        // cn(K-u), dn(K-u) when |sin| < |cos|
        
        if (Math.abs(sin_umu) < Math.abs(cos_umu))
        {
            t = sin_umu / cos_umu;
            
            c[n] = mu[n] * t;
            d[n] = 1.0;
            
            while (n > 0)
            {
                n = n - 1;
                c[n] = d[n+1] * c[n+1];
                r = (c[n+1] * c[n+1]) / mu[n+1];
                d[n] = (r + nu[n]) / (r + mu[n]);
            }
            
            dn.Double = Math.sqrt(1.0 - m) / d[n];
            cn.Double = dn.Double * GSL_SIGN(cos_umu) / gsl_sf_hypot(1.0, c[n]);
            sn.Double = cn.Double * c[n] / Math.sqrt(1.0 - m);
        }
        else
        {
            t = cos_umu / sin_umu;
            
            c[n] = mu[n] * t;
            d[n] = 1.0;
            
            while (n > 0)
            {
                n = n - 1;
                c[n] = d[n+1] * c[n+1];
                r = (c[n+1] * c[n+1]) / mu[n+1];
                d[n] = (r + nu[n]) / (r + mu[n]);
            }
            
            dn.Double = d[n];
            sn.Double = GSL_SIGN(sin_umu) / gsl_sf_hypot(1.0, c[n]);
            cn.Double = c[n] * sn.Double;
        }
    }

} // gsl_sf_elljac_e

// ----------------------------------------------------------------------------
// EOF SF-EllipticJacobi.mjs
