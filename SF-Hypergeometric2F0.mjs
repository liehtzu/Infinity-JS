// SF-Hypergeometric2F0.mjs
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
import { gsl_sf_hyperg_U_e }     from "./SF-HypergeometricU.mjs";

import { EVAL_RESULT_3D }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

export function gsl_sf_hyperg_2F0_e( a, b, x )
{
    var r = { val: 0.0, err: 0.0 }; // Result;

    if ( x < 0.0 )
    {
        // Use "definition" 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x).
        //
        var U   = { val: 0.0, err: 0.0 }; // Result;
        var pre = 0.0;

        pre = (-1.0 / x) ** a;
        U = gsl_sf_hyperg_U_e( a, 1.0 + a - b, -1.0 / x );
        r.val = pre * U.val;
        r.err = GSL_DBL_EPSILON * Math.abs( r.val ) + pre * U.err;
        return r;
    }
    else if ( x == 0.0 )
    {
        r.val = 1.0;
        r.err = 0.0;
        return r;
    }
    else
    {
        // Use asymptotic series. ??
        //
        // return hyperg_2F0_series(a, b, x, -1, result, &prec);
        throw "SF.DomainException";
    }

} // gsl_sf_hyperg_2F0_e

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_hyperg_2F0( a, b, x )
{ //- gsl_sf_hyperg_2F0
    return EVAL_RESULT_3D( gsl_sf_hyperg_2F0_e, { x: a, y: b, z: x }, "gsl_sf_hyperg_2F0" );
} // gsl_sf_hyperg_2F0

// ----------------------------------------------------------------------------
// EOF SF-Hypergeometric2F0.mjs
