// SF-BesselSequence.mjs
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

import { gsl_sf_bessel_Jnu_e }   from "./SF-BesselJnu.mjs";
import { GSL_MODE_PREC }         from "./SF-Mode.mjs";

// ----------------------------------------------------------------------------

function DYDX_p( p, u, x, nu )
{
    return (-p / x + ((nu * nu) / (x * x) - 1.0) * u);
} // DYDX_p

function DYDX_u( p, u, x )
{
    return p;
} // DYDX_u

function rk_step( nu, x, dx, Jp, J )
{
    var p_0 = Jp.Double;
    var u_0 = J.Double;
   
    var p_1 = dx * DYDX_p( p_0, u_0, x, nu );
    var u_1 = dx * DYDX_u( p_0, u_0, x );
   
    var p_2 = dx * DYDX_p( p_0 + 0.5 * p_1, u_0 + 0.5 * u_1, x + 0.5 * dx, nu );
    var u_2 = dx * DYDX_u( p_0 + 0.5 * p_1, u_0 + 0.5 * u_1, x + 0.5 * dx );
   
    var p_3 = dx * DYDX_p( p_0 + 0.5 * p_2, u_0 + 0.5 * u_2, x + 0.5 * dx, nu );
    var u_3 = dx * DYDX_u( p_0 + 0.5 * p_2, u_0 + 0.5 * u_2, x + 0.5 * dx );
   
    var p_4 = dx * DYDX_p( p_0 + p_3, u_0 + u_3, x + dx, nu );
    var u_4 = dx * DYDX_u( p_0 + p_3, u_0 + u_3, x + dx );

    Jp.Double = p_0 + p_1 / 6.0 + p_2 / 3.0 + p_3 / 3.0 + p_4 / 6.0;
    J.Double  = u_0 + u_1 / 6.0 + u_2 / 3.0 + u_3 / 3.0 + u_4 / 6.0;

} // rk_step

// ----------------------------------------------------------------------------

export function gsl_sf_bessel_sequence_Jnu_e( nu, mode, size, v )
{

    if ( nu < 0.0 )
    {
        throw "SF.DomainException";
    }
    else if ( size == 0 )
    {
        //GSL_ERROR("error", GSL_EINVAL);
        throw "SF.ArgumentException";
    }
    else
    {
        var goal = GSL_MODE_PREC( mode );
        const dx_array = [ 0.001, 0.03, 0.1 ]; // double, single, approx
        var dx_nominal = dx_array[goal];
        
        var cnu = 0;
        var nu13 = 0.0;
        const smalls =
            [
            0.01,
            0.02,
            0.4,
            0.7,
            1.3,
            2.0,
            2.5,
            3.2,
            3.5,
            4.5,
            6.0
            ];
        var x_small = 0.0;
        
        var J0 = { val: 0.0, err: 0.0 }; // Result;
        var J1 = { val: 0.0, err: 0.0 }; // Result;
        var Jp = { Double: 0.0 };
        var J  = { Double: 0.0 };
        var x  = 0.0;
        var i  = 0;

        var Nd = 0;
        var dv = 0.0;
        var dx = 0.0;
        var xj = 0.0;
           
        cnu = Math.trunc( Math.ceil( nu ) );
        nu13 = Math.pow( nu, 1.0 / 3.0 );
        if ( nu >= 10.0 )
        {
            x_small = nu - nu13;
        }
        else
        {
            x_small = smalls[cnu];
        }

        // Calculate the first point.
        x = v[0];
        try // Igor Izvarin 05-14-2020
        {
            J0 = gsl_sf_bessel_Jnu_e( nu, x );
            v[0] = J0.val;
        }
        catch (e)
        {
            v[0] = 1.0;
        }
        i = i + 1;
        
        // Step over the idiot case where the
        // first point was actually zero.
        //
        if ( x == 0.0 )
        {
            if ( v[1] <= x )
            {
                // Strict ordering failure.
                throw "SF.GenericException";
            }
            x = v[1];
            J0 = gsl_sf_bessel_Jnu_e( nu, x );
            v[1] = J0.val;
            i = i + 1;
        }
        
        // Calculate directly as long as the argument
        // is small. This is necessary because the
        // integration is not very good there.
        //
        while ( v[i] < x_small && i < size )
        {
            if ( v[i] <= x )
            {
                // Strict ordering failure.
                throw "SF.GenericException";
            }
            x = v[i];
            J0 = gsl_sf_bessel_Jnu_e( nu, x );
            v[i] = J0.val;
            i = i + 1;
        }
        
        // At this point we are ready to integrate.
        // The value of x is the last calculated
        // point, which has the value J0; v[i] is
        // the next point we need to calculate. We
        // calculate nu+1 at x as well to get
        // the derivative, then we go forward.
        //
        J1 = gsl_sf_bessel_Jnu_e( nu + 1.0, x );
        J.Double  = J0.val;
        Jp.Double = -J1.val + nu / x * J0.val;
        
        while ( i < size )
        {
            dv = v[i] - x;
            Nd = Math.trunc( Math.ceil( dv / dx_nominal ) );
            dx = dv / (Nd);
            
            if ( v[i] <= x )
            {
                // Strict ordering failure.
                throw "SF.GenericException";
            }
            
            // Integrate over interval up to next sample point.
            //
            xj = x;
            for ( let m = 0; m <= Nd - 1; m++ )
            {
                rk_step( nu, xj, dx, Jp, J );
                xj = xj + dx;
            }
            
            // Go to next interval.
            x = v[i];
            v[i] = J.Double;
            i = i + 1;
        }
    }

} // gsl_sf_bessel_sequence_Jnu_e

// ----------------------------------------------------------------------------
// EOF SF-BesselSequence.mjs
