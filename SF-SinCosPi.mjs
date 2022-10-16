// specfunc/sincos_pi.c
// 
// Copyright (C) 2017 Gerard Jungman, Konrad Griessinger (konradg@gmx.net)
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
//

// routines for computing sin(pi*x) and cos(pi*x), respectively, with argument reduction

import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";

// Any double precision number bigger than this is automatically an even integer.
const TWOBIG = (2.0 / GSL_DBL_EPSILON);


//Math.fmod = function (a,b) { return Number((a - (Math.floor(a / b) * b)).toPrecision(8)); };
Math.fmod = function (a,b)
{
    let mod;
    // Handling negative values
    if (a < 0)
        mod = -a;
    else
        mod =  a;
    if (b < 0)
        b = -b;
 
    // Finding mod by
    // repeated subtraction
     
    while (mod >= b)
        mod = mod - b;
 
    // Sign of result typically
    // depends on sign of a.
    if (a < 0)
        return -mod;
 
    return mod;
}

// routine computing sin(pi*x) using a Taylor expansion around the origin and
// otherwise the library function.
function sin_pi_taylor( x )
{
    var r = { val: 0.0, err: 0.0 };

    r.val = 0.0;
    r.err = 0.0;
    if ( 16.0 * Math.abs( x ) < 1.0 )
    {
        const y = M_PI * x;
        const a = y * y;
        r.val = y * (1.0 - a * (1.0 - a * (1.0 - a * (1.0 - a * (1.0 - a / 110.0) / 72.0) / 42.0) / 20.0) / 6.0);
    }
    else
    {
        r.val = Math.sin( M_PI * x );
    }
    
    r.err = GSL_DBL_EPSILON * Math.abs( r.val );
  
    return r;
}

// routine computing sin(pi*x) using a Taylor expansion around the origin
// and otherwise the library function.
function cos_pi_taylor( x )
{
    var r = { val: 0.0, err: 0.0 };

    r.val = 0.0;
    r.err = 0.0;
    if ( 20.0 * Math.abs( x ) < 1.0 )
    {
        const y = M_PI * x;
        const a = y * y;
        r.val = 1.0 - 0.5 * a * (1.0 - a * (1.0 - a * (1.0 - a * (1.0 - a / 90.0) / 56.0) / 30.0) / 12.0);
    }
    else
    {
        r.val = Math.cos( M_PI * x );
    }

    r.err = GSL_DBL_EPSILON * Math.abs(r.val);

    return r;
}

export function gsl_sf_sin_pi_e( x )
{
    var r = { val: 0.0, err: 0.0 };

    var intx = 0.0;
    var fracx = 0.0;
    var q;
    var sign = 1;

    r.val = 0.0;
    r.err = 0.0;
    //fracx = modf(x,&intx);
    intx = Math.trunc( x );
    fracx = x - intx;
    if ( fracx == 0.0 )
    {
        return r;
    }
    if ( Math.abs( intx ) >= TWOBIG )
    {
        return r; // to be sure. Actually should be covered by the line above
    }

    q = ( ( ( intx >= LONG_MIN ) && ( intx <= LONG_MAX ) ) ? intx : Math.fmod( intx, 2.0 ) );
    sign = ( q % 2 ? -1 : 1 );

    /* int sign = 1 - 2*((int)round(fmod(fabs(intx),2.0))); */
    if ( Math.abs( fracx ) == 0.5 )
    { // probably unnecessary
        if ( fracx < 0.0 ) sign = -sign;
        r.val = ( sign != 1 ? -1.0 : 1.0 );
        return r;
    }
    if ( Math.abs( fracx ) > 0.5)
    {
        sign = -sign;
        fracx = ( fracx > 0.0 ? fracx - 1.0 : fracx + 1.0 );
    }

    if ( fracx > 0.25 )
    {
        r = cos_pi_taylor( fracx - 0.5 );
    }
    else if ( fracx < -0.25 )
    {
        r = cos_pi_taylor( fracx + 0.5 );
        sign = -sign;
    }
    else
    {
        r = sin_pi_taylor( fracx );
    }
    if ( sign != 1 ) r.val = -r.val;
    return r;
}

export function gsl_sf_cos_pi_e( x )
{
    var r = { val: 0.0, err: 0.0 };

    var intx = 0.0;
    var fracx = 0.0;
    var q;
    var sign = 1;

    r.val = 0.0;
    r.err = 0.0;
    //fracx = modf(x,&intx);
    intx = Math.trunc( x );
    fracx = x - intx;
    if ( Math.abs( fracx ) == 0.5 )
    {
        return r;
    }
    
    if ( Math.abs( intx ) >= TWOBIG )
    {
        r.val = 1.0;
        return r;
    }

    q = ( ( (intx >= LONG_MIN) && (intx <= LONG_MAX) ) ? intx : Math.fmod( intx, 2.0 ) );
    sign = ( q % 2 ? -1 : 1 );

    // int sign = 1 - 2*((int)round(fmod(fabs(intx),2.0)));
    if ( fracx == 0.0 )
    { // probably unnecessary
        r.val = ( sign != 1 ? -1.0 : 1.0 );
        return r;
    }
    if ( Math.aabs( fracx ) > 0.5)
    {
        sign = -sign;
        fracx = ( fracx > 0.0 ? fracx - 1.0 : fracx + 1.0 );
    }

    if ( fracx > 0.25 )
    {
        r = sin_pi_taylor( fracx - 0.5 );
        sign = -sign;
    }
    else if ( fracx < -0.25 )
    {
        r = sin_pi_taylor( fracx + 0.5 );
    }
    else
    {
        r = cos_pi_taylor( fracx );
    }
    if ( sign != 1 ) r.val = -r.val;
    return r;
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

export function gsl_sf_sin_pi( x )
{
    return EVAL_RESULT_D( gsl_sf_sin_pi_e, x, "gsl_sf_sin_pi_e" );
}

export function gsl_sf_cos_pi( x )
{
    return EVAL_RESULT_D( gsl_sf_cos_pi_e, x, "gsl_sf_cos_pi_e" );
}

// ----------------------------------------------------------------------------
// EOF SF-SinCosPi.mjs
