// SF-Mode.mjs
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
// Translation to JavaScript: Igor Izvarin

import { GSL_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_FLT_EPSILON } from "./SF-Machine.mjs";
import { GSL_SFLT_EPSILON } from "./SF-Machine.mjs";

// Some functions can take a mode argument. This
// is a rough method to do things like control
// the precision of the algorithm. This mainly
// occurs in special functions, but we figured
// it was ok to have a general facility.
//
// The mode type is 32-bit field. Most of
// the fields are currently unused. Users
// '|' various predefined constants to get
// a desired mode.
//

// Here are the predefined constants.
// Note that the precision constants
// are special because they are used
// to index arrays, so do not change
// them. The precision information is
// in the low order 3 bits of gsl_mode_t
// (the third bit is currently unused).

// Note that "0" is double precision,
// so that you get that by default if
// you forget a flag.
//
export const GSL_PREC_DOUBLE = 0;
export const GSL_PREC_SINGLE = 1;
export const GSL_PREC_APPROX = 2;

// Here are some predefined generic modes.
export var GSL_MODE_DEFAULT = GSL_PREC_DOUBLE; // 0

export function GSL_MODE_PREC(mt)
{
    return mt & 7 != 0;
} // GSL_MODE_PREC

export const gsl_prec_eps       = [ GSL_DBL_EPSILON,       GSL_FLT_EPSILON,       GSL_SFLT_EPSILON ];
// export const gsl_prec_sqrt_eps  = [ GSL_SQRT_DBL_EPSILON,  GSL_SQRT_FLT_EPSILON,  GSL_SQRT_SFLT_EPSILON ];
// export const gsl_prec_root3_eps = [ GSL_ROOT3_DBL_EPSILON, GSL_ROOT3_FLT_EPSILON, GSL_ROOT3_SFLT_EPSILON ];
// export const gsl_prec_root4_eps = [ GSL_ROOT4_DBL_EPSILON, GSL_ROOT4_FLT_EPSILON, GSL_ROOT4_SFLT_EPSILON ];
// export const gsl_prec_root5_eps = [ GSL_ROOT5_DBL_EPSILON, GSL_ROOT5_FLT_EPSILON, GSL_ROOT5_SFLT_EPSILON ];
// export const gsl_prec_root6_eps = [ GSL_ROOT6_DBL_EPSILON, GSL_ROOT6_FLT_EPSILON, GSL_ROOT6_SFLT_EPSILON ];

// ----------------------------------------------------------------------------
// EOF SF-Mode.mjs
