// SF-Evaluate.mjs
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

// ----------------------------------------------------------------------------

export function EVAL_RESULT_D(func, x, funcname)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func(x);
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_D

export function EVAL_RESULT_DD(func, args, funcname)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func(args.x, args.y);
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_DD

export function EVAL_RESULT_3D( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.x, args.y, args.z );
        return r.val;
    }
    catch ( e )
    {
        throw ( e );
    }
} // EVAL_RESULT_3D

export function EVAL_RESULT_4D( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.x, args.dx, args.y, args.dy );
        return r.val;
    }
    catch ( e )
    {
        throw ( e );
    }
} // EVAL_RESULT_4D

export function EVAL_RESULT_ID(func, args, funcname)
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func(args.i, args.x);
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_ID

export function EVAL_RESULT_DI( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.x, args.i );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_DI

export function EVAL_RESULT_IDD( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.n, args.x, args.y );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_IDD

export function EVAL_RESULT_DM( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.x, args.m );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_DM

export function EVAL_RESULT_DDM( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.x, args.y, args.m );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_DDM

export function EVAL_RESULT_3DM( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.x, args.y, args.z, args.m );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_3DM

export function EVAL_RESULT_4DM( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.x, args.y, args.z, args.p, args.m );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_4DM

export function EVAL_RESULT_I( func, i, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( i );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_I

export function EVAL_RESULT_II( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.n, args.m );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_II

//     function EVAL_RESULT(func: FuncIIDD; args: ArgsIIDD; funcname: STRING) return LONG_FLOAT IS
//         r : Result;
//     BEGIN -- EVAL_RESULT
//         r = func(args.n, args.m, args.x, args.y);
//         return r.val;
//     EXCEPTION
//         WHEN OTHERS => RAISE;
//     END EVAL_RESULT;

export function EVAL_RESULT_IID( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.n, args.m, args.x );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_IID

export function EVAL_RESULT_IIDD( func, args, funcname )
{
    var r = { val: 0.0, err: 0.0 }; // Result;
    try
    {
        r = func( args.i1, args.i2, args.x, args.y );
        return r.val;
    }
    catch (e)
    {
        throw (e);
    }
} // EVAL_RESULT_IIDD

//     function EVAL_RESULT(func: Func6I; args: Args6I; funcname: STRING) return LONG_FLOAT IS
//         r : Result;
//     BEGIN -- EVAL_RESULT
//         r = func(args.i1, args.i2, args.i3, args.i4, args.i5, args.i6);
//         return r.val;
//     EXCEPTION
//         WHEN OTHERS => RAISE;
//     END EVAL_RESULT;

//     function EVAL_RESULT(func: Func9I; args: Args9I; funcname: STRING) return LONG_FLOAT IS
//         r : Result;
//     BEGIN -- EVAL_RESULT
//         r = func(args.i1, args.i2, args.i3, args.i4, args.i5, args.i6, args.i7, args.i8, args.i9);
//         return r.val;
//     EXCEPTION
//         WHEN OTHERS => RAISE;
//     END EVAL_RESULT;

// ----------------------------------------------------------------------------
// EOF SF-Evaluate.mjs
