// Test.mjs
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

import { gsl_test  }       from './SF-TestResult.mjs';
import { SetVerbose  }     from './SF-TestResult.mjs';
import { SetInternalLog  } from './SF-Test.mjs';
import { test_airy }       from './SF-Test.mjs';
import { test_gamma }      from './SF-Test.mjs';
import { test_erf    }     from './SF-Test.mjs';
import { test_expint }     from './SF-Test.mjs';
import { test_dawson }     from './SF-Test.mjs';
import { test_debye }      from './SF-Test.mjs';
import { test_psi }        from './SF-Test.mjs';
import { test_zeta }       from './SF-Test.mjs';
import { test_clausen }    from './SF-Test.mjs';
import { test_log }        from './SF-Test.mjs';
import { test_gegen }      from './SF-Test.mjs';
import { test_legendre }   from './SF-Test.mjs';
import { test_bessel }     from './SF-Test.mjs';
import { test_ellint }     from './SF-Test.mjs';
import { test_elementary } from './SF-Test.mjs';
import { test_exp }        from './SF-Test.mjs';
import { test_trig }       from "./SF-Test.mjs";
import { test_transport }  from "./SF-Test.mjs";
import { test_synch }      from "./SF-Test.mjs";
import { test_pow_int }    from "./SF-Test.mjs";
import { test_lambert }    from "./SF-Test.mjs";
import { test_laguerre }   from "./SF-Test.mjs";
import { test_jac }        from "./SF-Test.mjs";
import { test_dilog }      from "./SF-Test.mjs";
import { test_fermidirac } from "./SF-Test.mjs";
import { test_hyperg }     from "./SF-Test.mjs";
import { test_coupling }   from "./SF-Test.mjs";
import { test_coulomb }    from "./SF-Test.mjs";

console.log( "Special Functions Test Utility. Revision 1.0. February 2020." );
console.log( "Copyright (C) 2007, 2008 Brian Gough" );
console.log( "Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2004 Gerard Jungman" );

SetVerbose( true );
SetInternalLog( false );

// gsl_test( test_coupling( ),    "Coupling Coefficients" ); // NaN ???

// gsl_test( test_laguerre( ),    "Laguerre Polynomials" );        // OK. (ESLint verified.)
// gsl_test( test_erf( ),         "Error Functions" );             // OK. (ESLint verified.)
// gsl_test( test_exp( ),         "Exponential Functions" );       // OK. (ESLint verified.)
// gsl_test( test_log( ),         "Logarithm" );                   // OK. (ESLint verified.)
// gsl_test( test_pow_int( ),     "Integer Powers");               // OK. (ESLint verified.)
// gsl_test( test_dawson( ),      "Dawson Integral" );             // OK. (ESLint verified.)
// gsl_test( test_debye( ),       "Debye Functions" );             // OK. (ESLint verified.)
// gsl_test( test_clausen( ),     "Clausen Integral" );            // OK. (ESLint verified.)
// gsl_test( test_lambert( ),     "Lambert W Functions" );         // OK. (ESLint verified.)
// gsl_test( test_gegen( ),       "Gegenbauer Polynomials" );      // OK. (ESLint verified.)
// gsl_test( test_ellint( ),      "Elliptic Integrals" );          // OK. (ESLint verified.)
// gsl_test( test_synch( ),       "Synchrotron Functions" );       // OK. (ESLint verified.)
// gsl_test( test_expint( ),      "Exponential/Sine/Cosine Integrals" ); // OK. (ESLint verified.)
// gsl_test( test_trig( ),        "Trigonometric and Related Functions" );  // OK (ESLint verified.) (***)
// gsl_test( test_transport( ),   "Transport Functions" );         // OK. (ESLint verified.)
// gsl_test( test_psi( ),         "Psi Functions" );               // OK. (ESLint verified.)
// gsl_test( test_dilog( ),       "Dilogarithm" );                 // OK. (ESLint verified.)
// gsl_test( test_gamma( ),       "Gamma Functions" );             // OK. (ESLint verified.) (***)
// gsl_test( test_bessel( ),      "Bessel Functions" );            // OK. (ESLint verified.)
// gsl_test( test_jac( ),         "Elliptic Functions (Jacobi)" ); // OK. (ESLint verified.) (###)
// gsl_test( test_zeta( ),        "Zeta Functions" );              // OK. (ESLint verified.)
// gsl_test( test_elementary( ),  "Elementary Functions (Misc)" ); // OK. (ESLint verified.) (***)
// gsl_test( test_airy( ),        "Airy Functions" );              // OK. (ESLint verified.)
// gsl_test( test_legendre( ),    "Legendre Functions" );          // OK. (ESLint verified.)
// gsl_test( test_hyperg( ),      "Hypergeometric Functions" );    // OK. (ESLint verified.)
// gsl_test( test_coulomb( ),     "Coulomb Wave Functions" );      // OK. (ESLint verified.)
// gsl_test( test_fermidirac( ),  "Fermi-Dirac Functions" );       // OK. (ESLint verified.)

//gsl_test(test_mathieu,     "Mathieu Functions");
//gsl_test(test_psi_complex(), "Psi Function for complex argument");

// ----------------------------------------------------------------------------
// EOF Test.mjs
