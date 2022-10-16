// SF-Test-test_gamma.mjs
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

import { gsl_sf_lngamma_e } from "./SF-Gamma.mjs";

export function test_gamma()
{

    var s   = 0;
    var sgn = 0.0;

    console.log("Test Gamma Functions ...");

    console.log("    ... gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, -0.1,                     2.368961332728788655,     TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, -1.0 / 256.0,             5.547444766967471595,     TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, 1.0e-08,                 18.420680738180208905,     TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, 0.1,                      2.252712651734205,        TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, 1.0 + 1.0 / 256.0,       -0.0022422226599611501448, TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, 2.0 + 1.0 / 256.0,        0.0016564177556961728692, TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, 100.0,                  359.1342053695753,          TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, -1.0-1.0 / 65536.0,      11.090348438090047844,     TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, -1.0-1.0 / 268435456.0,  19.408121054103474300,     TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, -100.5,                -364.9009683094273518,       TEST_TOL0, "gsl_sf_lngamma_e");
//     TEST_SF_D(s,  gsl_sf_lngamma_e, -100.0-1.0 / 65536.0,  -352.6490910117097874,       TEST_TOL0, "gsl_sf_lngamma_e");

//     PUT_LINE("    ... gsl_sf_lngamma_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lngamma_sgn_e'Access, 0.7, 0.26086724653166651439, TEST_TOL1, 1.0, "gsl_sf_lngamma_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lngamma_sgn_e'Access, 0.1, 2.2527126517342059599, TEST_TOL0, 1.0, "gsl_sf_lngamma_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lngamma_sgn_e'Access, -0.1, 2.368961332728788655, TEST_TOL0, -1.0, "gsl_sf_lngamma_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lngamma_sgn_e'Access, -1.0-1.0/65536.0, 11.090348438090047844, TEST_TOL0, 1.0, "gsl_sf_lngamma_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lngamma_sgn_e'Access, -2.0-1.0/256.0, 4.848447725860607213, TEST_TOL0, -1.0, "gsl_sf_lngamma_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lngamma_sgn_e'Access, -2.0-1.0/65536.0, 10.397193628164674967, TEST_TOL0, -1.0, "gsl_sf_lngamma_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lngamma_sgn_e'Access, -3.0-1.0/8.0, 0.15431112768404182427, TEST_TOL2, 1.0, "gsl_sf_lngamma_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lngamma_sgn_e'Access, -100.5, -364.9009683094273518, TEST_TOL0, -1.0, "gsl_sf_lngamma_sgn_e");

//     PUT_LINE("    ... gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, 1.0 + 1.0/4096.0, 0.9998591371459403421 , TEST_TOL0, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, 1.0 + 1.0/32.0, 0.9829010992836269148 , TEST_TOL0, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, 2.0 + 1.0/256.0, 1.0016577903733583299 , TEST_TOL0, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, 9.0, 40320.0                   , TEST_TOL0, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, 10.0, 362880.0                  , TEST_TOL0, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, 100.0, 9.332621544394415268e+155 , TEST_TOL2, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, 170.0, 4.269068009004705275e+304 , TEST_TOL2, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, 171.0, 7.257415615307998967e+306 , TEST_TOL2, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, -10.5, -2.640121820547716316e-07  , TEST_TOL0, "gsl_sf_gamma_e");
//     TEST_SF(s,  gsl_sf_gamma_e'Access, -11.25, 6.027393816261931672e-08  , TEST_TOL0, "gsl_sf_gamma_e"); -- exp()... not my fault
//     TEST_SF(s,  gsl_sf_gamma_e'Access, -1.0+1.0/65536.0, -65536.42280587818970 , TEST_TOL0, "gsl_sf_gamma_e");

//     PUT_LINE("    ... gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 1.0e-08, 3989.423555759890865  , TEST_TOL1, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 1.0e-05, 126.17168469882690233 , TEST_TOL0, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 0.001, 12.708492464364073506 , TEST_TOL0, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 1.5, 1.0563442442685598666 , TEST_TOL0, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 3.0, 1.0280645179187893045 , TEST_TOL0, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 9.0, 1.0092984264218189715 , TEST_TOL0, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 11.0, 1.0076024283104962850 , TEST_TOL0, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 100.0, 1.0008336778720121418 , TEST_TOL0, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 1.0e+05, 1.0000008333336805529 , TEST_TOL0, "gsl_sf_gammastar_e");
//     TEST_SF(s,  gsl_sf_gammastar_e'Access, 1.0e+20, 1.0 , TEST_TOL0, "gsl_sf_gammastar_e");

//     PUT_LINE("    ... gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, 1.0, 1.0, TEST_TOL0, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, 2.0, 1.0, TEST_TOL0, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, 3.0, 0.5, TEST_TOL0, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, 4.0, 1.0/6.0, TEST_TOL0, "gsl_sf_gammainv_e");

//     TEST_SF(s,  gsl_sf_gammainv_e'Access, 10.0, 1.0/362880.0, TEST_TOL0, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, 100.0, 1.0715102881254669232e-156, TEST_TOL2, "gsl_sf_gammainv_e");

//     TEST_SF(s,  gsl_sf_gammainv_e'Access, 0.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, -1.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, -2.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, -3.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, -4.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");

//     TEST_SF(s,  gsl_sf_gammainv_e'Access, -10.5, -1.0/2.640121820547716316e-07, TEST_TOL2, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, -11.25, 1.0/6.027393816261931672e-08, TEST_TOL1, "gsl_sf_gammainv_e");
//     TEST_SF(s,  gsl_sf_gammainv_e'Access, -1.0+1.0/65536.0, -1.0/65536.42280587818970 , TEST_TOL1, "gsl_sf_gammainv_e");

//     PUT_LINE("    ... gsl_sf_lngamma_complex_e");
//     TEST_SF_2(s, gsl_sf_lngamma_complex_e'Access, (5.0, 2.0),
//             2.7487017561338026749, TEST_TOL0,
//             3.0738434100497007915, TEST_TOL0,
//             "gsl_sf_lngamma_complex_e");
            
//     TEST_SF_2(s, gsl_sf_lngamma_complex_e'Access, (100.0, 100.0),
//             315.07804459949331323, TEST_TOL1,
//             2.0821801804113110099, TEST_TOL3,
//             "gsl_sf_lngamma_complex_e");

//     TEST_SF_2(s, gsl_sf_lngamma_complex_e'Access, (100.0, -1000.0),
//             -882.3920483010362817000, TEST_TOL1,
//             -2.1169293725678813270, TEST_TOL3,
//             "gsl_sf_lngamma_complex_e");

//     TEST_SF_2(s, gsl_sf_lngamma_complex_e'Access, (-100.0, -1.0),
//             -365.0362469529239516000, TEST_TOL1,
//             -3.0393820262864361140, TEST_TOL1,
//             "gsl_sf_lngamma_complex_e");

//     PUT_LINE("    ... gsl_sf_taylorcoeff_e");
//     TEST_SF(s, gsl_sf_taylorcoeff_e'Access, (10,   1.0/1048576.0), 1.7148961854776073928e-67  , TEST_TOL0, "gsl_sf_taylorcoeff_e");
//     TEST_SF(s, gsl_sf_taylorcoeff_e'Access, (10,   1.0/1024.0), 2.1738891788497900281e-37  , TEST_TOL0, "gsl_sf_taylorcoeff_e");
//     TEST_SF(s, gsl_sf_taylorcoeff_e'Access, (10,   1.0), 2.7557319223985890653e-07  , TEST_TOL0, "gsl_sf_taylorcoeff_e");
//     TEST_SF(s, gsl_sf_taylorcoeff_e'Access, (10,   5.0), 2.6911444554673721340      , TEST_TOL0, "gsl_sf_taylorcoeff_e");
//     TEST_SF(s, gsl_sf_taylorcoeff_e'Access, (10,   500.0), 2.6911444554673721340e+20  , TEST_TOL0, "gsl_sf_taylorcoeff_e");
//     TEST_SF(s, gsl_sf_taylorcoeff_e'Access, (100,  100.0), 1.0715102881254669232e+42  , TEST_TOL1, "gsl_sf_taylorcoeff_e");
//     TEST_SF(s, gsl_sf_taylorcoeff_e'Access, (1000, 200.0), 2.6628790558154746898e-267 , TEST_TOL1, "gsl_sf_taylorcoeff_e");
//     TEST_SF(s, gsl_sf_taylorcoeff_e'Access, (1000, 500.0), 2.3193170139740855074e+131 , TEST_TOL1, "gsl_sf_taylorcoeff_e");

//     PUT_LINE("    ... gsl_sf_fact_e");
//     TEST_SF(s, gsl_sf_fact_e'Access, 0, 1.0 , TEST_TOL0, "gsl_sf_fact_e");
//     TEST_SF(s, gsl_sf_fact_e'Access, 1, 1.0 , TEST_TOL0, "gsl_sf_fact_e");
//     TEST_SF(s, gsl_sf_fact_e'Access, 7, 5040.0 , TEST_TOL0, "gsl_sf_fact_e");
//     TEST_SF(s, gsl_sf_fact_e'Access, 33, 8.683317618811886496e+36 , TEST_TOL0, "gsl_sf_fact_e");

//     PUT_LINE("    ... gsl_sf_doublefact_e");
//     TEST_SF(s, gsl_sf_doublefact_e'Access, 0, 1.0 , TEST_TOL0, "gsl_sf_doublefact_e");
//     TEST_SF(s, gsl_sf_doublefact_e'Access, 1, 1.0 , TEST_TOL0, "gsl_sf_doublefact_e");
//     TEST_SF(s, gsl_sf_doublefact_e'Access, 7, 105.0 , TEST_TOL0, "gsl_sf_doublefact_e");
//     TEST_SF(s, gsl_sf_doublefact_e'Access, 33, 6.332659870762850625e+18 , TEST_TOL0, "gsl_sf_doublefact_e");

//     PUT_LINE("    ... gsl_sf_lnfact_e");
//     TEST_SF(s, gsl_sf_lnfact_e'Access, 0, 0.0 , TEST_TOL0, "gsl_sf_lnfact_e");
//     TEST_SF(s, gsl_sf_lnfact_e'Access, 1, 0.0 , TEST_TOL0, "gsl_sf_lnfact_e");
//     TEST_SF(s, gsl_sf_lnfact_e'Access, 7, 8.525161361065414300 , TEST_TOL0, "gsl_sf_lnfact_e");
//     TEST_SF(s, gsl_sf_lnfact_e'Access, 33, 85.05446701758151741 , TEST_TOL0, "gsl_sf_lnfact_e");

//     PUT_LINE("    ... gsl_sf_lndoublefact_e");
//     TEST_SF(s, gsl_sf_lndoublefact_e'Access, 0, 0.0  , TEST_TOL0, "gsl_sf_lndoublefact_e");
//     TEST_SF(s, gsl_sf_lndoublefact_e'Access, 7, 4.653960350157523371  , TEST_TOL0, "gsl_sf_lndoublefact_e");
//     TEST_SF(s, gsl_sf_lndoublefact_e'Access, 33, 43.292252022541719660 , TEST_TOL0, "gsl_sf_lndoublefact_e");
//     TEST_SF(s, gsl_sf_lndoublefact_e'Access, 34, 45.288575519655959140 , TEST_TOL0, "gsl_sf_lndoublefact_e");
//     TEST_SF(s, gsl_sf_lndoublefact_e'Access, 1034, 3075.6383796271197707 , TEST_TOL0, "gsl_sf_lndoublefact_e");
//     TEST_SF(s, gsl_sf_lndoublefact_e'Access, 1035, 3078.8839081731809169 , TEST_TOL0, "gsl_sf_lndoublefact_e");

//     PUT_LINE("    ... gsl_sf_lnchoose_e");
//     TEST_SF(s, gsl_sf_lnchoose_e'Access, (7, 3), 3.555348061489413680 , TEST_TOL0, "gsl_sf_lnchoose_e");
//     TEST_SF(s, gsl_sf_lnchoose_e'Access, (5, 2), 2.302585092994045684 , TEST_TOL0, "gsl_sf_lnchoose_e");

//     PUT_LINE("    ... gsl_sf_choose_e");
//     TEST_SF(s, gsl_sf_choose_e'Access, (7, 3), 35.0 , TEST_TOL0, "gsl_sf_choose_e");
//     TEST_SF(s, gsl_sf_choose_e'Access, (7, 4), 35.0 , TEST_TOL0, "gsl_sf_choose_e");
//     TEST_SF(s, gsl_sf_choose_e'Access, (5, 2), 10.0 , TEST_TOL0, "gsl_sf_choose_e");
//     TEST_SF(s, gsl_sf_choose_e'Access, (5, 3), 10.0 , TEST_TOL0, "gsl_sf_choose_e");

//     TEST_SF(s, gsl_sf_choose_e'Access, (500, 495), 255244687600.0, TEST_TOL0, "gsl_sf_choose_e");
//     TEST_SF(s, gsl_sf_choose_e'Access, (500, 5), 255244687600.0, TEST_TOL0, "gsl_sf_choose_e");

//     TEST_SF(s, gsl_sf_choose_e'Access, (500, 200), 5.054949849935532221e+144 , TEST_TOL5, "gsl_sf_choose_e");
//     TEST_SF(s, gsl_sf_choose_e'Access, (500, 300), 5.054949849935532221e+144 , TEST_TOL5, "gsl_sf_choose_e");

//     PUT_LINE("    ... gsl_sf_lnpoch_e");
//     TEST_SF(s, gsl_sf_lnpoch_e'Access,  (5.0, 0.0), 0.0, TEST_TOL0, "gsl_sf_lnpoch_e");
//     TEST_SF(s, gsl_sf_lnpoch_e'Access,  (5.0, 1.0/65536.0), 0.000022981557571259389129, TEST_TOL0, "gsl_sf_lnpoch_e");
//     --TEST_SF(s, gsl_sf_lnpoch_e'Access,  (5.0, 1.0/256.0),   0.005884960217985189004,    TEST_TOL2, "gsl_sf_lnpoch_e");
//     TEST_SF(s, gsl_sf_lnpoch_e'Access,  (7.0, 3.0), 6.222576268071368616, TEST_TOL0, "gsl_sf_lnpoch_e");
//     TEST_SF(s, gsl_sf_lnpoch_e'Access,  (5.0, 2.0), 3.401197381662155375, TEST_TOL0, "gsl_sf_lnpoch_e");

//     PUT_LINE("    ... gsl_sf_lnpoch_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lnpoch_sgn_e'Access, (5.0, 0.0), 0.0, TEST_TOL1, 1.0, "gsl_sf_lnpoch_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lnpoch_sgn_e'Access, (-4.5, 0.25), 0.7430116475119920117, TEST_TOL1, 1.0, "gsl_sf_lnpoch_sgn_e");
//     TEST_SF_SGN(s, gsl_sf_lnpoch_sgn_e'Access, (-4.5, 1.25), 2.1899306304483174731, TEST_TOL1, -1.0, "gsl_sf_lnpoch_sgn_e");

//     PUT_LINE("    ... gsl_sf_poch_e");
//     TEST_SF(s,  gsl_sf_poch_e'Access, (5.0, 0.0), 1.0, TEST_TOL0, "gsl_sf_poch_e");
//     TEST_SF(s,  gsl_sf_poch_e'Access, (7.0, 3.0), 504.0 , TEST_TOL0, "gsl_sf_poch_e");
//     TEST_SF(s,  gsl_sf_poch_e'Access, (5.0, 2.0), 30.0  , TEST_TOL1, "gsl_sf_poch_e");
//     TEST_SF(s,  gsl_sf_poch_e'Access, (5.0, 1.0/256.0), 1.0059023106151364982 , TEST_TOL0, "gsl_sf_poch_e");

//     PUT_LINE("    ... gsl_sf_pochrel_e");
//     TEST_SF(s,  gsl_sf_pochrel_e'Access, (5.0, 0.0), 1.506117668431800472, TEST_TOL1, "gsl_sf_pochrel_e");
//     TEST_SF(s,  gsl_sf_pochrel_e'Access, (7.0, 3.0), 503.0/3.0, TEST_TOL0, "gsl_sf_pochrel_e");
//     TEST_SF(s,  gsl_sf_pochrel_e'Access, (5.0, 2.0), 29.0/2.0, TEST_TOL1, "gsl_sf_pochrel_e");
//     TEST_SF(s,  gsl_sf_pochrel_e'Access, (5.0, 0.01), 1.5186393661368275330, TEST_TOL2, "gsl_sf_pochrel_e");
//     TEST_SF(s,  gsl_sf_pochrel_e'Access, (-5.5, 0.01), 1.8584945633829063516, TEST_TOL1, "gsl_sf_pochrel_e");
//     TEST_SF(s,  gsl_sf_pochrel_e'Access, (-5.5, -1.0/8.0), 1.0883319303552135488, TEST_TOL1, "gsl_sf_pochrel_e");
//     --TEST_SF(s,  gsl_sf_pochrel_e'Access, (-5.5, -1.0/256.0), 1.7678268037726177453, TEST_TOL1, "gsl_sf_pochrel_e");
//     TEST_SF(s,  gsl_sf_pochrel_e'Access, (-5.5, -11.0), 0.09090909090939652475, TEST_TOL0, "gsl_sf_pochrel_e");

//     PUT_LINE("    ... gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (1.0e-100, 0.001), 1.0, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (0.001, 0.001), 0.9936876467088602902, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (0.001, 1.0), 0.9997803916424144436, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (0.001, 10.0), 0.9999999958306921828, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (1.0, 0.001), 0.0009995001666250083319, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (1.0, 1.01), 0.6357810204284766802, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (1.0, 10.0), 0.9999546000702375151, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (10.0, 10.01), 0.5433207586693410570, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (10.0, 20.0), 0.9950045876916924128, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (1000.0, 1000.1), 0.5054666401440661753, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (1000.0, 2000.0), 1.0, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
//     -- Test for failure of the Gautschi recurrence (now fixed) for x = a - 2 */
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (34.0, 32.0), 0.3849626436463866776322932129, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
//     -- and the next test is gamma_inc_P(37,35-20*eps) */
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (37.0, 3.499999999999999289e+01), 0.3898035054195570860969333039, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
  
//     -- Regression test Martin Jansche <jansche@ling.ohio-state.edu> BUG#12 */
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (10.0, 1.0e-16), 2.755731922398588814734648067e-167, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
  
//     -- Regression test for gsl_cdf_chisq_Pinv, (0.05, 1263131.0) */
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (1263131.0, 1261282.3637), 0.04994777516935182963821362168, TEST_TOL4, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (1263131.0, 1263131.0), 0.500118321758657770672882362502514254, TEST_TOL4, "gsl_sf_gamma_inc_P_e");

//     PUT_LINE("    ... gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (0.0, 0.001), 0.0, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (0.001, 0.001), 0.006312353291139709793, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (0.001, 1.0), 0.00021960835758555639171, TEST_TOL1, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (0.001, 2.0), 0.00004897691783098147880, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (0.001, 5.0), 1.1509813397308608541e-06, TEST_TOL1, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1.0, 0.001), 0.9990004998333749917, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1.0, 1.01), 0.3642189795715233198, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1.0, 10.0), 0.00004539992976248485154, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (10.0, 10.01), 0.4566792413306589430, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (10.0, 100.0), 1.1253473960842733885e-31, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1000.0, 1000.1), 0.4945333598559338247, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1000.0, 2000.0), 6.847349459614753180e-136, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");


//     -- designed to trap the a-x=1 problem
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (100.0,  99.0), 0.5266956696005394, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (200.0, 199.0), 0.5188414119121281, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (100.0,  99.0), 0.4733043303994607, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (200.0, 199.0), 0.4811585880878718, TEST_TOL2, "gsl_sf_gamma_inc_P_e");

//     -- Test for x86 cancellation problems
//     TEST_SF(s, gsl_sf_gamma_inc_P_e'Access, (5670.0, 4574.0),  3.063972328743934e-55, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (5670.0, 4574.0), 1.0000000000000000, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");

//     -- test suggested by Michel Lespinasse [gsl-discuss Sat, 13 Nov 2004]
//     --TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1.0e+06-1.0, 1.0e+06-2.0), 0.50026596175224547004, TEST_TOL3, "gsl_sf_gamma_inc_Q_e");

//     -- tests in asymptotic regime related to Lespinasse test
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1.0e+06+2.0, 1.0e+06+1.0), 0.50026596135330304336, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1.0e+06, 1.0e+06-2.0), 0.50066490399940144811, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
//     TEST_SF(s, gsl_sf_gamma_inc_Q_e'Access, (1.0e+07, 1.0e+07-2.0), 0.50021026104978614908, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");

//     -- non-normalized "Q" function
//     PUT_LINE("    ... gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-1.0/1048576.0, 1.0/1048576.0), 13.285819596290624271, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.001, 1.0/1048576.0), 13.381275128625328858, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-1.0,   1.0/1048576.0), 1.0485617142715768655e+06, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.00001,0.001), 6.3317681434563592142, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.0001,0.001), 6.3338276439767189385, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.001, 0.001), 6.3544709102510843793, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.5,   0.001), 59.763880515942196981, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-1.0,   0.001), 992.66896046923884234, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-3.5,   0.001), 9.0224404490639003706e+09, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-10.5,  0.001), 3.0083661558184815656e+30, TEST_TOL2, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.001, 0.1), 1.8249109609418620068, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.5,   0.1), 3.4017693366916154163, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-10.0,  0.1), 8.9490757483586989181e+08, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-10.5,  0.1), 2.6967403834226421766e+09, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.001, 1.0), 0.21928612679072766340, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.5,   1.0), 0.17814771178156069019, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-1.0,   1.0), 0.14849550677592204792, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-2.5,   1.0), 0.096556648631275160264, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-1.0,   10.0), 3.8302404656316087616e-07, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.001, 10.0), 4.1470562324807320961e-06, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-0.5,   10.0), 1.2609042613241570681e-06, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-1.0,   10.0), 3.8302404656316087616e-07, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-10.5,  10.0), 6.8404927328441566785e-17, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-100.0, 10.0), 4.1238327669858313997e-107, TEST_TOL2, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (-200.0, 10.0), 2.1614091830529343423e-207, TEST_TOL2, "gsl_sf_gamma_inc_e");

//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  0.0,     0.001), 6.3315393641361493320, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  0.001,   0.001), 6.3087159394864007261, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  1.0,     0.001), 0.99900049983337499167, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, ( 10.0,     0.001), 362880.0, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  0.0,     1.0), 0.21938393439552027368, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  0.001,   1.0), 0.21948181320730279613, TEST_TOL1, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  1.0,     1.0), 0.36787944117144232160, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, ( 10.0,     1.0), 362879.95956592242045, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (100.0,     1.0), 9.3326215443944152682e+155, TEST_TOL0, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  0.0,   100.0), 3.6835977616820321802e-46, TEST_TOL2, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  0.001, 100.0), 3.7006367674063550631e-46, TEST_TOL2, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (  1.0,   100.0), 3.7200759760208359630e-44, TEST_TOL2, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, ( 10.0,   100.0), 4.0836606309106112723e-26, TEST_TOL2, "gsl_sf_gamma_inc_e");
//     TEST_SF(s, gsl_sf_gamma_inc_e'Access, (100.0,   100.0), 4.5421981208626694294e+155, TEST_TOL1, "gsl_sf_gamma_inc_e");


//     PUT_LINE("    ... gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (1.0e-8, 1.0e-8),  19.113827924512310617 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (1.0e-8, 0.01),  18.420681743788563403 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (1.0e-8, 1.0),  18.420680743952365472 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (1.0e-8, 10.0),  18.420680715662683009 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (1.0e-8, 1000.0),  18.420680669107656949 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (0.1, 0.1), 2.9813614810376273949 , TEST_TOL1, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (0.1, 1.0),  2.3025850929940456840 , TEST_TOL1, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (0.1, 100.0),  1.7926462324527931217 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     --TEST_SF(s, gsl_sf_lnbeta_e'Access, (0.1, 1000.0),  1.5619821298353164928 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (1.0, 1.00025),  -0.0002499687552073570, TEST_TOL4, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (1.0, 1.01),  -0.009950330853168082848 , TEST_TOL3, "gsl_sf_lnbeta_e");
//     --TEST_SF(s, gsl_sf_lnbeta_e'Access, (1.0, 1000.0),  -6.907755278982137052 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (100.0, 100.0),  -139.66525908670663927 , TEST_TOL2, "gsl_sf_lnbeta_e");
//     TEST_SF(s, gsl_sf_lnbeta_e'Access, (100.0, 1000.0),  -336.4348576477366051 , TEST_TOL0, "gsl_sf_lnbeta_e");
//     --TEST_SF(s, gsl_sf_lnbeta_e'Access, (100.0, 1.0e+8),  -1482.9339185256447309 , TEST_TOL0, "gsl_sf_lnbeta_e");

//     PUT_LINE("    ... gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (1.0,   1.0), 1.0                   , TEST_TOL0, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (1.0, 1.001), 0.9990009990009990010 , TEST_TOL0, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (1.0,   5.0), 0.2                   , TEST_TOL1, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (1.0,  100.0), 0.01                  , TEST_TOL1, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (10.0, 100.0), 2.3455339739604649879e-15 , TEST_TOL2, "gsl_sf_beta_e");

//     -- Test negative arguments
//     TEST_SF(s, gsl_sf_beta_e'Access, (2.5, -0.1), -11.43621278354402041480, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (2.5, -1.1), 14.555179906328753255202, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (-0.25, -0.1), -13.238937960945229110, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (-1.25, -0.1), -14.298052997820847439, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (-100.1, -99.1), -1.005181917797644630375787297e60, TEST_TOL3, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (-100.1, 99.3), 0.0004474258199579694011200969001, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (100.1, -99.3), 1.328660939628876472028853747, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (-100.1, 1.2), 0.00365530364287960795444856281, TEST_TOL3, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (100.1, -1.2), 1203.895236907821059270698160, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (-100.1, -1.2), -3236.073671884748847700283841, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (-100.001, 0.0099), -853.946649365611147996495177, TEST_TOL4, "gsl_sf_beta_e");

//     -- Other test cases
//     TEST_SF(s, gsl_sf_beta_e'Access, (1.0e-32, 1.5), 1.0e32, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (1.0e-6, 0.5), 1000001.386293677092419390336, TEST_TOL2, "gsl_sf_beta_e");
//     TEST_SF(s, gsl_sf_beta_e'Access, (-1.5, 0.5), 0.0, TEST_TOL0, "gsl_sf_beta_e");

//     PUT_LINE("    ... gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (1.0, 1.0, 0.0), 0.0, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (1.0, 1.0, 1.0), 1.0, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (0.1, 0.1, 1.0), 1.0, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 1.0,  1.0, 0.5), 0.5, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 0.1,  1.0, 0.5), 0.9330329915368074160, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (10.0,  1.0, 0.5), 0.0009765625000000000000, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (50.0,  1.0, 0.5), 8.881784197001252323e-16, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 1.0,  0.1, 0.5), 0.06696700846319258402, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 1.0, 10.0, 0.5), 0.99902343750000000000, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 1.0, 50.0, 0.5), 0.99999999999999911180, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 1.0,  1.0, 0.1), 0.10, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 1.0,  2.0, 0.1), 0.19, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 1.0,  2.0, 0.9), 0.99, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (50.0, 60.0, 0.5), 0.8309072939016694143, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (90.0, 90.0, 0.5), 0.5, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, ( 500.0,  500.0, 0.6), 0.9999999999157549630, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (5000.0, 5000.0, 0.4), 4.518543727260666383e-91, TEST_TOL5, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (5000.0, 5000.0, 0.6), 1.0, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s, gsl_sf_beta_inc_e'Access, (5000.0, 2000.0, 0.6), 8.445388773903332659e-89, TEST_TOL5, "gsl_sf_beta_inc_e");

//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.1, -0.1, 1.0), 1.0, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.1, -0.2, 1.0), 1.0, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.2, -0.1, 1.0), 1.0, TEST_TOL2, "gsl_sf_beta_inc_e");

//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.1, -0.2, 0.5), 0.675252001958389971991335, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.2, -0.1, 0.5), 0.324747998041610028008665, TEST_TOL2, "gsl_sf_beta_inc_e");

//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.1, -0.1, 0.0), 0.0, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.1, -0.2, 0.0), 0.0, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.2, -0.1, 0.0), 0.0, TEST_TOL2, "gsl_sf_beta_inc_e");

//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.1, -0.2, 0.3), 0.7469186777964287252, TEST_TOL2, "gsl_sf_beta_inc_e");
//     TEST_SF(s,  gsl_sf_beta_inc_e'Access, (-0.2, -0.1, 0.3), 0.3995299653262016818, TEST_TOL2, "gsl_sf_beta_inc_e");

    return s;

} // test_gamma

// ----------------------------------------------------------------------------
// EOF SF-Test_test_gamma.mjs
