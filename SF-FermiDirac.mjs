// SF-FermiDirac.mjs
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

import { M_PI }                  from "./SF-Math.mjs";
import { GSL_DBL_EPSILON }       from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_EPSILON }  from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MIN }       from "./SF-Machine.mjs";
import { GSL_DBL_MAX }           from "./SF-Machine.mjs";
import { GSL_SQRT_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_ROOT3_DBL_MAX }     from "./SF-Machine.mjs";
import { gsl_sf_lngamma_e }      from "./SF-Gamma.mjs";
import { gsl_sf_fact_e }         from "./SF-Gamma.mjs";
import { gsl_sf_exp_err_e }      from "./SF-Exponential.mjs";
import { gsl_sf_exp_mult_err_e } from "./SF-Exponential.mjs";
import { cheb_eval_e }           from "./SF-Chebyshev.mjs";
import { gsl_sf_eta_int_e }      from "./SF-Zeta.mjs";
import { gsl_sf_hyperg_1F1_int_e } from "./SF-Hypergeometric1F1.mjs";
import { gsl_sf_hyperg_U_int_e } from "./SF-HypergeometricU.mjs";

import { EVAL_RESULT_D }         from "./SF-Evaluate.mjs";
import { EVAL_RESULT_ID }        from "./SF-Evaluate.mjs";
import { EVAL_RESULT_DD }        from "./SF-Evaluate.mjs";

// ----------------------------------------------------------------------------

// Chebyshev fit for F_{1}(t);  -1 < t < 1, -1 < x < 1
//
const fd_1_a_data =
    [
    1.8949340668482264365,
    0.7237719066890052793,
    0.1250000000000000000,
    0.0101065196435973942,
    0.0,
   -0.0000600615242174119,
    0.0,
    6.816528764623e-7,
    0.0,
   -9.5895779195e-9,
    0.0,
    1.515104135e-10,
    0.0,
   -2.5785616e-12,
    0.0,
    4.62270e-14,
    0.0,
   -8.612e-16,
    0.0,
    1.65e-17,
    0.0,
   -3.0e-19
    ];
const fd_1_a_cs = { length: 21, c: fd_1_a_data, order: 21, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{1}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
//
const fd_1_b_data =
    [
    10.409136795234611872,
    3.899445098225161947,
    0.513510935510521222,
    0.010618736770218426,
   -0.001584468020659694,
    0.000146139297161640,
   -1.408095734499e-6,
   -2.177993899484e-6,
    3.91423660640e-7,
   -2.3860262660e-8,
   -4.138309573e-9,
    1.283965236e-9,
   -1.39695990e-10,
   -4.907743e-12,
    4.399878e-12,
   -7.17291e-13,
    2.4320e-14,
    1.4230e-14,
   -3.446e-15,
    2.93e-16,
    3.7e-17,
   -1.6e-17
    ];
const fd_1_b_cs = { length: 21, c: fd_1_b_data, order: 21, a: -1.0, b: 1.0, order_sp: 11 };


// Chebyshev fit for F_{1}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
//
const fd_1_c_data =
    [
    56.78099449124299762,
    21.00718468237668011,
    2.24592457063193457,
    0.00173793640425994,
   -0.00058716468739423,
    0.00016306958492437,
   -0.00003817425583020,
    7.64527252009e-6,
   -1.31348500162e-6,
    1.9000646056e-7,
   -2.141328223e-8,
    1.23906372e-9,
    2.1848049e-10,
   -1.0134282e-10,
    2.484728e-11,
   -4.73067e-12,
    7.3555e-13,
   -8.740e-14,
    4.85e-15,
    1.23e-15,
   -5.6e-16,
    1.4e-16,
   -3.0e-17
    ];
const fd_1_c_cs = { length: 22, c: fd_1_c_data, order: 22, a: -1.0, b: 1.0, order_sp: 13 };


// Chebyshev fit for F_{1}(x) / x^2
// 10 < x < 30 
// -1 < t < 1
// t = 1/10 (x-10) - 1 = x/10 - 2
// x = 10(t+2)
//
const fd_1_d_data =
    [
    1.0126626021151374442,
   -0.0063312525536433793,
    0.0024837319237084326,
   -0.0008764333697726109,
    0.0002913344438921266,
   -0.0000931877907705692,
    0.0000290151342040275,
   -8.8548707259955e-6,
    2.6603474114517e-6,
   -7.891415690452e-7,
    2.315730237195e-7,
   -6.73179452963e-8,
    1.94048035606e-8,
   -5.5507129189e-9,
    1.5766090896e-9,
   -4.449310875e-10,
    1.248292745e-10,
   -3.48392894e-11,
    9.6791550e-12,
   -2.6786240e-12,
    7.388852e-13,
   -2.032828e-13,
    5.58115e-14,
   -1.52987e-14,
    4.1886e-15,
   -1.1458e-15,
    3.132e-16,
   -8.56e-17,
    2.33e-17,
   -5.9e-18
    ];
const fd_1_d_cs = { length: 29, c: fd_1_d_data, order: 29, a:-1.0, b: 1.0, order_sp: 14 };


// Chebyshev fit for F_{1}(x) / x^2
// 30 < x < Inf
// -1 < t < 1
// t = 60/x - 1
// x = 60/(t+1)
//
const fd_1_e_data =
    [
    1.0013707783890401683,
    0.0009138522593601060,
    0.0002284630648400133,
   -1.57e-17,
   -1.27e-17,
   -9.7e-18,
   -6.9e-18,
   -4.6e-18,
   -2.9e-18,
   -1.7e-18
    ];
const fd_1_e_cs = { length: 9, c: fd_1_e_data, order: 9, a: -1.0, b: 1.0, order_sp: 4 };


// Chebyshev fit for F_{2}(t);  -1 < t < 1, -1 < x < 1
//
const fd_2_a_data =
    [
    2.1573661917148458336,
    0.8849670334241132182,
    0.1784163467613519713,
    0.0208333333333333333,
    0.0012708226459768508,
    0.0,
   -5.0619314244895e-6,
    0.0,
    4.32026533989e-8,
    0.0,
   -4.870544166e-10,
    0.0,
    6.4203740e-12,
    0.0,
   -9.37424e-14,
    0.0,
    1.4715e-15,
    0.0,
   -2.44e-17,
    0.0,
    4.0e-19
    ];
const fd_2_a_cs = { length: 20, c: fd_2_a_data, order: 20, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{2}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
//
const fd_2_b_data =
    [
    16.508258811798623599,
    7.421719394793067988,
    1.458309885545603821,
    0.128773850882795229,
    0.001963612026198147,
   -0.000237458988738779,
    0.000018539661382641,
   -1.92805649479e-7,
   -2.01950028452e-7,
    3.2963497518e-8,
   -1.885817092e-9,
   -2.72632744e-10,
    8.0554561e-11,
   -8.313223e-12,
   -2.24489e-13,
    2.18778e-13,
   -3.4290e-14,
    1.225e-15,
    5.81e-16,
   -1.37e-16,
    1.2e-17,
    1.0e-18
    ];
const fd_2_b_cs = { length: 21, c: fd_2_b_data, order: 21, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{1}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
//
const fd_2_c_data =
    [
    168.87129776686440711,
    81.80260488091659458,
    15.75408505947931513,
    1.12325586765966440,
    0.00059057505725084,
   -0.00016469712946921,
    0.00003885607810107,
   -7.89873660613e-6,
    1.39786238616e-6,
   -2.1534528656e-7,
    2.831510953e-8,
   -2.94978583e-9,
    1.6755082e-10,
    2.234229e-11,
   -1.035130e-11,
    2.41117e-12,
   -4.3531e-13,
    6.447e-14,
   -7.39e-15,
    4.3e-16
    ];
const fd_2_c_cs = { length: 19, c: fd_2_c_data, order: 19, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{1}(x) / x^3
// 10 < x < 30 
// -1 < t < 1
// t = 1/10 (x-10) - 1 = x/10 - 2
// x = 10(t+2)
//
const fd_2_d_data =
    [
    0.3459960518965277589,
   -0.00633136397691958024,
    0.00248382959047594408,
   -0.00087651191884005114,
    0.00029139255351719932,
   -0.00009322746111846199,
    0.00002904021914564786,
   -8.86962264810663e-6,
    2.66844972574613e-6,
   -7.9331564996004e-7,
    2.3359868615516e-7,
   -6.824790880436e-8,
    1.981036528154e-8,
   -5.71940426300e-9,
    1.64379426579e-9,
   -4.7064937566e-10,
    1.3432614122e-10,
   -3.823400534e-11,
    1.085771994e-11,
   -3.07727465e-12,
    8.7064848e-13,
   -2.4595431e-13,
    6.938531e-14,
   -1.954939e-14,
    5.50162e-15,
   -1.54657e-15,
    4.3429e-16,
   -1.2178e-16,
    3.394e-17,
   -8.81e-18
    ];
const fd_2_d_cs = { length: 29, c: fd_2_d_data, order: 29, a: -1.0, b: 1.0, order_sp: 14 };


// Chebyshev fit for F_{2}(x) / x^3
// 30 < x < Inf
// -1 < t < 1
// t = 60/x - 1
// x = 60/(t+1)
//
const fd_2_e_data =
    [
    0.3347041117223735227,
    0.00091385225936012645,
    0.00022846306484003205,
    5.2e-19
    ];
const fd_2_e_cs = { length: 3, c: fd_2_e_data, order: 3, a: -1.0, b: 1.0, order_sp: 3 };


// Chebyshev fit for F_{-1/2}(t);  -1 < t < 1, -1 < x < 1
//
const fd_mhalf_a_data =//: CONSTANT Series(0..19) = --[20]
    [
    1.2663290042859741974,
    0.3697876251911153071,
    0.0278131011214405055,
   -0.0033332848565672007,
   -0.0004438108265412038,
    0.0000616495177243839,
    8.7589611449897e-6,
   -1.2622936986172e-6,
   -1.837464037221e-7,
    2.69495091400e-8,
    3.9760866257e-9,
   -5.894468795e-10,
   -8.77321638e-11,
    1.31016571e-11,
    1.9621619e-12,
   -2.945887e-13,
   -4.43234e-14,
    6.6816e-15,
    1.0084e-15,
   -1.561e-16
    ];
const fd_mhalf_a_cs = { length: 19, c: fd_mhalf_a_data, order: 19, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{-1/2}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
//
const fd_mhalf_b_data =
    [
    3.270796131942071484,
    0.5809004935853417887,
   -0.0299313438794694987,
   -0.0013287935412612198,
    0.0009910221228704198,
   -0.0001690954939688554,
    6.5955849946915e-6,
    3.5953966033618e-6,
   -9.430672023181e-7,
    8.75773958291e-8,
    1.06247652607e-8,
   -4.9587006215e-9,
    7.160432795e-10,
    4.5072219e-12,
   -2.3695425e-11,
    4.9122208e-12,
   -2.905277e-13,
   -9.59291e-14,
    3.00028e-14,
   -3.4970e-15
    ];
const fd_mhalf_b_cs = { length: 19, c: fd_mhalf_b_data, order: 19, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{-1/2}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
//
const fd_mhalf_c_data =
    [
    5.828283273430595507,
    0.677521118293264655,
   -0.043946248736481554,
    0.005825595781828244,
   -0.000864858907380668,
    0.000110017890076539,
   -6.973305225404e-6,
   -1.716267414672e-6,
    8.59811582041e-7,
   -2.33066786976e-7,
    4.8503191159e-8,
   -8.130620247e-9,
    1.021068250e-9,
   -5.3188423e-11,
   -1.9430559e-11,
    8.750506e-12,
   -2.324897e-12,
    4.83102e-13,
   -8.1207e-14,
    1.0132e-14,
   -4.64e-16,
   -2.24e-16,
    9.7e-17,
   -2.6e-17,
    5.0e-18
    ];
const fd_mhalf_c_cs = { length: 24, c: fd_mhalf_c_data, order: 24, a: -1.0, b: 1.0, order_sp: 13 };


// Chebyshev fit for F_{-1/2}(x) / x^(1/2)
// 10 < x < 30 
// -1 < t < 1
// t = 1/10 (x-10) - 1 = x/10 - 2
//
const fd_mhalf_d_data =
    [
    2.2530744202862438709,
    0.0018745152720114692,
   -0.0007550198497498903,
    0.0002759818676644382,
   -0.0000959406283465913,
    0.0000324056855537065,
   -0.0000107462396145761,
    3.5126865219224e-6,
   -1.1313072730092e-6,
    3.577454162766e-7,
   -1.104926666238e-7,
    3.31304165692e-8,
   -9.5837381008e-9,
    2.6575790141e-9,
   -7.015201447e-10,
    1.747111336e-10,
   -4.04909605e-11,
    8.5104999e-12,
   -1.5261885e-12,
    1.876851e-13,
    1.00574e-14,
   -1.82002e-14,
    8.6634e-15,
   -3.2058e-15,
    1.0572e-15,
   -3.259e-16,
    9.60e-17,
   -2.74e-17,
    7.6e-18,
   -1.9e-18
    ];
const fd_mhalf_d_cs = { length: 29, c: fd_mhalf_d_data, order: 29, a: -1.0, b: 1.0, order_sp: 15 };


// Chebyshev fit for F_{1/2}(t);  -1 < t < 1, -1 < x < 1
//
const fd_half_a_data =
    [
    1.7177138871306189157,
    0.6192579515822668460,
    0.0932802275119206269,
    0.0047094853246636182,
   -0.0004243667967864481,
   -0.0000452569787686193,
    5.2426509519168e-6,
    6.387648249080e-7,
   -8.05777004848e-8,
   -1.04290272415e-8,
    1.3769478010e-9,
    1.847190359e-10,
   -2.51061890e-11,
   -3.4497818e-12,
    4.784373e-13,
    6.68828e-14,
   -9.4147e-15,
   -1.3333e-15,
    1.898e-16,
    2.72e-17,
   -3.9e-18,
   -6.0e-19,
    1.0e-19
    ];
const fd_half_a_cs = { length: 22, c: fd_half_a_data, order: 22, a: -1.0, b: 1.0, order_sp:11 };


// Chebyshev fit for F_{1/2}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
//
const fd_half_b_data =
    [
    7.651013792074984027,
    2.475545606866155737,
    0.218335982672476128,
   -0.007730591500584980,
   -0.000217443383867318,
    0.000147663980681359,
   -0.000021586361321527,
    8.07712735394e-7,
    3.28858050706e-7,
   -7.9474330632e-8,
    6.940207234e-9,
    6.75594681e-10,
   -3.10200490e-10,
    4.2677233e-11,
   -2.1696e-14,
   -1.170245e-12,
    2.34757e-13,
   -1.4139e-14,
   -3.864e-15,
    1.202e-15
    ];
const fd_half_b_cs = { length: 19, c: fd_half_b_data, order: 19, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{1/2}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
//
const fd_half_c_data =
    [
    29.584339348839816528,
    8.808344283250615592,
    0.503771641883577308,
   -0.021540694914550443,
    0.002143341709406890,
   -0.000257365680646579,
    0.000027933539372803,
   -1.678525030167e-6,
   -2.78100117693e-7,
    1.35218065147e-7,
   -3.3740425009e-8,
    6.474834942e-9,
   -1.009678978e-9,
    1.20057555e-10,
   -6.636314e-12,
   -1.710566e-12,
    7.75069e-13,
   -1.97973e-13,
    3.9414e-14,
   -6.374e-15,
    7.77e-16,
   -4.0e-17,
   -1.4e-17
    ];
const fd_half_c_cs = { length: 22, c: fd_half_c_data, order: 22, a: -1.0, b: 1.0, order_sp: 13 };


// Chebyshev fit for F_{1/2}(x) / x^(3/2)
// 10 < x < 30 
// -1 < t < 1
// t = 1/10 (x-10) - 1 = x/10 - 2
//
const fd_half_d_data =
    [
    1.5116909434145508537,
   -0.0036043405371630468,
    0.0014207743256393359,
   -0.0005045399052400260,
    0.0001690758006957347,
   -0.0000546305872688307,
    0.0000172223228484571,
   -5.3352603788706e-6,
    1.6315287543662e-6,
   -4.939021084898e-7,
    1.482515450316e-7,
   -4.41552276226e-8,
    1.30503160961e-8,
   -3.8262599802e-9,
    1.1123226976e-9,
   -3.204765534e-10,
    9.14870489e-11,
   -2.58778946e-11,
    7.2550731e-12,
   -2.0172226e-12,
    5.566891e-13,
   -1.526247e-13,
    4.16121e-14,
   -1.12933e-14,
    3.0537e-15,
   -8.234e-16,
    2.215e-16,
   -5.95e-17,
    1.59e-17,
   -4.0e-18
    ];
const fd_half_d_cs = { length: 29, c: fd_half_d_data, order: 29, a: -1.0, b: 1.0, order_sp: 15 };



// Chebyshev fit for F_{3/2}(t);  -1 < t < 1, -1 < x < 1
//
const fd_3half_a_data =
    [
    2.0404775940601704976,
    0.8122168298093491444,
    0.1536371165644008069,
    0.0156174323847845125,
    0.0005943427879290297,
   -0.0000429609447738365,
   -3.8246452994606e-6,
    3.802306180287e-7,
    4.05746157593e-8,
   -4.5530360159e-9,
   -5.306873139e-10,
    6.37297268e-11,
    7.8403674e-12,
   -9.840241e-13,
   -1.255952e-13,
    1.62617e-14,
    2.1318e-15,
   -2.825e-16,
   -3.78e-17,
    5.1e-18
    ];
const fd_3half_a_cs = { length: 19, c: fd_3half_a_data, order: 19, a: -1.0, b: 1.0, order_sp: 11 };


// Chebyshev fit for F_{3/2}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
//
const fd_3half_b_data =
    [
    13.403206654624176674,
    5.574508357051880924,
    0.931228574387527769,
    0.054638356514085862,
   -0.001477172902737439,
   -0.000029378553381869,
    0.000018357033493246,
   -2.348059218454e-6,
    8.3173787440e-8,
    2.6826486956e-8,
   -6.011244398e-9,
    4.94345981e-10,
    3.9557340e-11,
   -1.7894930e-11,
    2.348972e-12,
   -1.2823e-14,
   -5.4192e-14,
    1.0527e-14,
   -6.39e-16,
   -1.47e-16,
    4.5e-17,
   -5.0e-18
    ];
const fd_3half_b_cs = { length: 21, c: fd_3half_b_data, order: 21, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{3/2}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
//
const fd_3half_c_data =
    [
    101.03685253378877642,
    43.62085156043435883,
    6.62241373362387453,
    0.25081415008708521,
   -0.00798124846271395,
    0.00063462245101023,
   -0.00006392178890410,
    6.04535131939e-6,
   -3.4007683037e-7,
   -4.072661545e-8,
    1.931148453e-8,
   -4.46328355e-9,
    7.9434717e-10,
   -1.1573569e-10,
    1.304658e-11,
   -7.4114e-13,
   -1.4181e-13,
    6.491e-14,
   -1.597e-14,
    3.05e-15,
   -4.8e-16
    ];
const fd_3half_c_cs = { length: 20, c: fd_3half_c_data, order: 20, a: -1.0, b: 1.0, order_sp: 12 };


// Chebyshev fit for F_{3/2}(x) / x^(5/2)
// 10 < x < 30 
// -1 < t < 1
// t = 1/10 (x-10) - 1 = x/10 - 2
//
const fd_3half_d_data =
    [
    0.6160645215171852381,
   -0.0071239478492671463,
    0.0027906866139659846,
   -0.0009829521424317718,
    0.0003260229808519545,
   -0.0001040160912910890,
    0.0000322931223232439,
   -9.8243506588102e-6,
    2.9420132351277e-6,
   -8.699154670418e-7,
    2.545460071999e-7,
   -7.38305056331e-8,
    2.12545670310e-8,
   -6.0796532462e-9,
    1.7294556741e-9,
   -4.896540687e-10,
    1.380786037e-10,
   -3.88057305e-11,
    1.08753212e-11,
   -3.0407308e-12,
    8.485626e-13,
   -2.364275e-13,
    6.57636e-14,
   -1.81807e-14,
    4.6884e-15
    ];
const fd_3half_d_cs = { length: 24, c: fd_3half_d_data, order: 24, a: -1.0, b: 1.0, order_sp: 16 };

// qsize : CONSTANT INTEGER = 100 + 1;
// TYPE ANum IS ARRAY(0 .. qsize - 1) OF LONG_FLOAT;
// TYPE ADen IS ARRAY(0 .. qsize - 1) OF LONG_FLOAT;

// Goano's modification of the Levin-u implementation.
// This is a simplification of the original WHIZ algorithm.
// See [Fessler et al., ACM Toms 9, 346 (1983)].
//
function fd_whiz(term, iterm, r0) //, qnum, qden, result, s)
{
    var factor = 0.0;
    var ratio  = 0.0;
    var c      = 0.0;
    var r      = { qnum: [], qden: [], result: 0.0, s: 0.0 };

    r.qnum = r0.qnum;
    r.qden = r0.qden;
    r.result = r0.result;
    r.s = r0.s;

    if (iterm == 0)
    {
        r.s = 0.0;
    }
  
    r.s = r.s + term;
  
    r.qden[iterm] = 1.0 / (term * (iterm + 1) * (iterm + 1));
    r.qnum[iterm] = r.s * r.qden[iterm];
  
    if (iterm > 0)
    {
        factor = 1.0;
        ratio  = (iterm) / (iterm + 1);
        for (let j = iterm - 1; j >= 0; j--)
        {
            c = factor * (j + 1) / (iterm + 1);
            factor = factor * ratio;
            r.qden[j] = r.qden[j+1] - c * r.qden[j];
            r.qnum[j] = r.qnum[j+1] - c * r.qnum[j];
        }
    }
  
    r.result = r.qnum[0] / r.qden[0];

    // r.qnum = qnum;
    // r.qden = qden;
    // r.result = result;
    // r.s = s;
    return r;

} // fd_whiz

// ----------------------------------------------------------------------------

// Handle case of integer j <= -2.
//
function fd_nint( j, x )
{
    const nsize  = 100 + 1;
    var qcoeff   = [];
    var r        = { val: 0.0, err: 0.0 };

    if ( j >= -1 )
    {
        throw "SF.SanityException";
    }
    else if ( j < -(nsize) )
    {
        throw "SF.NotImplementedException";
    }
    else
    {
        let a = 0.0;
        let p = 0.0;
        let f = 0.0;
        let n = -(j + 1);

        qcoeff[1] = 1.0;
        
        for ( let k = 2; k <= n; k++ )
        {
            qcoeff[k] = -qcoeff[k-1];
            for ( let i = k - 1; i >= 2; i-- )
            {
                qcoeff[i] = (i) * qcoeff[i] - (k - (i - 1)) * qcoeff[i - 1];
            }
        }
        
        if ( x >= 0.0 )
        {
            a = Math.exp( -x );
            f = qcoeff[1];
            for ( let i = 2; i <= n; i++ )
            {
                f = f * a + qcoeff[i];
            }
        }
        else
        {
            a = Math.exp( x );
            f = qcoeff[n];
            for ( let i = n - 1; i >= 1; i-- )
            {
                f = f * a + qcoeff[i];
            }
        }
        
        p = Math.pow_factor( 1.0 + a, j );
        r.val = f * a * p;
        r.err = 3.0 * GSL_DBL_EPSILON * Math.abs( f * a * p );
        return r;
    }

} // fd_nint

// ----------------------------------------------------------------------------

// x < 0
//
function fd_neg(j, x)
{
    const itmax = 100;
    var qnum  = []; //: ANum = (OTHERS => 0.0);
    var qden  = []; //: ADen = (OTHERS => 0.0);

    var r     = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        r.val = 0.0;
        r.err = 0.0;
    }
    else if (x < -1.0 && x < -Math.abs(j + 1.0))
    {
        // Simple series implementation. Avoid the
        // complexity and extra work of the series
        // acceleration method below.
        //
        let ex   = Math.exp(x);
        let term = ex;
        let sum  = term;
        let rat  = 0.0;
        let p    = 0.0;

        for (let n = 2; n <= 100 - 1; n++)
        {
            rat = (n - 1) / (n);
            p   = Math.pow(rat, j + 1.0);
            term = term * (-ex * p);
            sum  = sum + term;
            if (Math.abs(term / sum) < GSL_DBL_EPSILON)
            {
                break;
            }
        }
        r.val = sum;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(sum);
    }
    else
    {
        let s          = 0.0;
        let xn         = x;
        let ex         = -Math.exp(x);
        let enx        = -ex;
        let f          = 0.0;
        let f_previous = 0.0;
        let p          = 0.0;
        let term       = 0.0;
        let jterm      = 0;
        let rwhiz      = { qnum: [], qden: [], result: 0.0, s: 0.0 };

        jterm = 0;
        while ( jterm <= itmax )
        {
            p = Math.pow( jterm + 1, j + 1.0 );
            term = enx / p;
            f_previous = f;
            rwhiz = fd_whiz( term, jterm, rwhiz ); //, qnum, qden, f, s);
            qnum = rwhiz.qnum;
            qden = rwhiz.qden;
            f = rwhiz.result;
            s = rwhiz.s;
            xn = xn + x;
            if ( Math.abs( f - f_previous ) < Math.abs( f ) * 2.0 * GSL_DBL_EPSILON || xn < GSL_LOG_DBL_MIN )
            {
                break;
            }
            enx = enx * ex;
            jterm = jterm + 1;
        }
        
        r.val = f;
        r.err = Math.abs(f - f_previous);
        r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(f);
        
        if (jterm >= itmax)
        {
            throw "SF.MaxIterationsException";
        }
    }

    return r;

} // fd_neg

// ----------------------------------------------------------------------------

// asymptotic expansion
// j + 2.0 > 0.0
//
function fd_asymp(j, x)
{
    const j_integer = (Math.abs(j - Math.floor(j + 0.5)) < 100.0 * GSL_DBL_EPSILON);
    const itmax     = 200;
    var seqn_val  = 0.0;
    var seqn_err  = 0.0;
    var xm2       = 0.0;
    var xgam      = 0.0;
    var add       = 0.0;
    var cos_term  = 0.0;
    var ln_x      = 0.0;
    var ex_term_1 = 0.0;
    var ex_term_2 = 0.0;
    var add_previous = 0.0;
    var lg        = { val: 0.0, err: 0.0 }; // Result;
    var fneg      = { val: 0.0, err: 0.0 }; // Result;
    var ex_arg    = { val: 0.0, err: 0.0 }; // Result;
    var ex        = { val: 0.0, err: 0.0 }; // Result;
    var eta       = { val: 0.0, err: 0.0 }; // Result;
    var r         = { val: 0.0, err: 0.0 }; // Result;

    lg = gsl_sf_lngamma_e(j + 2.0);
    seqn_val = 0.5;
    seqn_err = 0.0;
    xm2  = (1.0 / x) / x;
    xgam = 1.0;
    add = GSL_DBL_MAX;

    for (let n = 1; n <= itmax; n++)
    {
        add_previous = add;
        eta = gsl_sf_eta_int_e(2 * n);
        xgam = xgam * xm2 * (j + 1.0 - (2 * n - 2)) * (j + 1.0 - (2 * n - 1));
        add  = eta.val * xgam;
        if (! j_integer && Math.abs(add) > Math.abs(add_previous))
        {
            break;
        }
        if (Math.abs(add / seqn_val) < GSL_DBL_EPSILON)
        {
            break;
        }
        seqn_val = seqn_val + add;
        seqn_err = seqn_err +2.0 * GSL_DBL_EPSILON * Math.abs(add);
    }
    seqn_err = seqn_err + Math.abs(add);
  
    fneg = fd_neg(j, -x);
    ln_x = Math.log(x);
    ex_term_1 = (j + 1.0) * ln_x;
    ex_term_2 = lg.val;
    ex_arg.val = ex_term_1 - ex_term_2; // (j+1.0)*ln_x - lg.val;
    ex_arg.err = GSL_DBL_EPSILON * (Math.abs(ex_term_1) + Math.abs(ex_term_2)) + lg.err;
    ex = gsl_sf_exp_err_e(ex_arg.val, ex_arg.err);
    cos_term  = Math.cos(j * M_PI);
    r.val = cos_term * fneg.val + 2.0 * seqn_val * ex.val;
    r.err = Math.abs(2.0 * ex.err * seqn_val);
    r.err = r.err + Math.abs(2.0 * ex.val * seqn_err);
    r.err = r.err + Math.abs(cos_term) * fneg.err;
    r.err = r.err + 4.0 * GSL_DBL_EPSILON * Math.abs(r.val);

    return r;

} // fd_asymp

// ----------------------------------------------------------------------------

// Series evaluation for small x > 0, integer j > 0; x < Pi.
// [Goano (8)]
//
function fd_series_int(j, x)
{
    var sum        = 0.0;
    var del        = 0.0;
    var pow_factor = 0.0;
    var sum2       = 0.0;
    var pre2       = 0.0;
    var eta_factor = { val: 0.0, err: 0.0 }; // Result;
    var jfact      = { val: 0.0, err: 0.0 }; // Result;
    var r          = { val: 0.0, err: 0.0 }; // Result;

    pow_factor = 1.0;
    eta_factor = gsl_sf_eta_int_e(j + 1);
    del = pow_factor * eta_factor.val;
    sum = sum + del;
  
    // Sum terms where the argument
    // of eta() is positive.
    //
    for (let n = 1; n <= j + 2; n++)
    {
        eta_factor = gsl_sf_eta_int_e(j + 1 - n);
        pow_factor = pow_factor * (x / (n));
        del = pow_factor * eta_factor.val;
        sum = sum + del;
        if (Math.abs(del / sum) < GSL_DBL_EPSILON)
        {
            break;
        }
    }
  
    // Now sum the terms where eta() is negative.
    // The argument of eta() must be odd as well,
    // so it is convenient to transform the series
    // as follows:
    //
    // Sum[ eta(j+1-n) x^n / n!, {n,j+4,Infinity}]
    // = x^j / j! Sum[ eta(1-2m) x^(2m) j! / (2m+j)! , {m,2,Infinity}]
    //
    // We do not need to do this sum if j is large enough.
    //
    if (j < 32)
    {
        jfact = gsl_sf_fact_e(j);
        pre2 = Math.pow(x, j) / jfact.val;
       
        eta_factor = gsl_sf_eta_int_e(-3);
        pow_factor = x * x * x * x / ((j + 4) * (j + 3) * (j + 2) * (j + 1));
        sum2 = eta_factor.val * pow_factor;
       
        for (let m = 3; m <= 24 - 1; m++)
        {
            eta_factor = gsl_sf_eta_int_e(1 - 2 * m);
            pow_factor = pow_factor * x * x / ((j + 2 * m) * (j + 2 * m - 1));
            sum2 = sum2 + eta_factor.val * pow_factor;
        }
       
        sum = sum + pre2 * sum2;
    }
  
    r.val = sum;
    r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(sum);
  
    return r;

} // fd_series_int

// ----------------------------------------------------------------------------

// series of hypergeometric functions for integer j > 0, x > 0
// [Goano (7)]
//
function fd_UMseries_int(j, x)
{
    const nmax         = 2000;
    var pre          = 0.0;
    var lnpre_val    = 0.0;
    var lnpre_err    = 0.0;
    var sum_even_val = 0.0;
    var sum_even_err = 0.0;
    var sum_odd_val  = 0.0;
    var sum_odd_err  = 0.0;
    var p            = 0.0;
    var lnx          = 0.0;
    var del_val      = 0.0;
    var del_err      = 0.0;
    var n            = 0;
    var g            = { val: 0.0, err: 0.0 }; // Result;
    var lg           = { val: 0.0, err: 0.0 }; // Result;
    var U            = { val: 0.0, err: 0.0 }; // Result;
    var M            = { val: 0.0, err: 0.0 }; // Result;
    var r            = { val: 0.0, err: 0.0 }; // Result;

    sum_even_val = 1.0;
  
    if (x < 500.0 && j < 80)
    {
        p = Math.pow(x, j + 1);
        g = gsl_sf_fact_e(j + 1); // Gamma(j+2)
        lnpre_val = 0.0;
        lnpre_err = 0.0;
        pre   = p / g.val;
    }
    else
    {
        lnx = Math.log(x);
        lg = gsl_sf_lngamma_e((j + 2));
        lnpre_val = (j + 1) * lnx - lg.val;
        lnpre_err = 2.0 * GSL_DBL_EPSILON * Math.abs((j + 1) * lnx) + lg.err;
        pre = 1.0;
    }
  
    // Add up the odd terms of the sum.
    //
    n = 1;
    for (let i = 1; i <= nmax / 2; i++)
    {
        U = gsl_sf_hyperg_U_int_e(1, j + 2, (n) * x);
        M = gsl_sf_hyperg_1F1_int_e(1, j + 2, -(n) * x);
        del_val = ((j + 1) * U.val - M.val);
        del_err = (Math.abs((j + 1)) * U.err + M.err);
        sum_odd_val = sum_odd_val + del_val;
        sum_odd_err = sum_odd_err + del_err;
        if (Math.abs(del_val / sum_odd_val) < GSL_DBL_EPSILON)
        {
            break;
        }
        n = n + 2;
    }
  
    // Add up the even terms of the sum.
    //
    n = 2;
    for (let i = 1; i <= nmax / 2; i++)
    {
        U = gsl_sf_hyperg_U_int_e(1, j + 2, (n) * x);
        M = gsl_sf_hyperg_1F1_int_e(1, j + 2, -(n) * x);
        del_val = ((j + 1) * U.val - M.val);
        del_err = (Math.abs((j + 1)) * U.err + M.err);
        sum_even_val = sum_even_val - del_val;
        sum_even_err = sum_even_err + del_err;
        if (Math.abs(del_val / sum_even_val) < GSL_DBL_EPSILON)
        {
            break;
        }
        n = n + 2;
    }
  
    if (n >= nmax)
    {
        throw "SF.MaxIterationsException";
    }
    r = gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err,
                          pre * (sum_even_val + sum_odd_val), pre * (sum_even_err + sum_odd_err));
    r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
  
    return r;

} // fd_UMseries_int

// //*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*

// [Goano (4)]
export function gsl_sf_fermi_dirac_m1_e(x)
{
    var ex = 0.0;

    var r  = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (x < 0.0)
    {
        ex = Math.exp(x);
        r.val = ex / (1.0 + ex);
        r.err = 2.0 * (Math.abs(x) + 1.0) * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        ex = Math.exp(-x);
        r.val = 1.0 / (1.0 + ex);
        r.err = 2.0 * GSL_DBL_EPSILON * (x + 1.0) * ex;
    }

    return r;

} // gsl_sf_fermi_dirac_m1_e

// ----------------------------------------------------------------------------

// [Goano (3)]
export function gsl_sf_fermi_dirac_0_e(x)
{
    var ex  = 0.0;
    var ser = 0.0;

    var r   = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (x < -5.0)
    {
        ex  = Math.exp(x);
        ser = 1.0 - ex * (0.5 - ex * (1.0 / 3.0 - ex * (1.0 / 4.0 - ex * (1.0 / 5.0 - ex / 6.0))));
        r.val = ex * ser;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < 10.0)
    {
        r.val = Math.log(1.0 + Math.exp(x));
        r.err = Math.abs(x * GSL_DBL_EPSILON);
    }
    else
    {
        ex = Math.exp(-x);
        r.val = x + ex * (1.0 - 0.5 * ex + ex * ex / 3.0 - ex * ex * ex / 4.0);
        r.err = (x + ex) * GSL_DBL_EPSILON;
    }

    return r;

} // gsl_sf_fermi_dirac_0_e

// ----------------------------------------------------------------------------

export function gsl_sf_fermi_dirac_1_e(x)
{
    var t    = 0.0;
    var ex   = 0.0;
    var term = 0.0;
    var sum  = 0.0;
    var rat  = 0.0;

    var c    = { val: 0.0, err: 0.0 }; // Result;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (x < -1.0)
    {
        // series [Goano (6)]
        //
        ex   = Math.exp(x);
        term = ex;
        sum  = term;
        for (let n = 2; n <= 100 - 1; n++)
        {
            rat = (n - 1) / (n);
            term = term * (-ex * rat * rat);
            sum  = sum + term;
            if (Math.abs(term / sum) < GSL_DBL_EPSILON)
            {
                break;
            }
        }
        r.val = sum;
        r.err = 2.0 * Math.abs(sum) * GSL_DBL_EPSILON;
    }
    else if (x < 1.0)
    {
        r = cheb_eval_e(fd_1_a_cs, x);
    }
    else if (x < 4.0)
    {
        t = 2.0 / 3.0 * (x - 1.0) - 1.0;
        r = cheb_eval_e(fd_1_b_cs, t);
    }
    else if (x < 10.0)
    {
        t = 1.0 / 3.0 * (x - 4.0) - 1.0;
        r = cheb_eval_e(fd_1_c_cs, t);
    }
    else if (x < 30.0)
    {
        t = 0.1 * x - 2.0;
        c = cheb_eval_e(fd_1_d_cs, t);
        r.val = c.val * x * x;
        r.err = c.err * x * x + GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < 1.0/GSL_SQRT_DBL_EPSILON)
    {
        t = 60.0 / x - 1.0;
        c = cheb_eval_e(fd_1_e_cs, t);
        r.val = c.val * x * x;
        r.err = c.err * x * x + GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < GSL_SQRT_DBL_MAX)
    {
        r.val = 0.5 * x * x;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

    return r;

} // gsl_sf_fermi_dirac_1_e

// ----------------------------------------------------------------------------

export function gsl_sf_fermi_dirac_2_e(x)
{
    var t    = 0.0;
    var ex   = 0.0;
    var term = 0.0;
    var sum  = 0.0;
    var rat  = 0.0;
    var c    = { val: 0.0, err: 0.0 }; // Result;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (x < -1.0)
    {
        // series [Goano (6)]
        //
        ex   = Math.exp(x);
        term = ex;
        sum  = term;
        for (let n = 2; n <= 100 - 1; n++)
        {
            rat = (n - 1) / (n);
            term = term * (-ex * rat * rat * rat);
            sum  = sum + term;
            if (Math.abs(term / sum) < GSL_DBL_EPSILON)
            {
                break;
            }
        }
        r.val = sum;
        r.err = 2.0 * GSL_DBL_EPSILON * Math.abs(sum);
    }
    else if (x < 1.0)
    {
        r = cheb_eval_e(fd_2_a_cs, x);
    }
    else if (x < 4.0)
    {
        t = 2.0 / 3.0 * (x - 1.0) - 1.0;
        r = cheb_eval_e(fd_2_b_cs, t);
    }
    else if (x < 10.0)
    {
        t = 1.0 / 3.0 * (x - 4.0) - 1.0;
        r = cheb_eval_e(fd_2_c_cs, t);
    }
    else if (x < 30.0)
    {
        t = 0.1 * x - 2.0;
        c = cheb_eval_e(fd_2_d_cs, t);
        r.val = c.val * x * x * x;
        r.err = c.err * x * x * x + 3.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < 1.0 / GSL_ROOT3_DBL_EPSILON)
    {
        t = 60.0 / x - 1.0;
        c = cheb_eval_e(fd_2_e_cs, t);
        r.val = c.val * x * x * x;
        r.err = c.err * x * x * x + 3.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else if (x < GSL_ROOT3_DBL_MAX)
    {
        r.val = 1.0 / 6.0 * x * x * x;
        r.err = 3.0 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        throw "SF.OverflowException";
    }

    return r;

} // gsl_sf_fermi_dirac_2_e

// ----------------------------------------------------------------------------

export function gsl_sf_fermi_dirac_int_e(j, x)
{
    var fasymp = { val: 0.0, err: 0.0 }; // Result;
    var r      = { val: 0.0, err: 0.0 }; // Result;

    if (j < -1)
    {
        return fd_nint(j, x);
    }
    else if (j == -1)
    {
        return gsl_sf_fermi_dirac_m1_e(x);
    }
    else if (j == 0)
    {
        return gsl_sf_fermi_dirac_0_e(x);
    }
    else if (j == 1)
    {
        return gsl_sf_fermi_dirac_1_e(x);
    }
    else if (j == 2)
    {
        return gsl_sf_fermi_dirac_2_e(x);
    }
    else if (x < 0.0)
    {
        return fd_neg((j), x);
    }
    else if (x == 0.0)
    {
        return gsl_sf_eta_int_e(j + 1);
    }
    else if (x < 1.5)
    {
        return fd_series_int(j, x);
    }
    else
    {
        try
        {
            fasymp = fd_asymp((j), x);
            r.val = fasymp.val;
            r.err = fasymp.err;
            r.err = r.err + 2.0 * GSL_DBL_EPSILON * Math.abs(r.val);
        }
        catch (e)
        {
            r = fd_UMseries_int(j, x);
        }
        return r;
    }

} // gsl_sf_fermi_dirac_int_e

// ----------------------------------------------------------------------------

export function gsl_sf_fermi_dirac_mhalf_e(x)
{
    var t    = 0.0;
    var ex   = 0.0;
    var term = 0.0;
    var sum  = 0.0;
    var rat  = 0.0;
    var rtx  = 0.0;
    var c    = { val: 0.0, err: 0.0 }; // Result;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (x < -1.0)
    {
        // series [Goano (6)]
        //
        ex   = Math.exp(x);
        term = ex;
        sum  = term;
        for (let n = 2; n <= 200 - 1; n++)
        {
            rat = (n - 1) / (n);
            term = term * (-ex * Math.sqrt(rat));
            sum  = sum + term;
            if (Math.abs(term / sum) < GSL_DBL_EPSILON)
            {
                break;
            }
        }
        r.val = sum;
        r.err = 2.0 * Math.abs(sum) * GSL_DBL_EPSILON;
    }
    else if (x < 1.0)
    {
        r = cheb_eval_e(fd_mhalf_a_cs, x);
    }
    else if (x < 4.0)
    {
        t = 2.0 / 3.0 * (x - 1.0) - 1.0;
        r = cheb_eval_e(fd_mhalf_b_cs, t);
    }
    else if (x < 10.0)
    {
        t = 1.0 / 3.0 * (x - 4.0) - 1.0;
        r = cheb_eval_e(fd_mhalf_c_cs, t);
    }
    else if (x < 30.0)
    {
        rtx = Math.sqrt(x);
        t = 0.1 * x - 2.0;
        c = cheb_eval_e(fd_mhalf_d_cs, t);
        r.val = c.val * rtx;
        r.err = c.err * rtx + 0.5 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        r = fd_asymp(-0.5, x);
    }

    return r;

} // gsl_sf_fermi_dirac_mhalf_e

// ----------------------------------------------------------------------------

export function gsl_sf_fermi_dirac_half_e(x)
{
    var t    = 0.0;
    var ex   = 0.0;
    var term = 0.0;
    var sum  = 0.0;
    var x32  = 0.0;
    var rat  = 0.0;
    var c    = { val: 0.0, err: 0.0 }; // Result;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (x < -1.0)
    {
        // series [Goano (6)]
        //
        ex   = Math.exp(x);
        term = ex;
        sum  = term;
        for (let n = 2; n <= 100 - 1; n++)
        {
            rat = (n - 1) / (n);
            term = term * (-ex * rat * Math.sqrt(rat));
            sum  = sum + term;
            if (Math.abs(term / sum) < GSL_DBL_EPSILON)
            {
                break;
            }
        }
        r.val = sum;
        r.err = 2.0 * Math.abs(sum) * GSL_DBL_EPSILON;
    }
    else if (x < 1.0)
    {
        r = cheb_eval_e(fd_half_a_cs, x);
    }
    else if (x < 4.0)
    {
        t = 2.0 / 3.0 * (x - 1.0) - 1.0;
        r = cheb_eval_e(fd_half_b_cs, t);
    }
    else if (x < 10.0)
    {
        t = 1.0 / 3.0 * (x - 4.0) - 1.0;
        r = cheb_eval_e(fd_half_c_cs, t);
    }
    else if (x < 30.0)
    {
        x32 = x * Math.sqrt(x);
        t = 0.1 * x - 2.0;
        c = cheb_eval_e(fd_half_d_cs, t);
        r.val = c.val * x32;
        r.err = c.err * x32 + 1.5 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        r = fd_asymp(0.5, x);
    }

    return r;

} // gsl_sf_fermi_dirac_half_e

// ----------------------------------------------------------------------------

export function gsl_sf_fermi_dirac_3half_e(x)
{
    var t    = 0.0;
    var ex   = 0.0;
    var term = 0.0;
    var sum  = 0.0;
    var x52  = 0.0;
    var rat  = 0.0;
    var c    = { val: 0.0, err: 0.0 }; // Result;
    var r    = { val: 0.0, err: 0.0 }; // Result;

    if (x < GSL_LOG_DBL_MIN)
    {
        throw "SF.UnderflowException";
    }
    else if (x < -1.0)
    {
        // series [Goano (6)]
        //
        ex   = Math.exp(x);
        term = ex;
        sum  = term;
        for (let n = 2; n <= 100 - 1; n++)
        {
            rat = (n - 1) / (n);
            term = term * (-ex * rat * rat * Math.sqrt(rat));
            sum  = sum + term;
            if (Math.abs(term / sum) < GSL_DBL_EPSILON)
            {
                break;
            }
        }
        r.val = sum;
        r.err = 2.0 * Math.abs(sum) * GSL_DBL_EPSILON;
    }
    else if (x < 1.0)
    {
        r = cheb_eval_e(fd_3half_a_cs, x);
    }
    else if (x < 4.0)
    {
        t = 2.0 / 3.0 * (x - 1.0) - 1.0;
        r = cheb_eval_e(fd_3half_b_cs, t);
    }
    else if (x < 10.0)
    {
        t = 1.0 / 3.0 * (x - 4.0) - 1.0;
        r = cheb_eval_e(fd_3half_c_cs, t);
    }
    else if (x < 30.0)
    {
        x52 = x * x * Math.sqrt(x);
        t = 0.1 * x - 2.0;
        c = cheb_eval_e(fd_3half_d_cs, t);
        r.val = c.val * x52;
        r.err = c.err * x52 + 2.5 * GSL_DBL_EPSILON * Math.abs(r.val);
    }
    else
    {
        r = fd_asymp(1.5, x);
    }

    return r;

} // gsl_sf_fermi_dirac_3half_e

// ----------------------------------------------------------------------------

// [Goano p. 222]
export function gsl_sf_fermi_dirac_inc_0_e( x, b )
{
    var r    = { val: 0.0, err: 0.0 }; // Result;

    if ( b < 0.0 )
    {
        throw "SF.DomainException"; //DOMAIN_ERROR(result);
    }
    else
    {
        var arg = b - x;
        var f0  = { val: 0.0, err: 0.0 }; // Result;
        f0 = gsl_sf_fermi_dirac_0_e( arg );
        r.val = f0.val - arg;
        r.err = f0.err + GSL_DBL_EPSILON * (Math.abs( x ) + Math.abs( b ));
    }

    return r;
}

//*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*

export function gsl_sf_fermi_dirac_m1( x )
{ // gsl_sf_fermi_dirac_m1
    return EVAL_RESULT_D( gsl_sf_fermi_dirac_m1_e, x, "gsl_sf_fermi_dirac_m1" );
} // gsl_sf_fermi_dirac_m1

export function gsl_sf_fermi_dirac_0( x )
{ // gsl_sf_fermi_dirac_0
    return EVAL_RESULT_D( gsl_sf_fermi_dirac_0_e, x, "gsl_sf_fermi_dirac_0" );
} // gsl_sf_fermi_dirac_0

export function gsl_sf_fermi_dirac_1( x )
{ // gsl_sf_fermi_dirac_1
    return EVAL_RESULT_D( gsl_sf_fermi_dirac_1_e, x, "gsl_sf_fermi_dirac_1" );
} // gsl_sf_fermi_dirac_1

export function gsl_sf_fermi_dirac_2( x )
{ // gsl_sf_fermi_dirac_2
    return EVAL_RESULT_D( gsl_sf_fermi_dirac_2_e, x, "gsl_sf_fermi_dirac_2" );
} // gsl_sf_fermi_dirac_2

export function gsl_sf_fermi_dirac_int( j, x )
{ // gsl_sf_fermi_dirac_int
    return EVAL_RESULT_ID( gsl_sf_fermi_dirac_int_e, { i: j, x: x }, "gsl_sf_fermi_dirac_int" );
} // gsl_sf_fermi_dirac_int

export function gsl_sf_fermi_dirac_mhalf( x )
{ // gsl_sf_fermi_dirac_mhalf
    return EVAL_RESULT_D( gsl_sf_fermi_dirac_mhalf_e, x, "gsl_sf_fermi_dirac_mhalf" );
} // gsl_sf_fermi_dirac_mhalf

export function gsl_sf_fermi_dirac_half( x )
{ // gsl_sf_fermi_dirac_half
    return EVAL_RESULT_D( gsl_sf_fermi_dirac_half_e, x, "gsl_sf_fermi_dirac_half" );
} // gsl_sf_fermi_dirac_half

export function gsl_sf_fermi_dirac_3half( x )
{ // gsl_sf_fermi_dirac_3half
    return EVAL_RESULT_D( gsl_sf_fermi_dirac_3half_e, x, "gsl_sf_fermi_dirac_3half" );
} // gsl_sf_fermi_dirac_3half

export function gsl_sf_fermi_dirac_inc_0( x, b )
{
    return EVAL_RESULT_DD( gsl_sf_fermi_dirac_inc_0_e, { x: x, y: b }, "gsl_sf_fermi_dirac_inc_0_e" );
}

// ----------------------------------------------------------------------------
// EOF SF-FermiDirac.mjs
