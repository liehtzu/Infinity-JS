// SF-Test.mjs
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

import { M_E }                   from "./SF-Math.mjs";
import { M_PI }                  from "./SF-Math.mjs";
import { M_PI_2 }                from "./SF-Math.mjs";
import { M_LN2 }                 from "./SF-Math.mjs";

import { gsl_test  }   from './SF-TestResult.mjs';

import { gsl_sf_expint_E1_e }        from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_E2_e }        from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_En_e }        from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_Ei_e }        from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_E1_scaled_e } from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_E2_scaled_e } from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_En_scaled_e } from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_Ei_scaled_e } from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_expint_3_e }         from "./SF-ExponentialIntegral.mjs";
import { gsl_sf_Shi_e }              from "./SF-HyperbolicIntegrals.mjs";
import { gsl_sf_Chi_e }              from "./SF-HyperbolicIntegrals.mjs";

import { gsl_sf_Si_e }               from "./SF-TrigonometricIntegrals.mjs";
import { gsl_sf_Ci_e }               from "./SF-TrigonometricIntegrals.mjs";
import { gsl_sf_atanint_e }          from "./SF-TrigonometricIntegrals.mjs";

import { gsl_sf_erfc_e }       from "./SF-Erfc.mjs";
import { gsl_sf_log_erfc_e }   from "./SF-Erfc.mjs";
import { gsl_sf_erf_e }        from "./SF-Erfc.mjs";
import { gsl_sf_erf_Z_e }      from "./SF-Erfc.mjs";
import { gsl_sf_erf_Q_e }      from "./SF-Erfc.mjs";
import { gsl_sf_hazard_e }     from "./SF-Erfc.mjs";

import { gsl_sf_dawson_e }     from "./SF-Dawson.mjs";

import { gsl_sf_debye_1_e }    from "./SF-Debye.mjs";
import { gsl_sf_debye_2_e }    from "./SF-Debye.mjs";
import { gsl_sf_debye_3_e }    from "./SF-Debye.mjs";
import { gsl_sf_debye_4_e }    from "./SF-Debye.mjs";
import { gsl_sf_debye_5_e }    from "./SF-Debye.mjs";
import { gsl_sf_debye_6_e }    from "./SF-Debye.mjs";

import { gsl_sf_psi_int_e }    from "./SF-Psi.mjs";
import { gsl_sf_psi_e }        from "./SF-Psi.mjs";
import { gsl_sf_psi_1piy_e }   from "./SF-Psi.mjs";
import { gsl_sf_psi_1_int_e }  from "./SF-Psi.mjs";
import { gsl_sf_psi_1_e }      from "./SF-Psi.mjs";
import { gsl_sf_psi_n_e }      from "./SF-Psi.mjs";

import { gsl_sf_lngamma_e }      from "./SF-Gamma.mjs";
import { gsl_sf_gamma_e }        from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_sgn_e }  from "./SF-Gamma.mjs";
import { gsl_sf_gammastar_e }    from "./SF-Gamma.mjs";
import { gsl_sf_gammainv_e }     from "./SF-Gamma.mjs";
import { gsl_sf_lnfact_e }       from "./SF-Gamma.mjs";
import { gsl_sf_fact_e }         from "./SF-Gamma.mjs";
import { gsl_sf_doublefact_e }   from "./SF-Gamma.mjs";
import { gsl_sf_lndoublefact_e } from "./SF-Gamma.mjs";
import { gsl_sf_lnchoose_e }     from "./SF-Gamma.mjs";
import { gsl_sf_choose_e }       from "./SF-Gamma.mjs";
import { gsl_sf_taylorcoeff_e }  from "./SF-Gamma.mjs";
import { gsl_sf_lngamma_complex_e } from "./SF-Gamma.mjs";

import { gsl_sf_hzeta_e}       from "./SF-Zeta.mjs";
import { gsl_sf_zeta_int_e }   from "./SF-Zeta.mjs";
import { gsl_sf_zetam1_int_e } from "./SF-Zeta.mjs";
import { gsl_sf_zeta_e }       from "./SF-Zeta.mjs";
import { gsl_sf_zetam1_e }     from "./SF-Zeta.mjs";
import { gsl_sf_eta_int_e }    from "./SF-Zeta.mjs";
import { gsl_sf_eta_e }        from "./SF-Zeta.mjs";

import { gsl_sf_lnpoch_e }      from "./SF-Pochhammer.mjs";
import { gsl_sf_lnpoch_sgn_e }  from "./SF-Pochhammer.mjs";
import { gsl_sf_poch_e }        from "./SF-Pochhammer.mjs";
import { gsl_sf_pochrel_e }     from "./SF-Pochhammer.mjs";

import { gsl_sf_gamma_inc_e }   from "./SF-GammaIncomplete.mjs";
import { gsl_sf_gamma_inc_P_e } from "./SF-GammaIncomplete.mjs";
import { gsl_sf_gamma_inc_Q_e } from "./SF-GammaIncomplete.mjs";

import { gsl_sf_lnbeta_e }      from "./SF-Beta.mjs";
import { gsl_sf_beta_e }        from "./SF-Beta.mjs";
import { gsl_sf_beta_inc_e }    from "./SF-BetaIncomplete.mjs";

import { gsl_sf_log_e }           from "./SF-Logarithmic.mjs";
import { gsl_sf_log_abs_e }       from "./SF-Logarithmic.mjs";
import { gsl_sf_log_1plusx_e }    from "./SF-Logarithmic.mjs";
import { gsl_sf_log_1plusx_mx_e } from "./SF-Logarithmic.mjs";
import { gsl_sf_complex_log_e }   from "./SF-Logarithmic.mjs";

import { gsl_sf_gegenpoly_1_e }   from "./SF-Gegenbauer.mjs";
import { gsl_sf_gegenpoly_2_e }   from "./SF-Gegenbauer.mjs";
import { gsl_sf_gegenpoly_3_e }   from "./SF-Gegenbauer.mjs";
import { gsl_sf_gegenpoly_n_e }   from "./SF-Gegenbauer.mjs";
import { gsl_sf_gegenpoly_array } from "./SF-Gegenbauer.mjs";

import { gsl_sf_clausen_e } from "./SF-Clausen.mjs";

import { GSL_MODE_DEFAULT } from "./SF-Mode.mjs";

import { gsl_sf_airy_Ai_e }              from "./SF-Airy.mjs";
import { gsl_sf_airy_Bi_e }              from "./SF-Airy.mjs";
import { gsl_sf_airy_Ai_scaled_e }       from "./SF-Airy.mjs";
import { gsl_sf_airy_Bi_scaled_e }       from "./SF-Airy.mjs";
import { gsl_sf_airy_Ai_deriv_e }        from "./SF-Airy.mjs";
import { gsl_sf_airy_Ai_deriv_scaled_e } from "./SF-Airy.mjs";
import { gsl_sf_airy_Bi_deriv_e }        from "./SF-Airy.mjs";
import { gsl_sf_airy_Bi_deriv_scaled_e } from "./SF-Airy.mjs";
import { gsl_sf_airy_zero_Ai_e }         from "./SF-AiryZeros.mjs";
import { gsl_sf_airy_zero_Bi_e }         from "./SF-AiryZeros.mjs";
import { gsl_sf_airy_zero_Ai_deriv_e }   from "./SF-AiryZeros.mjs";
import { gsl_sf_airy_zero_Bi_deriv_e }   from "./SF-AiryZeros.mjs";

import { gsl_sf_bessel_J0_e } from "./SF-BesselJ0.mjs";
import { gsl_sf_bessel_J1_e } from "./SF-BesselJ1.mjs";
import { gsl_sf_bessel_Jn_e } from "./SF-BesselJn.mjs";
import { gsl_sf_bessel_Y0_e } from "./SF-BesselY0.mjs";
import { gsl_sf_bessel_Y1_e } from "./SF-BesselY1.mjs";
import { gsl_sf_bessel_Yn_e } from "./SF-BesselYn.mjs";
import { gsl_sf_bessel_I0_scaled_e } from "./SF-BesselI0.mjs";
import { gsl_sf_bessel_I1_scaled_e } from "./SF-BesselI1.mjs";
import { gsl_sf_bessel_I0_e }  from "./SF-BesselI0.mjs";
import { gsl_sf_bessel_I1_e }  from "./SF-BesselI1.mjs";
import { gsl_sf_bessel_In_e }  from "./SF-BesselIn.mjs";
import { gsl_sf_bessel_In_scaled_e } from "./SF-BesselIn.mjs";
import { gsl_sf_bessel_K0_scaled_e } from "./SF-BesselK0.mjs";
import { gsl_sf_bessel_K1_scaled_e } from "./SF-BesselK1.mjs";
import { gsl_sf_bessel_Kn_scaled_e } from "./SF-BesselKn.mjs";
import { gsl_sf_bessel_Kn_e }  from "./SF-BesselKn.mjs";
import { gsl_sf_bessel_K0_e }  from "./SF-BesselK0.mjs";
import { gsl_sf_bessel_K1_e }  from "./SF-BesselK1.mjs";
import { gsl_sf_bessel_j0s_e } from "./SF-BesselJ.mjs";
import { gsl_sf_bessel_j1s_e } from "./SF-BesselJ.mjs";
import { gsl_sf_bessel_j2s_e } from "./SF-BesselJ.mjs";
import { gsl_sf_bessel_jl_e }  from "./SF-BesselJ.mjs";
import { gsl_sf_bessel_y0s_e } from "./SF-BesselY.mjs";
import { gsl_sf_bessel_y1s_e } from "./SF-BesselY.mjs";
import { gsl_sf_bessel_y2s_e } from "./SF-BesselY.mjs";
import { gsl_sf_bessel_yl_e }  from "./SF-BesselY.mjs";
import { gsl_sf_bessel_i0s_scaled_e } from "./SF-BesselI.mjs";
import { gsl_sf_bessel_i1s_scaled_e } from "./SF-BesselI.mjs";
import { gsl_sf_bessel_i2s_scaled_e } from "./SF-BesselI.mjs";
import { gsl_sf_bessel_il_scaled_e }  from "./SF-BesselI.mjs";
import { gsl_sf_bessel_k0s_scaled_e } from "./SF-BesselK.mjs";
import { gsl_sf_bessel_k1s_scaled_e } from "./SF-BesselK.mjs";
import { gsl_sf_bessel_k2s_scaled_e } from "./SF-BesselK.mjs";
import { gsl_sf_bessel_kl_scaled_e }  from "./SF-BesselK.mjs";
import { gsl_sf_bessel_Jnu_e } from "./SF-BesselJnu.mjs";
import { gsl_sf_bessel_Ynu_e } from "./SF-BesselYnu.mjs";
import { gsl_sf_bessel_Inu_scaled_e } from "./SF-BesselInu.mjs";
import { gsl_sf_bessel_Inu_e } from "./SF-BesselInu.mjs";
import { gsl_sf_bessel_Knu_scaled_e } from "./SF-BesselKnu.mjs";
import { gsl_sf_bessel_Knu_e } from "./SF-BesselKnu.mjs";
import { gsl_sf_bessel_lnKnu_e }      from "./SF-BesselKnu.mjs";
import { gsl_sf_bessel_Jn_array }     from "./SF-BesselJn.mjs";
import { gsl_sf_bessel_Yn_array }     from "./SF-BesselYn.mjs";
import { gsl_sf_bessel_In_scaled_array } from "./SF-BesselIn.mjs";
import { gsl_sf_bessel_In_array }     from "./SF-BesselIn.mjs";
import { gsl_sf_bessel_Kn_array }     from "./SF-BesselKn.mjs";
import { gsl_sf_bessel_Kn_scaled_array } from "./SF-BesselKn.mjs";
import { gsl_sf_bessel_jl_array }     from "./SF-BesselJ.mjs";
import { gsl_sf_bessel_jl_steed_array } from "./SF-BesselJ.mjs";
import { gsl_sf_bessel_yl_array }     from "./SF-BesselY.mjs";
import { gsl_sf_bessel_il_scaled_array } from "./SF-BesselI.mjs";
import { gsl_sf_bessel_kl_scaled_array } from "./SF-BesselK.mjs";
import { gsl_sf_bessel_zero_J0_e }    from "./SF-BesselZero.mjs";
import { gsl_sf_bessel_zero_J1_e }    from "./SF-BesselZero.mjs";
import { gsl_sf_bessel_zero_Jnu_e }   from "./SF-BesselZero.mjs";
import { gsl_sf_bessel_sequence_Jnu_e } from "./SF-BesselSequence.mjs";

import { gsl_sf_ellint_Kcomp_e } from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_Ecomp_e } from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_Pcomp_e } from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_Dcomp_e } from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_F_e }     from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_E_e }     from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_P_e }     from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_D_e }     from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_RF_e }    from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_RD_e }    from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_RC_e }    from "./SF-EllipticIntegrals.mjs";
import { gsl_sf_ellint_RJ_e }    from "./SF-EllipticIntegrals.mjs";

import { gsl_sf_multiply_e } from "./SF-Elementary.mjs";

import { gsl_sf_exp_e }              from "./SF-Exponential.mjs";
import { gsl_sf_exp_e10_e }          from "./SF-Exponential.mjs";
import { gsl_sf_exp_err_e }          from "./SF-Exponential.mjs";
import { gsl_sf_exp_err_e10_e }      from "./SF-Exponential.mjs";
import { gsl_sf_exp_mult_e }         from "./SF-Exponential.mjs";
import { gsl_sf_exp_mult_err_e }     from "./SF-Exponential.mjs";
import { gsl_sf_exp_mult_e10_e }     from "./SF-Exponential.mjs";
import { gsl_sf_exp_mult_err_e10_e } from "./SF-Exponential.mjs";
import { gsl_sf_expm1_e }            from "./SF-Exponential.mjs";
import { gsl_sf_exprel_e }           from "./SF-Exponential.mjs";
import { gsl_sf_exprel_2_e }         from "./SF-Exponential.mjs";
import { gsl_sf_exprel_n_e }         from "./SF-Exponential.mjs";
import { gsl_sf_exprel_n_CF_e }      from "./SF-Exponential.mjs";

import { gsl_sf_sin_e }                     from "./SF-Trigonometric.mjs";
import { gsl_sf_cos_e }                     from "./SF-Trigonometric.mjs";
import { gsl_sf_sinc_e }                    from "./SF-Trigonometric.mjs";
import { gsl_sf_lnsinh_e }                  from "./SF-Trigonometric.mjs";
import { gsl_sf_lncosh_e }                  from "./SF-Trigonometric.mjs";
import { gsl_sf_polar_to_rect }             from "./SF-Trigonometric.mjs";
import { gsl_sf_angle_restrict_pos_err_e }  from "./SF-Trigonometric.mjs";
import { gsl_sf_angle_restrict_symm_err_e } from "./SF-Trigonometric.mjs";
import { gsl_sf_angle_restrict_symm_e }     from "./SF-Trigonometric.mjs";
import { gsl_sf_angle_restrict_pos_e }      from "./SF-Trigonometric.mjs";
import { gsl_sf_complex_sin_e }             from "./SF-Trigonometric.mjs";
import { gsl_sf_complex_cos_e }             from "./SF-Trigonometric.mjs";
import { gsl_sf_complex_logsin_e }          from "./SF-Trigonometric.mjs";

import { gsl_sf_transport_2_e } from "./SF-Transport.mjs";
import { gsl_sf_transport_3_e } from "./SF-Transport.mjs";
import { gsl_sf_transport_4_e } from "./SF-Transport.mjs";
import { gsl_sf_transport_5_e } from "./SF-Transport.mjs";

import { gsl_sf_synchrotron_1_e } from "./SF-Synchrotron.mjs";
import { gsl_sf_synchrotron_2_e } from "./SF-Synchrotron.mjs";

import { gsl_sf_pow_int_e } from "./SF-Power.mjs";

import { gsl_sf_lambert_W0_e }  from "./SF-Lambert.mjs";
import { gsl_sf_lambert_Wm1_e } from "./SF-Lambert.mjs";

import { gsl_sf_laguerre_1_e } from "./SF-Laguerre.mjs";
import { gsl_sf_laguerre_2_e } from "./SF-Laguerre.mjs";
import { gsl_sf_laguerre_3_e } from "./SF-Laguerre.mjs";
import { gsl_sf_laguerre_n_e } from "./SF-Laguerre.mjs";

import { gsl_sf_legendre_P1_e } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_P2_e } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_P3_e } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_Pl_e } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_Pl_array } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_Pl_deriv_array } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_Plm_e } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_Plm_array } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_Plm_deriv_array } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_sphPlm_e } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_sphPlm_array } from "./SF-LegendrePolynomials.mjs";
import { gsl_sf_legendre_sphPlm_deriv_array } from "./SF-LegendrePolynomials.mjs";

import { gsl_sf_conicalP_half_e }    from "./SF-LegendreConical.mjs";
import { gsl_sf_conicalP_mhalf_e }   from "./SF-LegendreConical.mjs";
import { gsl_sf_conicalP_0_e }       from "./SF-LegendreConical.mjs";
import { gsl_sf_conicalP_1_e }       from "./SF-LegendreConical.mjs";
import { gsl_sf_conicalP_sph_reg_e } from "./SF-LegendreConical.mjs";
import { gsl_sf_conicalP_cyl_reg_e } from "./SF-LegendreConical.mjs";
import { gsl_sf_legendre_H3d_0_e }   from "./SF-LegendreH3D.mjs";
import { gsl_sf_legendre_H3d_1_e }   from "./SF-LegendreH3D.mjs";
import { gsl_sf_legendre_H3d_e }     from "./SF-LegendreH3D.mjs";
import { gsl_sf_legendre_H3d_array } from "./SF-LegendreH3D.mjs";
import { gsl_sf_legendre_Q0_e }      from "./SF-LegendreQn.mjs";
import { gsl_sf_legendre_Q1_e }      from "./SF-LegendreQn.mjs";
import { gsl_sf_legendre_Ql_e }      from "./SF-LegendreQn.mjs";

import { gsl_sf_elljac_e } from "./SF-EllipticJacobi.mjs";

import { gsl_sf_dilog_e }             from "./SF-Dilogarithm.mjs";
import { gsl_sf_complex_dilog_e }     from "./SF-Dilogarithm.mjs";
import { gsl_sf_complex_dilog_xy_e }  from "./SF-Dilogarithm.mjs";
import { gsl_sf_complex_spence_xy_e } from "./SF-Dilogarithm.mjs";

import { gsl_sf_fermi_dirac_m1_e }    from "./SF-FermiDirac.mjs";
import { gsl_sf_fermi_dirac_0_e }     from "./SF-FermiDirac.mjs";
import { gsl_sf_fermi_dirac_1_e }     from "./SF-FermiDirac.mjs";
import { gsl_sf_fermi_dirac_2_e }     from "./SF-FermiDirac.mjs";
import { gsl_sf_fermi_dirac_mhalf_e } from "./SF-FermiDirac.mjs";
import { gsl_sf_fermi_dirac_half_e }  from "./SF-FermiDirac.mjs";
import { gsl_sf_fermi_dirac_3half_e } from "./SF-FermiDirac.mjs";
import { gsl_sf_fermi_dirac_int_e }   from "./SF-FermiDirac.mjs";

import { gsl_sf_hyperg_0F1_e }      from "./SF-Hypergeometric0F1.mjs";
import { gsl_sf_hyperg_1F1_int_e }  from "./SF-Hypergeometric1F1.mjs";
import { gsl_sf_hyperg_1F1_e }      from "./SF-Hypergeometric1F1.mjs";
import { gsl_sf_hyperg_U_int_e }    from "./SF-HypergeometricU.mjs";
import { gsl_sf_hyperg_U_e }        from "./SF-HypergeometricU.mjs";
import { gsl_sf_hyperg_2F0_e }      from "./SF-Hypergeometric2F0.mjs";
import { gsl_sf_hyperg_2F1_e }      from "./SF-Hypergeometric2F1.mjs";
import { gsl_sf_hyperg_2F1_conj_e } from "./SF-Hypergeometric2F1.mjs";
import { gsl_sf_hyperg_2F1_renorm_e }      from "./SF-Hypergeometric2F1.mjs";
import { gsl_sf_hyperg_2F1_conj_renorm_e } from "./SF-Hypergeometric2F1.mjs";

import { GSL_SQRT_DBL_EPSILON } from "./SF-Machine.mjs";
import { GSL_DBL_MAX }          from "./SF-Machine.mjs";
import { GSL_DBL_MIN }          from "./SF-Machine.mjs";
import { GSL_LOG_DBL_MAX }      from "./SF-Machine.mjs";
import { GSL_DBL_EPSILON }      from "./SF-Machine.mjs";

import { gsl_sf_coupling_3j_e } from "./SF-Coupling.mjs";
import { gsl_sf_coupling_6j_e } from "./SF-Coupling.mjs";
import { gsl_sf_coupling_9j_e } from "./SF-Coupling.mjs";

import { gsl_sf_hydrogenicR_1_e }   from "./SF-CoulombBound.mjs";
import { gsl_sf_hydrogenicR_e }     from "./SF-CoulombBound.mjs";
import { gsl_sf_coulomb_wave_FG_e } from "./SF-Coulomb.mjs";

const TEST_FACTOR = 100.0;
const TEST_SIGMA  = 1.5;

export const TEST_TOL0 = 2.0 * Number.EPSILON;
export const TEST_TOL1 = 16.0 * Number.EPSILON;
export const TEST_TOL2 = 256.0 * Number.EPSILON;
export const TEST_TOL3 = 2048.0 * Number.EPSILON;
export const TEST_TOL4 = 16384.0 * Number.EPSILON;
export const TEST_TOL5 = 131072.0 * Number.EPSILON;
export const TEST_TOL6 = 1048576.0 * Number.EPSILON;
export const TEST_SQRT_TOL0 = 2.0 * GSL_SQRT_DBL_EPSILON;
export const TEST_SNGL = 1.0e-06;

const TEST_SF_INCONS = 1;
const TEST_SF_ERRNEG = 2;
const TEST_SF_TOLBAD = 4;
const TEST_SF_RETBAD = 8;
const TEST_SF_ERRBAD = 16;
const TEST_SF_ERRBIG = 32;
const TEST_SF_EXPBAD = 64;

// ----------------------------------------------------------------------------

function test_sf_frac_diff(x1, x2)
{ // test_sf_frac_diff
    //console.log("test_sf_frac_diff called");
    if (x1 == 0.0 && x2 == 0.0)
    {
        return 0.0;
    }
    else if (x1 == 0.0)
    {
        return Math.abs(x2);
    }
    else if (x1 <= Number.MAX_VALUE && x2 <= Number.MAX_VALUE && (x1 + x2 != 0.0))
    {
        return Math.abs((x1 - x2) / (x1 + x2));
    }
    else
    {
        return 1.0;
    }
} // test_sf_frac_diff

// ----------------------------------------------------------------------------
// Check a result against a given value and tolerance.
function test_sf_check_result(s, r, val, tol)
{
    var f = 0.0;
    var d = 0.0;

    s.Integer = 0;
    //console.log("test_sf_check_result called");

//  if (gsl_isnan(r.val) || gsl_isnan(val))
//    {
//      s = (gsl_isnan(r.val) != gsl_isnan(val)) ? TEST_SF_INCONS : s;
//    }
//  else if (gsl_isinf(r.val) || gsl_isinf(val))
//    {
//      s = (gsl_isinf(r.val) != gsl_isinf(val)) ? TEST_SF_INCONS : s;
//    }
//  else
//    {
    f = test_sf_frac_diff(val, r.val);
    d = Math.abs(val - r.val);

    if (d > 2.0 * TEST_SIGMA * r.err)
    {
        s.Integer = s.Integer | TEST_SF_INCONS;
    }
    if (r.err < 0.0)
    {
        s.Integer = s.Integer | TEST_SF_ERRNEG;
    }
        //IF (gsl_isinf(r.err)) THEN
        //    s = s OR TEST_SF_ERRBAD;
        //END IF;
//#if TEST_EXCESSIVE_ERROR
//      if(d > 0 && r.err > 1e4 * fabs(val)*tol)         s |= TEST_SF_ERRBIG;
//#endif
    if (f > TEST_FACTOR * tol)
    {
        s.Integer = s.Integer | TEST_SF_TOLBAD;
    }
//    }

    if ( s.Integer != 0 )
    {
        // process.stdout.write("  expected: %20.16f\n", val);
        // process.stdout.write("  obtained: %20.16f +/- %.16f (rel=%g)\n", r.val, r.err, r.err/(Math.abs(r.val) + r.err));
        // process.stdout.write("  fracdiff: %20.16f\n", f);
        // process.stdout.write(" tolerance: %20.16f\n", tol);
        process.stdout.write("  expected: " + val + "\n");
        process.stdout.write("  obtained: " + r.val + " +/- " + r.err + " (rel=" + r.err/(Math.abs(r.val) + r.err) + ")\n");
        process.stdout.write("  fracdiff: " + f + "\n");
        process.stdout.write(" tolerance: " + tol + "n");
    }

    if ( (s.Integer & TEST_SF_INCONS) != 0 )
    {
        process.stdout.write("  value/expected not consistent within reported error\n");
    }
    if ( (s.Integer & TEST_SF_ERRNEG) != 0 )
    {
        process.stdout.write("  reported error negative\n");
    }
    if ( (s.Integer & TEST_SF_ERRBAD) != 0 )
    {
        process.stdout.write("  reported error is bad\n");
    }
    if ( (s.Integer & TEST_SF_ERRBIG) != 0 )
    {
        process.stdout.write("  reported error is much too big\n");
    }
    if ( (s.Integer & TEST_SF_TOLBAD) != 0 )
    {
        process.stdout.write("  value not within tolerance of expected value\n");
    }

} // test_sf_check_result

// ----------------------------------------------------------------------------

// Check a result against a given value and tolerance.
function test_sf_check_e10(s, e10, e10_in)
{

    s.Integer = 0;

    if (e10 != e10_in)
    {
        s = TEST_SF_EXPBAD;
    }

    if ( s.Integer != 0 )
    {
        process.stdout.write("  expected exponent: 10^" + e10_in + "\n");
        process.stdout.write("  obtained exponent: 10^" + e10 + "\n");
    }

    if ( (s.Integer & TEST_SF_EXPBAD) != 0 )
    {
        process.stdout.write("  exponent is incorrect\n");
    }

} // test_sf_check_e10

// ----------------------------------------------------------------------------

function test_sf_check_val(s, rval, val, tol)
{
    var f = 0.0;

    f = test_sf_frac_diff(val, rval);

    if (f > TEST_FACTOR * tol)
    {
        s = s | TEST_SF_TOLBAD;
    }

    if (s != 0)
    {
        process.stdout.write("  expected: " + val + "\n");
        process.stdout.write("  obtained: " + rval + "\n");
        process.stdout.write("  fracdiff: " + f + "\n");
    }

    if ((s & TEST_SF_TOLBAD) != 0)
    {
        process.stdout.write("  value not within tolerance of expected value\n");
    }

} // test_sf_check_val

// ----------------------------------------------------------------------------

///* Relax the condition on the agreement and on the usefulness
// * of the error estimate.
// */
//int
//test_sf_check_result_relax(char * message_buff, Result r, double val, double tol)
//{
//  int    s = 0;
//  double f = test_sf_frac_diff(val, r.val);
//
//  if(f > GSL_MAX_DBL(TEST_SNGL, TEST_FACTOR * tol))   s |= TEST_SF_INCONS;
//  if(r.err < 0.0)     s |= TEST_SF_ERRNEG;
//  if(gsl_isinf(r.err))              s |= TEST_SF_ERRBAD;
//  if(f > TEST_FACTOR * tol)         s |= TEST_SF_TOLBAD;
//
//  if(s != 0) {
//    char buff[2048];
//    sprintf(buff, "  expected: %20.16e\n", val);
//    strcat(message_buff, buff);
//    sprintf(buff, "  obtained: %20.16e +/- %.16e  (rel=%g)\n", r.val, r.err, r.err/(fabs(r.val) + r.err));
//    strcat(message_buff, buff);
//    sprintf(buff, "  fracdiff: %20.16e\n", f);
//    strcat(message_buff, buff);
//  }
//
//  if(s & TEST_SF_INCONS) {
//    strcat(message_buff, "  value/expected not consistent MAX(tol,SINGLE_PREC)\n");
//  }
//  if(s & TEST_SF_ERRNEG) {
//    strcat(message_buff, "  reported error negative\n");
//  }
//  if(s & TEST_SF_ERRBAD) {
//    strcat(message_buff, "  reported error is bad\n");
//  }
//  if(s & TEST_SF_TOLBAD) {
//    strcat(message_buff, "  value not within tolerance of expected value\n");
//  }
//
//  return s;
//}
//

// ----------------------------------------------------------------------------
// Check a return value.
// PROCEDURE test_sf_check_return(s: IN OUT INTEGER; val_return: INTEGER; expected_return: INTEGER) IS
// BEGIN -- test_sf_check_return

//     IF (val_return /= expected_return) THEN
//         PUT("  unexpected return code: ");
//         PUT(val_return);
//         NEW_LINE;
//         s = TEST_SF_RETBAD;
//     ELSE
//         s = 0;
//     END IF;

// END test_sf_check_return;

// ----------------------------------------------------------------------------

function test_sf(r, val_in, tol, status, desc)
{
    var local_s = { Integer: 0 };
    var s = { Integer: 0 };
    //console.log("test_sf called");

    test_sf_check_result( s, r, val_in, tol );
    local_s.Integer = local_s.Integer | s.Integer;
    //test_sf_check_return(s, status, expect_return);
    //local_s = local_s OR s;

    gsl_test( local_s, desc );
    if (local_s.Integer != 0)
    {
        // /* printf("  %s %d\n", __FILE__, __LINE__); */
        process.stdout.write( "  " + r.val + "  " + r.err + "\n" );
    }

    return local_s.Integer;

} // test_sf

// ----------------------------------------------------------------------------

function test_sf_e10(re, val_in, e10_in, tol, status, desc)
{
    var local_s = { Integer: 0 };
    var s = { Integer: 0 };
    var r = { val: 0.0, err: 0.0 }; // Result;

    local_s.Integer = 0;
    r.val = re.val;
    r.err = re.err;

    test_sf_check_result( s, r, val_in, tol );
    local_s.Integer = local_s.Integer | s.Integer;
    test_sf_check_e10(s, re.e10, e10_in);
    local_s.Integer = local_s.Integer | s.Integer;

    gsl_test(local_s, desc);
    if (local_s.Integer != 0)
    {
        // printf("  %s %d\n", __FILE__, __LINE__);
        //process.stdout.write("  %24.18f  %24.18f  10^\n", re.val, re.err, re.e10);
        process.stdout.write("  " + re.val + "  " + re.err + "  10^" + re.e10 + "\n");
    }
    return local_s.Integer;

} // test_sf_e10

// ----------------------------------------------------------------------------

function test_sf_val(val, val_in, tol, desc)
{
    var local_s = { Integer: 0 };
    var s = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    test_sf_check_val( local_s, val, val_in, tol );

    gsl_test(local_s, desc);
    if ( local_s.Integer != 0 )
    {
        // printf("  %s %d\n", __FILE__, __LINE__); */
        //printf("  %22.18e\n", val);
        process.stdout.write( "  " + val );
    }

    return local_s.Integer;

} // test_sf_val

// ----------------------------------------------------------------------------

//int
//test_sf_rlx (Result r, double val_in, double tol, int status,
//             int expect_return, const char * desc)
//{
//  char message_buff[4096];
//  int local_s = 0;
//
//  message_buff[0] = '\0';
//
//  local_s |= test_sf_check_result_relax(message_buff, r, val_in, tol);
//  local_s |= test_sf_check_return(message_buff, status, expect_return);
//
//  gsl_test(local_s, desc);
//  if(local_s != 0) {
//    /* printf("  %s %d\n", __FILE__, __LINE__); */
//    printf("%s", message_buff);
//    printf("  %22.18e  %22.18e\n", r.val, r.err);
//  }
//  return local_s;
//}

// ----------------------------------------------------------------------------

function test_sf_2( r1, val1, tol1, r2, val2, tol2, status, desc )
{

    var local_s = { Integer: 0 };
    var s = { Integer: 0 };

    test_sf_check_result( s, r1, val1, tol1 );
    local_s.Integer = local_s.Integer | s.Integer;
    test_sf_check_result( s, r2, val2, tol2 );
    local_s.Integer = local_s.Integer | s.Integer;
    //test_sf_check_return(s, status, expect_return);
    //local_s = local_s OR s;

    gsl_test( local_s, desc );
    if ( local_s.Integer != 0 )
    {
        // printf("  %s %d\n", __FILE__, __LINE__);
        //printf("  %22.18e  %22.18e\n", r1.val, r1.err);
        process.stdout.write("  " + r.val + "  " + r.err + "\n");
        //printf("  %22.18e  %22.18e\n", r2.val, r2.err);
        process.stdout.write("  " + r2.val + "  " + r2.err + "\n");
    }

    return local_s.Integer;

} // test_sf_2

// ----------------------------------------------------------------------------

function test_sf_sgn(r, sgn, val_in, tol, expect_sgn, status, desc)
{
    var local_s = { Integer: 0 };
    var s       = { Integer: 0 };
    var local_r = { val: 0.0, err: 0.0 }; // Result;

    local_r.val = sgn;
    local_r.err = 0.0;
    test_sf_check_result(s, r, val_in, tol);
    local_s.Integer = local_s.Integer | s.Integer;
    test_sf_check_result(s, local_r, expect_sgn, 0.0);
    local_s.Integer = local_s.Integer | s.Integer;

    gsl_test(local_s, desc);
    if (local_s.Integer != 0)
    {
        // printf("  %s %d\n", __FILE__, __LINE__);
        // printf("  %22.18e  %22.18e\n", r.val, r.err);
        process.stdout.write( "  " + r.val + "  " + r.err + "\n" );
    }

    return local_s.Integer;

} // test_sf_sgn

// ----------------------------------------------------------------------------

// This variable is used only in TEST_SF_* routines to control the
// output to console the result of the function call
var sf_internalLog = false;

export function SetInternalLog( s )  { sf_internalLog = s; }

// Single LONG_FLOAT (D) parameter --------------------------------------------

function TEST_SF_D( stat, func, x, val_in, tol, funcName )
{
    var status = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( x );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_D

// Two LONG_FLOAT (D) parameters ----------------------------------------------

function TEST_SF_DD( stat, func, args, val_in, tol, funcName )
{
    var status = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.y );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_DD

// Three LONG_FLOAT (D) parameters --------------------------------------------

function TEST_SF_3D( stat, func, args, val_in, tol, funcName )
{
    var status = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.y, args.z );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf(r, val_in, tol, status, funcName);
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_3D

// Four LONG_FLOAT (D) parameters ---------------------------------------------

function TEST_SF_4D( stat, func, args, val_in, tol, funcName )
{
    var status = 0;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.dx, args.y, args.dy );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf(r, val_in, tol, status, funcName);
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_4D

// Single INTEGER (I) and single LONG_FLOAT (D) parameters --------------------

function TEST_SF_ID( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.i, args.x );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf(r, val_in, tol, status, funcName);
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_ID

// Single single LONG_FLOAT (D) and INTEGER (I) parameters --------------------

function TEST_SF_DI( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.i );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_DI

// Single INTEGER (I) and two LONG_FLOAT (D) parameters -----------------------

function TEST_SF_IDD( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.n, args.x, args.y );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_IDD

// Single Double (D) and single Mode (M) parameters -----------------------

export function TEST_SF_DM( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        //func(status, args.x, args.m, r);
        r = func( args.x, args.m );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_DM

// Two LONG_FLOAT (D) and single Mode (M) parameters --------------------------

export function TEST_SF_DDM( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.y, args.m );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_DDM

// Three LONG_FLOAT (D) and single Mode (M) parameters ------------------------

export function TEST_SF_3DM( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.y, args.z, args.m );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_3DM

// Four LONG_FLOAT (D) and single Mode (M) parameters -------------------------

export function TEST_SF_4DM( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.y, args.z, args.p, args.m );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_4DM

// Single INTEGER (I) parameter -----------------------------------------------

function TEST_SF_I( stat, func, n, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( n );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_I

// Two INTEGER (I) parameters -------------------------------------------------

function TEST_SF_II( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.n, args.m );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_II

// Two INTEGER (I) and two LONG_FLOAT (D) parameters --------------------------

// PROCEDURE TEST_SF(stat: IN OUT INTEGER; func: FuncIIDD; args: ArgsIIDD; val_in: LONG_FLOAT;
//     tol: LONG_FLOAT; funcName: STRING) IS
//     status : INTEGER = 0; --GSL_SUCCESS;
//     r : Result;
// BEGIN -- TEST_SF
//     r = func(args.n, args.m, args.x, args.y);
//     stat = stat + test_sf(r, val_in, tol, status, funcName);
// EXCEPTION
//     WHEN SF.DomainException =>
//         PUT_LINE("DOMAIN EXCEPTION");
//     WHEN SF.OverflowException =>
//         PUT_LINE("OVERFLOW EXCEPTION");
//     WHEN SF.MaxIterationsException =>
//         PUT_LINE("MAXIMUM ITERATIONS REACHED EXCEPTION");
//     WHEN OTHERS =>
//         PUT_LINE("GENERAL EXCEPTION");
// END TEST_SF;

// Two INTEGER (I) and single LONG_FLOAT (D) parameters -----------------------

function TEST_SF_IID( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.n, args.m, args.x );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_IID

// Two INTEGER (I) and two LONG_FLOAT (D) parameters -----------------------

function TEST_SF_IIDD( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.i1, args.i2, args.x, args.y );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_IID

// Six INTEGER (I) parameters -------------------------------------------------

function TEST_SF_6I( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.i1, args.i2, args.i3, args.i4, args.i5, args.i6 );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_6I

// Nine INTEGER (I) parameters ------------------------------------------------

function TEST_SF_9I( stat, func, args, val_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.i1, args.i2, args.i3, args.i4, args.i5, args.i6, args.i7, args.i8, args.i9 );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf( r, val_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_9I

// ----------------------------------------------------------------------------

// Single LONG_FLOAT (D) parameter --------------------------------------------

function TEST_SF_E10_D( stat, func, args, val_in, e10_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf_e10( r, val_in, e10_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_E10_D

// Two LONG_FLOAT (D) parameters ----------------------------------------------

function TEST_SF_E10_DD( stat, func, args, val_in, e10_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.y );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf_e10( r, val_in, e10_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_E10_DD

// Four LONG_FLOAT (D) parameters ---------------------------------------------

function TEST_SF_E10_4D( stat, func, args, val_in, e10_in, tol, funcName )
{
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.dx, args.y, args.dy );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf_e10( r, val_in, e10_in, tol, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_E10_4D

// ****************************************************************************

// Single LONG_FLOAT (D) parameter --------------------------------------------

function TEST_SF_SGN_D( stat, func, args, val_in, tol, expect_sgn, funcName )
{
    var sgn    = 0.0;
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;

    try
    {
        r = func( args );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf_sgn(r, r.sign, val_in, tol, expect_sgn, status, funcName);
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_SGN_D

// Two LONG_FLOAT (D) parameters ----------------------------------------------

function TEST_SF_SGN_DD(stat, func, args, val_in, tol, expect_sgn, funcName)
{
    var sgn    = 0.0;
    var status = 0; //GSL_SUCCESS;
    var r = { val: 0.0, err: 0.0, sign: 0.0 }; // Result;

    try
    {
        r = func( args.x, args.y );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf_sgn(r, r.sign, val_in, tol, expect_sgn, status, funcName);
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_SGN_DD

// ****************************************************************************

function TEST_SF_2(stat, func, args, val1, tol1, val2, tol2, funcName)
{
    var status = 0; //GSL_SUCCESS;
    var r1 = { val: 0.0, err: 0.0 }; // Result;
    var r2 = { val: 0.0, err: 0.0 }; // Result;

    try
    {
        func( args.x, args.y, r1, r2 );
        if ( sf_internalLog )
        {
            console.log( r1 );
            console.log( r2 );
        }
        stat.Integer = stat.Integer + test_sf_2( r1, val1, tol1, r2, val2, tol2, status, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_2

//#define TEST_SF_2(stat, func, args, val1, tol1, val2, tol2, expect_return)
// { int status = func args; stat += test_sf_2(r1, val1, tol1, r2, val2, tol2, status, expect_return, #func #args); }

// ****************************************************************************

// PROCEDURE TEST_SF_THETA(
//               stat     : IN OUT INTEGER;
//               func     : IN     FuncT;
//               args     : IN     LONG_FLOAT;
//               val_in   : IN     LONG_FLOAT;
//               tol      : IN     LONG_FLOAT;
//               funcName : IN     STRING) IS
//     status : INTEGER = 0; --GSL_SUCCESS;
//     theta : LONG_FLOAT = 0.0;
// BEGIN -- TEST_SF_THETA
//     theta = func(args);
//     stat = stat + test_sf_val(theta, val_in, tol, funcName);
// EXCEPTION
//     WHEN SF.DomainException =>
//         PUT_LINE("DOMAIN EXCEPTION");
//     WHEN SF.OverflowException =>
//         PUT_LINE("OVERFLOW EXCEPTION");
//     WHEN SF.MaxIterationsException =>
//         PUT_LINE("MAXIMUM ITERATIONS REACHED EXCEPTION");
//     WHEN OTHERS =>
//         PUT_LINE("GENERAL EXCEPTION");
// END TEST_SF_THETA;

// PROCEDURE TEST_SF_THETA(
//               stat     : IN OUT INTEGER;
//               func     : IN     ProcT;
//               args     : IN     LONG_FLOAT;
//               val_in   : IN     LONG_FLOAT;
//               tol      : IN     LONG_FLOAT;
//               funcName : IN     STRING) IS
//     status : INTEGER = 0; --GSL_SUCCESS;
//     theta : LONG_FLOAT = 0.0;
// BEGIN -- TEST_SF_THETA
//     func(theta);
//     stat = stat + test_sf_val(theta, val_in, tol, funcName);
// EXCEPTION
//     WHEN SF.DomainException =>
//         PUT_LINE("DOMAIN EXCEPTION");
//     WHEN SF.OverflowException =>
//         PUT_LINE("OVERFLOW EXCEPTION");
//     WHEN SF.MaxIterationsException =>
//         PUT_LINE("MAXIMUM ITERATIONS REACHED EXCEPTION");
//     WHEN OTHERS =>
//         PUT_LINE("GENERAL EXCEPTION");
// END TEST_SF_THETA;

// PROCEDURE TEST_SF_THETA_D(
//               stat     : IN OUT INTEGER;
//               func     : IN     FuncD;
//               args     : IN     LONG_FLOAT;
//               val_in   : IN     LONG_FLOAT;
//               tol      : IN     LONG_FLOAT;
//               funcName : IN     STRING) IS
//     status : INTEGER = 0; --GSL_SUCCESS;
//     r      : Result;
// BEGIN -- TEST_SF_THETA
//     r = func(args);
//     stat = stat + test_sf_val(r.val, val_in, tol, funcName);
// EXCEPTION
//     WHEN SF.DomainException =>
//         PUT_LINE("DOMAIN EXCEPTION");
//     WHEN SF.OverflowException =>
//         PUT_LINE("OVERFLOW EXCEPTION");
//     WHEN SF.MaxIterationsException =>
//         PUT_LINE("MAXIMUM ITERATIONS REACHED EXCEPTION");
//     WHEN OTHERS =>
//         PUT_LINE("GENERAL EXCEPTION");
// END TEST_SF_THETA_D;
function TEST_SF_THETA(stat, func, args, val_in, tol, funcName)
{
    var status = 0; //GSL_SUCCESS;
    var r = 0.0;

    try
    {
        r = func( args );
        if ( sf_internalLog )
        {
            console.log( r );
        }
        stat.Integer = stat.Integer + test_sf_val( r, val_in, tol, funcName );
    }
    catch ( e )
    {
        console.log( e );
        console.trace( );
    }
} // TEST_SF_THETA
//#define TEST_SF_VAL(stat, func, args, val_in, tol)
// { double val = func args; stat += test_sf_val(val, val_in, tol, #func #args); } 

// //#define TEST_SF_THETA(stat, func, args, val_in, tol)
// // { int status; theta=args; status = func (&theta);  stat += test_sf_val(theta, val_in, tol, #func #args); }

// ----------------------------------------------------------------------------

export function test_clausen( )
{
    var s = { Integer: 0 };

    console.log( "Test Clausen Function ..." );

    console.log( "    ... gsl_sf_clausen_e" );
    TEST_SF_D( s, gsl_sf_clausen_e, M_PI / 20.0,               0.4478882448133546, TEST_TOL0, "gsl_sf_clausen_e" );
    TEST_SF_D( s, gsl_sf_clausen_e, M_PI / 6.0,                0.8643791310538927, TEST_TOL0, "gsl_sf_clausen_e" );
    TEST_SF_D( s, gsl_sf_clausen_e, M_PI / 3.0,                1.0149416064096535, TEST_TOL0, "gsl_sf_clausen_e" );
    TEST_SF_D( s, gsl_sf_clausen_e, 2.0 * M_PI + M_PI / 3.0,   1.0149416064096535, TEST_TOL0, "gsl_sf_clausen_e" );
    TEST_SF_D( s, gsl_sf_clausen_e, 100.0 * M_PI + M_PI / 3.0, 1.0149416064096535, TEST_TOL0, "gsl_sf_clausen_e" );

    return s;

} // test_clausen

// ----------------------------------------------------------------------------

export function test_coupling( )
{

    var s = { Integer: 0 };

    console.log( "Test Coupling Functions ..." );

    // Test 3j

    console.log( "    ... gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 0, i2: 1, i3: 1, i4: 0, i5:  1, i6: -1 }, Math.sqrt( 1.0 / 2.0 ),   TEST_TOL0, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 1, i2: 1, i3: 2, i4: 1, i5: -1, i6:  0 }, Math.sqrt( 1.0 / 6.0 ),   TEST_TOL0, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 2, i2: 4, i3: 6, i4: 0, i5:  2, i6: -2 }, Math.sqrt( 8.0 / 105.0 ), TEST_TOL0, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 4, i2: 4, i3: 8, i4: 0, i5:  0, i6:  0 }, Math.sqrt( 2.0 / 35.0 ),  TEST_TOL0, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 4, i2: 4, i3: 8, i4: 2, i5: -2, i6:  0 }, 2.0 / 3.0 * Math.sqrt( 2.0 / 35.0 ), TEST_TOL2, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 4, i2: 4, i3: 8, i4: 4, i5: -4, i6:  0 }, 1.0 / (3.0 * Math.sqrt( 70.0 )), TEST_TOL2, "gsl_sf_coupling_3j_e" );

    // Test 3j error checking

//     --TEST_SF(s, gsl_sf_coupling_3j_e'Access, (-1, 1, 2, 1, -1, 0), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_3j_e");
//     --TEST_SF(s, gsl_sf_coupling_3j_e'Access, (1, -1, 2, 1, -1, 0), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_3j_e");
//     --TEST_SF(s, gsl_sf_coupling_3j_e'Access, (1, 1, -2, 1, -1, 0), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_3j_e");

    // Test |m_i|<=j_i

    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 1, i2: 1, i3: 2, i4: 2, i5: -1, i6: 0 }, 0.0, 0.0, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 1, i2: 1, i3: 2, i4: 1, i5: -2, i6: 0 }, 0.0, 0.0, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 1, i2: 1, i3: 2, i4: 1, i5: -1, i6: 3 }, 0.0, 0.0, "gsl_sf_coupling_3j_e" );

    // Test triangle condition j1 + j2 >= j, j >= j2 - j1, j>= j1 - j2

    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 1, i2: 1, i3: 3, i4: 1, i5: -1, i6: 0 }, 0.0, 0.0, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 1, i2: 4, i3: 2, i4: 1, i5: -1, i6: 0 }, 0.0, 0.0, "gsl_sf_coupling_3j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_3j_e, { i1: 4, i2: 1, i3: 2, i4: 1, i5: -1, i6: 0 }, 0.0, 0.0, "gsl_sf_coupling_3j_e" );

    // Test 6j

    console.log( "    ... gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 2, i2: 2, i3: 4, i4: 2, i5: 2, i6: 2 },  1.0 / 6.0, TEST_TOL0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 4, i2: 4, i3: 2, i4: 4, i5: 4, i6: 4 }, -1.0 / 10.0, TEST_TOL0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 4, i2: 4, i3: 2, i4: 4, i5: 4, i6: 2 },  1.0 / 6.0, TEST_TOL0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 4, i2: 4, i3: 2, i4: 2, i5: 2, i6: 2 }, -0.5 / Math.sqrt( 5.0 ), TEST_TOL0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 4, i2: 4, i3: 4, i4: 2, i5: 2, i6: 2 },  Math.sqrt( 7.0 / 3.0 ) / 10.0, TEST_TOL0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 6, i2: 6, i3: 6, i4: 4, i5: 4, i6: 4 }, -Math.sqrt( 3.0 / 5.0 ) / 14.0, TEST_TOL0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 6, i2: 6, i3: 6, i4: 4, i5: 4, i6: 2 }, -Math.sqrt( 3.0 / 5.0 ) / 7.0, TEST_TOL0, "gsl_sf_coupling_6j_e" );

    // Test 6j error checking

//     --TEST_SF(s, gsl_sf_coupling_6j_e'Access, (-2, 2, 4, 2, 2, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_6j_e");
//     --TEST_SF(s, gsl_sf_coupling_6j_e'Access, (2, -2, 4, 2, 2, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_6j_e");
//     --TEST_SF(s, gsl_sf_coupling_6j_e'Access, (2, 2, -4, 2, 2, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_6j_e");
//     --TEST_SF(s, gsl_sf_coupling_6j_e'Access, (2, 2, 4, -2, 2, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_6j_e");
//     --TEST_SF(s, gsl_sf_coupling_6j_e'Access, (2, 2, 4, 2, -2, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_6j_e");
//     --TEST_SF(s, gsl_sf_coupling_6j_e'Access, (2, 2, 4, 2, 2, -2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_6j_e");

    // Test 6j triangle conditions

    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 2, i2: 2, i3: 4, i4: 2, i5: 2, i6: 7}, 0.0, 0.0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 2, i2: 2, i3: 4, i4: 2, i5: 7, i6: 2}, 0.0, 0.0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 2, i2: 2, i3: 4, i4: 7, i5: 2, i6: 2}, 0.0, 0.0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 2, i2: 2, i3: 7, i4: 2, i5: 2, i6: 2}, 0.0, 0.0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 2, i2: 7, i3: 4, i4: 2, i5: 2, i6: 2}, 0.0, 0.0, "gsl_sf_coupling_6j_e" );
    TEST_SF_6I( s, gsl_sf_coupling_6j_e, { i1: 7, i2: 2, i3: 4, i4: 2, i5: 2, i6: 2}, 0.0, 0.0, "gsl_sf_coupling_6j_e" );

    // Test 9j

    console.log( "    ... gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 2, i3:  4, i4: 3, i5: 3, i6: 2, i7: 1, i8: 1, i9: 2 }, -Math.sqrt( 1.0 / 6.0 ) / 10.0, TEST_TOL2, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 8, i2: 4, i3: 10, i4: 7, i5: 3, i6: 8, i7: 1, i8: 1, i9: 2 },  Math.sqrt( 7.0 / 3.0 ) / 60.0, TEST_TOL2, "gsl_sf_coupling_9j_e" );

    // Test 9j error checking

//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (-4, 2, 4, 3, 3, 2, 1, 1, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");
//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (4, -2, 4, 3, 3, 2, 1, 1, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");
//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (4, 2, -4, 3, 3, 2, 1, 1, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");
//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (4, 2, 4, -3, 3, 2, 1, 1, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");
//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (4, 2, 4, 3, -3, 2, 1, 1, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");
//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (4, 2, 4, 3, 3, -2, 1, 1, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");
//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (4, 2, 4, 3, 3, 2, -1, 1, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");
//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (4, 2, 4, 3, 3, 2, 1, -1, 2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");
//     --TEST_SF(s, gsl_sf_coupling_9j_e'Access, (4, 2, 4, 3, 3, 2, 1, 1, -2), GSL_NAN, GSL_NAN, GSL_EDOM, "gsl_sf_coupling_9j_e");

    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 10, i2: 2, i3: 4, i4: 3, i5: 3, i6: 2, i7: 1, i8: 1, i9: 2 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 10, i3: 4, i4: 3, i5: 3, i6: 2, i7: 1, i8: 1, i9: 2 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 2, i3: 10, i4: 3, i5: 3, i6: 2, i7: 1, i8: 1, i9: 2 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 2, i3: 4, i4: 10, i5: 3, i6: 2, i7: 1, i8: 1, i9: 2 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 2, i3: 4, i4: 3, i5: 10, i6: 2, i7: 1, i8: 1, i9: 2 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 2, i3: 4, i4: 3, i5: 3, i6: 10, i7: 1, i8: 1, i9: 2 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 2, i3: 4, i4: 3, i5: 3, i6: 2, i7: 10, i8: 1, i9: 2 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 2, i3: 4, i4: 3, i5: 3, i6: 2, i7: 1, i8: 10, i9: 2 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );
    TEST_SF_9I( s, gsl_sf_coupling_9j_e, { i1: 4, i2: 2, i3: 4, i4: 3, i5: 3, i6: 2, i7: 1, i8: 1, i9: 10 }, 0.0, 0.0, "gsl_sf_coupling_9j_e" );

    return s;

} // test_coupling

// ----------------------------------------------------------------------------

export function test_dawson( )
{
    var s = { Integer: 0 };

    console.log( "Test Dawson Functions ..." );

    console.log( "    ... gsl_sf_dawson_e" );
    TEST_SF_D( s, gsl_sf_dawson_e, 1.0e-15, 1.0e-15,                  TEST_TOL0, "gsl_sf_dawson_e" );
    TEST_SF_D( s, gsl_sf_dawson_e, 0.5,     0.4244363835020222959,    TEST_TOL0, "gsl_sf_dawson_e" );
    TEST_SF_D( s, gsl_sf_dawson_e, 2.0,     0.30134038892379196603,   TEST_TOL0, "gsl_sf_dawson_e" );
    TEST_SF_D( s, gsl_sf_dawson_e, 1000.0,  0.0005000002500003750009, TEST_TOL0, "gsl_sf_dawson_e" );

    return s;

} // test_dawson

// ----------------------------------------------------------------------------

export function test_debye( )
{
    var s = { Integer: 0 };

    console.log( "Test Debye Functions ..." );

    console.log( "    ... gsl_sf_debye_1_e" );
    TEST_SF_D( s, gsl_sf_debye_1_e, 0.1,  0.975277750004723276, TEST_TOL0, "gsl_sf_debye_1_e" );
    TEST_SF_D( s, gsl_sf_debye_1_e, 1.0,  0.777504634112248239, TEST_TOL0, "gsl_sf_debye_1_e" );
    TEST_SF_D( s, gsl_sf_debye_1_e, 10.0, 0.164443465679946027, TEST_TOL0, "gsl_sf_debye_1_e" );

    console.log( "    ... gsl_sf_debye_2_e" );
    TEST_SF_D( s, gsl_sf_debye_2_e, 0.1,  0.967083287045302664,  TEST_TOL0, "gsl_sf_debye_2_e" );
    TEST_SF_D( s, gsl_sf_debye_2_e, 1.0,  0.70787847562782924,   TEST_TOL0, "gsl_sf_debye_2_e" );
    TEST_SF_D( s, gsl_sf_debye_2_e, 10.0, 0.0479714980201218708, TEST_TOL0, "gsl_sf_debye_2_e" );

    console.log( "    ... gsl_sf_debye_3_e" );
    TEST_SF_D( s, gsl_sf_debye_3_e, 0.1,  0.962999940487211048,  TEST_TOL0, "gsl_sf_debye_3_e" );
    TEST_SF_D( s, gsl_sf_debye_3_e, 1.0,  0.674415564077814667,  TEST_TOL0, "gsl_sf_debye_3_e" );
    TEST_SF_D( s, gsl_sf_debye_3_e, 10.0, 0.0192957656903454886, TEST_TOL0, "gsl_sf_debye_3_e" );

    console.log( "    ... gsl_sf_debye_4_e" );
    TEST_SF_D( s, gsl_sf_debye_4_e, 0.1,  0.960555486124335944,   TEST_TOL0, "gsl_sf_debye_4_e" );
    TEST_SF_D( s, gsl_sf_debye_4_e, 1.0,  0.654874068886737049,   TEST_TOL0, "gsl_sf_debye_4_e" );
    TEST_SF_D( s, gsl_sf_debye_4_e, 10.0, 0.00967367556027115896, TEST_TOL0, "gsl_sf_debye_4_e" );

    console.log( "    ... gsl_sf_debye_5_e" );
    TEST_SF_D( s, gsl_sf_debye_5_e, 0.1,  0.95892849428310568745,  TEST_TOL0, "gsl_sf_debye_5_e" );
    TEST_SF_D( s, gsl_sf_debye_5_e, 1.0,  0.6421002580217790246,   TEST_TOL0, "gsl_sf_debye_5_e" );
    TEST_SF_D( s, gsl_sf_debye_5_e, 10.0, 0.005701535852992908538, TEST_TOL0, "gsl_sf_debye_5_e" );

    console.log( "    ... gsl_sf_debye_6_e" );
    TEST_SF_D( s, gsl_sf_debye_6_e, 0.1,  0.95776777382605465878,  TEST_TOL0, "gsl_sf_debye_6_e" );
    TEST_SF_D( s, gsl_sf_debye_6_e, 1.0,  0.63311142583495107588,  TEST_TOL0, "gsl_sf_debye_6_e" );
    TEST_SF_D( s, gsl_sf_debye_6_e, 10.0, 3.7938493294615955279e-3, TEST_TOL0, "gsl_sf_debye_6_e ");

    return s;

} // test_debye

// ----------------------------------------------------------------------------

export function test_elementary( )
{
    var s = { Integer: 0 };
    var x = 0.2 * GSL_DBL_MAX;

    console.log( "Test Elementary Functions ..." );

    console.log( "    ... gsl_sf_multiply_e" );
    TEST_SF_DD( s,  gsl_sf_multiply_e, { x: -3.0, y: 2.0 },     -6.0,                TEST_TOL0, "gsl_sf_multiply_e" );
    TEST_SF_DD( s,  gsl_sf_multiply_e, { x: x,    y: 1.0 / x },  1.0,                TEST_TOL0, "gsl_sf_multiply_e" );
    TEST_SF_DD( s,  gsl_sf_multiply_e, { x: x,    y: 0.2 },      0.04 * GSL_DBL_MAX, TEST_TOL1, "gsl_sf_multiply_e" );
    TEST_SF_DD( s,  gsl_sf_multiply_e, { x: x,    y: 4.0 },      0.8 * GSL_DBL_MAX,  TEST_TOL1, "gsl_sf_multiply_e" );
    // try
    // {
    //     let r = { val: 0.0, err: 0.0 };
    //     r = gsl_sf_multiply_e( GSL_DBL_MAX, 1.1 );
    //     s.Integer = s.Integer + 1;
    // }
    // catch ( e )
    // {
    //     // Do nothing
    //     console.log( e );
    // }
    // try
    // {
    //     let r = { val: 0.0, err: 0.0 };
    //     r = gsl_sf_multiply_e( GSL_DBL_MIN, GSL_DBL_MIN );
    //     s.Integer = s.Integer + 1;
    // }
    // catch ( e )
    // {
    //     // Do nothing
    //     console.log( e );
    // }
    // try
    // {
    //     let r = { val: 0.0, err: 0.0 };
    //     r = gsl_sf_multiply_e( GSL_DBL_MIN, -GSL_DBL_MIN );
    //     s = s + 1;
    // }
    // catch ( e )
    // {
    //     // Do nothing
    //     console.log( e );
    // }

    return s;

} // test_elementary

// ----------------------------------------------------------------------------

export function test_ellint( )
{
    var s  = { Integer: 0 };
    var mode = GSL_MODE_DEFAULT;

    console.log("Test Elliptic Integrals ...");

    console.log("    ... gsl_sf_ellint_Kcomp_e");
    TEST_SF_DM(s, gsl_sf_ellint_Kcomp_e, { x:  0.99, m: mode }, 3.3566005233611923760, TEST_TOL0, "gsl_sf_ellint_Kcomp_e");
    TEST_SF_DM(s, gsl_sf_ellint_Kcomp_e, { x:  0.50, m: mode }, 1.6857503548125960429, TEST_TOL0, "gsl_sf_ellint_Kcomp_e");
    TEST_SF_DM(s, gsl_sf_ellint_Kcomp_e, { x: 0.010, m: mode }, 1.5708355989121522360, TEST_TOL0, "gsl_sf_ellint_Kcomp_e");

    // Bug report from Thies Heidecke
    TEST_SF_DM(s, gsl_sf_ellint_Kcomp_e, { x: 0.99999999906867742538, m: mode }, 11.4369284843320018031, TEST_SNGL, "gsl_sf_ellint_Kcomp_e");

    console.log("    ... gsl_sf_ellint_Ecomp_e");
    TEST_SF_DM(s,  gsl_sf_ellint_Ecomp_e, { x: 0.99, m: mode }, 1.0284758090288040010, TEST_TOL0, "gsl_sf_ellint_Ecomp_e");
    TEST_SF_DM(s,  gsl_sf_ellint_Ecomp_e, { x: 0.50, m: mode }, 1.4674622093394271555, TEST_TOL0, "gsl_sf_ellint_Ecomp_e");
    TEST_SF_DM(s,  gsl_sf_ellint_Ecomp_e, { x: 0.01, m: mode }, 1.5707570561503852873, TEST_TOL0, "gsl_sf_ellint_Ecomp_e");

    console.log("    ... gsl_sf_ellint_Pcomp_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_Pcomp_e, { x: 0.99, y: 0.1, m: mode }, 3.13792612351836506315593, TEST_TOL0, "gsl_sf_ellint_Pcomp_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_Pcomp_e, { x: 0.50, y: 0.1, m: mode }, 1.60455249360848890075108, TEST_TOL0, "gsl_sf_ellint_Pcomp_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_Pcomp_e, { x: 0.01, y: 0.1, m: mode }, 1.49773208536003801277453, TEST_TOL0, "gsl_sf_ellint_Pcomp_e");

    console.log("    ... gsl_sf_ellint_Dcomp_e");
    TEST_SF_DM(s,  gsl_sf_ellint_Dcomp_e, { x: 0.99, m: mode }, 2.375395076351788975665323192, TEST_TOL0, "gsl_sf_ellint_Dcomp_e");
    TEST_SF_DM(s,  gsl_sf_ellint_Dcomp_e, { x: 0.50, m: mode }, 0.8731525818926755496456335628, TEST_TOL0, "gsl_sf_ellint_Dcomp_e");
    TEST_SF_DM(s,  gsl_sf_ellint_Dcomp_e, { x: 0.01, m: mode }, 0.7854276176694868932799393751, TEST_TOL0, "gsl_sf_ellint_Dcomp_e");

    console.log("    ... gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI/3.0, y: 0.99, m: mode }, 1.3065333392738766762, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI/3.0, y: 0.50, m: mode }, 1.0895506700518854093, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI/3.0, y: 0.01, m: mode }, 1.0472129063770918952, TEST_TOL0, "gsl_sf_ellint_F_e");

    console.log("    ... gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI/3.0, y: 0.99, m: mode }, 0.8704819220377943536, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI/3.0, y: 0.50, m: mode }, 1.0075555551444720293, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI/3.0, y: 0.01, m: mode }, 1.0471821963889481104, TEST_TOL0, "gsl_sf_ellint_E_e");

    console.log("    ... gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI/3.0, y: 0.99, z: 0.5, m: mode }, 1.1288726598764099882, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI/3.0, y: 0.50, z: 0.5, m: mode }, 0.9570574331323584890, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI/3.0, y: 0.01, z: 0.5, m: mode }, 0.9228868127118118465, TEST_TOL0, "gsl_sf_ellint_P_e");

    console.log("    ... gsl_sf_ellint_RF_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_RF_e, { x: 5.0e-11, y: 1.0e-10, z: 1.0, m: mode }, 12.36441982979439, TEST_TOL0, "gsl_sf_ellint_RF_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_RF_e, { x: 1.0, y: 2.0, z: 3.0, m: mode }, 0.7269459354689082, TEST_TOL0, "gsl_sf_ellint_RF_e");

    console.log("    ... gsl_sf_ellint_RD_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_RD_e, { x: 5.0e-11, y: 1.0e-10, z: 1.0, m: mode }, 34.0932594919337362, TEST_TOL0, "gsl_sf_ellint_RD_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_RD_e, { x: 1.0, y: 2.0, z: 3.0, m: mode }, 0.2904602810289906, TEST_TOL0, "gsl_sf_ellint_RD_e");

    console.log("    ... gsl_sf_ellint_RC_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_RC_e, { x: 1.0, y: 2.0, m: mode }, 0.7853981633974482, TEST_TOL0, "gsl_sf_ellint_RC_e");

    console.log("    ... gsl_sf_ellint_RJ_e");
    TEST_SF_4DM(s,  gsl_sf_ellint_RJ_e, { x: 2.0, y: 3.0, z: 4.0, p: 5.0, m: mode }, 0.1429757966715675, TEST_TOL0, "gsl_sf_ellint_RJ_e");

    // E, argument phi > pi/2

    console.log("    ... gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI/2.0, y: 0.99, m: mode }, 1.02847580902880400098389, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI/2.0, y: 0.50, m: mode }, 1.46746220933942715545980, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI/2.0, y: 0.01, m: mode }, 1.57075705615038528733708, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 2.0*M_PI/3.0, y: 0.99, m: mode }, 1.18646969601981364833972, TEST_TOL1, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 2.0*M_PI/3.0, y: 0.50, m: mode }, 1.92736886353438228163734, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 2.0*M_PI/3.0, y: 0.01, m: mode }, 2.09433191591182246425715, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI, y: 0.99, m: mode }, 2.05695161805760800196777, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI, y: 0.50, m: mode }, 2.93492441867885431091959, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: M_PI, y: 0.01, m: mode }, 3.14151411230077057467416, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 4.0*M_PI/3.0, y: 0.99, m: mode }, 2.92743354009540235559582, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 4.0*M_PI/3.0, y: 0.50, m: mode }, 3.94247997382332634020184, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 4.0*M_PI/3.0, y: 0.01, m: mode }, 4.18869630868971868509117, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 3.0*M_PI/2.0, y: 0.99, m: mode }, 3.08542742708641200295166, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 3.0*M_PI/2.0, y: 0.50, m: mode }, 4.40238662801828146637939, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 3.0*M_PI/2.0, y: 0.01, m: mode }, 4.71227116845115586201123, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 5.0*M_PI/3.0, y: 0.99, m: mode }, 3.24342131407742165030750, TEST_TOL1, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 5.0*M_PI/3.0, y: 0.50, m: mode }, 4.86229328221323659255693, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 5.0*M_PI/3.0, y: 0.01, m: mode }, 5.23584602821259303893130, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 2.0*M_PI, y: 0.99, m: mode }, 4.11390323611521600393555, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 2.0*M_PI, y: 0.50, m: mode }, 5.86984883735770862183918, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 2.0*M_PI, y: 0.01, m: mode }, 6.28302822460154114934831, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 7.0*M_PI/3.0, y: 0.99, m: mode }, 4.98438515815301035756360, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 7.0*M_PI/3.0, y: 0.50, m: mode }, 6.87740439250218065112143, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: 7.0*M_PI/3.0, y: 0.01, m: mode }, 7.33021042099048925976532, TEST_TOL0, "gsl_sf_ellint_E_e");

    // Test some negative arguments, phi < 0

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -M_PI/2.0, y: 0.99, m: mode }, -1.02847580902880400098389, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -M_PI/2.0, y: 0.50, m: mode }, -1.46746220933942715545980, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -M_PI/2.0, y: 0.01, m: mode }, -1.57075705615038528733708, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -2.0*M_PI/3.0, y: 0.99, m: mode }, -1.18646969601981364833972, TEST_TOL1, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -2.0*M_PI/3.0, y: 0.50, m: mode }, -1.92736886353438228163734, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -2.0*M_PI/3.0, y: 0.01, m: mode }, -2.09433191591182246425715, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -M_PI, y: 0.99, m: mode }, -2.05695161805760800196777, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -M_PI, y: 0.50, m: mode }, -2.93492441867885431091959, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -M_PI, y: 0.01, m: mode }, -3.14151411230077057467416, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -4.0*M_PI/3.0, y: 0.99, m: mode }, -2.92743354009540235559582, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -4.0*M_PI/3.0, y: 0.50, m: mode }, -3.94247997382332634020184, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -4.0*M_PI/3.0, y: 0.01, m: mode }, -4.18869630868971868509117, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -3.0*M_PI/2.0, y: 0.99, m: mode }, -3.08542742708641200295166, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -3.0*M_PI/2.0, y: 0.50, m: mode }, -4.40238662801828146637939, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -3.0*M_PI/2.0, y: 0.01, m: mode }, -4.71227116845115586201123, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -5.0*M_PI/3.0, y: 0.99, m: mode }, -3.24342131407742165030750, TEST_TOL1, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -5.0*M_PI/3.0, y: 0.50, m: mode }, -4.86229328221323659255693, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -5.0*M_PI/3.0, y: 0.01, m: mode }, -5.23584602821259303893130, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -2.0*M_PI, y: 0.99, m: mode }, -4.11390323611521600393555, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -2.0*M_PI, y: 0.50, m: mode }, -5.86984883735770862183918, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -2.0*M_PI, y: 0.01, m: mode }, -6.28302822460154114934831, TEST_TOL0, "gsl_sf_ellint_E_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -7.0*M_PI/3.0, y: 0.99, m: mode }, -4.98438515815301035756360, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -7.0*M_PI/3.0, y: 0.50, m: mode }, -6.87740439250218065112143, TEST_TOL0, "gsl_sf_ellint_E_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_E_e, { x: -7.0*M_PI/3.0, y: 0.01, m: mode }, -7.33021042099048925976532, TEST_TOL0, "gsl_sf_ellint_E_e");

//     -- F, argument phi > pi/2

    console.log("    ... gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI/2.0, y: 0.99, m: mode }, 3.35660052336119237603347, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI/2.0, y: 0.50, m: mode }, 1.68575035481259604287120, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI/2.0, y: 0.01, m: mode }, 1.57083559891215223602641, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 2.0*M_PI/3.0, y: 0.99, m: mode }, 5.40666770744850807588478, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 2.0*M_PI/3.0, y: 0.50, m: mode }, 2.28195003957330667648585, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 2.0*M_PI/3.0, y: 0.01, m: mode }, 2.09445829144721257687207, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI, y: 0.99, m: mode }, 6.71320104672238475206694, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI, y: 0.50, m: mode }, 3.37150070962519208574241, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: M_PI, y: 0.01, m: mode }, 3.14167119782430447205281, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 4.0*M_PI/3.0, y: 0.99, m: mode }, 8.01973438599626142824910, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 4.0*M_PI/3.0, y: 0.50, m: mode }, 4.46105137967707749499897, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 4.0*M_PI/3.0, y: 0.01, m: mode }, 4.18888410420139636723356, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 3.0*M_PI/2.0, y: 0.99, m: mode }, 10.0698015700835771281004, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 3.0*M_PI/2.0, y: 0.50, m: mode }, 5.05725106443778812861361, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 3.0*M_PI/2.0, y: 0.01, m: mode }, 4.71250679673645670807922, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 5.0*M_PI/3.0, y: 0.99, m: mode }, 12.1198687541708928279517, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 5.0*M_PI/3.0, y: 0.50, m: mode }, 5.65345074919849876222825, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 5.0*M_PI/3.0, y: 0.01, m: mode }, 5.23612948927151704892488, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 2.0*M_PI, y: 0.99, m: mode }, 13.4264020934447695041339, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 2.0*M_PI, y: 0.50, m: mode }, 6.74300141925038417148481, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 2.0*M_PI, y: 0.01, m: mode }, 6.28334239564860894410562, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 7.0*M_PI/3.0, y: 0.99, m: mode }, 14.7329354327186461803160, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 7.0*M_PI/3.0, y: 0.50, m: mode }, 7.83255208930226958074138, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: 7.0*M_PI/3.0, y: 0.01, m: mode }, 7.33055530202570083928637, TEST_TOL0, "gsl_sf_ellint_F_e");

    // F, negative argument phi < 0

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -M_PI/2.0, y: 0.99, m: mode }, -3.35660052336119237603347, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -M_PI/2.0, y: 0.50, m: mode }, -1.68575035481259604287120, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -M_PI/2.0, y: 0.01, m: mode }, -1.57083559891215223602641, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -2.0*M_PI/3.0, y: 0.99, m: mode }, -5.40666770744850807588478, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -2.0*M_PI/3.0, y: 0.50, m: mode }, -2.28195003957330667648585, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -2.0*M_PI/3.0, y: 0.01, m: mode }, -2.09445829144721257687207, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -M_PI, y: 0.99, m: mode }, -6.71320104672238475206694, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -M_PI, y: 0.50, m: mode }, -3.37150070962519208574241, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -M_PI, y: 0.01, m: mode }, -3.14167119782430447205281, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -4.0*M_PI/3.0, y: 0.99, m: mode }, -8.01973438599626142824910, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -4.0*M_PI/3.0, y: 0.50, m: mode }, -4.46105137967707749499897, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -4.0*M_PI/3.0, y: 0.01, m: mode }, -4.18888410420139636723356, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -3.0*M_PI/2.0, y: 0.99, m: mode }, -10.0698015700835771281004, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -3.0*M_PI/2.0, y: 0.50, m: mode }, -5.05725106443778812861361, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -3.0*M_PI/2.0, y: 0.01, m: mode }, -4.71250679673645670807922, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -5.0*M_PI/3.0, y: 0.99, m: mode }, -12.1198687541708928279517, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -5.0*M_PI/3.0, y: 0.50, m: mode }, -5.65345074919849876222825, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -5.0*M_PI/3.0, y: 0.01, m: mode }, -5.23612948927151704892488, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -2.0*M_PI, y: 0.99, m: mode }, -13.4264020934447695041339, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -2.0*M_PI, y: 0.50, m: mode }, -6.74300141925038417148481, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -2.0*M_PI, y: 0.01, m: mode }, -6.28334239564860894410562, TEST_TOL0, "gsl_sf_ellint_F_e");

    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -7.0*M_PI/3.0, y: 0.99, m: mode }, -14.7329354327186461803160, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -7.0*M_PI/3.0, y: 0.50, m: mode }, -7.83255208930226958074138, TEST_TOL0, "gsl_sf_ellint_F_e");
    TEST_SF_DDM(s,  gsl_sf_ellint_F_e, { x: -7.0*M_PI/3.0, y: 0.01, m: mode }, -7.33055530202570083928637, TEST_TOL0, "gsl_sf_ellint_F_e");

    // P, argument phi > pi/2

    console.log("    ... gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI/2.0, y: 0.99, z: -0.1, m: mode }, 3.61678162163246646783050, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI/2.0, y: 0.50, z: -0.1, m: mode }, 1.78030349465454812629168, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI/2.0, y: 0.01, z: -0.1, m: mode }, 1.65580719756898353270922, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 2.0*M_PI/3.0, y: 0.99, z: -0.1, m: mode }, 5.88008918207571119911983, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 2.0*M_PI/3.0, y: 0.50, z: -0.1, m: mode }, 2.43655207300356008717867, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 2.0*M_PI/3.0, y: 0.01, z: -0.1, m: mode }, 2.23211110528200554950903, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI, y: 0.99, z: -0.1, m: mode }, 7.23356324326493293566099, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI, y: 0.50, z: -0.1, m: mode }, 3.56060698930909625258336, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: M_PI, y: 0.01, z: -0.1, m: mode }, 3.31161439513796706541844, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 4.0*M_PI/3.0, y: 0.99, z: -0.1, m: mode }, 8.58703730445415467220216, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 4.0*M_PI/3.0, y: 0.50, z: -0.1, m: mode }, 4.68466190561463241798805, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 4.0*M_PI/3.0, y: 0.01, z: -0.1, m: mode }, 4.39111768499392858132786, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 3.0*M_PI/2.0, y: 0.99, z: -0.1, m: mode }, 10.8503448648973994034915, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 3.0*M_PI/2.0, y: 0.50, z: -0.1, m: mode }, 5.34091048396364437887504, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 3.0*M_PI/2.0, y: 0.01, z: -0.1, m: mode }, 4.96742159270695059812767, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 5.0*M_PI/3.0, y: 0.99, z: -0.1, m: mode }, 13.1136524253406441347808, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 5.0*M_PI/3.0, y: 0.50, z: -0.1, m: mode }, 5.99715906231265633976204, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 5.0*M_PI/3.0, y: 0.01, z: -0.1, m: mode }, 5.54372550041997261492747, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 2.0*M_PI, y: 0.99, z: -0.1, m: mode }, 14.4671264865298658713220, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 2.0*M_PI, y: 0.50, z: -0.1, m: mode }, 7.12121397861819250516672, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 2.0*M_PI, y: 0.01, z: -0.1, m: mode }, 6.62322879027593413083689, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 7.0*M_PI/3.0, y: 0.99, z: -0.1, m: mode }, 15.8206005477190876078631, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 7.0*M_PI/3.0, y: 0.50, z: -0.1, m: mode }, 8.24526889492372867057141, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: 7.0*M_PI/3.0, y: 0.01, z: -0.1, m: mode }, 7.70273208013189564674630, TEST_TOL0, "gsl_sf_ellint_P_e");

    // P, negative argument phi < 0

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -M_PI/2.0, y: 0.99, z: -0.1, m: mode }, -3.61678162163246646783050, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -M_PI/2.0, y: 0.50, z: -0.1, m: mode }, -1.78030349465454812629168, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -M_PI/2.0, y: 0.01, z: -0.1, m: mode }, -1.65580719756898353270922, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -2.0*M_PI/3.0, y: 0.99, z: -0.1, m: mode }, -5.88008918207571119911983, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -2.0*M_PI/3.0, y: 0.50, z: -0.1, m: mode }, -2.43655207300356008717867, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -2.0*M_PI/3.0, y: 0.01, z: -0.1, m: mode }, -2.23211110528200554950903, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -M_PI, y: 0.99, z: -0.1, m: mode }, -7.23356324326493293566099, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -M_PI, y: 0.50, z: -0.1, m: mode }, -3.56060698930909625258336, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -M_PI, y: 0.01, z: -0.1, m: mode }, -3.31161439513796706541844, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -4.0*M_PI/3.0, y: 0.99, z: -0.1, m: mode }, -8.58703730445415467220216, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -4.0*M_PI/3.0, y: 0.50, z: -0.1, m: mode }, -4.68466190561463241798805, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -4.0*M_PI/3.0, y: 0.01, z: -0.1, m: mode }, -4.39111768499392858132786, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -3.0*M_PI/2.0, y: 0.99, z: -0.1, m: mode }, -10.8503448648973994034915, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -3.0*M_PI/2.0, y: 0.50, z: -0.1, m: mode }, -5.34091048396364437887504, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -3.0*M_PI/2.0, y: 0.01, z: -0.1, m: mode }, -4.96742159270695059812767, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -5.0*M_PI/3.0, y: 0.99, z: -0.1, m: mode }, -13.1136524253406441347808, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -5.0*M_PI/3.0, y: 0.50, z: -0.1, m: mode }, -5.99715906231265633976204, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -5.0*M_PI/3.0, y: 0.01, z: -0.1, m: mode }, -5.54372550041997261492747, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -2.0*M_PI, y: 0.99, z: -0.1, m: mode }, -14.4671264865298658713220, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -2.0*M_PI, y: 0.50, z: -0.1, m: mode }, -7.12121397861819250516672, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -2.0*M_PI, y: 0.01, z: -0.1, m: mode }, -6.62322879027593413083689, TEST_TOL0, "gsl_sf_ellint_P_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -7.0*M_PI/3.0, y: 0.99, z: -0.1, m: mode }, -15.8206005477190876078631, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -7.0*M_PI/3.0, y: 0.50, z: -0.1, m: mode }, -8.24526889492372867057141, TEST_TOL0, "gsl_sf_ellint_P_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_P_e, { x: -7.0*M_PI/3.0, y: 0.01, z: -0.1, m: mode }, -7.70273208013189564674630, TEST_TOL0, "gsl_sf_ellint_P_e");

    // D, argument phi > pi/2

    console.log("    ... gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: M_PI/2.0, y: 0.99, z: 0.0, m: mode }, 2.375395076351788975665323192, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: M_PI/2.0, y: 0.50, z: 0.0, m: mode }, 0.8731525818926755496456335628, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: M_PI/2.0, y: 0.01, z: 0.0, m: mode }, 0.7854276176694868932799393751, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 2.0*M_PI/3.0, y: 0.99, z: 0.0, m: mode }, 4.305885125424644860264320635, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 2.0*M_PI/3.0, y: 0.50, z: 0.0, m: mode }, 1.418324704155697579394036402, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 2.0*M_PI/3.0, y: 0.01, z: 0.0, m: mode }, 1.263755353901126149206022061, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: M_PI, y: 0.99, z: 0.0, m: mode }, 4.750790152703577951330646444, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: M_PI, y: 0.50, z: 0.0, m: mode }, 1.746305163785351099291267125, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: M_PI, y: 0.01, z: 0.0, m: mode }, 1.570855235338973786559878750, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 4.0*M_PI/3.0, y: 0.99, z: 0.0, m: mode }, 5.195695179982511042396972113, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 4.0*M_PI/3.0, y: 0.50, z: 0.0, m: mode }, 2.074285623415004619188497818, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 4.0*M_PI/3.0, y: 0.01, z: 0.0, m: mode }, 1.877955116776821423913735408, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 3.0*M_PI/2.0, y: 0.99, z: 0.0, m: mode }, 7.126185229055366926995969476, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 3.0*M_PI/2.0, y: 0.50, z: 0.0, m: mode }, 2.619457745678026648936900687, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 3.0*M_PI/2.0, y: 0.01, z: 0.0, m: mode }, 2.356282853008460679839818125, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 5.0*M_PI/3.0, y: 0.99, z: 0.0, m: mode }, 9.056675278128222811594967044, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 5.0*M_PI/3.0, y: 0.50, z: 0.0, m: mode }, 3.164629867941048678685303509, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 5.0*M_PI/3.0, y: 0.01, z: 0.0, m: mode }, 2.834610589240099935765900794, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 2.0*M_PI, y: 0.99, z: 0.0, m: mode }, 9.501580305407155902661292832, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 2.0*M_PI, y: 0.50, z: 0.0, m: mode }, 3.492610327570702198582534249, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 2.0*M_PI, y: 0.01, z: 0.0, m: mode }, 3.141710470677947573119757500, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 7.0*M_PI/3.0, y: 0.99, z: 0.0, m: mode }, 9.946485332686088993727618315, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 7.0*M_PI/3.0, y: 0.50, z: 0.0, m: mode }, 3.820590787200355718479764901, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: 7.0*M_PI/3.0, y: 0.01, z: 0.0, m: mode }, 3.448810352115795210473614120, TEST_TOL0, "gsl_sf_ellint_D_e");

    // P, negative argument phi < 0

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -M_PI/2.0, y: 0.99, z: 0.0, m: mode }, -2.375395076351788975665323192, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -M_PI/2.0, y: 0.50, z: 0.0, m: mode }, -0.8731525818926755496456335628, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -M_PI/2.0, y: 0.01, z: 0.0, m: mode }, -0.7854276176694868932799393751, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -2.0*M_PI/3.0, y: 0.99, z: 0.0, m: mode }, -4.305885125424644860264320635, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -2.0*M_PI/3.0, y: 0.50, z: 0.0, m: mode }, -1.418324704155697579394036402, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -2.0*M_PI/3.0, y: 0.01, z: 0.0, m: mode }, -1.263755353901126149206022061, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -M_PI, y: 0.99, z: 0.0, m: mode }, -4.750790152703577951330646444, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -M_PI, y: 0.50, z: 0.0, m: mode }, -1.746305163785351099291267125, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -M_PI, y: 0.01, z: 0.0, m: mode }, -1.570855235338973786559878750, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -4.0*M_PI/3.0, y: 0.99, z: 0.0, m: mode }, -5.195695179982511042396972113, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -4.0*M_PI/3.0, y: 0.50, z: 0.0, m: mode }, -2.074285623415004619188497818, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -4.0*M_PI/3.0, y: 0.01, z: 0.0, m: mode }, -1.877955116776821423913735408, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -3.0*M_PI/2.0, y: 0.99, z: 0.0, m: mode }, -7.126185229055366926995969476, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -3.0*M_PI/2.0, y: 0.50, z: 0.0, m: mode }, -2.619457745678026648936900687, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -3.0*M_PI/2.0, y: 0.01, z: 0.0, m: mode }, -2.356282853008460679839818125, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -5.0*M_PI/3.0, y: 0.99, z: 0.0, m: mode }, -9.056675278128222811594967044, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -5.0*M_PI/3.0, y: 0.50, z: 0.0, m: mode }, -3.164629867941048678685303509, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -5.0*M_PI/3.0, y: 0.01, z: 0.0, m: mode }, -2.834610589240099935765900794, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -2.0*M_PI, y: 0.99, z: 0.0, m: mode }, -9.501580305407155902661292832, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -2.0*M_PI, y: 0.50, z: 0.0, m: mode }, -3.492610327570702198582534249, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -2.0*M_PI, y: 0.01, z: 0.0, m: mode }, -3.141710470677947573119757500, TEST_TOL0, "gsl_sf_ellint_D_e");

    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -7.0*M_PI/3.0, y: 0.99, z: 0.0, m: mode }, -9.946485332686088993727618315, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -7.0*M_PI/3.0, y: 0.50, z: 0.0, m: mode }, -3.820590787200355718479764901, TEST_TOL0, "gsl_sf_ellint_D_e");
    TEST_SF_3DM(s,  gsl_sf_ellint_D_e, { x: -7.0*M_PI/3.0, y: 0.01, z: 0.0, m: mode }, -3.448810352115795210473614120, TEST_TOL0, "gsl_sf_ellint_D_e");

    return s;

} // test_ellint

// ----------------------------------------------------------------------------

export function test_erf( )
{
    var s = { Integer: 0 };

    console.log( "Test Error Functions ..." );

    console.log("    ... gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, -10.0, 2.0, TEST_TOL0, "gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, -5.0000002, 1.9999999999984625433, TEST_TOL0, "gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, -5.0, 1.9999999999984625402, TEST_TOL0, "gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, -1.0, 1.8427007929497148693, TEST_TOL0, "gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, -0.5, 1.5204998778130465377, TEST_TOL0, "gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, 1.0, 0.15729920705028513066, TEST_TOL0, "gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, 3.0, 0.000022090496998585441373, TEST_TOL1, "gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, 7.0, 4.183825607779414399e-23, TEST_TOL2, "gsl_sf_erfc_e");
    TEST_SF_D(s, gsl_sf_erfc_e, 10.0, 2.0884875837625447570e-45, TEST_TOL2, "gsl_sf_erfc_e");

    console.log("    ... gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, -1.0, Math.log(1.842700792949714869), TEST_TOL0, "gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, -0.1, 0.106576400586522485015, TEST_TOL0, "gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, -1.0e-10,  1.1283791670318505967e-10, TEST_TOL0, "gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, 0.0, Math.log(1.0), TEST_TOL0, "gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, 1.0e-10, -1.128379167159174551e-10, TEST_TOL0, "gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, 0.001, -0.0011290158896213548027, TEST_TOL0, "gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, 0.1, -0.119304973737395598329, TEST_TOL0, "gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, 1.0, Math.log(0.15729920705028513066), TEST_TOL0, "gsl_sf_log_erfc_e");
    TEST_SF_D(s, gsl_sf_log_erfc_e, 10.0, Math.log(2.0884875837625447570e-45), TEST_TOL0, "gsl_sf_log_erfc_e");

    console.log("    ... gsl_sf_erf_e");
    TEST_SF_D(s, gsl_sf_erf_e, -10.0, -1.0000000000000000000, TEST_TOL0, "gsl_sf_erf_e");
    TEST_SF_D(s, gsl_sf_erf_e, 0.5, 0.5204998778130465377, TEST_TOL0, "gsl_sf_erf_e");
    TEST_SF_D(s, gsl_sf_erf_e, 1.0, 0.8427007929497148693, TEST_TOL0, "gsl_sf_erf_e");
    TEST_SF_D(s, gsl_sf_erf_e, 10.0, 1.0000000000000000000, TEST_TOL0, "gsl_sf_erf_e");

    console.log("    ... gsl_sf_erf_Z_e");
    TEST_SF_D(s, gsl_sf_erf_Z_e, 1.0,  0.24197072451914334980,   TEST_TOL0, "gsl_sf_erf_Z_e");
    console.log("    ... gsl_sf_erf_Q_e");
    TEST_SF_D(s, gsl_sf_erf_Q_e, 10.0, 7.619853024160526066e-24, TEST_TOL2, "gsl_sf_erf_Q_e");

    console.log("    ... gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, -20.0, 5.5209483621597631896e-88, TEST_TOL2, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, -10.0, 7.6945986267064193463e-23, TEST_TOL2, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, -1.0, 0.28759997093917836123, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e,  0.0, 0.79788456080286535588, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e,  1.0, 1.5251352761609812091, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 10.0, 10.098093233962511963, TEST_TOL2, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 20.0, 20.049753068527850542, TEST_TOL2, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 30.0, 30.033259667433677037, TEST_TOL2, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 50.0, 50.019984031905639809, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 80.0, 80.012496096798234468, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 150.0, 150.00666607420571802, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 300.0, 300.00333325926337415, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 900.0, 900.00111110836764382, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 1001.0, 1001.0009989990049990, TEST_TOL0, "gsl_sf_hazard_e");
    TEST_SF_D(s, gsl_sf_hazard_e, 2000.0, 2000.0004999997500003, TEST_TOL0, "gsl_sf_hazard_e");

    return s;

} // test_erf

// ----------------------------------------------------------------------------

export function test_exp( )
{
    var s  = { Integer: 0 };
    var sa = { Integer: 0 };
    var x  = 0.0;
    var re = { val: 0.0, err: 0.0 }; // Result;

    console.log("Test Exponential Functions ...");

    console.log("    ... gsl_sf_exp_e");
    TEST_SF_D(s, gsl_sf_exp_e, -10.0, Math.exp(-10.0), TEST_TOL0, "gsl_sf_exp_e");
    TEST_SF_D(s, gsl_sf_exp_e,  10.0, Math.exp( 10.0), TEST_TOL0, "gsl_sf_exp_e");

    console.log("    ... gsl_sf_exp_e10_e");
    sa.Integer = 0;
    re = gsl_sf_exp_e10_e(1.0); console.log(re);
    if (test_sf_frac_diff(re.val, M_E ) > TEST_TOL0)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.err > TEST_TOL1)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.e10 != 0)
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test(sa, "  gsl_sf_exp_e10_e(1.0, re)");

    sa.Integer = 0;
    re = gsl_sf_exp_e10_e(2000.0); console.log(re);
    if (test_sf_frac_diff(re.val, 3.88118019428363725) > TEST_TOL3)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.err > TEST_TOL5)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.e10 /= 868)
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test(sa, "  gsl_sf_exp_e10_e(2000.0, re)");

    console.log("    ... gsl_sf_exp_err_e");
    TEST_SF_DD(s, gsl_sf_exp_err_e, { x: -10.0, y: TEST_TOL1 }, Math.exp(-10.0), TEST_TOL1, "gsl_sf_exp_err_e");
    TEST_SF_DD(s, gsl_sf_exp_err_e, { x:  10.0, y: TEST_TOL1 }, Math.exp( 10.0), TEST_TOL1, "gsl_sf_exp_err_e");

    console.log("    ... gsl_sf_exp_err_e10_e");
    sa.Integer = 0;
    re = gsl_sf_exp_err_e10_e(1.0, TEST_SQRT_TOL0); console.log(re);
    if (test_sf_frac_diff(re.val, M_E ) > TEST_TOL1)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.err > 32.0 * TEST_SQRT_TOL0)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.e10 != 0)
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test(sa, "  gsl_sf_exp_err_e10_e(1.0, TEST_SQRT_TOL0, re)");

    sa.Integer = 0;
    re = gsl_sf_exp_err_e10_e(2000.0, 1.0e-10); console.log(re);
    if (test_sf_frac_diff(re.val, 3.88118019428363725 ) > TEST_TOL3)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.err > 1.0e-07)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.e10 != 868)
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test(sa, "  gsl_sf_exp_err_e10_e(2000.0, 1.0e-10, re)");

    console.log("    ... gsl_sf_exp_mult_e");
    x = 0.8 * GSL_LOG_DBL_MAX;
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: -10.0, y: 1.0e-06 }, 1.0e-06*Math.exp(-10.0), TEST_TOL0, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: -10.0, y: 2.0 },     2.0*Math.exp(-10.0),     TEST_TOL0, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: -10.0, y:-2.0 },    -2.0*Math.exp(-10.0),     TEST_TOL0, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x:  10.0, y: 1.0e-06 }, 1.0e-06*Math.exp( 10.0), TEST_TOL0, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x:  10.0, y:-2.0 },    -2.0*Math.exp( 10.0),     TEST_TOL0, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: x, y: 1.00001 },      1.00001*Math.exp(x),     TEST_TOL3, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: x, y: 1.000001 },     1.000001*Math.exp(x),    TEST_TOL3, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: x, y: 1.000000001 },  1.000000001*Math.exp(x), TEST_TOL3, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: x, y: 100.0 },        100.0*Math.exp(x),       TEST_TOL3, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: x, y: 1.0e+20 },      1.0e+20*Math.exp(x),     TEST_TOL3, "gsl_sf_exp_mult_e");
    TEST_SF_DD(s, gsl_sf_exp_mult_e, { x: x, y: Math.exp(-x)*Math.exp(M_LN2) },  2.0, TEST_TOL4, "gsl_sf_exp_mult_e");

    console.log("    ... gsl_sf_exp_mult_err_e");
    TEST_SF_4D(s, gsl_sf_exp_mult_err_e, { x: -10.0, dx: TEST_SQRT_TOL0, y: 2.0, dy: TEST_SQRT_TOL0 }, 2.0 * Math.exp(-10.0), TEST_SQRT_TOL0, "gsl_sf_exp_mult_err_e");
    TEST_SF_4D(s, gsl_sf_exp_mult_err_e, { x: x, dx: TEST_SQRT_TOL0 * x, y: Math.exp(-x) * Math.exp(M_LN2), dy: TEST_SQRT_TOL0 * Math.exp(-x) * Math.exp(M_LN2) },  2.0, TEST_SQRT_TOL0, "gsl_sf_exp_mult_err_e");

    console.log("    ... gsl_sf_exp_mult_e10_e");
    sa.Integer = 0;
    re = gsl_sf_exp_mult_e10_e(1.0, 1.0); console.log(re);
    if (test_sf_frac_diff(re.val, M_E ) > TEST_TOL0)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.err > TEST_TOL2)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.e10 != 0)
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test(sa, "gsl_sf_exp_mult_e10_e(1.0, 1.0, re)");
    TEST_SF_E10_DD(s, gsl_sf_exp_mult_e10_e, { x: 1.0, y: 1.0 }, M_E, 0, TEST_TOL0, "gsl_sf_exp_mult_e10_e");

    sa.Integer = 0;
    re = gsl_sf_exp_mult_e10_e(1000.0, 1.0e+200); console.log(re);
    if (test_sf_frac_diff(re.val, 1.970071114017046993888879352) > TEST_TOL3)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.err > 1.0e-11)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.e10 != 634)
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test(sa, "gsl_sf_exp_mult_e10_e(1000.0, 1.0e+200, re)");
    TEST_SF_E10_DD(s, gsl_sf_exp_mult_e10_e, { x: 1000.0, y: 1.0e+200 }, 1.970071114017046993888879352, 634, TEST_TOL3, "gsl_sf_exp_mult_e10_e");

    console.log("    ... gsl_sf_exp_mult_err_e10_e");
    sa.Integer = 0;
    re = gsl_sf_exp_mult_err_e10_e(1.0, TEST_TOL0, 1.0, TEST_TOL0); console.log(re);
    if (test_sf_frac_diff(re.val, M_E ) > TEST_TOL0)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.err > TEST_TOL2)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.e10 != 0)
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test(sa, "gsl_sf_exp_mult_err_e10_e(1.0, TEST_TOL0, 1.0, TEST_TOL0, re)");
    TEST_SF_E10_4D(s, gsl_sf_exp_mult_err_e10_e, { x: 1.0, dx: TEST_TOL0, y: 1.0, dy: TEST_TOL0 }, M_E, 0, TEST_TOL0, "gsl_sf_exp_mult_err_e10_e");

    sa.Integer = 0;
    re = gsl_sf_exp_mult_err_e10_e(1000.0, 1.0e-12, 1.0e+200, 1.0e+190); console.log(re);
    if (test_sf_frac_diff(re.val, 1.9700711140165661 ) > TEST_TOL3)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.err > 1.0e-09)
    {
        sa.Integer = sa.Integer + 1;
    }
    if (re.e10 != 634)
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test(sa, "gsl_sf_exp_mult_err_e10_e(1.0e-12, 1.0e+200, re)");
    TEST_SF_E10_4D(s, gsl_sf_exp_mult_err_e10_e, { x: 1000.0, dx: 1.0e-12, y: 1.0e+200, dy: 1.0e+190 }, 1.9700711140165661,634, TEST_TOL3, "gsl_sf_exp_mult_err_e10_e");

    // Test cases from Szymon Jaroszewicz
    console.log("*** Test cases from Szymon Jaroszewicz ***");

    console.log("    ... gsl_sf_exp_mult_e10_e");
    TEST_SF_E10_DD(s, gsl_sf_exp_mult_e10_e, { x: 10000.0, y: 1.0 }, 8.806818225662921587261496007, 4342, TEST_TOL5, "gsl_sf_exp_mult_e10_e");
    TEST_SF_E10_DD(s, gsl_sf_exp_mult_e10_e, { x: 100.0,   y: 1.0 }, 2.688117141816135448412625551e43, 0, TEST_TOL2, "gsl_sf_exp_mult_e10_e");
    console.log("    ... gsl_sf_exp_e10_e");
    TEST_SF_E10_D(s, gsl_sf_exp_e10_e, 100.0, 2.688117141816135448412625551e43, 0, TEST_TOL2, "gsl_sf_exp_e10_e");
    TEST_SF_E10_D(s, gsl_sf_exp_e10_e, 1000.0, 1.970071114017046993888879352, 434, TEST_TOL3, "gsl_sf_exp_e10_e");
    TEST_SF_E10_D(s, gsl_sf_exp_e10_e, -100.0, 3.720075976020835962959695803e-44, 0, TEST_TOL2, "gsl_sf_exp_e10_e");
    TEST_SF_E10_D(s, gsl_sf_exp_e10_e, -1000.0, 5.075958897549456765291809479, -435, TEST_TOL3, "gsl_sf_exp_e10_e");

    console.log("    ... gsl_sf_expm1_e");
    TEST_SF_D(s, gsl_sf_expm1_e, -10.0, Math.exp(-10.0)-1.0, TEST_TOL0, "gsl_sf_expm1_e");
    TEST_SF_D(s, gsl_sf_expm1_e, -0.001, -0.00099950016662500845, TEST_TOL0, "gsl_sf_expm1_e");
    TEST_SF_D(s, gsl_sf_expm1_e, -1.0e-8, -1.0e-08 + 0.5e-16, TEST_TOL0, "gsl_sf_expm1_e");
    TEST_SF_D(s, gsl_sf_expm1_e,  1.0e-8, 1.0e-08 + 0.5e-16, TEST_TOL0, "gsl_sf_expm1_e");
    TEST_SF_D(s, gsl_sf_expm1_e,  0.001, 0.0010005001667083417, TEST_TOL0, "gsl_sf_expm1_e");
    TEST_SF_D(s, gsl_sf_expm1_e,  10.0, Math.exp(10.0)-1.0, TEST_TOL0, "gsl_sf_expm1_e");

    console.log("    ... gsl_sf_exprel_e");
    TEST_SF_D(s, gsl_sf_exprel_e, -10.0, 0.0999954600070237515, TEST_TOL0, "gsl_sf_exprel_e");
    TEST_SF_D(s, gsl_sf_exprel_e, -0.001, 0.9995001666250084, TEST_TOL0, "gsl_sf_exprel_e");
    TEST_SF_D(s, gsl_sf_exprel_e, -1.0e-8, 1.0 - 0.5e-08, TEST_TOL0, "gsl_sf_exprel_e");
    TEST_SF_D(s, gsl_sf_exprel_e,  1.0e-8, 1.0 + 0.5e-08, TEST_TOL0, "gsl_sf_exprel_e");
    TEST_SF_D(s, gsl_sf_exprel_e,  0.001, 1.0005001667083417, TEST_TOL0, "gsl_sf_exprel_e");
    TEST_SF_D(s, gsl_sf_exprel_e,  10.0, 2202.5465794806716517, TEST_TOL0, "gsl_sf_exprel_e");

    console.log("    ... gsl_sf_exprel_2_e");
    TEST_SF_D(s, gsl_sf_exprel_2_e, -10.0, 0.18000090799859524970, TEST_TOL0, "gsl_sf_exprel_2_e");
    TEST_SF_D(s, gsl_sf_exprel_2_e, -0.001, 0.9996667499833361107, TEST_TOL0, "gsl_sf_exprel_2_e");
    TEST_SF_D(s, gsl_sf_exprel_2_e, -1.0e-8, 0.9999999966666666750, TEST_TOL0, "gsl_sf_exprel_2_e");
    TEST_SF_D(s, gsl_sf_exprel_2_e,  1.0e-8, 1.0000000033333333417, TEST_TOL0, "gsl_sf_exprel_2_e");
    TEST_SF_D(s, gsl_sf_exprel_2_e,  0.001, 1.0003334166833361115, TEST_TOL0, "gsl_sf_exprel_2_e");
    TEST_SF_D(s, gsl_sf_exprel_2_e,  10.0, 440.3093158961343303, TEST_TOL0, "gsl_sf_exprel_2_e");

    console.log("    ... gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x: -1000.0 }, 0.00299400600000000000, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x: -100.0 }, 0.02940600000000000000, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x: -10.0}, 0.24599972760042142509, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x: -3.0 }, 0.5444917625849191238, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x: -0.001 }, 0.9997500499916678570, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x: -1.0e-8 }, 0.9999999975000000050, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x:  1.0e-8}, 1.0000000025000000050, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x:  0.001 }, 1.0002500500083345240, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x:  3.0 }, 2.5745637607083706091, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x:  3.1 }, 2.6772417068460206247, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x:  10.0 }, 131.79279476884029910, TEST_TOL1, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 3, x:  100.0 }, 1.6128702850896812690e+38, TEST_TOL2, "gsl_sf_exprel_n_e");

    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x: -1000.0 }, 0.04766231609253975959, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x: -100.0 }, 0.3348247572345889317, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x: -10.0 }, 0.8356287051853286482, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x: -3.0 }, 0.9443881609152163615, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x: -1.0 }, 0.980762245565660617, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x: -1.0e-8 }, 1.0 -1.0e-8/51.0, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x:  1.0e-8 }, 1.0 +1.0e-8/51.0, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x:  1.0 }, 1.01999216583666790, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x:  3.0 }, 1.0624205757460368307, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x:  48.0 }, 7.499573876877194416, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x:  50.1 }, 9.311803306230992272, TEST_TOL4, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x:  100.0 }, 8.175664432485807634e+07, TEST_TOL4, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 50, x:  500.0 }, 4.806352370663185330e+146, TEST_TOL3, "gsl_sf_exprel_n_e");

    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x: -1000.0 }, 0.3334815803127619256, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x: -100.0 }, 0.8335646217536183909, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x: -10.0 }, 0.9804297803131823066, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x: -3.0 }, 0.9940475488850672997, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x: -1.0 }, 0.9980079602383488808, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x: -1.0e-8 }, 1.0 -1.0e-8/501.0, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x:  1.0e-8 }, 1.0 +1.0e-8/501.0, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x:  1.0 }, 1.0019999920160634252, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x:  3.0 }, 1.0060240236632444934, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x:  48.0 }, 1.1059355517981272174, TEST_TOL0, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x:  100.0 }, 1.2492221464878287204, TEST_TOL1, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x:  500.0 }, 28.363019877927630858, TEST_TOL2, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x:  1000.0 }, 2.4037563160335300322e+68, TEST_TOL4, "gsl_sf_exprel_n_e");
    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 500, x:  1600.0 }, 7.899293535320607403e+226, TEST_TOL4, "gsl_sf_exprel_n_e");

    TEST_SF_ID(s, gsl_sf_exprel_n_e, { i: 1263131, x: 1261282.3637 }, 545.0113107238425900305428360, TEST_TOL4, "gsl_sf_exprel_n_e");

    console.log("    ... gsl_sf_exprel_n_CF_e");
    TEST_SF_DD(s, gsl_sf_exprel_n_CF_e, { x: 6.315655e+05, y: 6.302583168053568806e+05 }, 385.425369029433473098652465720, TEST_TOL4, "gsl_sf_exprel_n_CF_e");

    return s;

} // test_exp

// ----------------------------------------------------------------------------

export function test_expint( )
{
    var s = { Integer: 0 };

    console.log("Test Exponential Integral Functions ...");

    console.log("    ... gsl_sf_expint_E1_e");
    TEST_SF_D(s, gsl_sf_expint_E1_e, -1.0, -1.8951178163559367555, TEST_TOL0, "gsl_sf_expint_E1_e");
    TEST_SF_D(s, gsl_sf_expint_E1_e, 1.0e-10, 22.448635265138923980, TEST_TOL0, "gsl_sf_expint_E1_e");
    TEST_SF_D(s, gsl_sf_expint_E1_e, 1.0e-05, 10.935719800043695615, TEST_TOL0, "gsl_sf_expint_E1_e");
    TEST_SF_D(s, gsl_sf_expint_E1_e, 0.1, 1.82292395841939066610, TEST_TOL0, "gsl_sf_expint_E1_e");
    TEST_SF_D(s, gsl_sf_expint_E1_e, 1.0, 0.21938393439552027368, TEST_TOL0, "gsl_sf_expint_E1_e");
    TEST_SF_D(s, gsl_sf_expint_E1_e, 10.0, 4.156968929685324277e-06, TEST_TOL1, "gsl_sf_expint_E1_e");
    TEST_SF_D(s, gsl_sf_expint_E1_e, 50.0, 3.783264029550459019e-24, TEST_TOL2, "gsl_sf_expint_E1_e");
    TEST_SF_D(s, gsl_sf_expint_E1_e, 300.0, 1.710384276804510115e-133, TEST_TOL2, "gsl_sf_expint_E1_e");

    console.log("    ... gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, -1.0, 0.8231640121031084799, TEST_TOL1, "gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, 0.0, 1.0, TEST_TOL0, "gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, 1.0/4294967296.0, 0.9999999947372139168, TEST_TOL0, "gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, 1.0/65536.0, 0.9998243233207178845, TEST_TOL0, "gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, 0.1, 0.7225450221940205066, TEST_TOL0, "gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, 1.0, 0.14849550677592204792, TEST_TOL0, "gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, 10.0, 3.830240465631608762e-06, TEST_TOL1, "gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, 50.0, 3.711783318868827367e-24, TEST_TOL2, "gsl_sf_expint_E2_e");
    TEST_SF_D(s, gsl_sf_expint_E2_e, 300.0, 1.7047391998483433998e-133, TEST_TOL2, "gsl_sf_expint_E2_e");

    // Tests for E_n(x)

    console.log("    ... gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 1, x: -1.0 }, -1.8951178163559367555, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 1, x: 1.0e-10 }, 22.448635265138923980, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 1, x: 1.0e-05 }, 10.935719800043695615, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 1, x: 0.1 }, 1.82292395841939066610, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 1, x: 1.0 }, 0.21938393439552027368, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 1, x: 10.0 }, 4.156968929685324277e-06, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 1, x: 50.0 }, 3.783264029550459019e-24, TEST_TOL2, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 1, x: 300.0 }, 1.710384276804510115e-133, TEST_TOL2, "gsl_sf_expint_En_e");

    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: -1.0 }, 0.8231640121031084799, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: 0.0 }, 1.0, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: 1.0/4294967296.0 }, 0.9999999947372139168, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: 1.0/65536.0 }, 0.9998243233207178845, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: 0.1 }, 0.7225450221940205066, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: 1.0 }, 0.14849550677592204792, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: 10.0 }, 3.830240465631608762e-06, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: 50.0 }, 3.711783318868827367e-24, TEST_TOL2, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 2, x: 300.0 }, 1.7047391998483433998e-133, TEST_TOL2, "gsl_sf_expint_En_e");

    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 3, x: 0.0 }, 0.5, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 3, x: 1.0/4294967296.0 }, 0.499999999767169356972, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 3, x: 1.0/65536.0 }, 0.4999847426094515610, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 3, x: 0.1 }, 0.4162914579082787612543, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 3, x: 1.0 }, 0.10969196719776013683858, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 3, x: 10.0 },0.000003548762553084381959981, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 3, x: 50.0 }, 3.6429094264752049812e-24, TEST_TOL2, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 3, x: 300.0 },1.699131143349179084e-133, TEST_TOL2, "gsl_sf_expint_En_e");

    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 10, x: 0.0 }, 0.111111111111111111, TEST_TOL0, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 10, x: 1.0/4294967296.0 }, 0.111111111082007280658, TEST_TOL2, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 10, x: 1.0/65536.0 }, 0.11110920377910896018606, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 10, x: 0.1 }, 0.099298432000896813567905, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 10, x: 1.0 }, 0.036393994031416401634164534, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 10, x: 10.0 }, 0.00000232530265702821081778968, TEST_TOL1, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 10, x: 50.0 }, 3.223296586749110919572e-24, TEST_TOL2, "gsl_sf_expint_En_e");
    TEST_SF_ID(s, gsl_sf_expint_En_e, { i: 10, x: 300.0 }, 1.6608815083360041367294736e-133, TEST_TOL2, "gsl_sf_expint_En_e");

    // Tests for Ei(x)

    console.log("    ... gsl_sf_expint_Ei_e");
    TEST_SF_D(s, gsl_sf_expint_Ei_e, -1.0, -0.21938393439552027368, TEST_TOL0, "gsl_sf_expint_Ei_e");
    TEST_SF_D(s, gsl_sf_expint_Ei_e, 1.0/4294967296.0, -21.603494112783886397, TEST_TOL0, "gsl_sf_expint_Ei_e");
    TEST_SF_D(s, gsl_sf_expint_Ei_e, 1.0, 1.8951178163559367555, TEST_TOL0, "gsl_sf_expint_Ei_e");

    console.log("    ... gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, -10000.0, -0.00010001000200060024012, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, -1000.0, -0.0010010020060241207251, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, -10.0, -0.11314702047341077803, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, -1.0, -0.69717488323506606877, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 1.0e-10, 22.448635267383787506, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 1.0e-05, 10.935829157788483865, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 0.1, 2.0146425447084516791, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 1.0, 0.59634736232319407434, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 10.0, 0.091563333939788081876, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 50.0, 0.019615109930114870365, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 300.0, 0.0033222955652707070644, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 1000.0, 0.00099900199402388071500, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E1_scaled_e, 10000.0, 0.000099990001999400239880, TEST_TOL0, "gsl_sf_expint_E1_scaled_e");

    console.log("    ... gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, -10000.0, -0.00010002000600240120072, TEST_TOL3, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, -1000.0, -0.0010020060241207250807, TEST_TOL3, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, -10.0, -0.13147020473410778034, TEST_TOL1, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, -1.0, 0.30282511676493393123, TEST_TOL1, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 0.0, 1.0, TEST_TOL1, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 1.0/4294967296.0, 0.99999999497004455927, TEST_TOL0, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 1.0/65536.0, 0.99983957954556245453, TEST_TOL0, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 0.1, 0.79853574552915483209, TEST_TOL0, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 1.0, 0.40365263767680592566, TEST_TOL0, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 10.0, 0.084366660602119181239, TEST_TOL1, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 50.0, 0.019244503494256481735, TEST_TOL2, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 300.0, 0.0033113304187878806691, TEST_TOL0, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 1000.0, 0.00099800597611928500004, TEST_TOL0, "gsl_sf_expint_E2_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_E2_scaled_e, 10000.0, 0.000099980005997601199281, TEST_TOL0, "gsl_sf_expint_E2_scaled_e");

    // Tests for E_n(x)

    console.log("    ... gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: -10000.0 }, -0.00010001000200060024012, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: -1000.0 }, -0.0010010020060241207251, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: -10.0 }, -0.11314702047341077803, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: -1.0 }, -0.69717488323506606877, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 1.0e-10 }, 22.448635267383787506, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 1.0e-05 }, 10.935829157788483865, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 0.1 }, 2.0146425447084516791, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 1.0 }, 0.59634736232319407434, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 10.0 }, 0.091563333939788081876, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 50.0 }, 0.019615109930114870365, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 300.0 }, 0.0033222955652707070644, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 1000.0 }, 0.00099900199402388071500, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 1, x: 10000.0 }, 0.000099990001999400239880, TEST_TOL0, "gsl_sf_expint_En_scaled_e");

    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: -10000.0 }, -0.00010002000600240120072, TEST_TOL3, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: -1000.0 }, -0.0010020060241207250807, TEST_TOL3, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: -10.0 }, -0.13147020473410778034, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: -1.0 }, 0.30282511676493393123, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 0.0 }, 1.0, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 1.0/4294967296.0 }, 0.99999999497004455927, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 1.0/65536.0 }, 0.99983957954556245453, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 0.1 }, 0.79853574552915483209, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 1.0 }, 0.40365263767680592566, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 10.0 }, 0.084366660602119181239, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 50.0 }, 0.019244503494256481735, TEST_TOL2, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 300.0 }, 0.0033113304187878806691, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 1000.0 }, 0.00099800597611928500004, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 2, x: 10000.0 }, 0.000099980005997601199281, TEST_TOL0, "gsl_sf_expint_En_scaled_e");

    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 3, x: 0.0 }, 0.5, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 3, x: 1.0/4294967296.0 }, 0.4999999998835846787586, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 3, x: 1.0/65536.0 }, 0.4999923718293796877864492, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 3, x: 0.1 }, 0.4600732127235422583955, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 3, x: 1.0 }, 0.298173681161597037170539, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 3, x: 10.0 }, 0.07816669698940409380349, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 3, x: 50.0 }, 0.0188874126435879566345, TEST_TOL2, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 3, x: 300.0 }, 0.00330043718181789963028657675, TEST_TOL2, "gsl_sf_expint_En_scaled_e");

    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 10, x: 0.0 }, 0.111111111111111111, TEST_TOL0, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 10, x: 1.0/4294967296.0 }, 0.11111111110787735217158, TEST_TOL2, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 10, x: 1.0/65536.0 }, 0.1111108991839472074435, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 10, x: 0.1 }, 0.1097417392579033988025, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 10, x: 1.0 }, 0.09892913264064615521915, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 10, x: 10.0 }, 0.0512181994376050593314159875, TEST_TOL1, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 10, x: 50.0 }, 0.0167118436335939556034579, TEST_TOL2, "gsl_sf_expint_En_scaled_e");
    TEST_SF_ID(s, gsl_sf_expint_En_scaled_e, { i: 10, x: 300.0 }, 0.0032261400811599644878615, TEST_TOL2, "gsl_sf_expint_En_scaled_e");

    console.log("    ... gsl_sf_expint_Ei_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_Ei_scaled_e, -1000.0, -0.00099900199402388071500, TEST_TOL0, "gsl_sf_expint_Ei_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_Ei_scaled_e, -1.0, -0.59634736232319407434, TEST_TOL0, "gsl_sf_expint_Ei_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_Ei_scaled_e, 1.0/4294967296.0, -21.603494107753930958, TEST_TOL0, "gsl_sf_expint_Ei_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_Ei_scaled_e, 1.0, 0.69717488323506606877, TEST_TOL0, "gsl_sf_expint_Ei_scaled_e");
    TEST_SF_D(s, gsl_sf_expint_Ei_scaled_e, 1000.0, 0.0010010020060241207251, TEST_TOL0, "gsl_sf_expint_Ei_scaled_e");

    console.log("    ... gsl_sf_Shi_e");
    TEST_SF_D(s,  gsl_sf_Shi_e, -1.0, -1.0572508753757285146, TEST_TOL0, "gsl_sf_Shi_e");
    TEST_SF_D(s,  gsl_sf_Shi_e, 1.0/4294967296.0, 2.3283064365386962891e-10, TEST_TOL0, "gsl_sf_Shi_e");
    TEST_SF_D(s,  gsl_sf_Shi_e, 1.0/65536.0, 0.00001525878906269737298, TEST_TOL0, "gsl_sf_Shi_e");
    TEST_SF_D(s,  gsl_sf_Shi_e, 0.1, 0.1000555722250569955, TEST_TOL0, "gsl_sf_Shi_e");
    TEST_SF_D(s,  gsl_sf_Shi_e, 1.0, 1.0572508753757285146, TEST_TOL0, "gsl_sf_Shi_e");
    TEST_SF_D(s,  gsl_sf_Shi_e, 10.0, 1246.1144901994233444, TEST_TOL1, "gsl_sf_Shi_e");
    TEST_SF_D(s,  gsl_sf_Shi_e, 50.0, 5.292818448565845482e+19, TEST_TOL2, "gsl_sf_Shi_e");
    TEST_SF_D(s,  gsl_sf_Shi_e, 300.0, 3.248241254044332895e+127, TEST_TOL2, "gsl_sf_Shi_e");

    console.log("    ... gsl_sf_Chi_e");
    TEST_SF_D(s,  gsl_sf_Chi_e, -1.0, 0.8378669409802082409, TEST_TOL0, "gsl_sf_Chi_e");
    TEST_SF_D(s,  gsl_sf_Chi_e, 1.0/4294967296.0, -21.603494113016717041, TEST_TOL0, "gsl_sf_Chi_e");
    TEST_SF_D(s,  gsl_sf_Chi_e, 1.0/65536.0, -10.513139223999384429, TEST_TOL0, "gsl_sf_Chi_e");
    TEST_SF_D(s,  gsl_sf_Chi_e, 1.0/8.0, -1.4983170827635760646, TEST_TOL0, "gsl_sf_Chi_e");
    TEST_SF_D(s,  gsl_sf_Chi_e, 1.0, 0.8378669409802082409, TEST_TOL0, "gsl_sf_Chi_e");
    TEST_SF_D(s,  gsl_sf_Chi_e, 10.0, 1246.1144860424544147, TEST_TOL1, "gsl_sf_Chi_e");
    TEST_SF_D(s,  gsl_sf_Chi_e, 50.0, 5.292818448565845482e+19, TEST_TOL2, "gsl_sf_Chi_e");
    TEST_SF_D(s,  gsl_sf_Chi_e, 300.0, 3.248241254044332895e+127, TEST_TOL2, "gsl_sf_Chi_e");

    console.log("    ... gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 1.0e-10, 1.0e-10, TEST_TOL0, "gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 1.0e-05, 9.9999999999999975e-06, TEST_TOL0, "gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 0.1, 0.09997500714119079665122, TEST_TOL0, "gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 0.5, 0.48491714311363971332427, TEST_TOL0, "gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 1.0, 0.80751118213967145285833, TEST_TOL0, "gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 2.0, 0.89295351429387631138208, TEST_TOL0, "gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 5.0, 0.89297951156924921121856, TEST_TOL0, "gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 10.0, 0.89297951156924921121856, TEST_TOL0, "gsl_sf_expint_3_e");
    TEST_SF_D(s, gsl_sf_expint_3_e, 100.0, 0.89297951156924921121856, TEST_TOL0, "gsl_sf_expint_3_e");

    console.log("    ... gsl_sf_Si_e");
    TEST_SF_D(s, gsl_sf_Si_e, -1.0, -0.9460830703671830149, TEST_TOL0, "gsl_sf_Si_e");
    TEST_SF_D(s, gsl_sf_Si_e, 1.0e-10, 1.0e-10, TEST_TOL0, "gsl_sf_Si_e");
    TEST_SF_D(s, gsl_sf_Si_e, 1.0e-05, 9.999999999944444444e-06, TEST_TOL0, "gsl_sf_Si_e");
    TEST_SF_D(s, gsl_sf_Si_e, 0.1, 0.09994446110827695016, TEST_TOL0, "gsl_sf_Si_e");
    TEST_SF_D(s, gsl_sf_Si_e, 1.0, 0.9460830703671830149, TEST_TOL0, "gsl_sf_Si_e");
    TEST_SF_D(s, gsl_sf_Si_e, 10.0, 1.6583475942188740493, TEST_TOL0, "gsl_sf_Si_e");
    TEST_SF_D(s, gsl_sf_Si_e, 50.0, 1.5516170724859358947, TEST_TOL0, "gsl_sf_Si_e");
    TEST_SF_D(s, gsl_sf_Si_e, 300.0, 1.5708810882137495193, TEST_TOL0, "gsl_sf_Si_e");
    // *****
    TEST_SF_D(s, gsl_sf_Si_e, 1.0e+20, 1.5707963267948966192, TEST_TOL0, "gsl_sf_Si_e");
    // *****

    console.log("    ... gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 1.0/4294967296.0, -21.603494113016717041, TEST_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 1.0/65536.0, -10.513139224115799751, TEST_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 1.0/8.0, -1.5061295845296396649, TEST_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 1.0, 0.3374039229009681347, TEST_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 10.0, -0.04545643300445537263, TEST_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 50.0, -0.005628386324116305440, TEST_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 300.0, -0.003332199918592111780, TEST_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 65536.0, 0.000010560248837656279453, TEST_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 4294967296.0, -1.0756463261957757485e-10, TEST_SQRT_TOL0, "gsl_sf_Ci_e");
    TEST_SF_D(s, gsl_sf_Ci_e, 1099511627776.0, -3.689865584710764214e-13, 1024.0*TEST_SQRT_TOL0, "gsl_sf_Ci_e");

    console.log("    ... gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 1.0e-10, 1.0e-10, TEST_TOL0, "gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 1.0e-05, 9.99999999988888888889e-06, TEST_TOL0, "gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 0.1, 0.09988928686033618404, TEST_TOL0, "gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 1.0, 0.91596559417721901505, TEST_TOL0, "gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 2.0, 1.57601540344632342236, TEST_TOL0, "gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 10.0, 3.71678149306806859029, TEST_TOL0, "gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 50.0, 6.16499047850274874222, TEST_TOL0, "gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 300.0, 8.96281388924518959990, TEST_TOL0, "gsl_sf_atanint_e");
    TEST_SF_D(s, gsl_sf_atanint_e, 1.0e+5, 18.084471031038661920, TEST_TOL0, "gsl_sf_atanint_e");

    return s;

} // test_expint

// ----------------------------------------------------------------------------

export function test_fermidirac( )
{
    var s = { Integer: 0 };

    console.log("Test Fermi-Dirac Functions ...");

    console.log("    ... gsl_sf_fermi_dirac_m1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_m1_e, -10.0, 0.00004539786870243439450, TEST_TOL0, "gsl_sf_fermi_dirac_m1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_m1_e,  -1.0, 0.26894142136999512075, TEST_TOL0, "gsl_sf_fermi_dirac_m1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_m1_e,   1.0, 0.7310585786300048793, TEST_TOL0, "gsl_sf_fermi_dirac_m1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_m1_e,  10.0, 0.9999546021312975656, TEST_TOL0, "gsl_sf_fermi_dirac_m1_e");

    console.log("    ... gsl_sf_fermi_dirac_0_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_0_e, -10.0, 0.00004539889921686464677, TEST_TOL0, "gsl_sf_fermi_dirac_0_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_0_e,  -1.0, 0.31326168751822283405, TEST_TOL0, "gsl_sf_fermi_dirac_0_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_0_e,   1.0, 1.3132616875182228340, TEST_TOL0, "gsl_sf_fermi_dirac_0_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_0_e,  10.0, 10.000045398899216865, TEST_TOL0, "gsl_sf_fermi_dirac_0_e");

    console.log("    ... gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e, -10.0, 0.00004539941448447633524, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,  -2.0, 0.13101248471442377127, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,  -1.0, 0.3386479964034521798, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,  -0.4, 0.5825520806897909028, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,   0.4, 1.1423819861584355337, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,   1.0, 1.8062860704447742567, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,   1.5, 2.5581520872227806402, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,   2.5, 4.689474797599761667, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,  10.0, 51.64488866743374196, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,  12.0, 73.64492792264531092, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,  20.0, 201.64493406478707282, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_1_e,  50.0, 1251.6449340668482264, TEST_TOL0, "gsl_sf_fermi_dirac_1_e");

    console.log("    ... gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e, -10.0, 0.00004539967212174776662, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,  -2.0, 0.13313272938565030508, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,  -1.0, 0.3525648792978077590, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,  -0.4, 0.6229402647001272120, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,   0.4, 1.2915805581060844533, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,   1.0, 2.1641656128127008622, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,   1.5, 3.247184513920792475, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,   2.5, 6.797764392735056317, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,  10.0, 183.11605273482105278, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,  12.0, 307.73921494638635166, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,  20.0, 1366.2320146723590157, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e,  50.0, 20915.580036675744655, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_2_e, 200.0, 1.3336623201467029786e+06, TEST_TOL0, "gsl_sf_fermi_dirac_2_e");

    console.log("    ... gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e, -10.0, 0.00004539847236080549532, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,  -2.0, 0.12366562180120994266, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,  -1.0, 0.29402761761145122022, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,  -0.4, 0.4631755336886027800, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,   0.4, 0.7654084737661656915, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,   1.0, 1.0270571254743506890, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,   1.5, 1.2493233478527122008, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,   2.5, 1.6663128834358313625, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,  10.0, 3.552779239536617160, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,  12.0, 3.897268231925439359, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,  20.0, 5.041018507535328603, TEST_TOL0, "gsl_sf_fermi_dirac_mhalf_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_mhalf_e,  50.0, 7.977530858581869960, TEST_TOL1, "gsl_sf_fermi_dirac_mhalf_e");

    console.log("    ... gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e, -10.0, 0.00004539920105264132755, TEST_TOL1, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,  -2.0, 0.12929851332007559106, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,  -1.0, 0.3277951592607115477, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,  -0.4, 0.5522452153690688947, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,   0.4, 1.0386797503389389277, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,   1.0, 1.5756407761513002308, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,   1.5, 2.1448608775831140360, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,   2.5, 3.606975377950373251, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,  10.0, 24.084656964637653615, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,  12.0, 31.540203287044242593, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,  20.0, 67.49151222165892049, TEST_TOL0, "gsl_sf_fermi_dirac_half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_half_e,  50.0, 266.09281252136259343, TEST_TOL1, "gsl_sf_fermi_dirac_half_e");

    console.log("    ... gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e, -10.0, 0.00004539956540456176333, TEST_TOL0, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,  -2.0, 0.13224678225177236685, TEST_TOL0, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,  -1.0, 0.3466747947990574170, TEST_TOL0, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,  -0.4, 0.6056120213305040910, TEST_TOL0, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,   0.4, 1.2258236403963668282, TEST_TOL0, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,   1.0, 2.0022581487784644573, TEST_TOL0, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,   1.5, 2.9277494127932173068, TEST_TOL0, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,   2.5, 5.768879312210516582, TEST_TOL0, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,  10.0, 101.00510084332600020, TEST_TOL2, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,  12.0, 156.51518642795728036, TEST_TOL1, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,  20.0, 546.5630100657601959, TEST_TOL1, "gsl_sf_fermi_dirac_3half_e");
    TEST_SF_D(s, gsl_sf_fermi_dirac_3half_e,  50.0, 5332.353566687145552, TEST_TOL1, "gsl_sf_fermi_dirac_3half_e");

    console.log("    ... gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 3, x:  -2.0 }, 0.1342199155038680215, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 3, x:   0.0 }, 0.9470328294972459176, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 3, x:   0.1 }, 1.0414170610956165759, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 3, x:   1.0 }, 2.3982260822489407070, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 3, x:   3.0 }, 12.621635313399690724, TEST_TOL1, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 3, x: 100.0 }, 4.174893231066566793e+06, TEST_TOL1, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 3, x: 500.0 }, 2.604372285319088354e+09, TEST_TOL1, "gsl_sf_fermi_dirac_int_e");

    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 5, x:  -2.0 }, 0.13505242246823676478, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 5, x:   0.0 }, 0.9855510912974351041, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 5, x:   0.1 }, 1.0876519750101492782, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 5, x:   1.0 }, 2.6222337848692390539, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 5, x:   3.0 }, 17.008801618012113022, TEST_TOL1, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 5, x: 100.0 }, 1.3957522531334869874e+09, TEST_TOL1, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 5, x: 500.0 }, 2.1705672808114817955e+13, TEST_TOL2, "gsl_sf_fermi_dirac_int_e");

    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 7, x:  -2.0 }, 0.1352641105671255851, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 7, x:   0.0 }, 0.9962330018526478992, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 7, x:   0.1 }, 1.1005861815180315485, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 7, x:   1.0 }, 2.6918878172003129203, TEST_TOL0, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 7, x:   3.0 }, 19.033338976999367642, TEST_TOL2, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 7, x:  10.0 }, 5654.530932873610014, TEST_TOL1, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 7, x:  50.0 }, 1.005005069985066278e+09, TEST_TOL2, "gsl_sf_fermi_dirac_int_e");
    TEST_SF_ID(s, gsl_sf_fermi_dirac_int_e, { i: 7, x: 500.0 }, 9.691690268341569514e+16, TEST_TOL3, "gsl_sf_fermi_dirac_int_e");

    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 9, x:  -2.0 }, 0.1353174385330242691, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 9, x:   0.0 }, 0.9990395075982715656, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 9, x:   0.1 }, 1.1039997234712941212, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 9, x:   1.0 }, 2.7113648898129249947, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 9, x:   3.0 }, 19.768544008138602223, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 9, x:  10.0 }, 10388.990167312912478, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 9, x:  50.0 }, 2.85466960802601649e+10, TEST_TOL1, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 9, x: 500.0 }, 2.69273849842695876e+20, 2.0*TEST_TOL1, "gsl_sf_fermi_dirac_int_e" );

    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 10, x:  -2.0 }, 0.13532635396712288092, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 10, x:   0.0 }, 0.9995171434980607541, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 10, x:   0.1 }, 1.1045818238852612296, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 10, x:   1.0 }, 2.7147765350346120647, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 10, x:   3.0 }, 19.917151938411675171, TEST_TOL1, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 10, x:  10.0 }, 12790.918595516495955, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 10, x:  50.0 }, 1.3147703201869657654e+11, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 10, x: 500.0 }, 1.2241331244469204398e+22, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );

    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 11, x:  -2.0 }, 0.1353308162894847149, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 11, x:   0.0 }, 0.9997576851438581909, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 11, x:   0.1 }, 1.1048751811565850418, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 11, x:   1.0 }, 2.7165128749007313436, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 11, x:   3.0 }, 19.997483022044603065, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 11, x:  10.0 }, 14987.996005901818036, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 11, x:  50.0 }, 5.558322924078990628e+11, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 11, x: 500.0 }, 5.101293089606198280e+23, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );

    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 20, x:  -2.0 }, 0.13533527450327238373, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 20, x:   0.0 }, 0.9999995232582155428, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 20, x:   0.1 }, 1.1051703357941368203, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 20, x:   1.0 }, 2.7182783069905721654, TEST_TOL0, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 20, x:   3.0 }, 20.085345296028242734, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 20, x:  10.0 }, 21898.072920149606475, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 20, x:  50.0 }, 1.236873256595717618e+16, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );
    TEST_SF_ID( s, gsl_sf_fermi_dirac_int_e, { i: 20, x: 500.0 }, 9.358938204369557277e+36, TEST_TOL2, "gsl_sf_fermi_dirac_int_e" );

    return s;

} // test_fermidirac

// ----------------------------------------------------------------------------

export function test_gegen( )
{
    var s    = { Integer: 0 };
    var sa   = { Integer: 0 };
    var ga   = [];

    console.log( "Test Gegenbauer Functions ..." );

    console.log( "    ... gsl_sf_gegenpoly_1_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_1_e, { x:-0.2,   y: 1.0 }, -0.4,  TEST_TOL0, "gsl_sf_gegenpoly_1_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_1_e, { x: 0.0,   y: 1.0 }, 2.0,   TEST_TOL0, "gsl_sf_gegenpoly_1_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_1_e, { x: 1.0,   y: 1.0 }, 2.0,   TEST_TOL0, "gsl_sf_gegenpoly_1_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_1_e, { x: 1.0,   y: 0.5 }, 1.0,   TEST_TOL0, "gsl_sf_gegenpoly_1_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_1_e, { x: 5.0,   y: 1.0 }, 10.0,  TEST_TOL0, "gsl_sf_gegenpoly_1_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_1_e, { x: 100.0, y: 0.5 }, 100.0, TEST_TOL0, "gsl_sf_gegenpoly_1_e" );

    console.log( "    ... gs/l_sf_gegenpoly_2_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_2_e, { x:-0.2,   y: 0.5 }, 0.12,   TEST_TOL0, "gsl_sf_gegenpoly_2_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_2_e, { x: 0.0,   y: 1.0 }, 1.00,   TEST_TOL0, "gsl_sf_gegenpoly_2_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_2_e, { x: 1.0,   y: 1.0 }, 3.00,   TEST_TOL0, "gsl_sf_gegenpoly_2_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_2_e, { x: 1.0,   y: 0.1 }, -0.96,  TEST_TOL0, "gsl_sf_gegenpoly_2_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_2_e, { x: 5.0,   y: 1.0 }, 55.0,   TEST_TOL0, "gsl_sf_gegenpoly_2_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_2_e, { x: 100.0, y: 0.5 }, 4950.0, TEST_TOL0, "gsl_sf_gegenpoly_2_e" );

    console.log( "    ... gsl_sf_gegenpoly_3_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_3_e, { x:-0.2,   y: 0.5 }, 0.112,      TEST_TOL0, "gsl_sf_gegenpoly_3_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_3_e, { x: 0.0,   y: 1.0 }, -2.0/3.0,   TEST_TOL0, "gsl_sf_gegenpoly_3_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_3_e, { x: 1.0,   y: 1.0 }, 4.000,      TEST_TOL0, "gsl_sf_gegenpoly_3_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_3_e, { x: 1.0,   y: 0.1 }, -0.392,     TEST_TOL0, "gsl_sf_gegenpoly_3_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_3_e, { x: 5.0,   y: 1.0 }, 220.000,    TEST_TOL0, "gsl_sf_gegenpoly_3_e" );
    TEST_SF_DD( s, gsl_sf_gegenpoly_3_e, { x: 100.0, y: 0.5 }, 161600.000, TEST_TOL0, "gsl_sf_gegenpoly_3_e");
 
    console.log( "    ... gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 1,    x:    1.0, y: 1.0 }, 2.000,                   TEST_TOL0, "gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 10,   x:    1.0, y: 1.0 }, 11.000,                  TEST_TOL0, "gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 10,   x:    1.0, y: 0.1 }, -0.4542309376,           TEST_TOL0, "gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 10,   x:    5.0, y: 1.0 }, 9.23780e+4,              TEST_TOL0, "gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 10,   x:  100.0, y: 0.5 }, 1.5729338392690000e+13,  TEST_TOL0, "gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 1000, x:  100.0, y: 1.0 }, 3.3353666135627322e+232, TEST_TOL1, "gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 100,  x: 2000.0, y: 1.0 }, 5.8753432034937579e+202, TEST_TOL0, "gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 103,  x:  207.0, y: 2.0 }, 1.4210272202235983e+145, TEST_TOL0, "gsl_sf_gegenpoly_n_e" );
    TEST_SF_IDD( s, gsl_sf_gegenpoly_n_e, { n: 103,  x:   -0.4, y: 0.3 }, -1.64527498094522e-04,   TEST_TOL1, "gsl_sf_gegenpoly_n_e" );

    console.log( "    ... gsl_sf_gegenpoly_array" );
    sa.Integer = 0;
    gsl_sf_gegenpoly_array( 99, 5.0, 1.0, ga ); console.log( ga );
    if ( test_sf_frac_diff( ga[1], 10.0 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( ga[10], 9.23780e+4 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_gegenpoly_array" );
    s.Integer = s.Integer + sa.Integer;

    return s;

} // test_gegen

// ----------------------------------------------------------------------------

export function test_jac( )
{
    var s = { Integer: 0 };
    var sa = { Integer: 0 };
    var u  = 0.0;
    var m  = 0.0;
    var sn = { Double: 0.0 };
    var cn = { Double: 0.0 };
    var dn = { Double: 0.0 };

    var mc = 0.0;
    var K  = 0.0;
    var A  = 0.0;
    var B  = 0.0;
    var C  = 0.0;
    var C2 = 0.0;
    var eps = 0.0;

    console.log( "Test Elliptic Functions (Jacobi) ..." );

    console.log( "    ... gsl_sf_elljac_e" );
    u = 0.5;
    m = 0.5;
    sa.Integer = 0;
    gsl_sf_elljac_e( u, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val( sn.Double, 0.4707504736556572833, TEST_TOL0, "gsl_sf_elljac_e(0.5|0.5) sn" );
    sa.Integer = sa.Integer + test_sf_val( cn.Double, 0.8822663948904402865, TEST_TOL0, "gsl_sf_elljac_e(0.5|0.5) cn" );
    sa.Integer = sa.Integer + test_sf_val( dn.Double, 0.9429724257773856873, TEST_TOL0, "gsl_sf_elljac_e(0.5|0.5) dn" );
    gsl_test( s, "  gsl_sf_elljac_e(0.5|0.5)" );
    s.Integer = s.Integer + sa.Integer;

    u = 1.0;
    m = 0.3;
    sa.Integer = 0;
    gsl_sf_elljac_e(u, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 0.8187707145344889190, TEST_TOL0, "gsl_sf_elljac_e(1.0|0.3) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 0.5741206467465548795, TEST_TOL0, "gsl_sf_elljac_e(1.0|0.3) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 0.8938033089590823040, TEST_TOL0, "gsl_sf_elljac_e(1.0|0.3) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(1.0|0.3)");
    s.Integer = s.Integer + sa.Integer;

    u = 1.0;
    m = 0.6;
    sa.Integer = 0;
    gsl_sf_elljac_e(u, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 0.7949388393365780943, TEST_TOL0, "gsl_sf_elljac_e(1.0|0.6) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 0.6066895760718277578, TEST_TOL0, "gsl_sf_elljac_e(1.0|0.6) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 0.7879361300438814425, TEST_TOL0, "gsl_sf_elljac_e(1.0|0.6) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(1.0|0.6)");
    s.Integer = s.Integer + sa.Integer;

    u = 3.0;
    m = 0.6;
    sa.Integer = 0;
    gsl_sf_elljac_e(u, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double,  0.7432676860864044186, TEST_TOL0, " gsl_sf_elljac_e(3.0|0.6) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -0.6689941306317733154, TEST_TOL0, " gsl_sf_elljac_e(3.0|0.6) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double,  0.8176379933025723259, TEST_TOL0, " gsl_sf_elljac_e(3.0|0.6) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(3.0|0.6)");
    s.Integer = s.Integer + sa.Integer;

    u = 2.0;
    m = 0.999999;
    sa.Integer = 0;
    gsl_sf_elljac_e(u, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 0.96402778575700186570, TEST_TOL1, "gsl_sf_elljac_e(2.0|0.999999) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 0.26580148285600686381, TEST_TOL1, "gsl_sf_elljac_e(2.0|0.999999) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 0.26580323105264131136, TEST_TOL1, "gsl_sf_elljac_e(2.0|0.999999) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(2.0|0.999999)");
    s.Integer = s.Integer + sa.Integer;

    // test supplied by Ivan Panchenko
    u = 1.69695970624443;
    m = 0.270378013104138;
    sa.Integer = 0;
    gsl_sf_elljac_e(u, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(1.69..|0.27..) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(1.69..|0.27..) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 0.8541791304497336, TEST_TOL1, "gsl_sf_elljac_e(1.69..|0.27..) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(1.69695970624443|0.270378013104138)");
    s.Integer = s.Integer + sa.Integer;

    // Check known values from Abramowitz  Stegun, Table 16.5
    u = 0.0;
    m = 0.1;

    mc = 1.0 - m;
    // quarter period K is (pi/2)/agm(1,mc)
    K = (M_PI_2) / 0.9741726903999478375938128316;

    A = 1.0 / Math.sqrt( 1.0 + Math.sqrt( mc ) );
    B = (mc ** 0.25) / Math.sqrt( 1.0 + Math.sqrt( mc ) );
    C = (mc ** 0.25);
    C2 = Math.sqrt( mc );

    eps = 1.0e-10;

    sa.Integer = 0;
    gsl_sf_elljac_e(0.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 0.0, TEST_TOL0, "gsl_sf_elljac_e(0|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(0|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(0|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(0|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -eps, TEST_TOL0, "gsl_sf_elljac_e(-1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(-1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(-1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, eps, TEST_TOL0, "gsl_sf_elljac_e(1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(1.0e-30, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 1.0e-30, TEST_TOL0, "gsl_sf_elljac_e(1e-30|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(1e-30|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL0, "gsl_sf_elljac_e(1e-30|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(1e-30|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(K / 2.0 - eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, A - eps*B*C, TEST_TOL2, "gsl_sf_elljac_e(K/2-1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, B + eps*A*C, TEST_TOL2, "gsl_sf_elljac_e(K/2-1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C + m*eps*A*B, TEST_TOL2, "gsl_sf_elljac_e(K/2-1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(K/2-1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, A, TEST_TOL2, "gsl_sf_elljac_e(K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, B, TEST_TOL2, "gsl_sf_elljac_e(K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(K / 2.0 + eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, A + eps*B*C, TEST_TOL2, "gsl_sf_elljac_e(K/2+1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, B - eps*A*C, TEST_TOL2, "gsl_sf_elljac_e(K/2+1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C - m*eps*A*B, TEST_TOL2, "gsl_sf_elljac_e(K/2+1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(K/2+1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(K - eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(K-1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, eps * C2, 10.0 * TEST_SNGL, "gsl_sf_elljac_e(K-1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C2, TEST_TOL2, "gsl_sf_elljac_e(K-1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(K-1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(K, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(K|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(K|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C2, TEST_TOL2, "gsl_sf_elljac_e(K|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(K|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(K + eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(K+1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -eps * C2, 10.0 * TEST_SNGL, "gsl_sf_elljac_e(K+1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C2, TEST_TOL2, "gsl_sf_elljac_e(K+1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(K+1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(3.0 * K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, A, TEST_TOL2, "gsl_sf_elljac_e(3K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -B, TEST_TOL2, "gsl_sf_elljac_e(3K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(3K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(3K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(2.0 * K - eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, eps, 10.0 * TEST_SNGL, "gsl_sf_elljac_e(2K-1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -1.0, TEST_TOL1, "gsl_sf_elljac_e(2K-1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(2K-1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(2K-1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(2.0 * K, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(2K|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -1.0, TEST_TOL1, "gsl_sf_elljac_e(2K|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(2K|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(2K|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(2.0 * K + eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -eps, 10.0 * TEST_SNGL, "gsl_sf_elljac_e(2K+1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -1.0, TEST_TOL1, "gsl_sf_elljac_e(2K+1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(2K+1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(2K+1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(5.0 * K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -A, TEST_TOL2, "gsl_sf_elljac_e(5K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -B, TEST_TOL2, "gsl_sf_elljac_e(5K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(5K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(5K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(3.0 * K - eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -1.0, TEST_TOL1, "gsl_sf_elljac_e(3K-1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -C2 * eps, 10.0 * TEST_SNGL, "gsl_sf_elljac_e(3K-1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C2, TEST_TOL2, "gsl_sf_elljac_e(3K-1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(3K-1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(3.0 * K, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -1.0, TEST_TOL1, "gsl_sf_elljac_e(3K|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(3K|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C2, TEST_TOL2, "gsl_sf_elljac_e(3K|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(3K|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(3.0 * K + eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -1.0, TEST_TOL1, "gsl_sf_elljac_e(3K+1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, +C2 * eps, 10.0 * TEST_SNGL, "gsl_sf_elljac_e(3K+1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C2, TEST_TOL2, "gsl_sf_elljac_e(3K+1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(3K+1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(7.0 * K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -A, TEST_TOL2, "gsl_sf_elljac_e(7K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, B, TEST_TOL2, "gsl_sf_elljac_e(7K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(7K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(7K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(4.0 * K - eps, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -eps, 10.0 * TEST_SNGL, "gsl_sf_elljac_e(4K-1e-10|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(4K-1e-10|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(4K-1e-10|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(4K-1e-10|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(4.0 * K, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(4K|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(4K|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(4K|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(4K|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(9.0 * K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, A, TEST_TOL2, "gsl_sf_elljac_e(9K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, B, TEST_TOL2, "gsl_sf_elljac_e(9K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(9K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(9K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -A, TEST_TOL2, "gsl_sf_elljac_e(-K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, B, TEST_TOL2, "gsl_sf_elljac_e(-K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(-K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-K, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -1.0, TEST_TOL1, "gsl_sf_elljac_e(-K|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(-K|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C2, TEST_TOL2, "gsl_sf_elljac_e(-K|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-K|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-3.0 * K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, -A, TEST_TOL2, "gsl_sf_elljac_e(-3K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -B, TEST_TOL2, "gsl_sf_elljac_e(-3K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(-3K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-3K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-2.0 * K, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(-2K|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -1.0, TEST_TOL1, "gsl_sf_elljac_e(-2K|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(-2K|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-2K|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-5.0 * K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, A, TEST_TOL2, "gsl_sf_elljac_e(-5K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, -B, TEST_TOL2, "gsl_sf_elljac_e(-5K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(-5K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-5K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-3.0 * K, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(-3K|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(-3K|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C2, TEST_TOL2, "gsl_sf_elljac_e(-3K|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-3K|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-7.0 * K / 2.0, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, A, TEST_TOL2, "gsl_sf_elljac_e(-7K/2|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, B, TEST_TOL2, "gsl_sf_elljac_e(-7K/2|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, C, TEST_TOL2, "gsl_sf_elljac_e(-7K/2|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-7K/2|0.1)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_elljac_e(-4.0 * K, m, sn, cn, dn);
    sa.Integer = sa.Integer + test_sf_val(sn.Double, 0.0, TEST_TOL1, "gsl_sf_elljac_e(-4K|0.1) sn");
    sa.Integer = sa.Integer + test_sf_val(cn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(-4K|0.1) cn");
    sa.Integer = sa.Integer + test_sf_val(dn.Double, 1.0, TEST_TOL1, "gsl_sf_elljac_e(-4K|0.1) dn");
    gsl_test(sa, "  gsl_sf_elljac_e(-4K|0.1)");
    s.Integer = s.Integer + sa.Integer;

    return s;

} // test_jac

// ----------------------------------------------------------------------------

export function test_laguerre( )
{
    var s = { Integer: 0 };

    console.log("Test Laguerre Functions ...");

    console.log("    ... gsl_sf_laguerre_1_e");
    TEST_SF_DD(s, gsl_sf_laguerre_1_e, { x: 0.5, y:-1.0 }, 2.5, TEST_TOL0, "gsl_sf_laguerre_1_e");
    TEST_SF_DD(s, gsl_sf_laguerre_1_e, { x: 0.5, y: 1.0 }, 0.5, TEST_TOL0, "gsl_sf_laguerre_1_e");
    TEST_SF_DD(s, gsl_sf_laguerre_1_e, { x: 1.0, y: 1.0 }, 1.0, TEST_TOL0, "gsl_sf_laguerre_1_e");

    console.log("    ... gsl_sf_laguerre_2_e");
    TEST_SF_DD(s, gsl_sf_laguerre_2_e, { x:  0.5, y:-1.0 }, 4.875,  TEST_TOL0, "gsl_sf_laguerre_2_e");
    TEST_SF_DD(s, gsl_sf_laguerre_2_e, { x:  0.5, y: 1.0 }, -0.125, TEST_TOL0, "gsl_sf_laguerre_2_e");
    TEST_SF_DD(s, gsl_sf_laguerre_2_e, { x:  1.0, y: 1.0 },  0.5, TEST_TOL0, "gsl_sf_laguerre_2_e");
    TEST_SF_DD(s, gsl_sf_laguerre_2_e, { x: -1.0, y: 1.0 }, -0.5, TEST_TOL0, "gsl_sf_laguerre_2_e");
    TEST_SF_DD(s, gsl_sf_laguerre_2_e, { x: -2.0, y: 1.0 },  0.5, TEST_TOL0, "gsl_sf_laguerre_2_e");
    TEST_SF_DD(s, gsl_sf_laguerre_2_e, { x: -3.0, y: 1.0 },  2.5, TEST_TOL0, "gsl_sf_laguerre_2_e");

    console.log("    ... gsl_sf_laguerre_3_e");
    TEST_SF_DD(s, gsl_sf_laguerre_3_e, { x: 0.5,  y:-1.0 }, 8.479166666666666667,    TEST_TOL0, "gsl_sf_laguerre_3_e");
    TEST_SF_DD(s, gsl_sf_laguerre_3_e, { x: 0.5,  y: 1.0 }, -0.6041666666666666667,  TEST_TOL0, "gsl_sf_laguerre_3_e");
    TEST_SF_DD(s, gsl_sf_laguerre_3_e, { x: 1.0,  y: 1.0 }, -0.16666666666666666667, TEST_TOL1, "gsl_sf_laguerre_3_e");
    TEST_SF_DD(s, gsl_sf_laguerre_3_e, { x:  2.0, y: 1.0 }, 2.3333333333333333333,  TEST_TOL0, "gsl_sf_laguerre_3_e");
    TEST_SF_DD(s, gsl_sf_laguerre_3_e, { x: -2.0, y: 1.0 }, 1.0/3.0,  TEST_TOL0, "gsl_sf_laguerre_3_e");
    TEST_SF_DD(s, gsl_sf_laguerre_3_e, { x: -3.0, y: 1.0 }, -1.0/6.0, TEST_TOL0, "gsl_sf_laguerre_3_e");
    TEST_SF_DD(s, gsl_sf_laguerre_3_e, { x: -4.0, y: 1.0 }, -8.0/3.0, TEST_TOL0, "gsl_sf_laguerre_3_e");

    console.log("    ... gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 1, x: 0.5, y: 1.0 }, 0.5, TEST_TOL0, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 2, x: 1.0, y: 1.0 }, 0.5, TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 3, x: 2.0, y: 1.0 }, 2.3333333333333333333,   TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 4, x: 2.0, y: 0.5 }, 6.752604166666666667,    TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 90, x: 2.0, y:  0.5 }, -48.79047157201507897, TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 90, x: 2.0, y: -100.0 }, 2.5295879275042410902e+63, TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 90, x: 2.0, y:  100.0 }, -2.0929042259546928670e+20, TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: 2.0, y: -0.5 }, 2.2521795545919391405e+07,  TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: 2.0, y:  0.5 }, -28.764832945909097418,     TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 1000, x: 2.0, y: -0.5 }, 2.4399915170947549589e+21, TEST_TOL3, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 1000, x: 2.0, y:  0.5 }, -306.77440254315317525,    TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100000, x: 2.0, y: 1.0 }, 5107.73491348319,         TEST_TOL4, "gsl_sf_laguerre_n_e");

    // Compute these with the recurrence
    // L(0,alpha,x)=1;
    // L(1,alpha,x)=1+alpha-x;
    // L(n,alpha,x)=((2*n-1+alpha-x)*L(n-1,alpha,x)-(n+alpha-1)*L(n-2,alpha,x))/k

    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100000, x: 2.5, y: 2.5 },   -0.41491680394598644969113795e5, TEST_TOL4, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100001, x: 2.5, y: 2.5 }, -0.41629446949552321027514888e5, TEST_TOL4, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 1000001, x: 2.5, y: 2.5 }, -0.48017961545391273151977118e6, TEST_TOL4, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 5000001, x: 2.5, y: 2.5 }, -0.15174037401611122446089494e7, TEST_TOL6, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 8000001, x: 2.5, y: 2.5 },  0.63251509472091810994286362e6, TEST_SNGL, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 10000001, x: 2.5, y: 2.5 },  0.15299484685632983178033887e7, TEST_SNGL, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100000001, x: 2.5, y: 2.5 },  0.23645341644922756725290777e8, TEST_SNGL, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 1000000001, x: 2.5, y: 2.5 }, -0.17731002248958790286185878e8, 100.0*TEST_SNGL, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 1, x: -2.0, y:  1.0 },  -2.0,     TEST_TOL0, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 2, x: -2.0, y:  1.0 },   0.5,     TEST_TOL0, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 3, x: -2.0, y:  1.0 },   1.0/3.0, TEST_TOL0, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 10, x: -2.0, y:  1.0 }, -0.04654954805996472663,   TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 10, x: -5.0, y:  1.0 }, -0.0031385030864197530864, TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 10, x: -9.0, y:  1.0 }, -2.480158730158730159e-06, TEST_TOL5, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 10, x: -11.0,  y:  1.0 }, 2.7182818011463844797,    TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 10, x: -11.0,  y: -1.0 }, 0.3678794642857142857,    TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -2.0,  y:  1.0 },  -0.0027339992019526273866,  TEST_SQRT_TOL0, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -2.0,  y: -1.0 },   229923.09193402028290,     TEST_TOL5, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -10.0,  y:  1.0 },  3.25966665871244092e-11,   TEST_TOL6, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -10.0,  y: -1.0 },  0.00016484365618205810025, TEST_TOL6, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -20.0, y:  1.0 },  5.09567630343671251e-21,  TEST_TOL3, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -30.0, y:  1.0 },  3.46063150272466192e-34,  TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -50.0,  y:  1.0 },  1.20981872933162889e-65,  TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -50.0,  y: -1.0 },  8.60763477742332922e-65,  TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -50.5,  y:  1.0 },  4.84021010426688393e-31,  TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -50.5,  y: -1.0 },  8.49861345212160618e-33,  TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -101.0,  y:  1.0 }, 2.7182818284590452354,    TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -101.0,  y: -1.0 }, 0.3678794411714423216,    TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -102.0,  y:  1.0 }, 271.8281828459045235,    TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -102.0,  y: -1.0 }, 37.52370299948711680,    TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -110.0,  y:  1.0 }, 1.0666955248998831554e+13, TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -110.0,  y: -1.0 }, 1.7028306108058225871e+12, TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -200.0,  y:  1.0 }, 7.47851889721356628e+58,  TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -200.0,  y: -1.0 }, 2.73740299754732273e+58,  TEST_TOL1, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -50.0,   y: 10.0 }, 4.504712811317745591e-21,  TEST_SQRT_TOL0, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: -50.0, y: -10.0 }, 1.475165520610679937e-11,  TEST_TOL1, "gsl_sf_laguerre_n_e");

    // test cases for Ed Smith-Rowland

    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: 0.0, y: 0.5 }, 0.18682260367692278801, TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: 0.0, y: 10.5 }, 9.1796907354050059874, TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: 0.0, y: -10.5 }, 5.6329215744170606488e24, TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: 0.0, y: 100.5 }, -3.9844782875811907525e20, TEST_TOL2, "gsl_sf_laguerre_n_e");
    TEST_SF_IDD(s, gsl_sf_laguerre_n_e, { n: 100, x: 0.0, y: 150.0 }, -1.4463204337261709595e31, TEST_TOL2, "gsl_sf_laguerre_n_e");

    return s;

} // test_laguerre

// ----------------------------------------------------------------------------

export function test_lambert( )
{
    var s = { Integer: 0 };

    console.log( "Test Lambert W Functions ..." );

    console.log( "    ... gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, 0.0,      0.0,                        TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, 1.0,      0.567143290409783872999969, TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, 2.0,      0.852605502013725491346472, TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, 20.0,     2.205003278024059970493066, TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, 1000.0,   5.24960285240159622712606,  TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, 1.0e+6,   11.38335808614005262200016, TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, 1.0e+12,  24.43500440493491313826305, TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, 1.0e+308, 702.641362034106812081125,  TEST_TOL0, "gsl_sf_lambert_W0_e" );

    // Test case from Katrin Wolff <katrin_wolff@gmx.de> fails under
    // double-precision

    TEST_SF_D( s, gsl_sf_lambert_W0_e, (1.6849341956993852953416990), 0.775706963944252869680440,  TEST_TOL0, "gsl_sf_lambert_W0_e");

    TEST_SF_D( s, gsl_sf_lambert_W0_e, (-1.0/M_E - GSL_DBL_EPSILON),            -1.0,                         TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, (-1.0/M_E + 1.0/(1024.0*1024.0*1024.0)), -0.999928845560308370714970,  TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, (-1.0/M_E + 1.0/(1024.0*1024.0)),        -0.997724730359774141620354,  TEST_TOL0, "gsl_sf_lambert_W0_e");
    TEST_SF_D( s, gsl_sf_lambert_W0_e, (-1.0/M_E + 1.0/512.0),                  -0.900335676696088773044678,  TEST_TOL0, "gsl_sf_lambert_W0_e" );
    TEST_SF_D( s, gsl_sf_lambert_W0_e, (-1.0/M_E + 0.25),                       -0.1349044682661213545487599, TEST_TOL0, "gsl_sf_lambert_W0_e" );

    console.log( "    ... gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, 0.0,  0.0,                        TEST_TOL0, "gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, 1.0,  0.567143290409783872999969, TEST_TOL0, "gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, 2.0,  0.852605502013725491346472, TEST_TOL0, "gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, 20.0, 2.205003278024059970493066, TEST_TOL0, "gsl_sf_lambert_Wm1_e" );

    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, -1.0/M_E - GSL_DBL_EPSILON,            -1.0,                        TEST_TOL0, "gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, -1.0/M_E + 1.0/(1024.0*1024.0*1024.0), -1.000071157815154608049055, TEST_TOL1, "gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, -1.0/M_E + 1.0/(1024.0*1024.0),        -1.002278726118593023934693, TEST_TOL1, "gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, -1.0/M_E + 1.0/512.0,                  -1.106761200865743124599130, TEST_TOL1, "gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, -1.0/M_E + 1.0/64.0,                   -1.324240940341812125489772, TEST_TOL1, "gsl_sf_lambert_Wm1_e" );
    TEST_SF_D( s, gsl_sf_lambert_Wm1_e, -1.0/M_E + 0.25,                       -3.345798131120112,          TEST_TOL1, "gsl_sf_lambert_Wm1_e" );

    return s;

} // test_lambert

// ----------------------------------------------------------------------------

export function test_log( )
{
    var s = { Integer: 0 };

    console.log("Test Logarithmic Functions ...");

    console.log("    ... gsl_sf_log_e");
    TEST_SF_D(s, gsl_sf_log_e, 0.1, -2.3025850929940456840,  TEST_TOL0, "gsl_sf_log_e");
    TEST_SF_D(s, gsl_sf_log_e, 1.1, 0.09531017980432486004,  TEST_TOL1, "gsl_sf_log_e");
    TEST_SF_D(s, gsl_sf_log_e, 1000.0, 6.907755278982137052, TEST_TOL0, "gsl_sf_log_e");

    console.log("    ... gsl_sf_log_abs_e");
    TEST_SF_D(s, gsl_sf_log_abs_e, -0.1, -2.3025850929940456840,  TEST_TOL0, "gsl_sf_log_abs_e");
    TEST_SF_D(s, gsl_sf_log_abs_e, -1.1, 0.09531017980432486004,  TEST_TOL1, "gsl_sf_log_abs_e");
    TEST_SF_D(s, gsl_sf_log_abs_e, -1000.0, 6.907755278982137052, TEST_TOL0, "gsl_sf_log_abs_e");
    TEST_SF_D(s, gsl_sf_log_abs_e, 0.1, -2.3025850929940456840,  TEST_TOL0, "gsl_sf_log_abs_e");
    TEST_SF_D(s, gsl_sf_log_abs_e, 1.1, 0.09531017980432486004,  TEST_TOL1, "gsl_sf_log_abs_e");
    TEST_SF_D(s, gsl_sf_log_abs_e, 1000.0, 6.907755278982137052, TEST_TOL0, "gsl_sf_log_abs_e");

    console.log("    ... gsl_sf_complex_log_e");
    TEST_SF_2( s, gsl_sf_complex_log_e, { x: 1.0, y: 1.0 },
        0.3465735902799726547, TEST_TOL0,
        0.7853981633974483096, TEST_TOL0,
        "gsl_sf_complex_log_e" );

    TEST_SF_2( s, gsl_sf_complex_log_e, { x: 1.0, y: -1.0 },
         0.3465735902799726547, TEST_TOL0,
        -0.7853981633974483096, TEST_TOL0,
        "gsl_sf_complex_log_e" );

    TEST_SF_2( s, gsl_sf_complex_log_e, { x: 1.0, y: 100.0 },
        4.605220183488258022, TEST_TOL0,
        1.560796660108231381, TEST_TOL0,
        "gsl_sf_complex_log_e" );

    TEST_SF_2( s, gsl_sf_complex_log_e, { x: -1000.0, y: -1.0 },
        6.907755778981887052, TEST_TOL0,
        -3.1405926539231263718, TEST_TOL0,
        "gsl_sf_complex_log_e" );

    TEST_SF_2( s, gsl_sf_complex_log_e, { x: -1.0, y: 0.0 },
        0.0, TEST_TOL0,
        3.1415926535897932385, TEST_TOL0,
        "gsl_sf_complex_log_e" );

    console.log("    ... gsl_sf_log_1plusx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_e, 1.0e-10, 9.999999999500000000e-11, TEST_TOL0, "gsl_sf_log_1plusx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_e, 1.0e-8, 9.999999950000000333e-09, TEST_TOL0, "gsl_sf_log_1plusx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_e, 1.0e-4, 0.00009999500033330833533, TEST_TOL0, "gsl_sf_log_1plusx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_e, 0.1, 0.09531017980432486004, TEST_TOL0, "gsl_sf_log_1plusx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_e, 0.49, 0.3987761199573677730, TEST_TOL0, "gsl_sf_log_1plusx_e");

    TEST_SF_D(s, gsl_sf_log_1plusx_e, -0.49, -0.6733445532637655964, TEST_TOL0, "gsl_sf_log_1plusx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_e, 1.0, M_LN2, TEST_TOL0, "gsl_sf_log_1plusx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_e, -0.99, -4.605170185988091368, TEST_TOL0, "gsl_sf_log_1plusx_e");

    console.log("    ... gsl_sf_log_1plusx_mx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_mx_e, 1.0e-10, -4.999999999666666667e-21, TEST_TOL0, "gsl_sf_log_1plusx_mx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_mx_e, 1.0e-8, -4.999999966666666917e-17, TEST_TOL0, "gsl_sf_log_1plusx_mx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_mx_e, 1.0e-4, -4.999666691664666833e-09, TEST_TOL0, "gsl_sf_log_1plusx_mx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_mx_e, 0.1, -0.004689820195675139956, TEST_TOL0, "gsl_sf_log_1plusx_mx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_mx_e, 0.49, -0.09122388004263222704, TEST_TOL0, "gsl_sf_log_1plusx_mx_e");

    TEST_SF_D(s, gsl_sf_log_1plusx_mx_e, -0.49, -0.18334455326376559639, TEST_TOL0, "gsl_sf_log_1plusx_mx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_mx_e, 1.0, M_LN2-1.0, TEST_TOL0, "gsl_sf_log_1plusx_mx_e");
    TEST_SF_D(s, gsl_sf_log_1plusx_mx_e, -0.99, -3.615170185988091368, TEST_TOL0, "gsl_sf_log_1plusx_mx_e");

    return s;

} // test_log

// ----------------------------------------------------------------------------

export function test_pow_int( )
{
    var s = { Integer: 0 };

    console.log( "Test Integer Powers ..." );

    console.log( "    ... gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:  2.0, i:  3 }, 8.0,        TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x: -2.0, i:  3 }, -8.0,       TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:  2.0, i: -3 }, 1.0 / 8.0,  TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x: -2.0, i: -3 }, -1.0 / 8.0, TEST_TOL0, "gsl_sf_pow_int_e" );

    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:  10.0, i:  4 }, 1.0e+4, TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:  10.0, i: -4 }, 1.0e-4, TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x: -10.0, i:  4 }, 1.0e+4, TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x: -10.0, i: -4 }, 1.0e-4, TEST_TOL0, "gsl_sf_pow_int_e" );

    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:  10.0, i:  40 }, 1.0e+40,                  TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:   8.0, i: -40 }, 7.523163845262640051e-37, TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x: -10.0, i:  40 }, 1.0e+40,                  TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:  -8.0, i: -40 }, 7.523163845262640051e-37, TEST_TOL0, "gsl_sf_pow_int_e" );

    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:  10.0, i:  41 }, 1.0e+41,                   TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:   8.0, i: -41 }, 9.403954806578300064e-38,  TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x: -10.0, i:  41 }, -1.0e+41,                  TEST_TOL0, "gsl_sf_pow_int_e" );
    TEST_SF_DI( s,  gsl_sf_pow_int_e, { x:  -8.0, i: -41 }, -9.403954806578300064e-38, TEST_TOL0, "gsl_sf_pow_int_e" );

    return s;

} // test_pow_int

// ----------------------------------------------------------------------------

export function test_psi( )
{
    var s = { Integer: 0 };

    console.log( "Test Psi Functions ..." );

    // Test values taken 1-4 from gp-pari

    console.log( "    ... gsl_sf_psi_int_e" );
    TEST_SF_D( s, gsl_sf_psi_int_e, 1, -0.57721566490153286060, TEST_TOL0, "gsl_sf_psi_int_e" );
    TEST_SF_D( s, gsl_sf_psi_int_e, 2, 0.42278433509846713939,  TEST_TOL0, "gsl_sf_psi_int_e" );
    TEST_SF_D( s, gsl_sf_psi_int_e, 3, 0.92278433509846713939,  TEST_TOL0, "gsl_sf_psi_int_e" );
    TEST_SF_D( s, gsl_sf_psi_int_e, 4, 1.2561176684318004727,   TEST_TOL0, "gsl_sf_psi_int_e" );

    TEST_SF_D( s, gsl_sf_psi_int_e, 5,    1.5061176684318004727, TEST_TOL0, "gsl_sf_psi_int_e" );
    TEST_SF_D( s, gsl_sf_psi_int_e, 100,  4.600161852738087400,  TEST_TOL0, "gsl_sf_psi_int_e" );
    TEST_SF_D( s, gsl_sf_psi_int_e, 110,  4.695928024251535633,  TEST_TOL0, "gsl_sf_psi_int_e" );
    TEST_SF_D( s, gsl_sf_psi_int_e, 5000, 8.517093188082904107,  TEST_TOL0, "gsl_sf_psi_int_e" );

    console.log( "    ... gsl_sf_psi_e" );
    TEST_SF_D( s, gsl_sf_psi_e, 5000.0,        8.517093188082904107,  TEST_TOL0, "gsl_sf_psi_e" );
    TEST_SF_D( s, gsl_sf_psi_e, 5.0,           1.5061176684318004727, TEST_TOL0, "gsl_sf_psi_e" );
    TEST_SF_D( s, gsl_sf_psi_e, -10.5,         2.3982391295357816134, TEST_TOL0, "gsl_sf_psi_e" );
    TEST_SF_D( s, gsl_sf_psi_e, -100.5,        4.615124601338064117,  TEST_TOL2, "gsl_sf_psi_e" );
    TEST_SF_D( s, gsl_sf_psi_e, -1.0e+5-0.5,   11.512935464924395337, 4.0 * TEST_TOL4, "gsl_sf_psi_e" );
    TEST_SF_D( s, gsl_sf_psi_e, -262144.0-0.5, 12.476653064769611581, 4.0 * TEST_TOL4, "gsl_sf_psi_e" );

    console.log( "    ... gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, 0.8,   -0.07088340212750589223, TEST_TOL1, "gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, 1.0,    0.09465032062247697727, TEST_TOL0, "gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, 5.0,    1.6127848446157465854,  TEST_TOL2, "gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, 100.0,  4.605178519404762003,   TEST_TOL0, "gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, 2000.0, 7.600902480375416216,   TEST_TOL0, "gsl_sf_psi_1piy_e" );

    TEST_SF_D( s, gsl_sf_psi_1piy_e, -0.8,   -0.07088340212750589223, TEST_TOL1, "gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, -1.0,    0.09465032062247697727, TEST_TOL0, "gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, -5.0,    1.6127848446157465854,  TEST_TOL2, "gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, -100.0,  4.605178519404762003,   TEST_TOL0, "gsl_sf_psi_1piy_e" );
    TEST_SF_D( s, gsl_sf_psi_1piy_e, -2000.0, 7.600902480375416216,   TEST_TOL0, "gsl_sf_psi_1piy_e" );

    // Additional test values 1-4 computed using gp-pari and
    // Abramowitz & Stegun 6.4.6

    console.log( "    ... gsl_sf_psi_1_int_e" );
    TEST_SF_I( s, gsl_sf_psi_1_int_e, 1, 1.6449340668482264364,  TEST_TOL0, "gsl_sf_psi_1_int_e" );
    TEST_SF_I( s, gsl_sf_psi_1_int_e, 2, 0.64493406684822643647, TEST_TOL0, "gsl_sf_psi_1_int_e" );
    TEST_SF_I( s, gsl_sf_psi_1_int_e, 3, 0.39493406684822643647, TEST_TOL0, "gsl_sf_psi_1_int_e" );
    TEST_SF_I( s, gsl_sf_psi_1_int_e, 4, 0.28382295573711532536, TEST_TOL0, "gsl_sf_psi_1_int_e" );

    TEST_SF_I( s, gsl_sf_psi_1_int_e, 1, 1.6449340668482264365,      TEST_TOL0, "gsl_sf_psi_1_int_e" );
    TEST_SF_I( s, gsl_sf_psi_1_int_e, 5, 0.22132295573711532536,     TEST_TOL0, "gsl_sf_psi_1_int_e" );
    TEST_SF_I( s, gsl_sf_psi_1_int_e, 100, 0.010050166663333571395,  TEST_TOL0, "gsl_sf_psi_1_int_e" );
    TEST_SF_I( s, gsl_sf_psi_1_int_e, 110, 0.009132356622022545705,  TEST_TOL0, "gsl_sf_psi_1_int_e" );
    TEST_SF_I( s, gsl_sf_psi_1_int_e, 500, 0.0020020013333322666697, TEST_TOL0, "gsl_sf_psi_1_int_e" );

    console.log( "    ... gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, 1.0/32.0, 1025.5728544782377089,  TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, 1.0, 1.6449340668482264365,       TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, 5.0, 0.22132295573711532536,      TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, 100.0, 0.010050166663333571395,   TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, 110.0, 0.009132356622022545705,   TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, 500.0, 0.0020020013333322666697,  TEST_TOL0, "gsl_sf_psi_1_e" );

    TEST_SF_D( s, gsl_sf_psi_1_e, -1.0 - 1.0/128.0, 16386.648472598746587, TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, -1.50,   9.3792466449891237539, TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, -10.5,   9.7787577398148123845, TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, -15.5,   9.8071247184113896201, TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, -50.5,   9.8499971860824842274, TEST_TOL0, "gsl_sf_psi_1_e" );
    TEST_SF_D( s, gsl_sf_psi_1_e, -1000.5, 9.8686054001734414233, TEST_TOL0, "gsl_sf_psi_1_e" );

    console.log( "    ... gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: 1.0 }, 1.6449340668482264364,   TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: 2.0 }, 0.64493406684822643647,  TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: 3.0 }, 0.39493406684822643647,  TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: 4.0 }, 0.28382295573711532536,  TEST_TOL0, "gsl_sf_psi_n_e" );

    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: 5.0 },   0.22132295573711532536,   TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: 100.0 }, 0.010050166663333571395,  TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: 110.0 }, 0.009132356622022545705,  TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: 500.0 }, 0.0020020013333322666697, TEST_TOL0, "gsl_sf_psi_n_e" );

    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 3, x: 5.0   }, 0.021427828192755075022,   TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 3, x: 500.0 }, 1.6048063999872000683e-08, TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 10,x: 5.0   }, -0.08675107579196581317,   TEST_TOL1, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 10,x: 50.0  }, -4.101091112731268288e-12, TEST_TOL0, "gsl_sf_psi_n_e" );

    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 0, x: -1.5 }, 0.70315664064524318723,  TEST_TOL0, "gsl_sf_psi_n_e" );
    TEST_SF_ID( s, gsl_sf_psi_n_e, { i: 1, x: -1.5 }, 9.3792466449891237539,   TEST_TOL0, "gsl_sf_psi_n_e" );

    return s;

} // test_psi

// ----------------------------------------------------------------------------

// export function test_psi_complex( )
// {
//     var s = 0;

//     console.log( "Test Complex Psi Functions ..." );

//     console.log( "    ... gsl_sf_complex_psi_e" );
//    TEST_SF_2(s, gsl_sf_complex_psi_e'Access, (1.0e+07, 1.0e+06),
//              16.1230707668799525, TEST_TOL0,
//              0.09966865744165720, TEST_TOL0,
//              GSL_SUCCESS, "gsl_sf_complex_psi_e");
//
//    TEST_SF_2(s, gsl_sf_complex_psi_e'Access, (10.0, 50.0),
//              3.92973987174863660, TEST_TOL0,
//              1.38302847985210276, TEST_TOL0,
//              GSL_SUCCESS, "gsl_sf_complex_psi_e");
//
//    TEST_SF_2(s, gsl_sf_complex_psi_e'Access, (2.0, 21.0),
//              3.04697388853248195, TEST_TOL0,
//              1.49947549076817824, TEST_TOL0,
//              GSL_SUCCESS, "gsl_sf_complex_psi_e");
//
//    TEST_SF_2(s, gsl_sf_complex_psi_e'Access, (1.5, 0.0),
//              0.0364899739785765206, TEST_TOL2,
//              0.0, TEST_TOL1,
//              GSL_SUCCESS, "gsl_sf_complex_psi_e");
//
//    TEST_SF_2(s, gsl_sf_complex_psi_e'Access, (1.0, 5.0),
//              1.612784844615747, TEST_TOL1,
//              1.470796326794968, TEST_TOL1,
//              GSL_SUCCESS, "gsl_sf_complex_psi_e");
//
//    TEST_SF_2(s, gsl_sf_complex_psi_e'Access, (-1.5, 5.0),
//              1.68260717336484070, TEST_TOL0,
//              1.95230236730713338, TEST_TOL0,
//              GSL_SUCCESS, "gsl_sf_complex_psi_e");
//
//    TEST_SF_2(s, gsl_sf_complex_psi_e'Access, (-20.5, -20.5),
//              3.37919358657933066, TEST_TOL0,
//             -2.36829046481731091, TEST_TOL0,
//              GSL_SUCCESS, "gsl_sf_complex_psi_e");

//     return s;

// } // test_psi_complex

// ----------------------------------------------------------------------------

export function test_synch( )
{
    var s = { Integer: 0 };

    console.log( "Test Synchrotron Functions ..." );

    console.log( "    ... gsl_sf_synchrotron_1_e" );
    TEST_SF_D( s, gsl_sf_synchrotron_1_e, 0.01,  0.444972504114210632,    TEST_TOL0, "gsl_sf_synchrotron_1_e" );
    TEST_SF_D( s, gsl_sf_synchrotron_1_e, 1.0,   0.651422815355364504,    TEST_TOL1, "gsl_sf_synchrotron_1_e" );
    TEST_SF_D( s, gsl_sf_synchrotron_1_e, 10.0,  0.000192238264300868882, TEST_TOL1, "gsl_sf_synchrotron_1_e" );
    TEST_SF_D( s, gsl_sf_synchrotron_1_e, 100.0, 4.69759366592220221e-43, TEST_TOL1, "gsl_sf_synchrotron_1_e" );

    console.log( "    ... gsl_sf_synchrotron_2_e" );
    TEST_SF_D( s, gsl_sf_synchrotron_2_e, 0.01,  0.23098077342226277732,     TEST_TOL2, "gsl_sf_synchrotron_2_e" );
    TEST_SF_D( s, gsl_sf_synchrotron_2_e, 1.0,   0.4944750621042082670,      TEST_TOL1, "gsl_sf_synchrotron_2_e" );
    TEST_SF_D( s, gsl_sf_synchrotron_2_e, 10.0,  0.00018161187569530204281,  TEST_TOL1, "gsl_sf_synchrotron_2_e" );
    TEST_SF_D( s, gsl_sf_synchrotron_2_e, 256.0, 1.3272635474353774058e-110, TEST_TOL4, "gsl_sf_synchrotron_2_e" );  // exp()... not my fault

    return s;

} // test_synch

// ----------------------------------------------------------------------------

export function test_transport( )
{
    var s = { Integer: 0 };

    console.log( "Test Transport Functions ..." );

    console.log( "    ... gsl_sf_transport_2_e" );
    TEST_SF_D( s, gsl_sf_transport_2_e, 1.0e-10, 9.9999999999999999999e-11, TEST_TOL0, "gsl_sf_transport_2_e" );
    TEST_SF_D( s, gsl_sf_transport_2_e, 1.0,     0.97303256135517012845, TEST_TOL0, "gsl_sf_transport_2_e" );
    TEST_SF_D( s, gsl_sf_transport_2_e, 3.0,     2.41105004901695346199, TEST_TOL0, "gsl_sf_transport_2_e" );
    TEST_SF_D( s, gsl_sf_transport_2_e, 10.0,    3.28432911449795173575, TEST_TOL0, "gsl_sf_transport_2_e" );
    TEST_SF_D( s, gsl_sf_transport_2_e, 100.0,   3.28986813369645287294, TEST_TOL0, "gsl_sf_transport_2_e" );
    TEST_SF_D( s, gsl_sf_transport_2_e, 1.0e+05, 3.28986813369645287294, TEST_TOL0, "gsl_sf_transport_2_e" );

    console.log( "    ... gsl_sf_transport_3_e" );
    TEST_SF_D( s, gsl_sf_transport_3_e, 1.0e-10, 4.999999999999999999997e-21, TEST_TOL0, "gsl_sf_transport_3_e" );
    TEST_SF_D( s, gsl_sf_transport_3_e, 1.0,     0.479841006572417499939, TEST_TOL0, "gsl_sf_transport_3_e" );
    TEST_SF_D( s, gsl_sf_transport_3_e, 3.0,     3.210604662942246772338, TEST_TOL0, "gsl_sf_transport_3_e" );
    TEST_SF_D( s, gsl_sf_transport_3_e, 5.0,     5.614386613842273228585, TEST_TOL0, "gsl_sf_transport_3_e" );
    TEST_SF_D( s, gsl_sf_transport_3_e, 10.0,    7.150322712008592975030, TEST_TOL0, "gsl_sf_transport_3_e" );
    TEST_SF_D( s, gsl_sf_transport_3_e, 30.0,    7.212341416160946511930, TEST_TOL0, "gsl_sf_transport_3_e" );
    TEST_SF_D( s, gsl_sf_transport_3_e, 100.0,   7.212341418957565712398, TEST_TOL0, "gsl_sf_transport_3_e" );
    TEST_SF_D( s, gsl_sf_transport_3_e, 1.0e+05, 7.212341418957565712398, TEST_TOL0, "gsl_sf_transport_3_e" );

    console.log( "    ... gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 1.0e-10, 3.33333333333333333333e-31, TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 1.0e-07, 3.33333333333333166666e-22, TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 1.0e-04, 3.33333333166666666726e-13, TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 0.1,     0.000333166726172109903824, TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 1.0,     0.31724404523442648241,     TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 3.0,     5.96482239737147652446,     TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 5.0,     15.3597843168821829816,     TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 10.0,    25.2736676770304417334,     TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 30.0,    25.9757575220840937469,     TEST_TOL0, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 100.0,   25.9757576090673165963,     TEST_TOL1, "gsl_sf_transport_4_e" );
    TEST_SF_D( s, gsl_sf_transport_4_e, 1.0e+05, 25.9757576090673165963,     TEST_TOL2, "gsl_sf_transport_4_e" );

    console.log( "    ... gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 1.0e-10, 2.49999999999999999999e-41, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 1.0e-07, 2.49999999999999861111e-29, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 1.0e-04, 2.49999999861111111163e-17, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 0.1,     0.000024986116317791487410, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 1.0,     0.236615879239094789259153, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 3.0,     12.77055769104415951115760, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 5.0,     50.26309221817518778543615, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 10.0,    116.3807454024207107698556, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 30.0,    124.4313279083858954839911, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 100.0,   124.4313306172043911597639, TEST_TOL0, "gsl_sf_transport_5_e" );
    TEST_SF_D( s, gsl_sf_transport_5_e, 1.0e+05, 124.43133061720439115976,   TEST_TOL0, "gsl_sf_transport_5_e" );

    return s;

} // test_transport

// ----------------------------------------------------------------------------

export function test_trig()
{
    const D = 1.2246467991473531772e-16;
    var s = { Integer: 0 };
    var theta = 0.0;
    var sa    = { Integer: 0 };
    // var code  = 0;

    console.log("Test Trigonometric Functions ...");

    console.log("    ... gsl_sf_sin_e");
    TEST_SF_D(s, gsl_sf_sin_e, -10.0,       0.5440211108893698134,    TEST_TOL0, "gsl_sf_sin_e");
    TEST_SF_D(s, gsl_sf_sin_e, 1.0,         0.8414709848078965067,    TEST_TOL0, "gsl_sf_sin_e");
    TEST_SF_D(s, gsl_sf_sin_e, 1000.0,      0.8268795405320025603,    TEST_TOL0, "gsl_sf_sin_e");
    TEST_SF_D(s, gsl_sf_sin_e, 1048576.75,  0.8851545351115651914,    TEST_TOL1, "gsl_sf_sin_e");
    TEST_SF_D(s, gsl_sf_sin_e, 62831853.75, 0.6273955953485000827,    TEST_TOL3, "gsl_sf_sin_e");
    TEST_SF_D(s, gsl_sf_sin_e, 1073741822.5, -0.8284043541754465988,  TEST_SQRT_TOL0, "gsl_sf_sin_e");
    TEST_SF_D(s, gsl_sf_sin_e, 1073741824.0, -0.6173264150460421708,  TEST_SQRT_TOL0, "gsl_sf_sin_e");
    TEST_SF_D(s, gsl_sf_sin_e, 1073741825.5,  0.7410684679436226926,  TEST_SQRT_TOL0, "gsl_sf_sin_e");
    //
    TEST_SF_D(s, gsl_sf_sin_e, 1099511627776.0, -0.4057050115328287198, 32.0 * TEST_SQRT_TOL0, "gsl_sf_sin_e");
    //

    console.log("    ... gsl_sf_cos_e");
    TEST_SF_D(s, gsl_sf_cos_e, -10.0,      -0.8390715290764524523,    TEST_TOL0, "gsl_sf_cos_e");
    TEST_SF_D(s, gsl_sf_cos_e, 1.0,         0.5403023058681397174,    TEST_TOL0, "gsl_sf_cos_e");
    TEST_SF_D(s, gsl_sf_cos_e, 1000.0,      0.5623790762907029911,    TEST_TOL1, "gsl_sf_cos_e");
    TEST_SF_D(s, gsl_sf_cos_e, 1048576.75,  0.4652971620066351799,    TEST_TOL2, "gsl_sf_cos_e");
    TEST_SF_D(s, gsl_sf_cos_e, 62831853.75, 0.7787006914966116436,    TEST_TOL2, "gsl_sf_cos_e");
    TEST_SF_D(s, gsl_sf_cos_e, 1073741822.5,   -0.5601305436977716102,  TEST_SQRT_TOL0, "gsl_sf_cos_e");
    TEST_SF_D(s, gsl_sf_cos_e, 1073741824.0,    0.7867071229411881196,  TEST_SQRT_TOL0, "gsl_sf_cos_e");
    
    TEST_SF_D(s, gsl_sf_cos_e, 1099511627776.0, -0.9140040719915570023, 128.0 * TEST_SQRT_TOL0, "gsl_sf_cos_e");
    

    console.log("    ... gsl_sf_sinc_e");
    TEST_SF_D(s, gsl_sf_sinc_e, 1.0/1024.0, 0.9999984312693665404, TEST_TOL0, "gsl_sf_sinc_e");
    TEST_SF_D(s, gsl_sf_sinc_e, 1.0/2.0,    2.0/M_PI,              TEST_TOL0, "gsl_sf_sinc_e");
    TEST_SF_D(s, gsl_sf_sinc_e, 80.5,       0.0039541600768172754, TEST_TOL0, "gsl_sf_sinc_e");
    TEST_SF_D(s, gsl_sf_sinc_e, 100.5,      0.0031672625490924445, TEST_TOL0, "gsl_sf_sinc_e");
    TEST_SF_D(s, gsl_sf_sinc_e, 1.0e+06 + 0.5, 3.18309727028927157e-07, TEST_TOL0, "gsl_sf_sinc_e");

    console.log("    ... gsl_sf_complex_sin_e");
    TEST_SF_2(s, gsl_sf_complex_sin_e, { x: 1.0, y: 5.0 },
              62.44551846769653403, TEST_TOL0,
              40.09216577799840254, TEST_TOL0,
              "gsl_sf_complex_sin_e");

    console.log("    ... gsl_sf_complex_cos_e");
    TEST_SF_2(s, gsl_sf_complex_cos_e, { x: 1.0, y: 5.0 },
               40.09580630629882573, TEST_TOL0,
              -62.43984868079963017, TEST_TOL0,
              "gsl_sf_complex_cos_e");

    console.log("    ... gsl_sf_complex_logsin_e");
    TEST_SF_2(s, gsl_sf_complex_logsin_e, { x: 1.0, y: 100.0 },
              99.3068528194400546900, TEST_TOL0,
              0.5707963267948966192, TEST_TOL0,
              "gsl_sf_complex_logsin_e");

    TEST_SF_2(s, gsl_sf_complex_logsin_e, { x: 1.0, y: -100.0 },
               99.3068528194400546900, TEST_TOL1,
              -0.5707963267948966192, TEST_TOL1,
              "gsl_sf_complex_logsin_e");

    TEST_SF_2(s, gsl_sf_complex_logsin_e, { x: 5.0, y: 5.0 },
              4.3068909128079757420, TEST_TOL0,
              2.8540063315538773952, TEST_TOL0,
              "gsl_sf_complex_logsin_e");

    console.log("    ... gsl_sf_lnsinh_e");
    TEST_SF_D(s,  gsl_sf_lnsinh_e, 0.1,  -2.3009189815304652235,  TEST_TOL0, "gsl_sf_lnsinh_e");
    TEST_SF_D(s,  gsl_sf_lnsinh_e, 1.0,   0.16143936157119563361, TEST_TOL0, "gsl_sf_lnsinh_e");
    TEST_SF_D(s,  gsl_sf_lnsinh_e, 5.0,   4.306807418479684201,   TEST_TOL0, "gsl_sf_lnsinh_e");
    TEST_SF_D(s,  gsl_sf_lnsinh_e, 100.0, 99.30685281944005469,   TEST_TOL0, "gsl_sf_lnsinh_e");

    console.log("    ... gsl_sf_lncosh_e");
    TEST_SF_D(s,  gsl_sf_lncosh_e, 0.125, 0.007792239318898252791, TEST_TOL0, "gsl_sf_lncosh_e");
    TEST_SF_D(s,  gsl_sf_lncosh_e, 1.0,   0.4337808304830271870,   TEST_TOL0, "gsl_sf_lncosh_e");
    TEST_SF_D(s,  gsl_sf_lncosh_e, 5.0,   4.306898218339271555, TEST_TOL0, "gsl_sf_lncosh_e");
    TEST_SF_D(s,  gsl_sf_lncosh_e, 100.0, 99.30685281944005469, TEST_TOL0, "gsl_sf_lncosh_e");

    console.log("    ... gsl_sf_polar_to_rect");
    TEST_SF_2(s, gsl_sf_polar_to_rect, { x: 10.0, y: M_PI / 6.0 },
              (10.0 * Math.sqrt(3.0) / 2.0), TEST_TOL0,
              (10.0 * 0.5), TEST_TOL0,
              "gsl_sf_polar_to_rect");

    TEST_SF_2(s, gsl_sf_polar_to_rect, { x: 10.0, y: -2.0 / 3.0 * M_PI },
              (10.0 * (-0.5)), TEST_TOL1,
              (10.0 * (-Math.sqrt(3.0) / 2.0)), TEST_TOL1,
              "gsl_sf_polar_to_rect");

    // In double precision M_PI = \pi - 1.2246467991473531772e-16,
    // i.e. the nearest machine number is slightly below the exact value
    // of \pi.  The true value of \pi satisfies
    //
    //     M_PI < \pi < nextafter(M_PI,+Inf)
    //
    // where nextafter(M_PI,+Inf) = M_PI + 2*DBL_EPSILON
    //
    // This also means that 2*M_PI is less than \pi by 2.449e-16. The
    // true value of 2\pi satisfies
    //
    //     2*M_PI < 2\pi < nextafter(2*M_PI,+Inf)
    //
    // where nextafter(2*M_PI,+Inf) = 2*M_PI + 4*DBL_EPSILON
    //
    // BJG 25/9/06

    console.log("    ... gsl_sf_angle_restrict_pos_e");
    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (2.0 * M_PI), 2.0*M_PI, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (-2.0*M_PI), 2.0*D, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (2.0*M_PI+4.0*GSL_DBL_EPSILON), 4.0*GSL_DBL_EPSILON-2.0*D, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (-2.0*M_PI-4.0*GSL_DBL_EPSILON), 2.0*M_PI-4.0*GSL_DBL_EPSILON+2.0*D, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (4.0*M_PI+8.0*GSL_DBL_EPSILON), 8.0*GSL_DBL_EPSILON-4.0*D, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (-4.0*M_PI-8.0*GSL_DBL_EPSILON), 2.0*M_PI-8.0*GSL_DBL_EPSILON+4.0*D, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (1.0e9), 0.5773954235013851694, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (1.0e12), 5.625560548042800009446, TEST_SNGL, "gsl_sf_angle_restrict_pos_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (-1.0e9), 5.7057898836782013075, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (-1.0e12), 0.6576247591367864674792517289, 100.0*TEST_SNGL, "gsl_sf_angle_restrict_pos_e" );

//  #ifdef EXTENDED -- SF.AccuracyLossException Generated
    // TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (1.0e15), 2.1096981170701125979, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );
    // TEST_SF_THETA( s, gsl_sf_angle_restrict_pos_e, (-1.0e15), 4.1734871901094738790, TEST_TOL1, "gsl_sf_angle_restrict_pos_e" );
//  #endif

    console.log("    ... gsl_sf_angle_restrict_pos_err_e");
    TEST_SF_D(s, gsl_sf_angle_restrict_pos_err_e, 2.0*M_PI, 2.0*M_PI, TEST_TOL1, "gsl_sf_angle_restrict_pos_err_e");
    TEST_SF_D(s, gsl_sf_angle_restrict_pos_err_e, -2.0*M_PI, 2.0*D, TEST_TOL1, "gsl_sf_angle_restrict_pos_err_e");

    TEST_SF_D(s, gsl_sf_angle_restrict_pos_err_e, 1.0e9, 0.5773954235013851694, TEST_TOL1, "gsl_sf_angle_restrict_pos_err_e");
    TEST_SF_D(s, gsl_sf_angle_restrict_pos_err_e, 1.0e12, 5.625560548042800009446, TEST_SNGL, "gsl_sf_angle_restrict_pos_err_e");

    TEST_SF_D(s, gsl_sf_angle_restrict_pos_err_e, -1.0e9, 5.7057898836782013075, TEST_TOL1, "gsl_sf_angle_restrict_pos_err_e");
    TEST_SF_D(s, gsl_sf_angle_restrict_pos_err_e, -1.0e12, 0.6576247591367864674792517289, 100.0*TEST_SNGL, "gsl_sf_angle_restrict_pos_err_e");

    // TEST_SF_D(s, gsl_sf_angle_restrict_pos_err_e, 1.0e15, GSL_NAN, TEST_TOL1, GSL_ELOSS, "gsl_sf_angle_restrict_pos_err_e");
    // TEST_SF_D(s, gsl_sf_angle_restrict_pos_err_e, -1.0e15, GSL_NAN, TEST_TOL1, GSL_ELOSS, "gsl_sf_angle_restrict_pos_err_e");

    console.log( "    ... gsl_sf_angle_restrict_symm_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (2.0*M_PI), -2.0*D, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (-2.0*M_PI), 2.0*D, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (M_PI), M_PI, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (-M_PI), -M_PI, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (M_PI+2.0*GSL_DBL_EPSILON), -M_PI+2.0*(GSL_DBL_EPSILON-D), TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (-M_PI-2.0*GSL_DBL_EPSILON), M_PI-2.0*(GSL_DBL_EPSILON-D), TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (3.0*M_PI+6.0*GSL_DBL_EPSILON), -M_PI+6.0*GSL_DBL_EPSILON-4.0*D, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (-3.0*M_PI-6.0*GSL_DBL_EPSILON), M_PI-6.0*GSL_DBL_EPSILON+4.0*D, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (1.0e9), 0.5773954235013851694, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (1.0e12), -0.6576247591367864674792517289, 100.0*TEST_SNGL, "gsl_sf_angle_restrict_symm_e" );

    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (-1.0e9), -0.5773954235013851694, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );
    TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (-1.0e12), 0.6576247591367864674792517289, 100.0*TEST_SNGL, "gsl_sf_angle_restrict_symm_e" );

//  #ifdef EXTENDED  Accuracy Loss Exception generated
    // TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (1e15), 2.1096981170701125979, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );
    // TEST_SF_THETA( s, gsl_sf_angle_restrict_symm_e, (-1e15), -2.1096981170701125979, TEST_TOL1, "gsl_sf_angle_restrict_symm_e" );
//  #endif

    console.log("    ... gsl_sf_angle_restrict_symm_err_e");
    TEST_SF_D(s, gsl_sf_angle_restrict_symm_err_e, 2.0*M_PI, -2.0*D, TEST_TOL1, "gsl_sf_angle_restrict_symm_err_e");
    TEST_SF_D(s, gsl_sf_angle_restrict_symm_err_e, -2.0*M_PI, 2.0*D, TEST_TOL1, "gsl_sf_angle_restrict_symm_err_e");

    TEST_SF_D(s, gsl_sf_angle_restrict_symm_err_e, 1.0e9, 0.5773954235013851694, TEST_TOL1, "gsl_sf_angle_restrict_symm_err_e");
    TEST_SF_D(s, gsl_sf_angle_restrict_symm_err_e, 1.0e12, -0.6576247591367864674792517289, 100.0*TEST_SNGL, "gsl_sf_angle_restrict_symm_err_e");

    TEST_SF_D(s, gsl_sf_angle_restrict_symm_err_e, -1.0e9, -0.5773954235013851694, TEST_TOL1, "gsl_sf_angle_restrict_symm_err_e");
    TEST_SF_D(s, gsl_sf_angle_restrict_symm_err_e, -1.0e12, 0.6576247591367864674792517289, 100.0*TEST_SNGL, "gsl_sf_angle_restrict_symm_err_e");

//    TEST_SF (s, gsl_sf_angle_restrict_symm_err_e'Access, 1e15, GSL_NAN, TEST_TOL1, GSL_ELOSS, "gsl_sf_angle_restrict_symm_err_e");
//    TEST_SF (s, gsl_sf_angle_restrict_symm_err_e'Access, -1e15, GSL_NAN, TEST_TOL1, GSL_ELOSS, "gsl_sf_angle_restrict_symm_err_e");

    console.log( "    ... gsl_sf_angle_restrict_pos_e" );
    var r = 0.0;
    theta = 5.0 * M_PI + 5.0 * D + M_PI / 2.0;
    r = gsl_sf_angle_restrict_pos_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, 3.0 / 2.0 * M_PI ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_pos_e: theta =  11/2 Pi" );
    s.Integer = s.Integer + sa.Integer;

    theta = -5.0 * M_PI - 5.0 * D - M_PI / 2.0;
    r = gsl_sf_angle_restrict_pos_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, M_PI/2.0 ) > 2.0*TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_pos_e: theta = -11/2 Pi" );
    s.Integer = s.Integer + sa.Integer;

    theta = 50000.0 + 1.0/65536.0;
    r = gsl_sf_angle_restrict_pos_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, 4.6945260308194656055 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_pos_e: theta = 50000.0 + 1.0/65536.0" );
    s.Integer = s.Integer + sa.Integer;

    theta = 5000000.0 + 1.0/65536.0;
    r = gsl_sf_angle_restrict_pos_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, 4.49537973053997376 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_pos_e: theta = 5000000.0 + 1.0/65536.0" );
    s.Integer = s.Integer + sa.Integer;

    // theta = 140737488355328.0;
    // r = gsl_sf_angle_restrict_pos_e( theta );
    // console.log( r );
    // sa.Integer = 0;
    // if ( test_sf_frac_diff( r, 3.20652300406795792638 ) > TEST_TOL0 )
    // {
    //     sa.Integer = sa.Integer + 1;
    // }
    // gsl_test( sa, "  gsl_angle_restrict_pos_e: theta = 2^47" );
    // s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_angle_restrict_symm_e" );
    theta = 5.0*M_PI + (5.5*D + M_PI/2.0);
    r = gsl_sf_angle_restrict_symm_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, -M_PI/2.0 ) > 2.0*TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_symm_e: theta =  11/2 Pi" );
    s.Integer = s.Integer + sa.Integer;

    theta = -5.0*M_PI - (5.5*D + M_PI/2.0);
    r = gsl_sf_angle_restrict_symm_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, M_PI/2.0 ) > 2.0*TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_symm_e: theta = -11/2 Pi" );
    s.Integer = s.Integer + sa.Integer;

    theta = 5.0 * M_PI + 5.0 * D - M_PI / 2.0;
    r = gsl_sf_angle_restrict_symm_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, M_PI/2.0 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_symm_e: theta = -9/2 Pi" );
    s.Integer = s.Integer + sa.Integer;

    theta = 3.0 / 2.0 * M_PI + 3.0 / 2.0 * D;
    r = gsl_sf_angle_restrict_symm_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, -M_PI/2.0 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_symm_e: theta =  3/2 Pi" );
    s.Integer = s.Integer + sa.Integer;

    theta = -3.0 / 2.0 * M_PI;
    r = gsl_sf_angle_restrict_symm_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, M_PI/2.0 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_symm_e: theta = -3/2 Pi" );
    s.Integer = s.Integer + sa.Integer;

    theta = 50000.0 + 1.0 / 65536.0;
    r = gsl_sf_angle_restrict_symm_e( theta );
    console.log( r );
    sa.Integer = 0;
    if ( test_sf_frac_diff( r, -1.5886592763601208714 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_angle_restrict_symm_e: theta = 50000.0 + 1.0/65536.0" );
    s.Integer = s.Integer + sa.Integer;

    return s;

} // test_trig

// ----------------------------------------------------------------------------

// I computed the values of zeta for s = -1e-10, 0, 1e-10 using the
// Jensen formula,
//
// zeta(s) = -1/2 + 1/(1-s)
//   + integ(sin(s arctan(t))/((1+t^2)^(s/2)(exp(2pi*t)-1)), t, 0, inf)
//
// transforming the integral from a semi-infinite range to the range
// [0,pi/2] using the substitution t = tan(u).  After Taylor expansion
// in s and numerical evaluation of the integrals this gave,
//
// zeta(s) = 1/2 + 1/(1-s)
//           + (0.0810614667944862 +/- 2e-16) s
//           + (-3.17822795429232e-3 +/- 2e-17) s^2
//           + ....
//
// for an expansion about s = 0  [BJG 7/01]
//
//
export function test_zeta( )
{
    var s = { Integer: 0 };

    console.log("Test Zeta Functions ...");

    console.log("    ... gsl_sf_zeta_int_e");
    TEST_SF_I(s, gsl_sf_zeta_int_e, (-61), -3.30660898765775767257e+34, TEST_TOL0, "gsl_sf_zeta_int_e");

    TEST_SF_I(s, gsl_sf_zeta_int_e, (-8), 0.0, TEST_TOL0, "gsl_sf_zeta_int_e");
    TEST_SF_I(s, gsl_sf_zeta_int_e, (-6), 0.0, TEST_TOL0, "gsl_sf_zeta_int_e");
    TEST_SF_I(s, gsl_sf_zeta_int_e, (-5),  -0.003968253968253968253968, TEST_TOL0, "gsl_sf_zeta_int_e");

    TEST_SF_I(s, gsl_sf_zeta_int_e, (-4), 0.0, TEST_TOL0, "gsl_sf_zeta_int_e");
    TEST_SF_I(s, gsl_sf_zeta_int_e, (-3), 1.0/120.0, TEST_TOL0, "gsl_sf_zeta_int_e");
    TEST_SF_I(s, gsl_sf_zeta_int_e, (-2), 0.0, TEST_TOL0, "gsl_sf_zeta_int_e");
    TEST_SF_I(s, gsl_sf_zeta_int_e, (-1), -1.0/12.0, TEST_TOL0, "gsl_sf_zeta_int_e");

    TEST_SF_I(s, gsl_sf_zeta_int_e, ( 5), 1.0369277551433699263313655, TEST_TOL0, "gsl_sf_zeta_int_e");
    TEST_SF_I(s, gsl_sf_zeta_int_e, (31), 1.0000000004656629065033784, TEST_TOL0, "gsl_sf_zeta_int_e");

    console.log("    ... gsl_sf_zetam1_int_e");
    TEST_SF_I(s, gsl_sf_zetam1_int_e, (-61), -3.30660898765775767257e+34, TEST_TOL0, "gsl_sf_zetam1_int_e");
    TEST_SF_I(s, gsl_sf_zetam1_int_e, (-5),  -1.003968253968253968253968, TEST_TOL0, "gsl_sf_zetam1_int_e");

    TEST_SF_I(s, gsl_sf_zetam1_int_e, (-8), -1.0, TEST_TOL0, "gsl_sf_zetam1_int_e");
    TEST_SF_I(s, gsl_sf_zetam1_int_e, (-6), -1.0, TEST_TOL0, "gsl_sf_zetam1_int_e");
    TEST_SF_I(s, gsl_sf_zetam1_int_e, (-4), -1.0, TEST_TOL0, "gsl_sf_zetam1_int_e");
    TEST_SF_I(s, gsl_sf_zetam1_int_e, (-3), -119.0/120.0, TEST_TOL0, "gsl_sf_zetam1_int_e");
    TEST_SF_I(s, gsl_sf_zetam1_int_e, (-2), -1.0, TEST_TOL0, "gsl_sf_zetam1_int_e");
    TEST_SF_I(s, gsl_sf_zetam1_int_e, (-1), -13.0/12.0, TEST_TOL0, "gsl_sf_zetam1_int_e");

    TEST_SF_I(s, gsl_sf_zetam1_int_e, ( 5), 0.0369277551433699263313655, TEST_TOL0, "gsl_sf_zetam1_int_e");
    TEST_SF_I(s, gsl_sf_zetam1_int_e, (31), 0.0000000004656629065033784, TEST_TOL0, "gsl_sf_zetam1_int_e");

    console.log("    ... gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (-151.0), 8.195215221831378294e+143, TEST_TOL2, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (-51.0), 9.68995788746359406565e+24, TEST_TOL1, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (-5.0), -0.003968253968253968253968, TEST_TOL1, "gsl_sf_zeta_e");

    TEST_SF_D(s, gsl_sf_zeta_e, (-8.0), 0.0, TEST_TOL1, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (-6.0), 0.0, TEST_TOL1, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (-4.0), 0.0, TEST_TOL1, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (-3.0), 1.0/120.0, TEST_TOL1, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (-2.0), 0.0, TEST_TOL1, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (-1.0), -1.0/12.0, TEST_TOL1, "gsl_sf_zeta_e");

    TEST_SF_D(s, gsl_sf_zeta_e, (-0.5), -0.207886224977354566017307, TEST_TOL1, "gsl_sf_zeta_e");

    TEST_SF_D(s, gsl_sf_zeta_e, (-1.0e-10), -0.49999999990810614668948, TEST_TOL1, "gsl_sf_zeta_e"); //!!!
    TEST_SF_D(s, gsl_sf_zeta_e, (0.0),    -0.5, TEST_TOL0, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (1.0e-10),  -0.50000000009189385333058, TEST_TOL0, "gsl_sf_zeta_e");

    TEST_SF_D(s, gsl_sf_zeta_e, (0.5), -1.460354508809586812889499, TEST_TOL0, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (1.0-1.0/1024.0), -1023.4228554489429787, TEST_TOL0, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (1.0+1.0/1048576.0), 1.0485765772157343441e+06, TEST_TOL0, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (5.0), 1.036927755143369926331365, TEST_TOL0, "gsl_sf_zeta_e");
    TEST_SF_D(s, gsl_sf_zeta_e, (25.5), 1.000000021074106110269959, TEST_TOL0, "gsl_sf_zeta_e");

    console.log("    ... gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (-8.0), -1.0, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (-6.0), -1.0, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (-4.0), -1.0, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (-3.0), -119.0/120.0, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (-2.0), -1.0, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (-1.0), -13.0/12.0, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (-0.5), -1.207886224977354566017307, TEST_TOL1, "gsl_sf_zetam1_e");
    // TEST_SF_D(s, gsl_sf_zetam1_e, (-1.0e-10), -1.49999999990810614668948, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (0.0),    -1.5, TEST_TOL0, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (1.0e-10),  -1.50000000009189385333058, TEST_TOL0, "gsl_sf_zetam1_e");

    TEST_SF_D(s, gsl_sf_zetam1_e, (0.5), -2.460354508809586812889499, TEST_TOL0, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (2.0),  0.64493406684822643647,     TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (3.0),  0.20205690315959428540,     TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (5.0),  0.0369277551433699263314,   TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (9.5),  0.0014125906121736622712,   TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (10.5), 0.000700842641736155219500, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (12.5), 0.000173751733643178193390, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (13.5), 0.000086686727462338155188, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (15.5), 0.000021619904246069108133, TEST_TOL1, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (16.5), 0.000010803124900178547671, TEST_TOL0, "gsl_sf_zetam1_e");
    TEST_SF_D(s, gsl_sf_zetam1_e, (25.5), 0.000000021074106110269959, TEST_TOL0, "gsl_sf_zetam1_e");

    console.log("    ... gsl_sf_hzeta_e");
    TEST_SF_DD(s, gsl_sf_hzeta_e, { x: 2.0,  y: 1.0 },  1.6449340668482264365, TEST_TOL0, "gsl_sf_hzeta_e");
    TEST_SF_DD(s, gsl_sf_hzeta_e, { x: 2.0,  y: 10.0 },  0.1051663356816857461, TEST_TOL0, "gsl_sf_hzeta_e");
    TEST_SF_DD(s, gsl_sf_hzeta_e, { x: 5.0,  y:  1.0 },  1.0369277551433699263, TEST_TOL0, "gsl_sf_hzeta_e");
    TEST_SF_DD(s, gsl_sf_hzeta_e, { x: 5.0,  y: 10.0 },  0.000030413798676470276, TEST_TOL0, "gsl_sf_hzeta_e");
    TEST_SF_DD(s, gsl_sf_hzeta_e, { x: 9.0,  y:  0.1 },  1.0000000004253980e+09, TEST_TOL0, "gsl_sf_hzeta_e");
    TEST_SF_DD(s, gsl_sf_hzeta_e, { x: 30.0, y:  0.5 },  1.0737418240000053e+09, TEST_TOL0, "gsl_sf_hzeta_e");
    TEST_SF_DD(s, gsl_sf_hzeta_e, { x: 30.0, y:  0.9 },  2.3589824880264765e+01, TEST_TOL1, "gsl_sf_hzeta_e");
    TEST_SF_DD(s, gsl_sf_hzeta_e, { x: 75.0, y:  0.25}, 1.4272476927059599e+45, TEST_TOL1, "gsl_sf_hzeta_e");

    console.log("    ... gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, (-91), -4.945598888750002040e+94, TEST_TOL0, "gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, (-51), -4.363969073121683116e+40, TEST_TOL0, "gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, (-5), 0.25, TEST_TOL0, "gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, (-1), 0.25, TEST_TOL0, "gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, ( 0), 0.5, TEST_TOL0, "gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, ( 5), 0.9721197704469093059, TEST_TOL0, "gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, ( 6), 0.9855510912974351041, TEST_TOL0, "gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, ( 20), 0.9999990466115815221, TEST_TOL0, "gsl_sf_eta_int_e");
    TEST_SF_I(s, gsl_sf_eta_int_e, ( 1000), 1.0, TEST_TOL0, "gsl_sf_eta_int_e");

    console.log("    ... gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, (-51.5), -1.2524184036924703656e+41, TEST_TOL2, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, (-5.0), 0.25, TEST_TOL0, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, (0.5), 0.6048986434216303702, TEST_TOL0, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, (0.999), 0.6929872789683383574, TEST_TOL0, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, (1.0), 0.6931471805599453094, TEST_TOL0, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, (1.0+1.0e-10), 0.6931471805759321998, TEST_TOL0, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, ( 5.0), 0.9721197704469093059, TEST_TOL0, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, ( 5.2), 0.9755278712546684682, TEST_TOL0, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, ( 6.0), 0.9855510912974351041, TEST_TOL0, "gsl_sf_eta_e");
    TEST_SF_D(s, gsl_sf_eta_e, ( 20.0), 0.9999990466115815221, TEST_TOL0, "gsl_sf_eta_e");

    return s;

} // test_zeta

// ----------------------------------------------------------------------------

// //int test_results(void)
// //{
// //  int s = 0;
// //
// //  ResultE10 re;
// //  Result r;
// //
// //  re.val = -1.0;
// //  re.err = 0.5;
// //  re.e10 = 0;
// //  Result_smash_e(&re, &r);
// //  s += ( test_sf_frac_diff(r.val, -1.0) > TEST_TOL0 );
// //  s += ( test_sf_frac_diff(r.err,  0.5) > TEST_TOL0 );
// //
// //  re.val = -1.0;
// //  re.err = 0.5;
// //  re.e10 = 10;
// //  Result_smash_e(&re, &r);
// //  s += ( test_sf_frac_diff(r.val, -1.0e+10) > TEST_TOL1 );
// //  s += ( test_sf_frac_diff(r.err,  0.5e+10) > TEST_TOL1 );
// //
// //  re.val = 1.0;
// //  re.err = 0.5;
// //  re.e10 = 10000;
// //  s += ( Result_smash_e(&re, &r) != GSL_EOVRFLW );
// //
// //  re.val = 1.0;
// //  re.err = 0.5;
// //  re.e10 = -10000;
// //  s += ( Result_smash_e(&re, &r) != GSL_EUNDRFLW );
// //
// //  return s;
// //}

// ----------------------------------------------------------------------------

export function test_gamma( )
{

    var s   = { Integer: 0 };
    var sgn = 0.0;

    console.log("Test Gamma Functions ...");

    console.log("    ... gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, -0.1,                     2.368961332728788655,     TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, -1.0 / 256.0,             5.547444766967471595,     TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, 1.0e-08,                 18.420680738180208905,     TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, 0.1,                      2.252712651734205,        TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, 1.0 + 1.0 / 256.0,       -0.0022422226599611501448, TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, 2.0 + 1.0 / 256.0,        0.0016564177556961728692, TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, 100.0,                  359.1342053695753,          TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, -1.0-1.0 / 65536.0,      11.090348438090047844,     TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, -1.0-1.0 / 268435456.0,  19.408121054103474300,     TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, -100.5,                -364.9009683094273518,       TEST_TOL0, "gsl_sf_lngamma_e");
    TEST_SF_D(s,  gsl_sf_lngamma_e, -100.0-1.0 / 65536.0,  -352.6490910117097874,       TEST_TOL0, "gsl_sf_lngamma_e");

    TEST_SF_D(s,  gsl_sf_lngamma_e, 1.5, -0.120782237635245222345518445781647212251852727902599468363868473, TEST_TOL0, "gsl_sf_lngamma_e");

    console.log("    ... gsl_sf_lngamma_sgn_e");
    TEST_SF_SGN_D(s, gsl_sf_lngamma_sgn_e, 0.7, 0.26086724653166651439, TEST_TOL1, 1.0, "gsl_sf_lngamma_sgn_e");
    TEST_SF_SGN_D(s, gsl_sf_lngamma_sgn_e, 0.1, 2.2527126517342059599, TEST_TOL0, 1.0, "gsl_sf_lngamma_sgn_e");
    TEST_SF_SGN_D(s, gsl_sf_lngamma_sgn_e, -0.1, 2.368961332728788655, TEST_TOL0, -1.0, "gsl_sf_lngamma_sgn_e");
    TEST_SF_SGN_D(s, gsl_sf_lngamma_sgn_e, -1.0-1.0/65536.0, 11.090348438090047844, TEST_TOL0, 1.0, "gsl_sf_lngamma_sgn_e");
    TEST_SF_SGN_D(s, gsl_sf_lngamma_sgn_e, -2.0-1.0/256.0, 4.848447725860607213, TEST_TOL0, -1.0, "gsl_sf_lngamma_sgn_e");
    TEST_SF_SGN_D(s, gsl_sf_lngamma_sgn_e, -2.0-1.0/65536.0, 10.397193628164674967, TEST_TOL0, -1.0, "gsl_sf_lngamma_sgn_e");
    TEST_SF_SGN_D(s, gsl_sf_lngamma_sgn_e, -3.0-1.0/8.0, 0.15431112768404182427, TEST_TOL2, 1.0, "gsl_sf_lngamma_sgn_e");
    TEST_SF_SGN_D(s, gsl_sf_lngamma_sgn_e, -100.5, -364.9009683094273518, TEST_TOL0, -1.0, "gsl_sf_lngamma_sgn_e");

    console.log("    ... gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, 1.0 + 1.0/4096.0, 0.9998591371459403421 , TEST_TOL0, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, 1.0 + 1.0/32.0, 0.9829010992836269148 , TEST_TOL0, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, 2.0 + 1.0/256.0, 1.0016577903733583299 , TEST_TOL0, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, 9.0, 40320.0                   , TEST_TOL0, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, 10.0, 362880.0                  , TEST_TOL0, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, 100.0, 9.332621544394415268e+155 , TEST_TOL2, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, 170.0, 4.269068009004705275e+304 , TEST_TOL2, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, 171.0, 7.257415615307998967e+306 , TEST_TOL2, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, -10.5, -2.640121820547716316e-07  , TEST_TOL0, "gsl_sf_gamma_e");
    TEST_SF_D(s,  gsl_sf_gamma_e, -11.25, 6.027393816261931672e-08  , TEST_TOL0, "gsl_sf_gamma_e"); // exp()... not my fault
    TEST_SF_D(s,  gsl_sf_gamma_e, -1.0+1.0/65536.0, -65536.42280587818970 , TEST_TOL0, "gsl_sf_gamma_e");

    console.log("    ... gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 1.0e-08, 3989.423555759890865  , TEST_TOL1, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 1.0e-05, 126.17168469882690233 , TEST_TOL0, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 0.001, 12.708492464364073506 , TEST_TOL0, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 1.5, 1.0563442442685598666 , TEST_TOL0, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 3.0, 1.0280645179187893045 , TEST_TOL0, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 9.0, 1.0092984264218189715 , TEST_TOL0, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 11.0, 1.0076024283104962850 , TEST_TOL0, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 100.0, 1.0008336778720121418 , TEST_TOL0, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 1.0e+05, 1.0000008333336805529 , TEST_TOL0, "gsl_sf_gammastar_e");
    TEST_SF_D(s,  gsl_sf_gammastar_e, 1.0e+20, 1.0 , TEST_TOL0, "gsl_sf_gammastar_e");

    console.log("    ... gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, 1.0, 1.0, TEST_TOL0, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, 2.0, 1.0, TEST_TOL0, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, 3.0, 0.5, TEST_TOL0, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, 4.0, 1.0/6.0, TEST_TOL0, "gsl_sf_gammainv_e");

    TEST_SF_D(s,  gsl_sf_gammainv_e, 10.0, 1.0/362880.0, TEST_TOL0, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, 100.0, 1.0715102881254669232e-156, TEST_TOL2, "gsl_sf_gammainv_e");

    TEST_SF_D(s,  gsl_sf_gammainv_e, 0.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, -1.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, -2.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, -3.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, -4.0, 0.0, TEST_TOL0, "gsl_sf_gammainv_e");

    TEST_SF_D(s,  gsl_sf_gammainv_e, -10.5, -1.0/2.640121820547716316e-07, TEST_TOL2, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, -11.25, 1.0/6.027393816261931672e-08, TEST_TOL1, "gsl_sf_gammainv_e");
    TEST_SF_D(s,  gsl_sf_gammainv_e, -1.0+1.0/65536.0, -1.0/65536.42280587818970 , TEST_TOL1, "gsl_sf_gammainv_e");

    console.log( "    ... gsl_sf_lngamma_complex_e" );
    TEST_SF_2( s, gsl_sf_lngamma_complex_e, { x: 5.0, y: 2.0 },
            2.7487017561338026749, TEST_TOL0,
            3.0738434100497007915, TEST_TOL0,
            "gsl_sf_lngamma_complex_e");

    TEST_SF_2( s, gsl_sf_lngamma_complex_e, { x: 100.0, y: 100.0 },
            315.07804459949331323, TEST_TOL1,
            2.0821801804113110099, TEST_TOL3,
            "gsl_sf_lngamma_complex_e");

    TEST_SF_2( s, gsl_sf_lngamma_complex_e, { x: 100.0, y: -1000.0 },
            -882.3920483010362817000, TEST_TOL1,
            -2.1169293725678813270, TEST_TOL3,
            "gsl_sf_lngamma_complex_e");

    TEST_SF_2( s, gsl_sf_lngamma_complex_e, { x: -100.0, y: -1.0 },
            -365.0362469529239516000, TEST_TOL1,
            -3.0393820262864361140, TEST_TOL1,
            "gsl_sf_lngamma_complex_e");

    console.log("    ... gsl_sf_taylorcoeff_e");
    TEST_SF_ID(s, gsl_sf_taylorcoeff_e, { i: 10,   x: 1.0/1048576.0 }, 1.7148961854776073928e-67  , TEST_TOL0, "gsl_sf_taylorcoeff_e");
    TEST_SF_ID(s, gsl_sf_taylorcoeff_e, { i: 10,   x: 1.0/1024.0 }, 2.1738891788497900281e-37  , TEST_TOL0, "gsl_sf_taylorcoeff_e");
    TEST_SF_ID(s, gsl_sf_taylorcoeff_e, { i: 10,   x: 1.0 }, 2.7557319223985890653e-07  , TEST_TOL0, "gsl_sf_taylorcoeff_e");
    TEST_SF_ID(s, gsl_sf_taylorcoeff_e, { i: 10,   x: 5.0 }, 2.6911444554673721340      , TEST_TOL0, "gsl_sf_taylorcoeff_e");
    TEST_SF_ID(s, gsl_sf_taylorcoeff_e, { i: 10,   x: 500.0 }, 2.6911444554673721340e+20  , TEST_TOL0, "gsl_sf_taylorcoeff_e");
    TEST_SF_ID(s, gsl_sf_taylorcoeff_e, { i: 100,  x: 100.0 }, 1.0715102881254669232e+42  , TEST_TOL1, "gsl_sf_taylorcoeff_e");
    TEST_SF_ID(s, gsl_sf_taylorcoeff_e, { i: 1000, x: 200.0 }, 2.6628790558154746898e-267 , TEST_TOL1, "gsl_sf_taylorcoeff_e");
    TEST_SF_ID(s, gsl_sf_taylorcoeff_e, { i: 1000, x: 500.0 }, 2.3193170139740855074e+131 , TEST_TOL1, "gsl_sf_taylorcoeff_e");

    console.log("    ... gsl_sf_fact_e");
    TEST_SF_D(s, gsl_sf_fact_e, 0, 1.0 , TEST_TOL0, "gsl_sf_fact_e");
    TEST_SF_D(s, gsl_sf_fact_e, 1, 1.0 , TEST_TOL0, "gsl_sf_fact_e");
    TEST_SF_D(s, gsl_sf_fact_e, 7, 5040.0 , TEST_TOL0, "gsl_sf_fact_e");
    TEST_SF_D(s, gsl_sf_fact_e, 33, 8.683317618811886496e+36 , TEST_TOL0, "gsl_sf_fact_e");

    console.log("    ... gsl_sf_doublefact_e");
    TEST_SF_D(s, gsl_sf_doublefact_e, 0, 1.0 , TEST_TOL0, "gsl_sf_doublefact_e");
    TEST_SF_D(s, gsl_sf_doublefact_e, 1, 1.0 , TEST_TOL0, "gsl_sf_doublefact_e");
    TEST_SF_D(s, gsl_sf_doublefact_e, 7, 105.0 , TEST_TOL0, "gsl_sf_doublefact_e");
    TEST_SF_D(s, gsl_sf_doublefact_e, 33, 6.332659870762850625e+18 , TEST_TOL0, "gsl_sf_doublefact_e");

    console.log("    ... gsl_sf_lnfact_e");
    TEST_SF_I(s, gsl_sf_lnfact_e, 0, 0.0 , TEST_TOL0, "gsl_sf_lnfact_e");
    TEST_SF_I(s, gsl_sf_lnfact_e, 1, 0.0 , TEST_TOL0, "gsl_sf_lnfact_e");
    TEST_SF_I(s, gsl_sf_lnfact_e, 7, 8.525161361065414300 , TEST_TOL0, "gsl_sf_lnfact_e");
    TEST_SF_I(s, gsl_sf_lnfact_e, 33, 85.05446701758151741 , TEST_TOL0, "gsl_sf_lnfact_e");

    console.log("    ... gsl_sf_lndoublefact_e");
    TEST_SF_D(s, gsl_sf_lndoublefact_e, 0, 0.0  , TEST_TOL0, "gsl_sf_lndoublefact_e");
    TEST_SF_D(s, gsl_sf_lndoublefact_e, 7, 4.653960350157523371  , TEST_TOL0, "gsl_sf_lndoublefact_e");
    TEST_SF_D(s, gsl_sf_lndoublefact_e, 33, 43.292252022541719660 , TEST_TOL0, "gsl_sf_lndoublefact_e");
    TEST_SF_D(s, gsl_sf_lndoublefact_e, 34, 45.288575519655959140 , TEST_TOL0, "gsl_sf_lndoublefact_e");
    TEST_SF_D(s, gsl_sf_lndoublefact_e, 1034, 3075.6383796271197707 , TEST_TOL0, "gsl_sf_lndoublefact_e");
    TEST_SF_D(s, gsl_sf_lndoublefact_e, 1035, 3078.8839081731809169 , TEST_TOL0, "gsl_sf_lndoublefact_e");

    console.log("    ... gsl_sf_lnchoose_e");
    TEST_SF_II(s, gsl_sf_lnchoose_e, { n: 7, m: 3 }, 3.555348061489413680 , TEST_TOL0, "gsl_sf_lnchoose_e");
    TEST_SF_II(s, gsl_sf_lnchoose_e, { n: 5, m: 2 }, 2.302585092994045684 , TEST_TOL0, "gsl_sf_lnchoose_e");

    console.log("    ... gsl_sf_choose_e");
    TEST_SF_II(s, gsl_sf_choose_e, { n: 7, m: 3 }, 35.0 , TEST_TOL0, "gsl_sf_choose_e");
    TEST_SF_II(s, gsl_sf_choose_e, { n: 7, m: 4 }, 35.0 , TEST_TOL0, "gsl_sf_choose_e");
    TEST_SF_II(s, gsl_sf_choose_e, { n: 5, m: 2 }, 10.0 , TEST_TOL0, "gsl_sf_choose_e");
    TEST_SF_II(s, gsl_sf_choose_e, { n: 5, m: 3 }, 10.0 , TEST_TOL0, "gsl_sf_choose_e");

    TEST_SF_II(s, gsl_sf_choose_e, { n: 500, m: 495 }, 255244687600.0, TEST_TOL0, "gsl_sf_choose_e");
    TEST_SF_II(s, gsl_sf_choose_e, { n: 500, m: 5 }, 255244687600.0, TEST_TOL0, "gsl_sf_choose_e");

    TEST_SF_II(s, gsl_sf_choose_e, { n: 500, m: 200 }, 5.054949849935532221e+144 , TEST_TOL5, "gsl_sf_choose_e");
    TEST_SF_II(s, gsl_sf_choose_e, { n: 500, m: 300 }, 5.054949849935532221e+144 , TEST_TOL5, "gsl_sf_choose_e");

    console.log("    ... gsl_sf_lnpoch_e");
    TEST_SF_DD(s, gsl_sf_lnpoch_e,  { x: 5.0, y: 0.0 }, 0.0, TEST_TOL0, "gsl_sf_lnpoch_e");
    TEST_SF_DD(s, gsl_sf_lnpoch_e,  { x: 5.0, y: 1.0/65536.0 }, 0.000022981557571259389129, TEST_TOL0, "gsl_sf_lnpoch_e");
    TEST_SF_DD(s, gsl_sf_lnpoch_e,  { x: 5.0, y: 1.0/256.0 },   0.005884960217985189004,    TEST_TOL2, "gsl_sf_lnpoch_e");
    TEST_SF_DD(s, gsl_sf_lnpoch_e,  { x: 7.0, y: 3.0 }, 6.222576268071368616, TEST_TOL0, "gsl_sf_lnpoch_e");
    TEST_SF_DD(s, gsl_sf_lnpoch_e,  { x: 5.0, y: 2.0 }, 3.401197381662155375, TEST_TOL0, "gsl_sf_lnpoch_e");

    console.log("    ... gsl_sf_lnpoch_sgn_e");
    TEST_SF_SGN_DD(s, gsl_sf_lnpoch_sgn_e, { x: 5.0,  y: 0.0 }, 0.0, TEST_TOL1, 1.0, "gsl_sf_lnpoch_sgn_e");
    TEST_SF_SGN_DD(s, gsl_sf_lnpoch_sgn_e, { x: -4.5, y: 0.25 }, 0.7430116475119920117, TEST_TOL1, 1.0, "gsl_sf_lnpoch_sgn_e");
    TEST_SF_SGN_DD(s, gsl_sf_lnpoch_sgn_e, { x: -4.5, y: 1.25 }, 2.1899306304483174731, TEST_TOL1, -1.0, "gsl_sf_lnpoch_sgn_e");

    console.log("    ... gsl_sf_poch_e");
    TEST_SF_DD(s,  gsl_sf_poch_e, { x: 5.0, y: 0.0 }, 1.0, TEST_TOL0, "gsl_sf_poch_e");
    TEST_SF_DD(s,  gsl_sf_poch_e, { x: 7.0, y: 3.0 }, 504.0 , TEST_TOL0, "gsl_sf_poch_e");
    TEST_SF_DD(s,  gsl_sf_poch_e, { x: 5.0, y: 2.0 }, 30.0  , TEST_TOL1, "gsl_sf_poch_e");
    TEST_SF_DD(s,  gsl_sf_poch_e, { x: 5.0, y: 1.0/256.0 }, 1.0059023106151364982 , TEST_TOL0, "gsl_sf_poch_e");

    console.log("    ... gsl_sf_pochrel_e");
    TEST_SF_DD(s,  gsl_sf_pochrel_e, { x: 5.0,  y: 0.0 }, 1.506117668431800472, TEST_TOL1, "gsl_sf_pochrel_e");
    TEST_SF_DD(s,  gsl_sf_pochrel_e, { x: 7.0,  y: 3.0 }, 503.0/3.0, TEST_TOL0, "gsl_sf_pochrel_e");
    TEST_SF_DD(s,  gsl_sf_pochrel_e, { x: 5.0,  y: 2.0 }, 29.0/2.0, TEST_TOL1, "gsl_sf_pochrel_e");
    TEST_SF_DD(s,  gsl_sf_pochrel_e, { x: 5.0,  y: 0.01 }, 1.5186393661368275330, TEST_TOL2, "gsl_sf_pochrel_e");
    TEST_SF_DD(s,  gsl_sf_pochrel_e, { x: -5.5, y: 0.01 }, 1.8584945633829063516, TEST_TOL1, "gsl_sf_pochrel_e");
    TEST_SF_DD(s,  gsl_sf_pochrel_e, { x: -5.5, y: -1.0/8.0 }, 1.0883319303552135488, TEST_TOL1, "gsl_sf_pochrel_e");
    TEST_SF_DD(s,  gsl_sf_pochrel_e, { x: -5.5, y: -1.0/256.0 }, 1.7678268037726177453, TEST_TOL1, "gsl_sf_pochrel_e");
    TEST_SF_DD(s,  gsl_sf_pochrel_e, { x: -5.5, y: -11.0 }, 0.09090909090939652475, TEST_TOL0, "gsl_sf_pochrel_e");

    console.log("    ... gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 1.0e-100, y: 0.001 }, 1.0, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 0.001, y: 0.001 }, 0.9936876467088602902, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 0.001, y: 1.0 }, 0.9997803916424144436, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 0.001, y: 10.0 }, 0.9999999958306921828, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 1.0, y: 0.001 }, 0.0009995001666250083319, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 1.0, y: 1.01 }, 0.6357810204284766802, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 1.0, y: 10.0 }, 0.9999546000702375151, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 10.0, y: 10.01 }, 0.5433207586693410570, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 10.0, y: 20.0 }, 0.9950045876916924128, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 1000.0, y: 1000.1 }, 0.5054666401440661753, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 1000.0, y: 2000.0 }, 1.0, TEST_TOL0, "gsl_sf_gamma_inc_P_e");
    // Test for failure of the Gautschi recurrence (now fixed) for x = a - 2 */
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 34.0, y: 32.0 }, 0.3849626436463866776322932129, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
    // and the next test is gamma_inc_P(37,35-20*eps) */
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 37.0, y: 3.499999999999999289e+01 }, 0.3898035054195570860969333039, TEST_TOL2, "gsl_sf_gamma_inc_P_e");

    // Regression test Martin Jansche <jansche@ling.ohio-state.edu> BUG#12 */
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 10.0, y: 1.0e-16 }, 2.755731922398588814734648067e-167, TEST_TOL2, "gsl_sf_gamma_inc_P_e");

    // Regression test for gsl_cdf_chisq_Pinv, (0.05, 1263131.0) */
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 1263131.0, y: 1261282.3637 }, 0.04994777516935182963821362168, TEST_TOL4, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 1263131.0, y: 1263131.0 }, 0.500118321758657770672882362502514254, TEST_TOL4, "gsl_sf_gamma_inc_P_e");

    console.log("    ... gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 0.0, y: 0.001 }, 0.0, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 0.001, y: 0.001 }, 0.006312353291139709793, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 0.001, y: 1.0 }, 0.00021960835758555639171, TEST_TOL1, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 0.001, y: 2.0 }, 0.00004897691783098147880, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 0.001, y: 5.0 }, 1.1509813397308608541e-06, TEST_TOL1, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1.0, y: 0.001 }, 0.9990004998333749917, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1.0, y: 1.01 }, 0.3642189795715233198, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1.0, y: 10.0 }, 0.00004539992976248485154, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 10.0, y: 10.01 }, 0.4566792413306589430, TEST_TOL0, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 10.0, y: 100.0 }, 1.1253473960842733885e-31, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1000.0, y: 1000.1 }, 0.4945333598559338247, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1000.0, y: 2000.0 }, 6.847349459614753180e-136, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");


    // designed to trap the a-x=1 problem
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 100.0, y:  99.0 }, 0.5266956696005394, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 200.0, y: 199.0 }, 0.5188414119121281, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 100.0, y:  99.0 }, 0.4733043303994607, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 200.0, y: 199.0 }, 0.4811585880878718, TEST_TOL2, "gsl_sf_gamma_inc_P_e");

    // Test for x86 cancellation problems
    TEST_SF_DD(s, gsl_sf_gamma_inc_P_e, { x: 5670.0, y: 4574.0 },  3.063972328743934e-55, TEST_TOL2, "gsl_sf_gamma_inc_P_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 5670.0, y: 4574.0 }, 1.0000000000000000, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");

    // test suggested by Michel Lespinasse [gsl-discuss Sat, 13 Nov 2004]
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1.0e+06-1.0, y: 1.0e+06-2.0 }, 0.50026596175224547004, TEST_TOL3, "gsl_sf_gamma_inc_Q_e");

    // tests in asymptotic regime related to Lespinasse test
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1.0e+06+2.0, y: 1.0e+06+1.0 }, 0.50026596135330304336, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1.0e+06, y: 1.0e+06-2.0 }, 0.50066490399940144811, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_Q_e, { x: 1.0e+07, y: 1.0e+07-2.0 }, 0.50021026104978614908, TEST_TOL2, "gsl_sf_gamma_inc_Q_e");

    // non-normalized "Q" function
    console.log("    ... gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -1.0/1048576.0, y: 1.0/1048576.0 }, 13.285819596290624271, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.001, y: 1.0/1048576.0 }, 13.381275128625328858, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -1.0,   y: 1.0/1048576.0 }, 1.0485617142715768655e+06, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.00001, y: 0.001 }, 6.3317681434563592142, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.0001, y: 0.001 }, 6.3338276439767189385, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.001, y: 0.001 }, 6.3544709102510843793, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.5,   y: 0.001 }, 59.763880515942196981, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -1.0,   y: 0.001 }, 992.66896046923884234, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -3.5,   y: 0.001 }, 9.0224404490639003706e+09, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -10.5,  y: 0.001 }, 3.0083661558184815656e+30, TEST_TOL2, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.001, y: 0.1 }, 1.8249109609418620068, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.5,   y: 0.1 }, 3.4017693366916154163, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -10.0,  y: 0.1 }, 8.9490757483586989181e+08, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -10.5,  y: 0.1 }, 2.6967403834226421766e+09, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.001, y: 1.0 }, 0.21928612679072766340, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.5,   y: 1.0 }, 0.17814771178156069019, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -1.0,   y: 1.0 }, 0.14849550677592204792, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -2.5,   y: 1.0 }, 0.096556648631275160264, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -1.0,   y: 10.0 }, 3.8302404656316087616e-07, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.001, y: 10.0 }, 4.1470562324807320961e-06, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -0.5,   y: 10.0 }, 1.2609042613241570681e-06, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -1.0,   y: 10.0 }, 3.8302404656316087616e-07, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -10.5,  y: 10.0 }, 6.8404927328441566785e-17, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -100.0, y: 10.0 }, 4.1238327669858313997e-107, TEST_TOL2, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: -200.0, y: 10.0 }, 2.1614091830529343423e-207, TEST_TOL2, "gsl_sf_gamma_inc_e");

    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   0.0,   y:   0.001 }, 6.3315393641361493320, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   0.001, y:   0.001 }, 6.3087159394864007261, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   1.0,   y:   0.001 }, 0.99900049983337499167, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:  10.0,   y:   0.001 }, 362880.0, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   0.0,   y:   1.0 }, 0.21938393439552027368, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   0.001, y:   1.0 }, 0.21948181320730279613, TEST_TOL1, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   1.0,   y:   1.0 }, 0.36787944117144232160, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:  10.0,   y:   1.0 }, 362879.95956592242045, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: 100.0,   y:   1.0 }, 9.3326215443944152682e+155, TEST_TOL0, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   0.0,   y: 100.0 }, 3.6835977616820321802e-46, TEST_TOL2, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   0.001, y: 100.0 }, 3.7006367674063550631e-46, TEST_TOL2, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:   1.0,   y: 100.0 }, 3.7200759760208359630e-44, TEST_TOL2, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x:  10.0,   y: 100.0 }, 4.0836606309106112723e-26, TEST_TOL2, "gsl_sf_gamma_inc_e");
    TEST_SF_DD(s, gsl_sf_gamma_inc_e, { x: 100.0,   y: 100.0 }, 4.5421981208626694294e+155, TEST_TOL1, "gsl_sf_gamma_inc_e");

    console.log("    ... gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 1.0e-8, y: 1.0e-8 },  19.113827924512310617 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 1.0e-8, y: 0.01 },  18.420681743788563403 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 1.0e-8, y: 1.0 },  18.420680743952365472 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 1.0e-8, y: 10.0 },  18.420680715662683009 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 1.0e-8, y: 1000.0 },  18.420680669107656949 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 0.1, y: 0.1 }, 2.9813614810376273949 , TEST_TOL1, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 0.1, y: 1.0 },  2.3025850929940456840 , TEST_TOL1, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 0.1, y: 100.0 },  1.7926462324527931217 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 0.1, y: 1000.0 },  1.5619821298353164928 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 1.0, y: 1.00025 },  -0.0002499687552073570, TEST_TOL4, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 1.0, y: 1.01 },  -0.009950330853168082848 , TEST_TOL3, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 1.0, y: 1000.0 },  -6.907755278982137052 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 100.0, y: 100.0 },  -139.66525908670663927 , TEST_TOL2, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 100.0, y: 1000.0 },  -336.4348576477366051 , TEST_TOL0, "gsl_sf_lnbeta_e");
    TEST_SF_DD(s, gsl_sf_lnbeta_e, { x: 100.0, y: 1.0e+8 },  -1482.9339185256447309 , TEST_TOL0, "gsl_sf_lnbeta_e");

    console.log("    ... gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 1.0,  y:   1.0 }, 1.0                   , TEST_TOL0, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 1.0,  y: 1.001 }, 0.9990009990009990010 , TEST_TOL0, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 1.0,  y:   5.0 }, 0.2                   , TEST_TOL1, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 1.0,  y:  100.0 }, 0.01                  , TEST_TOL1, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 10.0, y:  100.0 }, 2.3455339739604649879e-15 , TEST_TOL2, "gsl_sf_beta_e");

    // Test negative arguments
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 2.5, y: -0.1 }, -11.43621278354402041480, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 2.5, y: -1.1 }, 14.555179906328753255202, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: -0.25, y: -0.1 }, -13.238937960945229110, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: -1.25, y: -0.1 }, -14.298052997820847439, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: -100.1, y: -99.1 }, -1.005181917797644630375787297e60, TEST_TOL3, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: -100.1, y: 99.3 }, 0.0004474258199579694011200969001, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 100.1, y: -99.3 }, 1.328660939628876472028853747, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: -100.1, y: 1.2 }, 0.00365530364287960795444856281, TEST_TOL3, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 100.1, y: -1.2 }, 1203.895236907821059270698160, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: -100.1, y: -1.2 }, -3236.073671884748847700283841, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: -100.001, y: 0.0099 }, -853.946649365611147996495177, TEST_TOL4, "gsl_sf_beta_e");

    // Other test cases
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 1.0e-32, y: 1.5 }, 1.0e32, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: 1.0e-6, y: 0.5 }, 1000001.386293677092419390336, TEST_TOL2, "gsl_sf_beta_e");
    TEST_SF_DD(s, gsl_sf_beta_e, { x: -1.5, y: 0.5 }, 0.0, TEST_TOL0, "gsl_sf_beta_e");

    console.log( "    ... gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 1.0, y: 1.0, z: 0.0 }, 0.0, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 1.0, y: 1.0, z: 1.0 }, 1.0, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 0.1, y: 0.1, z: 1.0 }, 1.0, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  1.0,  y: 1.0, z: 0.5 }, 0.5, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  0.1,  y: 1.0, z: 0.5 }, 0.9330329915368074160, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 10.0,  y: 1.0, z: 0.5 }, 0.0009765625000000000000, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 50.0,  y: 1.0, z: 0.5 }, 8.881784197001252323e-16, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  1.0,  y: 0.1, z: 0.5 }, 0.06696700846319258402, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  1.0, y: 10.0, z: 0.5 }, 0.99902343750000000000, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  1.0, y: 50.0, z: 0.5 }, 0.99999999999999911180, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  1.0,  y: 1.0, z: 0.1 }, 0.10, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  1.0,  y: 2.0, z: 0.1 }, 0.19, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  1.0,  y: 2.0, z: 0.9 }, 0.99, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 50.0, y: 60.0, z: 0.5 }, 0.8309072939016694143, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 90.0, y: 90.0, z: 0.5 }, 0.5, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x:  500.0,  y: 500.0, z: 0.6 }, 0.9999999999157549630, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 5000.0, y: 5000.0, z: 0.4 }, 4.518543727260666383e-91, TEST_TOL5, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 5000.0, y: 5000.0, z: 0.6 }, 1.0, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s, gsl_sf_beta_inc_e, { x: 5000.0, y: 2000.0, z: 0.6 }, 8.445388773903332659e-89, TEST_TOL5, "gsl_sf_beta_inc_e" );

    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.1, y: -0.1, z: 1.0 }, 1.0, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.1, y: -0.2, z: 1.0 }, 1.0, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.2, y: -0.1, z: 1.0 }, 1.0, TEST_TOL2, "gsl_sf_beta_inc_e" );

    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.1, y: -0.2, z: 0.5 }, 0.675252001958389971991335, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.2, y: -0.1, z: 0.5 }, 0.324747998041610028008665, TEST_TOL2, "gsl_sf_beta_inc_e" );

    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.1, y: -0.1, z: 0.0 }, 0.0, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.1, y: -0.2, z: 0.0 }, 0.0, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.2, y: -0.1, z: 0.0 }, 0.0, TEST_TOL2, "gsl_sf_beta_inc_e" );

    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.1, y: -0.2, z: 0.3 }, 0.7469186777964287252, TEST_TOL2, "gsl_sf_beta_inc_e" );
    TEST_SF_3D( s,  gsl_sf_beta_inc_e, { x: -0.2, y: -0.1, z: 0.3 }, 0.3995299653262016818, TEST_TOL2, "gsl_sf_beta_inc_e" );

    return s;

} // test_gamma

// ****************************************************************************

export function test_dilog( )
{

    var s = { Integer: 0 };

    console.log( "Test Dilogarithm Functions ..." );

    // real dilog

    console.log( "    ... gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, -3.0,   -1.9393754207667089531,     TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, -0.5,   -0.4484142069236462024,     TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, -0.001, -0.0009997501110486510834,  TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 0.1,     0.1026177910993911,        TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 0.7,     0.8893776242860387386,     TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 1.0,     1.6449340668482260,        TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 1.5,     2.3743952702724802007,     TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 2.0,     2.4674011002723397,        TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e,  5.0,    1.7837191612666306277,     TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e,  11.0,   0.3218540439999117111,     TEST_TOL1, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 12.59,   0.0010060918167266208634,  TEST_TOL3, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 12.595,  0.00003314826006436236810, TEST_TOL5, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 13.0,   -0.07806971248458575855,    TEST_TOL2, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 20.0,   -1.2479770861745251168,     TEST_TOL2, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 150.0,  -9.270042702348657270,      TEST_TOL0, "gsl_sf_dilog_e" );
    TEST_SF_D( s, gsl_sf_dilog_e, 1100.0, -21.232504073931749553,     TEST_TOL0, "gsl_sf_dilog_e" );

    // complex dilog

    console.log( "    ... gsl_sf_complex_dilog_e" );
    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99999, y: M_PI / 2.0 },
              -0.20561329262779687646, TEST_TOL0,
               0.91595774018131512060, TEST_TOL0,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.991, y: M_PI / 2.0 },
              -0.20250384721077806127, TEST_TOL0,
               0.90888544355846447810, TEST_TOL0,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.98, y: M_PI / 2.0 },
              -0.19871638377785918403, TEST_TOL2,
               0.90020045882981847610, TEST_TOL2,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.98, y: -M_PI / 2.0 },
              -0.19871638377785918403, TEST_TOL2,
              -0.90020045882981847610, TEST_TOL2,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.95, y: M_PI / 2.0 },
              -0.18848636456893572091, TEST_TOL1,
               0.87633754133420277830, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.8, y: M_PI / 2.0 },
              -0.13980800855429037810, TEST_TOL0,
               0.75310609092419884460, TEST_TOL0,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.8, y: -M_PI / 2.0 },
              -0.13980800855429037810, TEST_TOL0,
              -0.75310609092419884460, TEST_TOL0,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.5, y: M_PI / 2.0 },
              -0.05897507442156586346, TEST_TOL1,
               0.48722235829452235710, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.5, y: -M_PI / 2.0 },
              -0.05897507442156586346, TEST_TOL1,
              -0.48722235829452235710, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.01, y: M_PI / 2.0 },
              -0.000024999375027776215378, TEST_TOL3,
               0.009999888892888684820, TEST_TOL3,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.01, y: -M_PI / 2.0 },
              -0.000024999375027776215378, TEST_TOL3,
              -0.009999888892888684820, TEST_TOL3,
               "gsl_sf_complex_dilog_e" );


    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99, y: M_PI / 4.0 },
              0.56273366219795547757, TEST_TOL3,
              0.97009284079274560384, TEST_TOL3,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99, y: -M_PI / 4.0 },
              0.56273366219795547757, TEST_TOL3,
             -0.97009284079274560384, TEST_TOL3,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99, y: 3.0 * M_PI / 4.0 },
              -0.66210902664245926235, TEST_TOL1,
               0.51995305609998319025, TEST_TOL1,
              "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99, y: 5.0 * M_PI / 4.0 },
              -0.66210902664245926235, TEST_TOL1,
              -0.51995305609998319025, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99, y: 3.0 * M_PI / 2.0 },
              -0.20215874509123277909, TEST_TOL1,
              -0.90809733095648731408, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.25, y: 3.0 * M_PI / 2.0 },
              -0.01538741178141053563, TEST_TOL1,
              -0.24830175098230686908, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.25, y: 15.0 / 8.0 * M_PI },
              0.24266162342377302235, TEST_TOL1,
             -0.10860883369274445067, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );


    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99, y: M_PI / 8.0 },
              1.0571539648820244720, TEST_TOL0,
              0.7469145254610851318, TEST_TOL0,
              "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99, y: M_PI / 64.0 },
              1.5381800285902999666, TEST_TOL0,
              0.1825271634987756651, TEST_TOL0,
              "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 0.99, y: -M_PI / 8.0 },
              1.05715396488202447202, TEST_TOL1,
             -0.74691452546108513176, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.00001, y: M_PI / 2.0 },
              -0.20562022409960237363, TEST_TOL1,
               0.91597344814458309320, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 10.0, y: M_PI / 2.0 },
              -3.0596887943287347304, TEST_TOL0,
               3.7167814930680685900, TEST_TOL0,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 100.0, y: M_PI / 2.0 },
              -11.015004738293824854, TEST_TOL0,
               7.2437843013083534970, TEST_TOL0,
               "gsl_sf_complex_dilog_e" );


    // tests brought up by Jim McElwaine bug report

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.1, y: -M_PI / 2.0 },
              -0.24099184177382733037, TEST_TOL1,
              -0.99309132538137822631, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.1, y: 3.0 * M_PI / 2.0 },
              -0.24099184177382733037, TEST_TOL1,
              -0.99309132538137822631, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.1, y: -3.0 * M_PI / 2.0 },
              -0.24099184177382733037, TEST_TOL1,
               0.99309132538137822631, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.1, y: -M_PI - 0.25 * M_PI },
              -0.72908565537087935118, TEST_TOL1,
               0.56225783937234862649, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.1, y: M_PI + 0.25 * M_PI },
              -0.72908565537087935118, TEST_TOL1,
              -0.56225783937234862649, TEST_TOL1,
               "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.1, y: -M_PI / 128.0 },
              1.8881719454909716580, TEST_TOL1,
             -0.3556738764969238976, TEST_TOL1,
              "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.1, y:  M_PI / 128.0 },
              1.8881719454909716580, TEST_TOL1,
              0.3556738764969238976, TEST_TOL1,
              "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.5, y:  M_PI / 8.0 },
              1.3498525763442498343, TEST_TOL1,
              1.4976532712229749493, TEST_TOL1,
              "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.5, y: -M_PI / 8.0 },
              1.3498525763442498343, TEST_TOL1,
             -1.4976532712229749493, TEST_TOL1,
              "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.5, y: 2.0 * M_PI + M_PI / 8.0 },
              1.3498525763442498343, TEST_TOL1,
              1.4976532712229749493, TEST_TOL1,
              "gsl_sf_complex_dilog_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_e, { x: 1.5, y: 2.0 * M_PI - M_PI / 8.0 },
              1.3498525763442498343, TEST_TOL1,
             -1.4976532712229749493, TEST_TOL1,
              "gsl_sf_complex_dilog_e" );

    // tests of the (x,y) function, which is now the underlying implementation

    console.log( "    ... gsl_sf_complex_dilog_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_xy_e, { x: 0.0, y: 0.5 },
              -0.05897507442156586346, TEST_TOL1,
               0.48722235829452235710, TEST_TOL1,
               "gsl_sf_complex_dilog_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_xy_e, { x: 0.0, y: -0.5 },
              -0.05897507442156586346, TEST_TOL1,
              -0.48722235829452235710, TEST_TOL1,
               "gsl_sf_complex_dilog_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_xy_e, { x: 0.91464073718617389108, y: 0.37885659804143889673 },
              1.0571539648820244720, TEST_TOL0,
              0.7469145254610851318, TEST_TOL0,
              "gsl_sf_complex_dilog_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_xy_e, { x: 0.91464073718617389108, y: -0.37885659804143889673 },
              1.05715396488202447202, TEST_TOL1,
             -0.74691452546108513176, TEST_TOL1,
              "gsl_sf_complex_dilog_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_xy_e, { x: -1.5, y: 0.0 },
             -1.1473806603755707541, TEST_TOL1,
              0.0, TEST_TOL1,
              "gsl_sf_complex_dilog_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_xy_e, { x: 0.5, y: 0.0 },
              0.58224052646501250590, TEST_TOL1,
              0.0, TEST_TOL1,
              "gsl_sf_complex_dilog_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_dilog_xy_e, { x: 1.5, y: 0.0 },
              2.3743952702724802007, TEST_TOL1,
             -1.2738062049196005309, TEST_TOL1,
              "gsl_sf_complex_dilog_xy_e" );

    // small set of spence tests, mostly to check the value on the cut

    console.log( "    ... gsl_sf_complex_spence_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_spence_xy_e, { x: 1.5, y: 0.0 },
             -0.44841420692364620244, TEST_TOL1,
              0.0, TEST_TOL1,
              "gsl_sf_complex_spence_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_spence_xy_e, { x: 0.5, y: 0.0 },
              0.58224052646501250590, TEST_TOL1,
              0.0, TEST_TOL1,
              "gsl_sf_complex_spence_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_spence_xy_e, { x: 0.0, y: 0.0 },
              1.6449340668482264365, TEST_TOL1,
              0.0, TEST_TOL1,
              "gsl_sf_complex_spence_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_spence_xy_e, { x: -0.5, y: 0.0 },
              2.3743952702724802007, TEST_TOL1,
             -1.2738062049196005309, TEST_TOL1,
              "gsl_sf_complex_spence_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_spence_xy_e, { x: -0.5, y: 1.0 / 1024.0 },
              2.3723507455234125018, TEST_TOL1,
             -1.2742581376517839070, TEST_TOL1,
              "gsl_sf_complex_spence_xy_e" );

    TEST_SF_2( s, gsl_sf_complex_spence_xy_e, { x: -0.5, y: -1.0 / 1024.0 },
              2.3723507455234125018, TEST_TOL1,
              1.2742581376517839070, TEST_TOL1,
              "gsl_sf_complex_spence_xy_e" );

    return s;

} // test_dilog

// ****************************************************************************

export function test_legendre( )
{

    var s  = { Integer: 0 };
    var sa = { Integer: 0 };

    var L  = []; // ArrayResult(0..255);
    var DL = []; // ArrayResult(0..255);

    console.log("Test Legendre Functions ...");

    console.log("    ... gsl_sf_legendre_P1_e");
    TEST_SF_D(s,  gsl_sf_legendre_P1_e, -0.5, -0.5, TEST_TOL0, "gsl_sf_legendre_P1_e");
    TEST_SF_D(s,  gsl_sf_legendre_P1_e,  0.5, 0.5, TEST_TOL0, "gsl_sf_legendre_P1_e");

    console.log("    ... gsl_sf_legendre_P2_e");
    TEST_SF_D(s,  gsl_sf_legendre_P2_e, 0.0, -0.5  , TEST_TOL0, "gsl_sf_legendre_P2_e");
    TEST_SF_D(s,  gsl_sf_legendre_P2_e, 0.5, -0.125, TEST_TOL0, "gsl_sf_legendre_P2_e");
    TEST_SF_D(s,  gsl_sf_legendre_P2_e, 1.0, 1.0  , TEST_TOL0, "gsl_sf_legendre_P2_e");
    TEST_SF_D(s,  gsl_sf_legendre_P2_e, 100.0, 14999.5  , TEST_TOL0, "gsl_sf_legendre_P2_e");

    console.log("    ... gsl_sf_legendre_P3_e");
    TEST_SF_D(s,  gsl_sf_legendre_P3_e,  -0.5, 0.4375, TEST_TOL0, "gsl_sf_legendre_P3_e");
    TEST_SF_D(s,  gsl_sf_legendre_P3_e,   0.5, -0.4375, TEST_TOL0, "gsl_sf_legendre_P3_e");
    TEST_SF_D(s,  gsl_sf_legendre_P3_e,   1.0, 1.0        , TEST_TOL0, "gsl_sf_legendre_P3_e");
    TEST_SF_D(s,  gsl_sf_legendre_P3_e, 100.0, 2.49985e+06, TEST_TOL0, "gsl_sf_legendre_P3_e");

    console.log("    ... gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 1,   x: -0.5 }, -0.5, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 1,   x: 1.0e-8 }, 1.0e-08, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 1,   x: 0.5 }, 0.5, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 1,   x: 1.0 }, 1.0, TEST_TOL0, "gsl_sf_legendre_Pl_e");

    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 10,  x: -0.5 }, -0.1882286071777345, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 10,  x:  1.0e-8 }, -0.24609374999999864648, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 10,  x:  0.5 }, -0.18822860717773437500, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 10,  x:  1.0 }, 1.0, TEST_TOL0, "gsl_sf_legendre_Pl_e");

    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 99,  x: -0.5 }, 0.08300778172138770477, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 99,  x:  1.0e-8 }, -7.958923738716563193e-08, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 99,  x:  0.5 }, -0.08300778172138770477, TEST_TOL0, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 99,  x:  0.999 }, -0.3317727359254778874, TEST_TOL2, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 99,  x:  1.0 }, 1.0, TEST_TOL0, "gsl_sf_legendre_Pl_e");

    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 1000, x: -0.5 },   -0.019168251091650277878, TEST_TOL2, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 1000, x:  1.0e-8 }, 0.0252250181770982897470252620,  TEST_TOL2, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 1000, x:  0.5 },   -0.019168251091650277878, TEST_TOL2, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 1000, x:  1.0 },    1.0,                     TEST_TOL0, "gsl_sf_legendre_Pl_e");

    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 4000, x: -0.5 }, -0.009585404456573080972, TEST_TOL2, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 4000, x:  0.5 }, -0.009585404456573080972, TEST_TOL2, "gsl_sf_legendre_Pl_e");
    TEST_SF_ID(s, gsl_sf_legendre_Pl_e, { i: 4000, x:  1.0 }, 1.0, TEST_TOL0, "gsl_sf_legendre_Pl_e");

    console.log("    ... gsl_sf_legendre_Pl_array");
    sa.Integer = 0;
    gsl_sf_legendre_Pl_array(100, 0.5, L);
    sa.Integer = sa.Integer + test_sf_val(L[0],    1.0, TEST_TOL1, "gsl_sf_legendre_Pl_array");
    sa.Integer = sa.Integer + test_sf_val(L[10],  -0.18822860717773437500, TEST_TOL1, "gsl_sf_legendre_Pl_array");
    sa.Integer = sa.Integer + test_sf_val(L[100], -0.06051802596186118687, TEST_TOL1, "gsl_sf_legendre_Pl_array");
    gsl_test(sa, "gsl_sf_legendre_Pl_array(100, 0.5)");
    s.Integer = s.Integer + sa.Integer;

    console.log("    ... gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = 0;
    gsl_sf_legendre_Pl_deriv_array(100, 0.5, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],    0.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],    1.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10],  -2.3171234130859375000, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[100], -7.0331691653942815112, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(100, 0.5)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_Pl_deriv_array(10, 1.0, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],   0.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],   1.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], 55.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(10, 1.0)");
    s.Integer = s.Integer + sa.Integer;

    gsl_sf_legendre_Pl_deriv_array(10, 1.0 - 1.0e-11, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],   0.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],   1.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], 54.999999985150000001, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(10, 1.0 - 1.0e-11)");
    s.Integer = s.Integer + sa.Integer;

    gsl_sf_legendre_Pl_deriv_array(10, -1.0, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],    0.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],    1.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], -55.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(10, -1.0)");
    s.Integer = s.Integer + sa.Integer;

    gsl_sf_legendre_Pl_deriv_array(10, -1.0 + 1.0e-11, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],    0.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],    1.0, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], -54.999999985150000001, TEST_TOL1, "gsl_sf_legendre_Pl_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(10, -1.0 + 1.0e-11)");
    s.Integer = s.Integer + sa.Integer;

    console.log("    ... gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 0, x: -0.5 }, -0.18822860717773437500, TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 0, x: 1.0e-08 }, -0.24609374999999864648, TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 0, x: 0.5 }, -0.18822860717773437500, TEST_TOL0, "gsl_sf_legendre_Plm_e");

    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 1, x: -0.5 }, -2.0066877394361256516, TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 1, x: 1.0e-08 }, -2.7070312499999951725e-07, TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 1, x: 0.5 }, 2.0066877394361256516, TEST_TOL0, "gsl_sf_legendre_Plm_e");

    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 5, x: -0.5 },    -30086.169706116174977,    TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 5, x: 1.0e-08 }, -0.0025337812499999964949, TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 5, x: 0.5 },      30086.169706116174977,    TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 10,  m: 5, x: 0.999 },   -0.5036411489013270406,    TEST_TOL1, "gsl_sf_legendre_Plm_e");

    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 100, m: 5, x: -0.5 }, -6.617107444248382171e+08, TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 100, m: 5, x: 1.0e-08 }, 817.8987598063712851, TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 100, m: 5, x: 0.5 }, 6.617107444248382171e+08, TEST_TOL0, "gsl_sf_legendre_Plm_e");
    TEST_SF_IID(s, gsl_sf_legendre_Plm_e, { n: 100, m: 5, x: 0.999 }, -1.9831610803806212189e+09, TEST_TOL2, "gsl_sf_legendre_Plm_e");

    console.log("    ... gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = 0;
    gsl_sf_legendre_Plm_deriv_array(100, 2, -1.0 + 1.0/1125899906842624.0, L, DL);
    sa.Integer = sa.Integer + test_sf_val(L[0],   5.3290705182007490275e-15, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(L[1],  -2.6645352591003721471e-14, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(L[98],  2.2646284847349109694e-08, TEST_TOL2, "gsl_sf_legendre_Plm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, -1.0 + 2^(-50)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_Plm_deriv_array(100, 2, 1.0 - 1.0/1125899906842624.0, L, DL);
    sa.Integer = sa.Integer + test_sf_val(L[0],  5.3290705182007490275e-15, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(L[1],  2.6645352591003721471e-14, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(L[10], 5.3343995887188313290e-12, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(L[98], 2.2646284847349109694e-08, TEST_TOL2, "gsl_sf_legendre_Plm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, 1.0 - 2^(-50)");
    s.Integer = s.Integer + sa.Integer;

    console.log("    ... gsl_sf_legendre_Plm_array");
    sa.Integer = 0;
    gsl_sf_legendre_Plm_array(100, 5, 0.5, L);
    sa.Integer = sa.Integer + test_sf_val(L[0],  -460.3466286991656682, TEST_TOL1, "gsl_sf_legendre_Plm_array");
    sa.Integer = sa.Integer + test_sf_val(L[10],  38852.51334152290535, TEST_TOL1, "gsl_sf_legendre_Plm_array");
    sa.Integer = sa.Integer + test_sf_val(L[95],  6.617107444248382171e+08, TEST_TOL1, "gsl_sf_legendre_Plm_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_array(100, 5, 0.5)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_Plm_array(100, 5, 0.999, L);
    sa.Integer = sa.Integer + test_sf_val(L[0],  -0.00016883550990916552255, TEST_TOL2, "gsl_sf_legendre_Plm_array");
    sa.Integer = sa.Integer + test_sf_val(L[10], -30.651334850159821525, TEST_TOL2, "gsl_sf_legendre_Plm_array");
    sa.Integer = sa.Integer + test_sf_val(L[95], -1.9831610803806212189e+09, TEST_TOL2, "gsl_sf_legendre_Plm_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_array(100, 5, 0.999)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_Plm_array(100, 5, -0.999, L);
    sa.Integer = sa.Integer + test_sf_val(L[0],  -0.00016883550990916552255, TEST_TOL2, "gsl_sf_legendre_Plm_array");
    sa.Integer = sa.Integer + test_sf_val(L[10], -30.651334850159821525, TEST_TOL2, "gsl_sf_legendre_Plm_array");
    sa.Integer = sa.Integer + test_sf_val(L[95],  1.9831610803806212189e+09, TEST_TOL2, "gsl_sf_legendre_Plm_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_array(100, 5, -0.999)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_Plm_deriv_array(100, 2, 0.999, L, DL);
    sa.Integer = sa.Integer + test_sf_val(L[0],   0.00599700000000000000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(L[1],   0.02995501500000000000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[0],  -5.9940000000000000000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],  -29.910045000000000000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[2],  -89.490629790000000000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], -5703.9461633355291972, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[95], 6.4518473603456858414E+06, TEST_TOL3, "gsl_sf_legendre_Plm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, 0.999)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_Plm_deriv_array(100, 2, 1.0 - 1.0e-15, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],  -5.9999999999999940000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],  -29.999999999999910000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[2],  -89.999999999999490000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], -6005.9999999996936940, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[95], -2.2586255999928454270e+07, TEST_TOL3, "gsl_sf_legendre_Plm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, 1.0 - 1.0e-15)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_Plm_deriv_array(100, 2, -1.0 + 1.0e-15, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],   5.9999999999999940000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],  -29.999999999999910000, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[95], -2.2586255999928454270e+07, TEST_TOL3, "gsl_sf_legendre_Plm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, -1.0 + 1.0e-15)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_Plm_deriv_array(100, 5, 0.999, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],  0.42187762481054616565, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],  4.6341560284340909936, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[2],  27.759505566959219127, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], 76051.795860179545484, TEST_TOL1, "gsl_sf_legendre_Plm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[95], 3.0344503083851936814e+12, TEST_TOL3, "gsl_sf_legendre_Plm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 5, 0.999)");
    s.Integer = s.Integer + sa.Integer;

    console.log("    ... gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 0,  x: -0.5 }, -0.24332702369300133776, TEST_TOL0, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 0,  x: 0.5 }, -0.24332702369300133776, TEST_TOL0, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 0,  x: 0.999 }, 1.2225754122797385990, TEST_TOL1, "gsl_sf_legendre_sphPlm_e");

    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 5,  x: -0.5 },    -0.3725739049803293972,     TEST_TOL0, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 5,  x: 1.0e-08 }, -3.1377233589376792243e-08, TEST_TOL0, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 5,  x: 0.5 },      0.3725739049803293972,     TEST_TOL0, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 5,  x: 0.999 },   -6.236870674727370094e-06,  TEST_TOL2, "gsl_sf_legendre_sphPlm_e");

    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 10, x:  -0.5 }, 0.12876871185785724117, TEST_TOL1, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 10, x:  0.5 }, 0.12876871185785724117,  TEST_TOL1, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 10,  m: 10, x:  0.999 }, 1.7320802307583118647e-14, TEST_TOL2, "gsl_sf_legendre_sphPlm_e");

    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 200, m:  1, x:  -0.5 },   0.3302975570099492931, TEST_TOL1, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 200, m:  1, x:  0.5 },   -0.3302975570099492931, TEST_TOL1, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 200, m:  1, x:  0.999 }, -1.4069792055546256912, TEST_TOL2, "gsl_sf_legendre_sphPlm_e");

    // Test case from alberto@physik.fu-berlin.de

    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 3, m: 1, x: 0.0 }, 0.323180184114150653007, TEST_TOL2, "gsl_sf_legendre_sphPlm_e");

    // Other test cases

    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 200, m: 1, x: -0.5 }, 0.3302975570099492931418227583, TEST_TOL2, "gsl_sf_legendre_sphPlm_e");
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 140, m: 135, x: 1.0 }, 0.0, TEST_TOL2, "gsl_sf_legendre_sphPlm_e");

// --#ifdef EXTENDED
    TEST_SF_IID(s, gsl_sf_legendre_sphPlm_e, { n: 140, m: 135, x: 0.99998689456491752 }, -6.54265253269093276310395668335e-305, TEST_TOL6, "gsl_sf_legendre_sphPlm_e");
// --#endif

    console.log("    ... gsl_sf_legendre_sphPlm_array");
    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_array(100, 5, 0.5, L);
    sa.Integer = sa.Integer + test_sf_val(L[0],  -0.22609703187800460722, TEST_TOL1, "gsl_sf_legendre_sphPlm_array");
    sa.Integer = sa.Integer + test_sf_val(L[10],  0.07452710323813558940, TEST_TOL1, "gsl_sf_legendre_sphPlm_array");
    sa.Integer = sa.Integer + test_sf_val(L[95],  0.25865355990880161717, TEST_TOL1, "gsl_sf_legendre_sphPlm_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_array(100, 5, 0.5)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_array(100, 2, 1.0 - 1.0/1125899906842624.0, L);
    sa.Integer = sa.Integer + test_sf_val(L[0],  6.8616082064776657177e-16, TEST_TOL2, "gsl_sf_legendre_sphPlm_array");
    sa.Integer = sa.Integer + test_sf_val(L[10], 4.8543150313086787324e-14, TEST_TOL2, "gsl_sf_legendre_sphPlm_array");
    sa.Integer = sa.Integer + test_sf_val(L[95], 8.3138984963650838973e-12, TEST_TOL2, "gsl_sf_legendre_sphPlm_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_array(100, 2, 1.0 - 2^(-50))");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_array(100, 2, -1.0 + 1.0/1125899906842624.0, L);
    sa.Integer = sa.Integer + test_sf_val(L[0],   6.8616082064776657177e-16, TEST_TOL2, "gsl_sf_legendre_sphPlm_array");
    sa.Integer = sa.Integer + test_sf_val(L[95], -8.3138984963650838973e-12, TEST_TOL2, "gsl_sf_legendre_sphPlm_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_array(100, 2, -1.0 + 2^(-50))");
    s.Integer = s.Integer + sa.Integer;

    console.log("    ... gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_deriv_array(100, 0, 0.5, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],   0.0, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], -2.9953934850252897591, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[95], -36.411811015111761007, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 0, 0.5)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_deriv_array(100, 1, 0.5, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],   0.19947114020071633897, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],  -0.44603102903819277863, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10],  1.3658895325030216565, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[99], -27.925571865639037118, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 1, 0.5)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_deriv_array(100, 1, 1.0 - 1.0/1125899906842624.0, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],  8.1973898803378530946e+06, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],  1.8329921010504257405e+07, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], 1.8439572562895384115e+08, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[99], 4.7682463136232210552e+09, TEST_TOL3, "gsl_sf_legendre_sphPlm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 1, 1.0 - 2^(-50))");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_deriv_array(100, 2, 0.5, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],  -0.38627420202318958034, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],   0.25549636910832059085, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[2],   1.5053547230039006279, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10],  0.73576559668648243477, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[98], 28.444589950264378407, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 2, 0.5)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_deriv_array(100, 5, 0.5, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],   0.75365677292668202407, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],   0.54346962777757450534, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[2],  -0.98309969029001383773, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[3],  -2.7728270988954534293, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[10], -5.7407133315443482193, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[95], -25.893934624747394561, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 5, 0.5)");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_deriv_array(100, 5, 1.0 - 1.0/1125899906842624.0, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],  1.7374288379067753301e-22, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],  6.2643887625426827113e-22, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[2],  1.6482697200734667281e-21, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[95], 3.9890549466071349506e-15, TEST_TOL2, "gsl_sf_legendre_sphPlm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 5, 1.0 - 2^(-50))");
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_legendre_sphPlm_deriv_array(100, 5, -1.0 + 1.0/1125899906842624.0, L, DL);
    sa.Integer = sa.Integer + test_sf_val(DL[0],  -1.7374288379067753301e-22, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[1],   6.2643887625426827113e-22, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[2],  -1.6482697200734667281e-21, TEST_TOL1, "gsl_sf_legendre_sphPlm_deriv_array");
    sa.Integer = sa.Integer + test_sf_val(DL[95],  3.9890549466071349506e-15, TEST_TOL3, "gsl_sf_legendre_sphPlm_deriv_array");
    gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 5, -1.0 + 2^(-50))");
    s.Integer = s.Integer + sa.Integer;

    console.log("    ... gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 0.0,   y: -0.5 },   0.8573827581049917129, TEST_TOL0, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 0.0,   y: 0.5 },   0.8573827581049917129, TEST_TOL0, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 0.0,   y: 2.0 },   0.6062611623284649811, TEST_TOL0, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 0.0,   y: 100.0 }, 0.07979045091636735635, TEST_TOL0, "gsl_sf_conicalP_half_e");

    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 10.0,  y: -0.5 },    5.345484922591867188e+08, TEST_TOL1, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 10.0,  y:  0.5 },    15137.910380385258370, TEST_TOL1, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 10.0,  y:  2.0 },    0.4992680691891618544, TEST_TOL1, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 10.0,  y:  100.0 }, -0.07272008163718195685, TEST_TOL2, "gsl_sf_conicalP_half_e");

    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 200.0, y: -1.0e-3 },  1.3347639529084185010e+136, TEST_TOL2, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 200.0, y:  1.0e-8 },  1.0928098010940058507e+136, TEST_TOL2, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 200.0, y:  0.5 },     3.895546021611205442e+90,   TEST_TOL2, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 200.0, y:  10.0 },   -0.04308567180833581268,     TEST_TOL3, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 200.0, y:  100.0 },  -0.04694669186576399194,     TEST_TOL3, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 200.0, y:  1000.0 },  0.023698140704121273277,    TEST_TOL3, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 200.0, y:  1.0e+8 }, -0.00006790983312124277891,  TEST_TOL3, "gsl_sf_conicalP_half_e");

    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 1.0e+8, y: 1.1 },   1.1599311133054742944,  TEST_SQRT_TOL0, "gsl_sf_conicalP_half_e");
    TEST_SF_DD(s, gsl_sf_conicalP_half_e, { x: 1.0e+8, y: 100.0 }, 0.07971967557381557875, TEST_SQRT_TOL0, "gsl_sf_conicalP_half_e");

    console.log("    ... gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 0.0,   y: -0.5 },  1.7956982494514644808, TEST_TOL0, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 0.0,   y: 0.5 },  0.8978491247257322404, TEST_TOL0, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 0.0,   y: 2.0 },  0.7984204253272901551, TEST_TOL0, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 0.0,   y: 100.0 },  0.4227531369388072584, TEST_TOL0, "gsl_sf_conicalP_mhalf_e");

    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 10.0,  y: -0.5 },  5.345484922591867181e+07, TEST_TOL1, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 10.0,  y:  0.5 },  1513.7910356104985334, TEST_TOL1, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 10.0,  y:  2.0 },  0.03439243987215615642, TEST_TOL1, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 10.0,  y:  100.0 },  0.003283756665952609624, TEST_TOL2, "gsl_sf_conicalP_mhalf_e");

    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 200.0, y: -0.5 },  1.7699538115312304280e+179, TEST_TOL2, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 200.0, y:  1.0e-8 },  5.464049005470029253e+133, TEST_TOL2, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 200.0, y:  0.5 },  1.9477730108056027211e+88, TEST_TOL2, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 200.0, y:  10.0 },  0.0012462575917716355362, TEST_TOL2, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 200.0, y:  100.0 },  -0.0003225881344802625149, TEST_TOL2, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 200.0, y:  1000.0 }, -0.00004330652890886567623, TEST_TOL3, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 200.0, y:  1.0e+8 },  2.0943091278037078483e-07, TEST_TOL3, "gsl_sf_conicalP_mhalf_e");

    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 1.0e+8, y: 1.1 }, 2.092320445620989618e-09, 16.0*TEST_SQRT_TOL0, "gsl_sf_conicalP_mhalf_e");
    TEST_SF_DD(s, gsl_sf_conicalP_mhalf_e, { x: 1.0e+8, y: 100.0 },  -3.359967833599016923e-11, 256.0*TEST_SQRT_TOL0, "gsl_sf_conicalP_mhalf_e");

    console.log("    ... gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 0.0, y:-0.5 },  1.3728805006183501647, TEST_TOL0, "gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 0.0, y: 0.5 },  1.0731820071493643751, TEST_TOL0, "gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 0.0, y: 2.0 },  0.9012862993604472987, TEST_TOL0, "gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 0.0, y: 100.0 },  0.30091748588199264556, TEST_TOL0, "gsl_sf_conicalP_0_e");

    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 10.0, y:-0.5 },  1.6795592815421804669e+08, TEST_TOL1, "gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 10.0, y: 0.5 },  4826.034132009618240,      TEST_TOL1, "gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 10.0, y: 2.0 },  0.18798468917758716146,    TEST_TOL2, "gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 10.0, y: 100.0 }, -0.008622130749987962529, TEST_TOL2, "gsl_sf_conicalP_0_e");

    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 200.0, y: -0.5 }, 2.502194818646823e+180, TEST_TOL4, "gsl_sf_conicalP_0_e");

    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 1000.0, y: 100.0 },   0.0017908817653497715844, TEST_TOL3, "gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 1000.0, y: 1000.0 }, -0.0006566893804926284301, TEST_TOL3, "gsl_sf_conicalP_0_e");
    TEST_SF_DD(s, gsl_sf_conicalP_0_e, { x: 1000.0, y: 1.0e+8 },  2.3167213561756390068e-06, TEST_TOL4, "gsl_sf_conicalP_0_e");

    console.log( "    ... gsl_sf_conicalP_1_e" );
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 0.0, y: -0.5 },    0.4939371126656998499,  TEST_TOL1, "gsl_sf_conicalP_1_e" );
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 0.0, y:  0.5 },    0.14933621085538265636, TEST_TOL1, "gsl_sf_conicalP_1_e" );
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 0.0, y:  2.0 },   -0.13666874968871549533, TEST_TOL1, "gsl_sf_conicalP_1_e" );
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 0.0, y:  100.0 }, -0.10544528203156629098, TEST_TOL2, "gsl_sf_conicalP_1_e" );

    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 10.0, y: -0.5 },    1.7253802958788312520e+09, TEST_TOL2, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 10.0, y:  0.5 },    46781.02294059967988,      TEST_TOL1, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 10.0, y:  2.0 },    0.26613342643657444400,    TEST_TOL2, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 10.0, y:  100.0 }, -0.23281959695501029796,    TEST_TOL2, "gsl_sf_conicalP_1_e");

    // FIXME: Mathematica gets some brain-damaged numbers for
    // these x < 0 points. I have checked what I am doing in detail,
    // and it must be right because you can do it by summing
    // manifestly positive definite quantities.

    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 200.0, y: -0.999 }, 2.71635193199341135e+270, TEST_TOL2, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 200.0, y: -0.9 },   4.2952493176812905e+234,  TEST_TOL2, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 200.0, y: -0.5 },   5.01159205956053439e+182, TEST_TOL3, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 200.0, y:  0.999 }, 195733.0396081538,        TEST_TOL2, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 200.0, y:  10.0 }, -2.9272610662414349553,    TEST_TOL2, "gsl_sf_conicalP_1_e");

    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 1000.0, y: 100.0  },  -1.7783258105862399857,    TEST_TOL6, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 1000.0, y: 1000.0 },  0.4535161075156427179,    TEST_TOL4, "gsl_sf_conicalP_1_e");
    TEST_SF_DD( s, gsl_sf_conicalP_1_e, { x: 1000.0, y: 1.0e+8 },  0.0009983414549874888478, TEST_SQRT_TOL0, "gsl_sf_conicalP_1_e");

    console.log( "    ... gsl_sf_conicalP_sph_reg_e" );
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 2,  x: 1.0,   y: -0.5 },  1.6406279287008789526,      TEST_TOL0, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 10, x: 1.0,   y: -0.5 },  0.000029315266725049129448, TEST_TOL1, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 20, x: 1.0,   y: -0.5 },  7.335769429462034431e-15,   TEST_TOL1, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 30, x: 1.0,   y: -0.5 },  1.3235612394267378871e-26,  TEST_TOL2, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 10, x: 1.0,   y: 0.5 },  2.7016087199857873954e-10, TEST_TOL1, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 20, x: 1.0,   y: 0.5 },  1.1782569701435933399e-24, TEST_TOL1, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 30, x: 1.0,   y: 0.5 },  3.636240588303797919e-41,  TEST_TOL1, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 10, x: 1.0,   y: 2.0 },  2.4934929626284934483e-10, TEST_TOL1, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 20, x: 1.0,   y: 2.0 },  1.1284762488012616191e-24, TEST_TOL2, "gsl_sf_conicalP_sph_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_sph_reg_e, { n: 30, x: 100.0, y: 100.0 },  -1.6757772087159526048e-64, TEST_TOL6, "gsl_sf_conicalP_sph_reg_e");

    console.log( "    ... gsl_sf_conicalP_cyl_reg_e" );
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 2,  x: 1.0,    y: -0.5 },   2.2048510472375258708,       TEST_TOL0, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 10, x:  1.0,   y: -0.5 },  0.00007335034531618655690,   TEST_TOL1, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 20, x:  1.0,   y: -0.5 },  2.5419860619212164696e-14,   TEST_TOL1, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 30, x:  1.0,   y: -0.5 },  5.579714972260536827e-26,    TEST_TOL2, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 10, x:  1.0,   y:  0.5 },  1.1674078819646475282e-09,    TEST_TOL0, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 20, x:  1.0,   y:  0.5 },  7.066408031229072207e-24,     TEST_TOL1, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 30, x:  1.0,   y:  0.5 },  2.6541973286862588488e-40,    TEST_TOL1, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 10, x:  1.0,   y:  2.0 },  1.0736109751890863051e-09,    TEST_TOL2, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 20, x:  1.0,   y:  2.0 },  6.760965304863386741e-24,     TEST_TOL2, "gsl_sf_conicalP_cyl_reg_e");
    TEST_SF_IDD( s, gsl_sf_conicalP_cyl_reg_e, { n: 30, x:  100.0, y: 100.0 }, -4.268753482520651007e-63, TEST_TOL4, "gsl_sf_conicalP_cyl_reg_e");

    console.log( "    ... gsl_sf_legendre_H3d_0_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 1.0e-06, y: 1.0e-06 }, 0.9999999999998333333    , TEST_TOL0, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 1.0, y: 0.0 }, 1.0                      , TEST_TOL0, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 1.0, y: 1.0 }, 0.7160229153604338713    , TEST_TOL0, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 1.0, y: 100.0 }, -3.767437313149604566e-44 , TEST_TOL2, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 1.0, y: 500.0 }, -6.665351935878582205e-218, TEST_TOL2, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 100.0, y: 1.0 }, -0.004308757035378200029  , TEST_TOL0, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 100.0, y: 10.0 }, 7.508054627912986427e-07 , TEST_TOL0, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 1000.0, y: 1.0 }, 0.0007036067909088818319 , TEST_TOL0, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 1.0e+08, y: 1.0 }, 7.927485371429105968e-09 , TEST_TOL3, "gsl_sf_legendre_H3d_0_e");
    TEST_SF_DD( s, gsl_sf_legendre_H3d_0_e, { x: 1.0e+08, y: 100.0 }, -3.627118904186918957e-52 , 32.0*TEST_SQRT_TOL0, "gsl_sf_legendre_H3d_0_e");

    console.log( "    ... gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1.0e-06, y: 1.0e-06 },  3.333333333334222222e-07,  TEST_TOL0, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1.0,     y: 1.0e-10 },  4.714045207910316829e-11,  TEST_TOL0, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1.0,     y: 1.0     },  0.3397013994799344639,     TEST_TOL0, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1.0,     y: 100.0   }, -7.200624449531811272e-44,  TEST_TOL2, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1.0,     y: 500.0   },  4.192260336821728677e-218, TEST_TOL2, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 100.0,   y: 0.01    },  0.30117664944267412324   , TEST_TOL1, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 100.0,   y: 1.0     }, -0.007393833425336299309  , TEST_TOL0, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 100.0,   y: 10.0    }, -5.031062029821254982e-07 , TEST_TOL0, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1000.0,  y: 0.001   },  0.30116875865090396421   , TEST_TOL0, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1000.0,  y: 1.0     }, -0.0004776144516074971885 , TEST_TOL0, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1.0e+08, y: 1.0e-08 },  0.30116867893975679722   , TEST_TOL1, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1.0e+08, y: 1.0     },  3.0921097047369081582e-09, TEST_TOL4, "gsl_sf_legendre_H3d_1_e" );
    TEST_SF_DD( s, gsl_sf_legendre_H3d_1_e, { x: 1.0e+08, y: 100.0   }, -6.496142701296286936e-52 , 32.0 * TEST_SQRT_TOL0, "gsl_sf_legendre_H3d_1_e" );

    console.log( "    ... gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0e-06, y: 1.0e-06 },  1.1544011544013627977e-32,  TEST_TOL2, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0,     y: 1.0e-10 },  2.0224912016958766992e-52,  TEST_TOL2, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0,     y: 1.0     },  0.011498635037491577728,    TEST_TOL1, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0,     y: 5.0     },  0.0020696945662545205776,   TEST_TOL4, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0,     y: 7.0     }, -0.0017555303787488993676,   TEST_TOL4, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0,     y: 10.0    },  0.00008999979724504887101,  TEST_TOL2, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0,     y: 100.0   }, -4.185397793298567945e-44,   TEST_TOL2, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0,     y: 500.0   },  1.4235113901091961263e-217, TEST_TOL3, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 100.0,   y: 0.001   },  9.642762597222417946e-10,   TEST_TOL2, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 100.0,   y: 0.002   },  3.0821201254308036109e-08,  TEST_TOL2, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 100.0,   y: 0.01    },  0.00009281069019005840532,  TEST_TOL1, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 100.0,   y: 1.0     }, -0.008043100696178624653,    TEST_TOL2, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 100.0,   y: 10.0    }, -3.927678432813974207e-07,   TEST_TOL3, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1000.0,  y: 0.001   },  0.00009256365284253254503,  TEST_TOL1, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1000.0,  y: 0.01    }, -0.05553733815473079983,     TEST_TOL0, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0e+08, y: 1.0e-08 },  0.00009256115861125841299,  TEST_TOL2, "gsl_sf_legendre_H3d_e" );
    TEST_SF_IDD( s, gsl_sf_legendre_H3d_e, { n: 5, x: 1.0e+08, y: 100.0   }, -6.496143209092860765e-52 ,  128.0 * TEST_SQRT_TOL0, "gsl_sf_legendre_H3d_e" );

    //#if FIXME
    // sa.Integer = 0;
    // gsl_sf_legendre_H3d_array( 100, 1.0, 3.0, L );
    // TEST_SF_VAL( sa, L[0],   +0.0, gsl_sf_legendre_H3d( 0, 1.0, 3.0 ),   1.0e-12 );
    // TEST_SF_VAL( sa, L[1],   +0.0, gsl_sf_legendre_H3d( 1, 1.0, 3.0 ),   1.0e-12 );
    // TEST_SF_VAL( sa, L[10],  +0.0, gsl_sf_legendre_H3d( 10, 1.0, 3.0 ),  1.0e-12 );
    // TEST_SF_VAL( sa, L[100], +0.0, gsl_sf_legendre_H3d( 100, 1.0, 3.0 ), 1.0e-12 );
    // gsl_test( sa, "  gsl_sf_legendre_H3d_array(100, 1.0, 3.0)" );
    // s.Integer = s.Integer + sa.Integer;
    //#endif

    console.log( "    ... gsl_sf_legendre_Q0_e" );
    // x = -1 + 2^-16
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, -0.9999847412109375, -5.8917472200477175158028143531855, TEST_TOL4, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, -0.5,                -0.5493061443340548457,             TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, -1.0e-10,            -1.000000000000000000e-10,          TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 0.0,                  0.0,                               TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 1.0e-10,              1.000000000000000000e-10,          TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    // x = 1 - 2^-16
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 0.9999847412109375, 5.8917472200477175158028143531855, TEST_TOL4, "gsl_sf_legendre_Q0_e" );
    // x = 1 + 2^-16
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 1.0000152587890625, 5.8917548494422489138325509750429, TEST_TOL4, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 1.5,                0.8047189562170501873,             TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 9.99,               0.1004364599660005447,             TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 10.0,               0.1003353477310755806,             TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 10.01,              0.1002344395571710243,             TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 100.0,              0.010000333353334762015,           TEST_TOL0, "gsl_sf_legendre_Q0_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q0_e, 1.0e10,             1.000000000000000000e-10,          TEST_TOL0, "gsl_sf_legendre_Q0_e" );

    console.log( "    ... gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e, -0.9999847412109375,  4.8916573191196772369,    TEST_TOL4, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e, -0.5,                -0.7253469278329725772,    TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e, -0.01,               -0.9998999966664666524,    TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e, -1.0e-10,            -0.999999999999999999,     TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  0.0,                -1.0,                      TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  1.0e-10,            -0.999999999999999999,     TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  0.0001,             -0.9999999899999999667,    TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  0.01,               -0.9998999966664666524,    TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  0.5,                -0.7253469278329725772,    TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  0.9999847412109375,  4.8916573191196772369,    TEST_TOL4, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  1.0000152587890625,  4.8918447504867045145,    TEST_TOL4, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  1.5,                 0.20707843432557528095,   TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  9.99,                3.360235060345441639e-3,  TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  10.0,                3.353477310755806357e-3,  TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  10.01,               3.346739967281953346e-3,  TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  100.0,               3.333533347620158821e-5,  TEST_TOL0, "gsl_sf_legendre_Q1_e" );
    TEST_SF_D( s, gsl_sf_legendre_Q1_e,  1.0e10,              3.333333333333333333e-21, TEST_TOL0, "gsl_sf_legendre_Q1_e" );

    console.log( "    ... gsl_sf_legendre_Ql_e" );
    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 10,   x: -0.5 }, -0.29165813966586752393,     TEST_TOL0, "gsl_sf_legendre_Ql_e" );
    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 10,   x:  0.5 },  0.29165813966586752393,     TEST_TOL0, "gsl_sf_legendre_Ql_e" );
    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 10,   x:  1.5 },  0.000014714232718207477406, TEST_TOL0, "gsl_sf_legendre_Ql_e" );

    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 100,  x: -0.5 }, -0.09492507395207282096,    TEST_TOL1, "gsl_sf_legendre_Ql_e" );
    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 100,  x:  0.5 },  0.09492507395207282096,    TEST_TOL1, "gsl_sf_legendre_Ql_e" );
    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 100,  x:  1.5 },  1.1628163435044121988e-43, TEST_TOL2, "gsl_sf_legendre_Ql_e" );

    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 1000, x:  -0.5 }, -0.030105074974005303500,    TEST_TOL1, "gsl_sf_legendre_Ql_e" );
    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 1000, x:   0.5 },  0.030105074974005303500,    TEST_TOL1, "gsl_sf_legendre_Ql_e" );
    TEST_SF_ID( s, gsl_sf_legendre_Ql_e, { i: 1000, x:   1.1 },  1.0757258447825356443e-194, TEST_TOL3, "gsl_sf_legendre_Ql_e" );

    return s;

} // test_legendre

// ****************************************************************************

export function test_airy( )
{

    var s = { Integer: 0 };
    var m = GSL_MODE_DEFAULT;

    console.log("Test Airy Functions ...");

    // functions

    console.log("    ... gsl_sf_airy_Ai_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: -500.0, m: m },              0.0725901201040411396, TEST_TOL4, "gsl_sf_airy_Ai_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: -5.0, m: m },                0.3507610090241142,    TEST_TOL0, "gsl_sf_airy_Ai_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: -0.3000000000000094, m: m }, 0.4309030952855831,    TEST_TOL0, "gsl_sf_airy_Ai_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: 0.6999999999999907, m: m },  0.1891624003981519,    TEST_TOL0, "gsl_sf_airy_Ai_e");

    // This original value seemed to be slightly inaccurate in the last place.
    // I recomputed it with pari to get the new value which end in 885
    // instead of 882

    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: 1.649999999999991, m:  m },   0.05831058618720882,   TEST_TOL0, "gsl_sf_airy_Ai_e");


    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: 1.649999999999991, m:  m },   0.0583105861872088521,   TEST_TOL0, "gsl_sf_airy_Ai_e");

    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: 2.54999999999999,  m: m },    0.01446149513295428,   TEST_TOL0, "gsl_sf_airy_Ai_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: 3.499999999999987, m:  m },   0.002584098786989702,  TEST_TOL1, "gsl_sf_airy_Ai_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_e, { x: 5.39999999999998,  m: m },    4.272986169411866e-05, TEST_TOL0, "gsl_sf_airy_Ai_e");

    console.log("    ... gsl_sf_airy_Ai_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_scaled_e, { x: -5.0, m: m },                  0.3507610090241142, TEST_TOL0, "gsl_sf_airy_Ai_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_scaled_e, { x: 0.6999999999999907, m: m }, 0.2795125667681217, TEST_TOL0, "gsl_sf_airy_Ai_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_scaled_e, { x: 1.649999999999991, m: m },  0.2395493001442741, TEST_TOL0, "gsl_sf_airy_Ai_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_scaled_e, { x: 2.54999999999999, m: m },   0.2183658595899388, TEST_TOL0, "gsl_sf_airy_Ai_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_scaled_e, { x: 3.499999999999987, m: m },  0.2032920808163519, TEST_TOL0, "gsl_sf_airy_Ai_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_scaled_e, { x: 5.39999999999998, m: m },   0.1836050093282229, TEST_TOL0, "gsl_sf_airy_Ai_scaled_e");

    console.log("    ... gsl_sf_airy_Bi_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_e, { x: -500.0, m: m },             -0.094688570132991028, TEST_TOL4, "gsl_sf_airy_Bi_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_e, { x: -5.0, m: m },               -0.1383691349016005,   TEST_TOL1, "gsl_sf_airy_Bi_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_e, { x: 0.6999999999999907, m: m },  0.9733286558781599,   TEST_TOL0, "gsl_sf_airy_Bi_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_e, { x: 1.649999999999991, m: m },   2.196407956850028,    TEST_TOL0, "gsl_sf_airy_Bi_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_e, { x: 2.54999999999999, m: m },    6.973628612493443,    TEST_TOL0, "gsl_sf_airy_Bi_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_e, { x: 3.499999999999987, m: m },   33.05550675461069,    TEST_TOL1, "gsl_sf_airy_Bi_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_e, { x: 5.39999999999998, m: m },    1604.476078241272,    TEST_TOL1, "gsl_sf_airy_Bi_e");

    console.log("    ... gsl_sf_airy_Bi_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_scaled_e, { x: -5.0, m: m },                  -0.1383691349016005, TEST_TOL1, "gsl_sf_airy_Bi_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_scaled_e, { x: 0.6999999999999907, m: m },  0.6587080754582302, TEST_TOL0, "gsl_sf_airy_Bi_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_scaled_e, { x: 1.649999999999991, m: m },   0.5346449995597539, TEST_TOL0, "gsl_sf_airy_Bi_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_scaled_e, { x: 2.54999999999999, m: m },    0.461835455542297,  TEST_TOL0, "gsl_sf_airy_Bi_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_scaled_e, { x: 3.499999999999987, m: m },   0.4201771882353061, TEST_TOL1, "gsl_sf_airy_Bi_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_scaled_e, { x: 5.39999999999998, m: m },    0.3734050675720473, TEST_TOL0, "gsl_sf_airy_Bi_scaled_e");

    // derivatives

    console.log("    ... gsl_sf_airy_Ai_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_e, { x: -5.0, m: m },                 0.3271928185544435,       TEST_TOL1, "gsl_sf_airy_Ai_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_e, { x: -0.5500000000000094, m: m }, -0.1914604987143629,    TEST_TOL0, "gsl_sf_airy_Ai_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_e, { x: 0.4999999999999906, m: m },  -0.2249105326646850,    TEST_TOL0, "gsl_sf_airy_Ai_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_e, { x: 1.899999999999992, m: m },   -0.06043678178575718,   TEST_TOL0, "gsl_sf_airy_Ai_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_e, { x: 3.249999999999988, m: m },   -0.007792687926790889,  TEST_TOL0, "gsl_sf_airy_Ai_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_e, { x: 5.199999999999981, m: m },   -0.0001589434526459543, TEST_TOL1, "gsl_sf_airy_Ai_deriv_e");

    console.log("    ... gsl_sf_airy_Ai_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_scaled_e, { x: -5.0, m: m },                0.3271928185544435, TEST_TOL1, "gsl_sf_airy_Ai_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_scaled_e, { x: 0.5499999999999906, m: m }, -0.2874057279170166, TEST_TOL0, "gsl_sf_airy_Ai_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_scaled_e, { x: 1.499999999999991, m: m },  -0.3314199796863637, TEST_TOL0, "gsl_sf_airy_Ai_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_scaled_e, { x: 2.49999999999999, m: m },   -0.3661089384751620, TEST_TOL0, "gsl_sf_airy_Ai_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_scaled_e, { x: 3.649999999999986, m: m },  -0.3974033831453963, TEST_TOL0, "gsl_sf_airy_Ai_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Ai_deriv_scaled_e, { x: 6.299999999999977, m: m },  -0.4508799189585947, TEST_TOL0, "gsl_sf_airy_Ai_deriv_scaled_e");

    console.log("    ... gsl_sf_airy_Bi_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_e, { x: -5.0, m: m },                0.778411773001899,  TEST_TOL0, "gsl_sf_airy_Bi_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_e, { x: -0.5500000000000094, m: m }, 0.5155785358765014, TEST_TOL0, "gsl_sf_airy_Bi_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_e, { x: 0.4999999999999906, m: m },  0.5445725641405883, TEST_TOL0, "gsl_sf_airy_Bi_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_e, { x: 1.899999999999992, m: m },   3.495165862891568,  TEST_TOL0, "gsl_sf_airy_Bi_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_e, { x: 3.249999999999988, m: m },   36.55485149250338,  TEST_TOL0, "gsl_sf_airy_Bi_deriv_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_e, { x: 5.199999999999981, m: m },   2279.748293583233,  TEST_TOL1, "gsl_sf_airy_Bi_deriv_e");

    console.log("    ... gsl_sf_airy_Bi_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_scaled_e, { x: -5.0, m: m },               0.778411773001899,  TEST_TOL0, "gsl_sf_airy_Bi_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_scaled_e, { x: 0.5499999999999906, m: m }, 0.4322811281817566, TEST_TOL0, "gsl_sf_airy_Bi_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_scaled_e, { x: 1.499999999999991, m: m },  0.5542307563918037, TEST_TOL0, "gsl_sf_airy_Bi_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_scaled_e, { x: 2.49999999999999, m: m },   0.6755384441644985, TEST_TOL0, "gsl_sf_airy_Bi_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_scaled_e, { x: 3.649999999999986, m: m },  0.7613959373000228, TEST_TOL0, "gsl_sf_airy_Bi_deriv_scaled_e");
    TEST_SF_DM(s, gsl_sf_airy_Bi_deriv_scaled_e, { x: 6.299999999999977, m: m },  0.8852064139737571, TEST_TOL0, "gsl_sf_airy_Bi_deriv_scaled_e");

    console.log("    ... gsl_sf_airy_zero_Ai_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_e, 2,    -4.087949444130970617, TEST_TOL0, "gsl_sf_airy_zero_Ai_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_e, 50,   -38.02100867725525443, TEST_TOL0, "gsl_sf_airy_zero_Ai_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_e, 100,  -60.45555727411669871, TEST_TOL0, "gsl_sf_airy_zero_Ai_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_e, 110,  -64.43135670991324811, TEST_TOL0, "gsl_sf_airy_zero_Ai_e");

    console.log("    ... gsl_sf_airy_zero_Bi_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_e, 2,   -3.271093302836352716, TEST_TOL0, "gsl_sf_airy_zero_Bi_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_e, 50,  -37.76583438165180116, TEST_TOL0, "gsl_sf_airy_zero_Bi_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_e, 100, -60.25336482580837088, TEST_TOL0, "gsl_sf_airy_zero_Bi_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_e, 110, -64.2355167606561537,  TEST_TOL0, "gsl_sf_airy_zero_Bi_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_e, 111, -64.6268994819519378,  TEST_TOL0, "gsl_sf_airy_zero_Bi_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_e, 200, -95.88699147356682665, TEST_TOL0, "gsl_sf_airy_zero_Bi_e");

    console.log("    ... gsl_sf_airy_zero_Ai_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_deriv_e, 2,    -3.248197582179836561, TEST_TOL0, "gsl_sf_airy_zero_Ai_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_deriv_e, 50,   -37.76565910053887108, TEST_TOL0, "gsl_sf_airy_zero_Ai_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_deriv_e, 100,  -60.25329596442479317, TEST_TOL0, "gsl_sf_airy_zero_Ai_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_deriv_e, 110,  -64.23545617243546956, TEST_TOL0, "gsl_sf_airy_zero_Ai_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Ai_deriv_e, 1000, -280.9378080358935071, TEST_TOL0, "gsl_sf_airy_zero_Ai_deriv_e");

    console.log("    ... gsl_sf_airy_zero_Bi_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_deriv_e, 2,    -4.073155089071828216, TEST_TOL0, "gsl_sf_airy_zero_Bi_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_deriv_e, 50,   -38.02083574095788210, TEST_TOL0, "gsl_sf_airy_zero_Bi_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_deriv_e, 100,  -60.45548887257140819, TEST_TOL0, "gsl_sf_airy_zero_Bi_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_deriv_e, 110,  -64.43129648944845060, TEST_TOL0, "gsl_sf_airy_zero_Bi_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_deriv_e, 111,  -64.82208737584206093, TEST_TOL0, "gsl_sf_airy_zero_Bi_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_deriv_e, 200,  -96.04731050310324450, TEST_TOL0, "gsl_sf_airy_zero_Bi_deriv_e");
    TEST_SF_I(s, gsl_sf_airy_zero_Bi_deriv_e, 1000, -281.0315164471118527, TEST_TOL0, "gsl_sf_airy_zero_Bi_deriv_e");

    return s;

} // test_airy

// ****************************************************************************

export function test_bessel( )
{

    var s    = { Integer: 0 };
    var sa   = { Integer: 0 };

    var J   = [];
    var Y   = [];
    var I   = [];
    var K   = [];

    console.log("Test Bessel Functions ...");

    console.log("    ... gsl_sf_bessel_J0_e");
    TEST_SF_D(s, gsl_sf_bessel_J0_e, 0.1,     0.99750156206604003230,    TEST_TOL0, "gsl_sf_bessel_J0_e");
    TEST_SF_D(s, gsl_sf_bessel_J0_e, 2.0,     0.22389077914123566805,    TEST_TOL0, "gsl_sf_bessel_J0_e");
    TEST_SF_D(s, gsl_sf_bessel_J0_e, 100.0,   0.019985850304223122424,   TEST_TOL0, "gsl_sf_bessel_J0_e");
    TEST_SF_D(s, gsl_sf_bessel_J0_e, 1.0e+10, 2.1755917502468917269e-06, TEST_SQRT_TOL0, "gsl_sf_bessel_J0_e");

    console.log("    ... gsl_sf_bessel_J1_e");
    TEST_SF_D(s, gsl_sf_bessel_J1_e, (0.1),      0.04993752603624199756,   TEST_TOL0, "gsl_sf_bessel_J1_e");
    TEST_SF_D(s, gsl_sf_bessel_J1_e, (2.0),      0.57672480775687338720,   TEST_TOL0, "gsl_sf_bessel_J1_e");
    TEST_SF_D(s, gsl_sf_bessel_J1_e, (100.0),   -0.07714535201411215803,   TEST_TOL0, "gsl_sf_bessel_J1_e");
    TEST_SF_D(s, gsl_sf_bessel_J1_e, (1.0e+10), -7.676508175684157103e-06, TEST_TOL4, "gsl_sf_bessel_J1_e");

    console.log("    ... gsl_sf_bessel_Jn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 4,   x: 0.1 },     2.6028648545684032338e-07, TEST_TOL0, "gsl_sf_bessel_Jn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 5,   x: 2.0 },     0.007039629755871685484,   TEST_TOL0, "gsl_sf_bessel_Jn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 10,  x: 20.0 },   0.18648255802394508321,    TEST_TOL0, "gsl_sf_bessel_Jn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 100, x: 100.0 }, 0.09636667329586155967,    TEST_TOL0, "gsl_sf_bessel_Jn_e");

    // exercise the BUG#3 problem
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 2, x: 900.0 }, -0.019974345269680646400, TEST_TOL4, "gsl_sf_bessel_Jn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 2, x: 15000.0}, -0.0020455820181216382666, TEST_TOL4, "gsl_sf_bessel_Jn_e");

//#ifdef TEST_LARGE
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 0, x: 1.0e+10 }, 2.1755917502468917269e-06, TEST_SQRT_TOL0, "gsl_sf_bessel_Jn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 1, x: 1.0e+10 }, -7.676508175684157103e-06, TEST_TOL4, "gsl_sf_bessel_Jn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Jn_e, { i: 0, x: 20000.0 }, 0.00556597490495494615709982972, TEST_TOL4, "gsl_sf_bessel_Jn_e");
//#endif

    console.log("    ... gsl_sf_bessel_Y0_e");
    TEST_SF_D(s, gsl_sf_bessel_Y0_e, 0.1,         -1.5342386513503668441,    TEST_TOL0, "gsl_sf_bessel_Y0_e");
    TEST_SF_D(s, gsl_sf_bessel_Y0_e, 2.0,          0.5103756726497451196,    TEST_TOL0, "gsl_sf_bessel_Y0_e");
    TEST_SF_D(s, gsl_sf_bessel_Y0_e, 256.0,       -0.03381290171792454909 ,  TEST_TOL0, "gsl_sf_bessel_Y0_e");
    TEST_SF_D(s, gsl_sf_bessel_Y0_e, 4294967296.0, 3.657903190017678681e-06, TEST_SQRT_TOL0, "gsl_sf_bessel_Y0_e");

    console.log("    ... gsl_sf_bessel_Y1_e");
    TEST_SF_D(s, gsl_sf_bessel_Y1_e, 0.1,         -6.45895109470202698800,     TEST_TOL0, "gsl_sf_bessel_Y1_e");
    TEST_SF_D(s, gsl_sf_bessel_Y1_e, 2.0,           -0.10703243154093754689,     TEST_TOL0, "gsl_sf_bessel_Y1_e");
    TEST_SF_D(s, gsl_sf_bessel_Y1_e, 100.0,       -0.020372312002759793305,    TEST_TOL0, "gsl_sf_bessel_Y1_e");
    TEST_SF_D(s, gsl_sf_bessel_Y1_e, 4294967296.0, 0.000011612249378370766284, TEST_TOL4, "gsl_sf_bessel_Y1_e");

    console.log("    ... gsl_sf_bessel_Yn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Yn_e, { i: 4, x: 0.1 },            -305832.29793353160319,    TEST_TOL1, "gsl_sf_bessel_Yn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Yn_e, { i: 5, x: 2.0 },              -9.935989128481974981,     TEST_TOL0, "gsl_sf_bessel_Yn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Yn_e, { i: 100, x: 100.0 },        -0.16692141141757650654,   TEST_TOL0, "gsl_sf_bessel_Yn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Yn_e, { i: 100, x: 4294967296.0 },  3.657889671577715808e-06, TEST_SQRT_TOL0, "gsl_sf_bessel_Yn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Yn_e, { i: 1000, x: 4294967296.0 }, 3.656551321485397501e-06, 2.0e-05, "gsl_sf_bessel_Yn_e");

    TEST_SF_ID(s, gsl_sf_bessel_Yn_e, { i: 2, x: 15000.0 }, -0.006185217273358617849, TEST_TOL4, "gsl_sf_bessel_Yn_e");

    console.log("    ... gsl_sf_bessel_I0_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I0_scaled_e, 1.0e-10, 0.99999999990000000001,   TEST_TOL0, "gsl_sf_bessel_I0_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I0_scaled_e, 0.1,     0.90710092578230109640,   TEST_TOL0, "gsl_sf_bessel_I0_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I0_scaled_e, 2.0,     0.30850832255367103953,   TEST_TOL0, "gsl_sf_bessel_I0_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I0_scaled_e, 100.0,   0.03994437929909668265,   TEST_TOL0, "gsl_sf_bessel_I0_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I0_scaled_e, 65536.0, 0.0015583712551952223537, TEST_TOL0, "gsl_sf_bessel_I0_scaled_e");

    console.log("    ... gsl_sf_bessel_I1_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I1_scaled_e, 0.1,     0.04529844680880932501,   TEST_TOL0, "gsl_sf_bessel_I1_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I1_scaled_e, 2.0,     0.21526928924893765916,   TEST_TOL0, "gsl_sf_bessel_I1_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I1_scaled_e, 100.0,   0.03974415302513025267,   TEST_TOL0, "gsl_sf_bessel_I1_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_I1_scaled_e, 65536.0, 0.0015583593657207350452, TEST_TOL0, "gsl_sf_bessel_I1_scaled_e");

    console.log("    ... gsl_sf_bessel_In_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_In_scaled_e, { i:  -4, x:   0.1 }, 2.3575258620054605307e-07, TEST_TOL0, "gsl_sf_bessel_In_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_In_scaled_e, { i:   4, x:   0.1 }, 2.3575258620054605307e-07, TEST_TOL0, "gsl_sf_bessel_In_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_In_scaled_e, { i:   5, x:   2.0 }, 0.0013297610941881578142, TEST_TOL0, "gsl_sf_bessel_In_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_In_scaled_e, { i: 100, x: 100.0 }, 1.7266862628167695785e-22, TEST_TOL0, "gsl_sf_bessel_In_scaled_e");

    // BJG: the "exact" values in the following two tests were originally computed from the
    // taylor series for I_nu using "long double" and rescaling.  The last few digits
    // were inaccurate due to cumulative roundoff.
    //
    // BJG: 2006/05 I have now replaced these with the term asymptotic
    // expansion from A&S 9.7.1 which should be fully accurate.

    TEST_SF_ID(s, gsl_sf_bessel_In_scaled_e, { i: 2, x: 1.0e7 }, 1.261566024466416433e-4, TEST_TOL2, "gsl_sf_bessel_In_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_In_scaled_e, { i: 2, x: 1.0e8 }, 3.989422729212649531e-5, TEST_TOL2, "gsl_sf_bessel_In_scaled_e");

    console.log("    ... gsl_sf_bessel_I0_e");
    TEST_SF_D(s, gsl_sf_bessel_I0_e, (0.1), 1.0025015629340956014, TEST_TOL0, "gsl_sf_bessel_I0_e");
    TEST_SF_D(s, gsl_sf_bessel_I0_e, (2.0), 2.2795853023360672674, TEST_TOL0, "gsl_sf_bessel_I0_e");
    TEST_SF_D(s, gsl_sf_bessel_I0_e, (100.0), 1.0737517071310738235e+42, TEST_TOL2, "gsl_sf_bessel_I0_e");

    console.log("    ... gsl_sf_bessel_I1_e");
    TEST_SF_D(s, gsl_sf_bessel_I1_e, (0.1), 0.05006252604709269211,      TEST_TOL0, "gsl_sf_bessel_I1_e");
    TEST_SF_D(s, gsl_sf_bessel_I1_e, (2.0), 1.59063685463732906340,      TEST_TOL0, "gsl_sf_bessel_I1_e");
    TEST_SF_D(s, gsl_sf_bessel_I1_e, (100.0), 1.0683693903381624812e+42, TEST_TOL2, "gsl_sf_bessel_I1_e");

    console.log("    ... gsl_sf_bessel_In_e");
    TEST_SF_ID(s, gsl_sf_bessel_In_e, { i:   4, x:   0.1 }, 2.6054690212996573677e-07, TEST_TOL0, "gsl_sf_bessel_In_e");
    TEST_SF_ID(s, gsl_sf_bessel_In_e, { i:   5, x:   2.0 }, 0.009825679323131702321,   TEST_TOL0, "gsl_sf_bessel_In_e");
    TEST_SF_ID(s, gsl_sf_bessel_In_e, { i: 100, x: 100.0 }, 4.641534941616199114e+21,  TEST_TOL2, "gsl_sf_bessel_In_e");

    console.log("    ... gsl_sf_bessel_K0_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_K0_scaled_e, (0.1), 2.6823261022628943831, TEST_TOL0, "gsl_sf_bessel_K0_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_K0_scaled_e, (2.0), 0.8415682150707714179, TEST_TOL0, "gsl_sf_bessel_K0_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_K0_scaled_e, (100.0), 0.1251756216591265789, TEST_TOL0, "gsl_sf_bessel_K0_scaled_e");

    console.log("    ... gsl_sf_bessel_K1_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_K1_scaled_e, (0.1), 10.890182683049696574, TEST_TOL0, "gsl_sf_bessel_K1_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_K1_scaled_e, (2.0), 1.0334768470686885732, TEST_TOL0, "gsl_sf_bessel_K1_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_K1_scaled_e, (100.0), 0.1257999504795785293, TEST_TOL0, "gsl_sf_bessel_K1_scaled_e");

    console.log("    ... gsl_sf_bessel_Kn_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_Kn_scaled_e, { i:   4,  x:   0.1 }, 530040.2483725626207, TEST_TOL1, "gsl_sf_bessel_Kn_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_Kn_scaled_e, { i:   5,  x:   2.0 }, 69.68655087607675118, TEST_TOL0, "gsl_sf_bessel_Kn_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_Kn_scaled_e, { i: 100,  x: 100.0 }, 2.0475736731166756813e+19, TEST_TOL1, "gsl_sf_bessel_Kn_scaled_e");

    console.log("    ... gsl_sf_bessel_K0_e");
    TEST_SF_D(s, gsl_sf_bessel_K0_e, 0.1, 2.4270690247020166125, TEST_TOL0, "gsl_sf_bessel_K0_e");
    TEST_SF_D(s, gsl_sf_bessel_K0_e, 2.0, 0.11389387274953343565, TEST_TOL0, "gsl_sf_bessel_K0_e");
    TEST_SF_D(s, gsl_sf_bessel_K0_e, 100.0, 4.656628229175902019e-45, TEST_TOL2, "gsl_sf_bessel_K0_e");

    console.log("    ... gsl_sf_bessel_K1_e");
    TEST_SF_D(s, gsl_sf_bessel_K1_e, 0.1, 9.853844780870606135,       TEST_TOL0, "gsl_sf_bessel_K1_e");
    TEST_SF_D(s, gsl_sf_bessel_K1_e, 2.0, 0.13986588181652242728,     TEST_TOL0, "gsl_sf_bessel_K1_e");
    TEST_SF_D(s, gsl_sf_bessel_K1_e, 100.0, 4.679853735636909287e-45, TEST_TOL2, "gsl_sf_bessel_K1_e");

    console.log("    ... gsl_sf_bessel_Kn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Kn_e, { i:   4, x:   0.1 }, 479600.2497925682849,     TEST_TOL1, "gsl_sf_bessel_Kn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Kn_e, { i:   5, x:   2.0 }, 9.431049100596467443,     TEST_TOL0, "gsl_sf_bessel_Kn_e");
    TEST_SF_ID(s, gsl_sf_bessel_Kn_e, { i: 100, x: 100.0 }, 7.617129630494085416e-25, TEST_TOL2, "gsl_sf_bessel_Kn_e");

    console.log("    ... gsl_sf_bessel_j0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j0s_e, -10.0, -0.05440211108893698134, TEST_TOL0, "gsl_sf_bessel_j0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j0s_e, 0.001, 0.9999998333333416667, TEST_TOL0, "gsl_sf_bessel_j0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j0s_e,   1.0, 0.84147098480789650670, TEST_TOL0, "gsl_sf_bessel_j0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j0s_e,  10.0, -0.05440211108893698134, TEST_TOL0, "gsl_sf_bessel_j0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j0s_e, 100.0, -0.005063656411097587937, TEST_TOL1, "gsl_sf_bessel_j0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j0s_e, 1048576.0, 3.1518281938718287624e-07, TEST_TOL2, "gsl_sf_bessel_j0s_e");

    console.log("    ... gsl_sf_bessel_j1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j1s_e, (-10.0), -0.07846694179875154709, TEST_TOL0, "gsl_sf_bessel_j1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j1s_e, (0.01), 0.003333300000119047399, TEST_TOL0, "gsl_sf_bessel_j1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j1s_e, (  1.0), 0.30116867893975678925, TEST_TOL0, "gsl_sf_bessel_j1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j1s_e, ( 10.0), 0.07846694179875154709, TEST_TOL0, "gsl_sf_bessel_j1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j1s_e, (100.0), -0.008673825286987815220, TEST_TOL0, "gsl_sf_bessel_j1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j1s_e, (1048576.0), -9.000855242905546158e-07, TEST_TOL0, "gsl_sf_bessel_j1s_e");

    console.log("    ... gsl_sf_bessel_j2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j2s_e, (-10.0), 0.07794219362856244547, TEST_TOL0, "gsl_sf_bessel_j2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j2s_e, (0.01), 6.666619047751322551e-06, TEST_TOL0, "gsl_sf_bessel_j2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j2s_e, (  1.0), 0.06203505201137386110, TEST_TOL0, "gsl_sf_bessel_j2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j2s_e, ( 10.0), 0.07794219362856244547, TEST_TOL0, "gsl_sf_bessel_j2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j2s_e, (100.0), 0.004803441652487953480, TEST_TOL1, "gsl_sf_bessel_j2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_j2s_e, (1048576.0), -3.1518539455252413111e-07, TEST_TOL2, "gsl_sf_bessel_j2s_e");

    console.log("    ... gsl_sf_bessel_jl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 0, x: 0.0 }, 1.0, TEST_TOL0, "gsl_sf_bessel_jl_e");

    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 1,    x:    10.0 },   0.07846694179875154709000, TEST_TOL0, "gsl_sf_bessel_jl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 5,    x:     1.0 },   0.00009256115861125816357, TEST_TOL0, "gsl_sf_bessel_jl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 10,   x:    10.0 },   0.06460515449256426427,    TEST_TOL0, "gsl_sf_bessel_jl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 100,  x:   100.0 },   0.010880477011438336539,   TEST_TOL1, "gsl_sf_bessel_jl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 2000, x: 1048576.0 }, 7.449384239168568534e-07,  TEST_SQRT_TOL0, "gsl_sf_bessel_jl_e");

    // related to BUG#3 problem
    TEST_SF_ID(s, gsl_sf_bessel_jl_e, { i: 2, x: 900.0 },   -0.0011089115568832940086,  TEST_TOL4, "gsl_sf_bessel_jl_e");
    TEST_SF_ID(s, gsl_sf_bessel_jl_e, { i: 2, x: 15000.0 }, -0.00005955592033075750554, TEST_TOL4, "gsl_sf_bessel_jl_e");

    // Bug report by Mario Santos, value computed from AS 10.1.8
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 100, x: 1000.0 }, -0.00025326311230945818285, TEST_TOL4, "gsl_sf_bessel_jl_e");

    // Bug reported by Koichi Takahashi <ktakahashi@molsci.org>,
    // computed from Pari besseljh(n,x) and AS 10.1.1

    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 30, x: 3878.62 }, -0.00023285567034330878410434732790, TEST_TOL4, "gsl_sf_bessel_jl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 49, x: 9912.63 }, 5.2043354544842669214485107019E-5 , TEST_TOL4, "gsl_sf_bessel_jl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 49, x: 9950.35 }, 5.0077368819565969286578715503E-5 , TEST_TOL4, "gsl_sf_bessel_jl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_jl_e, { i: 52, x: 9930.51 }, -7.4838588266727718650124475651E-6 , TEST_TOL4, "gsl_sf_bessel_jl_e");

    console.log("    ... gsl_sf_bessel_y0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y0s_e, 0.001, -999.99950000004166670, TEST_TOL0, "gsl_sf_bessel_y0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y0s_e,   1.0, -0.5403023058681397174, TEST_TOL0, "gsl_sf_bessel_y0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y0s_e,  10.0, 0.08390715290764524523, TEST_TOL0, "gsl_sf_bessel_y0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y0s_e, 100.0, -0.008623188722876839341, TEST_TOL0, "gsl_sf_bessel_y0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y0s_e, 65536.0, 0.000011014324202158573930, TEST_TOL0, "gsl_sf_bessel_y0s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y0s_e, 4294967296.0, 2.0649445131370357007e-10, TEST_SQRT_TOL0, "gsl_sf_bessel_y0s_e");

    console.log("    ... gsl_sf_bessel_y1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y1s_e,  0.01, -10000.499987500069444, TEST_TOL0, "gsl_sf_bessel_y1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y1s_e,   1.0, -1.3817732906760362241, TEST_TOL0, "gsl_sf_bessel_y1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y1s_e,  10.0, 0.06279282637970150586, TEST_TOL0, "gsl_sf_bessel_y1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y1s_e, 100.0, 0.004977424523868819543, TEST_TOL0, "gsl_sf_bessel_y1s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y1s_e, 4294967296.0, 1.0756463271573404688e-10, TEST_SQRT_TOL0, "gsl_sf_bessel_y1s_e");

    console.log("    ... gsl_sf_bessel_y2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y2s_e,  0.01, -3.0000500012499791668e+06, TEST_TOL0, "gsl_sf_bessel_y2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y2s_e,   1.0, -3.605017566159968955, TEST_TOL0, "gsl_sf_bessel_y2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y2s_e,  10.0, -0.06506930499373479347, TEST_TOL0, "gsl_sf_bessel_y2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y2s_e, 100.0, 0.008772511458592903927, TEST_TOL0, "gsl_sf_bessel_y2s_e");
    TEST_SF_D(s,  gsl_sf_bessel_y2s_e, 4294967296.0, -2.0649445123857054207e-10, TEST_SQRT_TOL0, "gsl_sf_bessel_y2s_e");

    console.log("    ... gsl_sf_bessel_yl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 0,    x:     0.01 }, -99.995000041666528,    TEST_TOL0, "gsl_sf_bessel_yl_e"); 
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 0,    x:     1.0 },  -0.54030230586813972,   TEST_TOL0, "gsl_sf_bessel_yl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 1,    x:    10.0 },   0.062792826379701506,   TEST_TOL0, "gsl_sf_bessel_yl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 5,    x:     1.0 },  -999.44034339223641,     TEST_TOL0, "gsl_sf_bessel_yl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 10,   x:     0.01 }, -6.5473079797378378e+30, TEST_TOL0, "gsl_sf_bessel_yl_e"); 
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 10,   x:    10.0 },  -0.172453672088057849,    TEST_TOL0, "gsl_sf_bessel_yl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 100,  x:     1.0 },  -6.6830794632586775e+186, TEST_TOL1, "gsl_sf_bessel_yl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 100,  x:   100.0 },  -0.0229838504915622811,   TEST_TOL1, "gsl_sf_bessel_yl_e");
    TEST_SF_ID(s,  gsl_sf_bessel_yl_e, { i: 2000, x: 1048576.0 }, 5.9545201447146155e-07,  TEST_SQRT_TOL0, "gsl_sf_bessel_yl_e");

    console.log("    ... gsl_sf_bessel_i0s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i0s_scaled_e, 0.0, 1.0, TEST_TOL0, "gsl_sf_bessel_i0s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i0s_scaled_e, 0.1, 0.9063462346100907067, TEST_TOL0, "gsl_sf_bessel_i0s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i0s_scaled_e, 2.0, 0.24542109027781645493, TEST_TOL0, "gsl_sf_bessel_i0s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i0s_scaled_e, 100.0, 0.005000000000000000000, TEST_TOL0, "gsl_sf_bessel_i0s_scaled_e");

    console.log("    ... gsl_sf_bessel_i1s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i1s_scaled_e, 0.0, 0.0, TEST_TOL0, "gsl_sf_bessel_i1s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i1s_scaled_e, 0.1, 0.030191419289002226846, TEST_TOL0, "gsl_sf_bessel_i1s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i1s_scaled_e, 2.0, 0.131868364583275317610, TEST_TOL0, "gsl_sf_bessel_i1s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i1s_scaled_e, 100.0, 0.004950000000000000000, TEST_TOL0, "gsl_sf_bessel_i1s_scaled_e");

    console.log("    ... gsl_sf_bessel_i2s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i2s_scaled_e, 0.0, 0.0, TEST_TOL0, "gsl_sf_bessel_i2s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i2s_scaled_e, 0.1, 0.0006036559400239012567, TEST_TOL0, "gsl_sf_bessel_i2s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i2s_scaled_e, 2.0, 0.0476185434029034785100, TEST_TOL0, "gsl_sf_bessel_i2s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_i2s_scaled_e, 100.0, 0.0048515000000000000000, TEST_TOL0, "gsl_sf_bessel_i2s_scaled_e");

    console.log("    ... gsl_sf_bessel_il_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_il_scaled_e, { i:   0, x: 0.0 }, 1.0, TEST_TOL0, "gsl_sf_bessel_il_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_il_scaled_e, { i:   1, x: 0.0 }, 0.0, TEST_TOL0, "gsl_sf_bessel_il_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_il_scaled_e, { i:   4, x: 0.001 }, 1.0571434341190365013e-15, TEST_TOL0, "gsl_sf_bessel_il_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_il_scaled_e, { i:   4, x:   0.1 }, 9.579352242057134927e-08,  TEST_TOL1, "gsl_sf_bessel_il_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_il_scaled_e, { i:   5, x:   2.0 }, 0.0004851564602127540059,  TEST_TOL0, "gsl_sf_bessel_il_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_il_scaled_e, { i:   5, x: 100.0 }, 0.004300446777500000000,   TEST_TOL1, "gsl_sf_bessel_il_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_il_scaled_e, { i: 100, x: 100.0 }, 1.3898161964299132789e-23, TEST_TOL0, "gsl_sf_bessel_il_scaled_e");

    console.log("    ... gsl_sf_bessel_k0s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k0s_scaled_e, (0.1), 15.707963267948966192, TEST_TOL0, "gsl_sf_bessel_k0s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k0s_scaled_e, (2.0), 0.7853981633974483096, TEST_TOL0, "gsl_sf_bessel_k0s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k0s_scaled_e, (100.0), 0.015707963267948966192, TEST_TOL0, "gsl_sf_bessel_k0s_scaled_e");

    console.log("    ... gsl_sf_bessel_k1s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k1s_scaled_e, (0.1), 172.78759594743862812, TEST_TOL0, "gsl_sf_bessel_k1s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k1s_scaled_e, (2.0), 1.1780972450961724644, TEST_TOL0, "gsl_sf_bessel_k1s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k1s_scaled_e, (100.0), 0.015865042900628455854, TEST_TOL0, "gsl_sf_bessel_k1s_scaled_e");

    console.log("    ... gsl_sf_bessel_k2s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k2s_scaled_e, (0.1), 5199.335841691107810, TEST_TOL0, "gsl_sf_bessel_k2s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k2s_scaled_e, (2.0), 2.5525440310417070063, TEST_TOL0, "gsl_sf_bessel_k2s_scaled_e");
    TEST_SF_D(s, gsl_sf_bessel_k2s_scaled_e, (100.0), 0.016183914554967819868, TEST_TOL0, "gsl_sf_bessel_k2s_scaled_e");

    console.log("    ... gsl_sf_bessel_kl_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_kl_scaled_e, { i:   4, x: 1.0/256.0 }, 1.8205599816961954439e+14, TEST_TOL0, "gsl_sf_bessel_kl_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_kl_scaled_e, { i:   4, x: 1.0/8.0 },   6.1173217814406597530e+06, TEST_TOL0, "gsl_sf_bessel_kl_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_kl_scaled_e, { i:   5, x:   2.0 },     138.10735829492005119,     TEST_TOL0, "gsl_sf_bessel_kl_scaled_e");
    TEST_SF_ID(s, gsl_sf_bessel_kl_scaled_e, { i: 100, x: 100.0 },     3.985930768060258219e+18,  TEST_TOL1, "gsl_sf_bessel_kl_scaled_e");

    console.log("    ... gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 0.0001, y: 1.0 },         0.7652115411876708497,  TEST_TOL2, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 0.0001, y: 10.0 },       -0.2459270166445205,     TEST_TOL2, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 0.0009765625, y: 10.0 }, -0.2458500798634692,     TEST_TOL2, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 0.75, y: 1.0 },           0.5586524932048917478,  TEST_TOL2, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 0.75, y: 10.0 },         -0.04968928974751508135, TEST_TOL2, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x:  1.0, y: 0.001 }, 0.0004999999375000026,     TEST_TOL0, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x:  1.0,   y: 1.0 }, 0.4400505857449335160,     TEST_TOL0, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x:  1.75,  y: 1.0 }, 0.168593922545763170103,     TEST_TOL1, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 30.0,   y: 1.0 }, 3.482869794251482902e-42,  TEST_TOL0, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 30.0, y: 100.0 }, 0.08146012958117222297,    TEST_TOL1, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 10.0,   y: 1.0 }, 2.6306151236874532070e-10, TEST_TOL0, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 10.0, y: 100.0 }, -0.05473217693547201474,   TEST_TOL2, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 10.2, y: 100.0 }, -0.03548919161046526864,   TEST_TOL2, "gsl_sf_bessel_Jnu_e");

    // related to BUG#3 problem
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 2.0, y: 900.0 },   -0.019974345269680646400,  TEST_TOL4, "gsl_sf_bessel_Jnu_e");
    TEST_SF_DD(s, gsl_sf_bessel_Jnu_e, { x: 2.0, y: 15000.0 }, -0.0020455820181216382666, TEST_TOL4, "gsl_sf_bessel_Jnu_e");

    console.log("    ... gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x: 0.0001, y: 1.0   },  0.08813676933044478439,    TEST_TOL2, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x: 0.0001, y: 10.0  },  0.05570979797521875261,    TEST_TOL2, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x:  0.75,  y: 1.0   }, -0.6218694174429746383,     TEST_TOL0, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x:  0.75,  y: 10.0  },  0.24757063446760384953,    TEST_TOL0, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x:  1.0,   y: 0.001 }, -636.6221672311394281,      TEST_TOL0, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x:  1.0,   y: 1.0   }, -0.7812128213002887165,     TEST_TOL0, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x: 30.0,   y: 1.0   }, -3.0481287832256432162e+39, TEST_TOL0, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x: 30.0,   y: 100.0 },  0.006138839212010033452,   TEST_TOL2, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x: 10.0,   y: 1.0   }, -1.2161801427868918929e+08, TEST_TOL0, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x: 10.0,   y: 100.0 },  0.05833157423641492875,    TEST_TOL2, "gsl_sf_bessel_Ynu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Ynu_e, { x: 10.2,   y: 100.0 },  0.07169383985546287091,    TEST_TOL2, "gsl_sf_bessel_Ynu_e");

    console.log("    ... gsl_sf_bessel_Inu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_scaled_e, { x: 0.0001,y: 10.0  }, 0.12783333709581669672,    TEST_TOL0, "gsl_sf_bessel_Inu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_scaled_e, { x:  1.0,  y: 0.001 }, 0.0004995003123542213370,  TEST_TOL0, "gsl_sf_bessel_Inu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_scaled_e, { x:  1.0,  y:   1.0 }, 0.20791041534970844887,    TEST_TOL0, "gsl_sf_bessel_Inu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_scaled_e, { x: 30.0,  y:   1.0 }, 1.3021094983785914437e-42, TEST_TOL0, "gsl_sf_bessel_Inu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_scaled_e, { x: 30.0,  y: 100.0 }, 0.0004486987756920986146,  TEST_TOL3, "gsl_sf_bessel_Inu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_scaled_e, { x: 10.0,  y:   1.0 }, 1.0127529864692066036e-10, TEST_TOL0, "gsl_sf_bessel_Inu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_scaled_e, { x: 10.0,  y: 100.0 }, 0.024176682718258828365,   TEST_TOL3, "gsl_sf_bessel_Inu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_scaled_e, { x: 10.2,  y: 100.0 }, 0.023691628843913810043,   TEST_TOL3, "gsl_sf_bessel_Inu_scaled_e");

    console.log("    ... gsl_sf_bessel_Inu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_e, { x: 0.0001,y:  10.0 }, 2815.7166269770030352,     TEST_TOL0, "gsl_sf_bessel_Inu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_e, { x:  1.0,  y: 0.001 }, 0.0005000000625000026042,  TEST_TOL0, "gsl_sf_bessel_Inu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_e, { x:  1.0,  y:   1.0 }, 0.5651591039924850272,     TEST_TOL0, "gsl_sf_bessel_Inu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_e, { x: 30.0,  y:   1.0 }, 3.539500588106447747e-42,  TEST_TOL0, "gsl_sf_bessel_Inu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_e, { x: 30.0,  y: 100.0 }, 1.2061548704498434006e+40, TEST_TOL2, "gsl_sf_bessel_Inu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_e, { x: 10.0,  y:   1.0 }, 2.7529480398368736252e-10, TEST_TOL0, "gsl_sf_bessel_Inu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_e, { x: 10.0,  y: 100.0 }, 6.498975524720147799e+41,  TEST_TOL2, "gsl_sf_bessel_Inu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Inu_e, { x: 10.2,  y: 100.0 }, 6.368587361287030443e+41,  TEST_TOL2, "gsl_sf_bessel_Inu_e");

    console.log("    ... gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x: 0.0001,y:   10.0 }, 0.3916319346235421817, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x:  1.0,  y:  0.001 }, 1000.9967345590684524, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x:  1.0,  y:    1.0 }, 1.6361534862632582465, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x: 30.0,  y:    1.0 }, 1.2792629867539753925e+40, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x: 30.0,  y:  100.0 }, 10.673443449954850040, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x: 10.0,  y:    1.0 }, 4.912296520990198599e+08, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x: 10.0,  y:  100.0 }, 0.20578687173955779807, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x: 10.0,  y: 1000.0 }, 0.04165905142800565788, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x: 10.0,  y: 1.0e+8 }, 0.00012533147624060789938, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_scaled_e, { x: 10.2,  y: 100.0  }, 0.20995808355244385075, TEST_TOL0, "gsl_sf_bessel_Knu_scaled_e");

    console.log("    ... gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x: 0.0001,y: 0.001 }, 7.023689431812884141,      TEST_TOL0, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x: 0.0001,y:  10.0 }, 0.000017780062324654874306, TEST_TOL0, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x:  1.0,  y: 0.001 }, 999.9962381560855743,      TEST_TOL0, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x:  1.0,  y:   1.0 }, 0.6019072301972345747,     TEST_TOL0, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x: 10.0,  y: 0.001 }, 1.8579455483904008064e+38, TEST_TOL0, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x: 10.0,  y:   1.0 }, 1.8071328990102945469e+08, TEST_TOL0, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x: 10.0,  y: 100.0 }, 7.655427977388100611e-45,  TEST_TOL2, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x: 10.2,  y: 100.0 }, 7.810600225948217841e-45,  TEST_TOL2, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x: 30.0,  y:   1.0 }, 4.706145526783626883e+39,  TEST_TOL1, "gsl_sf_bessel_Knu_e");
    TEST_SF_DD( s, gsl_sf_bessel_Knu_e, { x: 30.0,  y: 100.0 }, 3.970602055959398739e-43,  TEST_TOL2, "gsl_sf_bessel_Knu_e");

    console.log("    ... gsl_sf_bessel_lnKnu_e");
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 0.0001,  y: 1.0e-100 }, 5.439794449319847,       TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 0.0001,  y: 0.0001   }, 2.232835507214331,       TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 0.0001,  y: 10.0     }, -10.93743282256098,      TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x:  1.0,    y: 1.0e-100 }, 230.2585092994045,       TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x:  1.0,    y: 1.0e-10  }, 23.025850929940456840,   TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x:  1.0,    y: 0.001    }, 6.907751517131146,       TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x:  1.0,    y: 1.0      }, -0.5076519482107523309,  TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 30.0,    y: 1.0e-100 }, 6999.113586185543475,    TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 30.0,    y: 1.0      }, 91.34968784026325464,    TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 30.0,    y: 100.0    }, -97.63224126416760932,   TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 100.0,   y: 1.0e-100 }, 23453.606706185466825,   TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 100.0,   y: 1.0      }, 427.7532510250188083,    TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 100.0,   y: 100.0    }, -55.53422771502921431,   TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 1000.0,  y: 1.0e-100 }, 236856.183755993135,     TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );
    TEST_SF_DD( s, gsl_sf_bessel_lnKnu_e, { x: 10000.0, y: 1.0e-100 }, 2.39161558914890695e+06, TEST_TOL0, "gsl_sf_bessel_lnKnu_e" );

    console.log("    ... gsl_sf_bessel_Jn_array");
    sa.Integer = 0;
    gsl_sf_bessel_Jn_array( 3, 38, 1.0, J );
    if ( test_sf_frac_diff( J[0],  0.0195633539826684059190  ) > TEST_TOL1 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[1],  0.0024766389641099550438  ) > TEST_TOL1 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[10], 1.9256167644801728904e-14 ) > TEST_TOL1 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[35], 6.911232970971626272e-57  ) > TEST_TOL1 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_Jn_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log("    ... gsl_sf_bessel_Yn_array");
    sa.Integer = 0;
    gsl_sf_bessel_Yn_array( 3, 38, 1.0, Y );
    if ( test_sf_frac_diff( Y[0],  -5.821517605964728848      ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[1],  -33.27842302897211870      ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[10], -1.2753618701519837595e+12 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[35], -1.2124435663593357154e+54 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_Yn_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log("    ... gsl_sf_bessel_In_scaled_array");
    sa.Integer = 0;
    gsl_sf_bessel_In_scaled_array(3, 38, 1.0, I);
    if ( test_sf_frac_diff( I[0],  0.0081553077728142938170  ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[1],  0.0010069302573377758637  ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[10], 7.341518665628926244e-15  ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[35], 2.5753065298357542893e-57 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_In_scaled_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_In_array" );
    sa.Integer = 0;
    gsl_sf_bessel_In_array( 3, 38, 1.0, Y );
    if ( test_sf_frac_diff( Y[0],   0.0221684249243319024760  ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[1],   0.0027371202210468663251  ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[10],  1.9956316782072007564e-14 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[35],  7.000408942764452901e-57  ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_In_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_Kn_array" );
    sa.Integer = 0;
    gsl_sf_bessel_Kn_array( 3, 38, 1.0, K );
    if ( test_sf_frac_diff( K[0],  7.101262824737944506 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[1],  44.23241584706284452 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[10], 1.9215763927929940846e+12 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[35], 1.8789385023806051223e+54 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_Kn_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_Kn_scaled_array" );
    sa.Integer = 0;
    gsl_sf_bessel_Kn_scaled_array( 3, 38, 1.0, K );
    if ( test_sf_frac_diff( K[0],  19.303233695596904277 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[1],  120.23617222591483717 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[10], 5.223386190525076473e+12 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[35], 5.107484387813251411e+54 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_Kn_scaled_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_jl_array" );
    sa.Integer = 0;
    gsl_sf_bessel_jl_array( 50, 1.0, J );
    if ( test_sf_frac_diff( J[0],  0.84147098480789650670   ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[1],  0.30116867893975678925   ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[10], 7.116552640047313024e-11 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[50], 3.615274717489787311e-81 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_jl_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_jl_steed_array" );
    sa.Integer = 0;
    gsl_sf_bessel_jl_steed_array( 99, 1.0, J );
    if ( test_sf_frac_diff( J[0],  0.84147098480789650670   ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[1],  0.30116867893975678925   ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[10], 7.116552640047313024e-11 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[50], 3.615274717489787311e-81 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( J[80], 1.136352423414503264e-144 ) > TEST_TOL1 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_jl_steed_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_yl_array" );
    sa.Integer = 0;
    gsl_sf_bessel_yl_array( 50, 1.0, Y );
    if ( test_sf_frac_diff( Y[0],  -0.5403023058681397174 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[1],  -1.3817732906760362241 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[10], -6.722150082562084436e+08  ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( Y[50], -2.7391922846297571576e+78 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_yl_array" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_yl_array");
    var Y0 = [];
    sa.Integer = 0;
    gsl_sf_bessel_yl_array( 0, 1.0, Y0 );
    if ( test_sf_frac_diff( Y0[0],  -0.5403023058681397174 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_yl_array (lmax=0)" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_il_scaled_array" );
    sa.Integer = 0;
    gsl_sf_bessel_il_scaled_array( 50, 1.0, I );
    if ( test_sf_frac_diff( I[0],  0.43233235838169365410 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[1],  0.13533528323661269189 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[10], 2.7343719371837065460e-11 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[50], 1.3429606061892023653e-81 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_il_scaled_array" );
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    gsl_sf_bessel_il_scaled_array( 50, 0.0, I );
    if ( test_sf_frac_diff( I[0],  1.0 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[1],  0.0 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[10], 0.0 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( I[50], 0.0 ) > TEST_TOL2 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_il_scaled_array (L=0)" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_kl_scaled_array" );
    sa.Integer = 0;
    gsl_sf_bessel_kl_scaled_array(50, 1.0, K);
    if ( test_sf_frac_diff( K[0],  1.5707963267948966192     ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[1],  3.1415926535897932385     ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[10], 2.7231075458948147010e+09 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    if ( test_sf_frac_diff( K[50], 1.1578440432804522544e+79 ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_kl_scaled_array" );
    s.Integer = s.Integer + sa.Integer;

    var K0 = [];
    sa.Integer = 0;
    gsl_sf_bessel_kl_scaled_array( 0, 1.0, K0 );
    if ( test_sf_frac_diff( K[0],  1.5707963267948966192     ) > TEST_TOL0 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_bessel_kl_scaled_array (lmax=0)" );
    s.Integer = s.Integer + sa.Integer;

    console.log( "    ... gsl_sf_bessel_zero_J0_e" );
    sa.Integer = 0;
    //sa += ( gsl_sf_bessel_zero_J0_e(0, &r) != GSL_EINVAL );
    //sa += ( r.val != 0.0 );
    s.Integer = s.Integer + sa.Integer;
    TEST_SF_I( s, gsl_sf_bessel_zero_J0_e,   1,  2.404825557695771, TEST_TOL1, "gsl_sf_bessel_zero_J0_e" );
    TEST_SF_I( s, gsl_sf_bessel_zero_J0_e,   2,  5.520078110286304, TEST_TOL1, "gsl_sf_bessel_zero_J0_e" );
    TEST_SF_I( s, gsl_sf_bessel_zero_J0_e,  20, 62.048469190227081, TEST_TOL1, "gsl_sf_bessel_zero_J0_e" );
    TEST_SF_I( s, gsl_sf_bessel_zero_J0_e,  25, 77.756025630388058, TEST_TOL1, "gsl_sf_bessel_zero_J0_e" );
    TEST_SF_I( s, gsl_sf_bessel_zero_J0_e, 100, 313.37426607752784, TEST_TOL1, "gsl_sf_bessel_zero_J0_e" );

    console.log( "    ... gsl_sf_bessel_zero_J1_e" );
    sa.Integer = 0;
    //sa += ( gsl_sf_bessel_zero_J1_e(0, &r) != GSL_SUCCESS );
    //sa += ( r.val != 0.0 );
    s.Integer = s.Integer + sa.Integer;
    TEST_SF_I( s, gsl_sf_bessel_zero_J1_e,   1, 3.831705970207512, TEST_TOL2, "gsl_sf_bessel_zero_J1_e" );
    TEST_SF_I( s, gsl_sf_bessel_zero_J1_e,   2, 7.015586669815619, TEST_TOL2, "gsl_sf_bessel_zero_J1_e" );
    TEST_SF_I( s, gsl_sf_bessel_zero_J1_e,  20, 63.61135669848124, TEST_TOL2, "gsl_sf_bessel_zero_J1_e" );
    TEST_SF_I( s, gsl_sf_bessel_zero_J1_e,  25, 79.32048717547630, TEST_TOL2, "gsl_sf_bessel_zero_J1_e" );
    TEST_SF_I( s, gsl_sf_bessel_zero_J1_e, 100, 314.9434728377672, TEST_TOL2, "gsl_sf_bessel_zero_J1_e" );

    console.log( "    ... gsl_sf_bessel_zero_Jnu_e");
    sa.Integer = 0;
    //sa += ( gsl_sf_bessel_zero_Jnu_e(0.0, 0, &r) != GSL_EINVAL );
    //sa += ( r.val != 0.0 );
    s.Integer = s.Integer + sa.Integer;
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 0.0, i:   1 },  2.404825557695771, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 0.0, i:   2 },  5.520078110286304, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 0.0, i:  20 }, 62.048469190227081, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 0.0, i:  25 }, 77.756025630388058, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 0.0, i: 100 }, 313.37426607752784, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );

    sa.Integer = 0;
    //sa += ( gsl_sf_bessel_zero_Jnu_e(1.0, 0, &r) != GSL_SUCCESS );
    //sa += (r.val != 0.0);
    s.Integer = s.Integer + sa.Integer;
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 1.5, i: 1 },  4.4934094579090641, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 5.0, i: 1 },  8.7714838159599540, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 1.5, i: 2 },  7.7252518369377072, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 5.0, i: 2 },  12.338604197466944, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 1.5, i: 3 },  10.904121659428900, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 5.0, i: 3 },  15.700174079711671, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 1.5, i: 4 },  14.066193912831473, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 5.0, i: 4 },  18.980133875179921, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 1.5, i: 5 },  17.220755271930768, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );

    // Something wrong with the tolerances on these
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 5.0, i: 5 },  22.217799896561268, TEST_SQRT_TOL0, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 8.0, i: 5 },  26.266814641176644, TEST_SQRT_TOL0, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:20.0, i: 5 },  41.413065513892636, TEST_SQRT_TOL0, "gsl_sf_bessel_zero_Jnu_e" );

    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  1.5, i: 6 },  20.371302959287563, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 6 },  25.430341154222704, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  8.0, i: 6 },  29.545659670998550, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  1.5, i: 7 },  23.519452498689007, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 7 },  28.626618307291138, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  8.0, i: 7 },  32.795800037341462, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  1.5, i: 8 },  26.666054258812674, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 8 },  31.811716724047763, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 10.0, i: 8 },  38.761807017881651, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  1.5, i: 9 },  29.811598790892959, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 9 },  34.988781294559295, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 10.0, i: 9 },  42.004190236671805, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  1.5, i: 10 },  32.956389039822477, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 10 },  38.159868561967132, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 10 },  52.017241278881633, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );

    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 11 }, 41.326383254047406, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 11 }, 55.289204146560061, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 12 }, 44.4893191232197314, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 12 }, 58.5458289043850856, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 13 }, 47.6493998066970948, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 13 }, 61.7897598959450550, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 14 }, 50.8071652030063595, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 14 }, 65.0230502510422545, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 15 }, 53.9630265583781707, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 15 }, 68.2473219964207837, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 16 }, 57.1173027815042647, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 16 }, 71.4638758850226630, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 17 }, 60.2702450729428077, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 17 }, 74.6737687121404241, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 18 }, 63.4220540458757799, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 18 }, 77.8778689734863729, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 19 }, 66.5728918871182703, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 19 }, 81.0768977206328326, TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  5.0, i: 20 }, 69.722891161716742,  TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 15.0, i: 20 }, 84.271459069716442,  TEST_TOL1, "gsl_sf_bessel_zero_Jnu_e" );

    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  23.0, i: 11 }, 65.843393469524653, TEST_TOL6, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  30.0, i: 11 }, 74.797306585175426, TEST_TOL6, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  32.0, i: 15 }, 90.913637691861741, TEST_TOL6, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  50.0, i: 15 }, 113.69747988073942, TEST_TOL6, "gsl_sf_bessel_zero_Jnu_e" );

    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:   5.0, i: 22 }, 76.020793430591605, TEST_TOL2, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  10.0, i: 22 }, 83.439189796105756, TEST_TOL3, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x:  12.0, i: 22 }, 86.345496520534055, TEST_TOL6, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 100.0, i: 22 }, 199.82150220122519, TEST_TOL4, "gsl_sf_bessel_zero_Jnu_e" );
    TEST_SF_DI( s, gsl_sf_bessel_zero_Jnu_e, { x: 500.0, i: 22 }, 649.34132440891735, TEST_TOL2, "gsl_sf_bessel_zero_Jnu_e" );

    console.log( "    ... gsl_sf_bessel_sequence_Jnu_e" );
    sa.Integer = 0;
    for ( let i = 0; i <= 100 - 1; i++ )
    {
        J[i] = (i) / 10.0;
    }
    gsl_sf_bessel_sequence_Jnu_e( 2.0, GSL_MODE_DEFAULT, 100, J );
    console.log( "J[1] = ", J[1] );
    if ( test_sf_frac_diff( J[1],  0.0012489586587999188454 ) > TEST_TOL1 )
    {
        console.log( "expected: ", 0.0012489586587999188454 );
        console.log( "obtained: ", J[1] );
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[20] = ", J[20] );
    if ( test_sf_frac_diff( J[20], 0.3528340286156377192 ) > TEST_TOL4 )
    {
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[50] = ", J[50] );
    if ( test_sf_frac_diff( J[50], 0.0465651162777522155 ) > TEST_TOL4 )
    {
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_sequence_Jnu_e(2)" );
    s.Integer = s.Integer + sa.Integer;

    sa.Integer = 0;
    for ( let i = 0; i <= 100 - 1; i++ )
    {
        J[i] = i;
    }
    gsl_sf_bessel_sequence_Jnu_e( 12.0, GSL_MODE_DEFAULT, 100, J );
    console.log( "J[1] = ", J[1] );
    if ( test_sf_frac_diff( J[1],  4.999718179448405289e-13 ) > TEST_TOL1 )
    {
        console.log( "expected: ", 4.999718179448405289e-13 );
        console.log( "obtained: ", J[1] );
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[5] = ", J[5] );
    if ( test_sf_frac_diff( J[5],  7.627813166084551355e-05 ) > TEST_TOL1 )
    {
        console.log( "expected: ", 7.627813166084551355e-05 );
        console.log( "obtained: ", J[5] );
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[7] = ", J[7] );
    if ( test_sf_frac_diff( J[7],  2.655620035894568062e-03 ) > TEST_TOL3 )
    {
        console.log( "expected: ", 2.655620035894568062e-03 );
        console.log( "obtained: ", J[7] );
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[10] = ", J[10] );
    if ( test_sf_frac_diff( J[10],  6.337025497015601509e-02 ) > TEST_TOL3 )
    {
        console.log( "expected: ", 6.337025497015601509e-02 );
        console.log( "obtained: ", J[10] );
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[15] = ", J[15] );
    if ( test_sf_frac_diff( J[15],  0.23666584405476805591 ) > TEST_TOL3 )
    {
        console.log( "expected: ", 0.23666584405476805591 );
        console.log( "obtained: ", J[15] );
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[30] = ", J[30] );
    if ( test_sf_frac_diff( J[30],  0.14825335109966010021 ) > TEST_TOL3 )
    {
        console.log( "expected: ", 0.14825335109966010021 );
        console.log( "obtained: ", J[30] );
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[70] = ", J[70] );
    if ( test_sf_frac_diff( J[70],  0.04109913716555364262 ) > TEST_TOL4 )
    {
        console.log( "expected: ", 0.04109913716555364262 );
        console.log( "obtained: ", J[70] );
        sa.Integer = sa.Integer + 1;
    }
    console.log( "J[99] = ", J[99] );
    if ( test_sf_frac_diff( J[99],  -0.0015052760501176728 ) > TEST_TOL5 )
    {
        console.log( "expected: ", -0.0015052760501176728 );
        console.log( "obtained: ", J[99] );
        sa.Integer = sa.Integer + 1;
    }
    gsl_test( sa, "  gsl_sf_sequence_Jnu_e(12)" );
    s.Integer = s.Integer + sa.Integer;

    // sa.Integer = 0;
    // for ( let i = 0; i <= 100 - 1; i++ )
    // {
    //     J[i] = i * 20.0;
    // }
    // gsl_sf_bessel_sequence_Jnu_e( 1000.0, GSL_MODE_DEFAULT, 100, J );
    // console.log( "J[50] = ", J[50] );
    // if ( test_sf_frac_diff( J[50],  0.04473067294796404088 ) > TEST_TOL5 )
    // {
    //     console.log( "expected: ", 0.04473067294796404088 );
    //     console.log( "obtained: ", J[50] );
    //     sa.Integer = sa.Integer + 1;
    // }
    // console.log( "J[99] = ", J[99] );
    // if ( test_sf_frac_diff( J[99],  0.01400619760349853902 ) > TEST_TOL6 )
    // {
    //     console.log( "expected: ", 0.01400619760349853902 );
    //     console.log( "obtained: ", J[99] );
    //     sa.Integer = sa.Integer + 1;
    // }
    // gsl_test( sa, "  gsl_sf_sequence_Jnu_e(1000)" );
    // s.Integer = s.Integer + sa.Integer;

    return s;

} // test_bessel

// ****************************************************************************

export function test_hyperg()
{

    var s = { Integer: 0 };
    var sgn = 0.0;

    console.log("Test Hypergeometric Functions ...");

    // 0F1

    console.log("    ... gsl_sf_hyperg_0F1_e");
    TEST_SF_DD(s, gsl_sf_hyperg_0F1_e, { x:    1.0, y:  0.5 }, 1.5660829297563505373, TEST_TOL0, "gsl_sf_hyperg_0F1_e");
    TEST_SF_DD(s, gsl_sf_hyperg_0F1_e, { x:    5.0, y:  0.5 }, 1.1042674404828684574, TEST_TOL1, "gsl_sf_hyperg_0F1_e");
    TEST_SF_DD(s, gsl_sf_hyperg_0F1_e, { x:  100.0, y: 30.0 }, 1.3492598639485110176, TEST_TOL2, "gsl_sf_hyperg_0F1_e");
    TEST_SF_DD(s, gsl_sf_hyperg_0F1_e, { x:   -0.5, y:  3.0 }, -39.29137997543434276,  TEST_TOL1, "gsl_sf_hyperg_0F1_e");
    TEST_SF_DD(s, gsl_sf_hyperg_0F1_e, { x: -100.5, y: 50.0 }, 0.6087930289227538496, TEST_TOL3, "gsl_sf_hyperg_0F1_e");
    TEST_SF_DD(s, gsl_sf_hyperg_0F1_e, { x:    1.0, y: -5.0 }, -0.3268752818235339109, TEST_TOL0, "gsl_sf_hyperg_0F1_e");
    TEST_SF_DD(s, gsl_sf_hyperg_0F1_e, { x:   -0.5, y: -5.0 }, -4.581634759005381184,  TEST_TOL1, "gsl_sf_hyperg_0F1_e");

    // 1F1 for integer parameters

    console.log( "    ... gsl_sf_hyperg_1F1_int_e" );
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 1, m: 1, x: 0.5 }, 1.6487212707001281468, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e" );
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 1, m: 2, x: 500.0 }, 2.8071844357056748215e+214, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e" );
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 1, m: 2, x: -500.0 }, 0.002, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e" );
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 8, m: 1, x: 0.5 }, 13.108875178030540372, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e" );

    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 1, x: 1.0 },  131.63017574352619931, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 1, x: 10.0 }, 8.514625476546280796e+09, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 1, x: 100.0 },  1.5671363646800353320e+56, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 20, x: 1.0 },  1.6585618002669675465, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 20, x: 10.0 },  265.26686430340188871, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 20, x: 100.0 }, 3.640477355063227129e+34, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 100, x: 1.0 },  1.1056660194025527099, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 100, x: 10.0 },  2.8491063634727594206, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 100, x: 40.0 },  133.85880835831230986, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 100, x: 80.0 },  310361.16228011433406, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 100, x: 100.0 },  8.032171336754168282e+07, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10, m: 100, x: 500.0 },  7.633961202528731426e+123, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 1, x: 1.0 },  6.892842729046469965e+07, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 1, x: 10.0 },  2.4175917112200409098e+28, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 1, x: 100.0 },  1.9303216896309102993e+110, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 200, x: 1.0 },  1.6497469106162459226, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 200, x: 10.0 },  157.93286197349321981, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 200, x: 100.0 },  2.1819577501255075240e+24, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 200, x: 400.0 },  3.728975529926573300e+119, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 400, x: 10.0 },  12.473087623658878813, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 400, x: 100.0 },  9.071230376818550241e+11, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 400, x: 150.0 },  7.160949515742170775e+18, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 400, x: 200.0 },   2.7406690412731576823e+26, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 400, x:300.0 },  6.175110613473276193e+43, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 400, x: 400.0 },  1.1807417662711371440e+64, TEST_TOL3, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 400, x: 600.0 },  2.4076076354888886030e+112, TEST_TOL3, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10,  m: 1, x: -1.0 },      0.11394854824644542810,   TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10,  m: 1, x: -10.0 },     0.0006715506365396127863, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10,  m: 1, x: -100.0 },   -4.208138537480269868e-32, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10,  m: 50, x: -1.0 },     0.820006196079380,        TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10,  m: 100, x: -10.0 },   0.38378859043466243,      TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10,  m: 100, x: -100.0 },  0.0008460143401464189061, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10,  m: 100, x: -500.0 },  1.1090822141973655929e-08, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 10,  m: 100, x: -10000.0 }, 5.173783508088272292e-21, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 50,  m: 1, x: -90.0 },    -1.6624258547648311554e-21, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 50,  m: 1, x: -100.0 },    4.069661775122048204e-24, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 50,  m: 1, x: -110.0 },    1.0072444993946236025e-25, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 10, x: -100.0 }, -2.7819353611733941962e-37, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 1, x: -90.0 },  7.501705041159802854e-22, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 1, x: -100.0 },  6.305128893152291187e-25, TEST_TOL3, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 1, x: -110.0 },  -7.007122115422439755e-26, TEST_TOL3, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 100, m: 10, x: -100.0 },  -2.7819353611733941962e-37, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 200, m: 50, x: -1.0 },  0.016087060191732290813, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 200, m: 50, x: -300.0 },  -4.294975979706421471e-121, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 200, m: 100, x: -1.0 },  0.13397521083325179687, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 200, m: 100, x: -10.0 },  5.835134393749807387e-10, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 200, m: 100, x: -100.0 },  4.888460453078914804e-74, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: 200, m: 100, x: -500.0 },  -1.4478509059582015053e-195, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");

    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -1, m: 1, x: 2.0 },  -1.0, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -1, m: -2, x: 2.0 },  2.0, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -2, m: -3, x: 2.0 },  3.0, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: 1, x: 1.0 },  0.4189459325396825397, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: 1, x: 10.0 },  27.984126984126984127, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: 1, x: 100.0 },  9.051283795429571429e+12, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -100, m: 20, x: 1.0 },  0.0020203016320697069566, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -20, x: 1.0 },  1.6379141878548080173, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -20, x: 10.0 },  78.65202404521289970, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -20, x: 100.0 },  4.416169713262624315e+08, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -100, x: 1.0 },  1.1046713999681950919, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -100, x: 10.0 },  2.6035952191039006838, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -100, x: 100.0 },  1151.6852040836932392, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -100, m: -200, x: 1.0 },  1.6476859702535324743, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -100, m: -200, x: 10.0 },  139.38026829540687270, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -100, m: -200, x: 100.0 },  1.1669433576237933752e+19, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -20, x: -1.0 },  0.6025549561148035735, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -20, x: -10.0 },  0.00357079636732993491, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -20, x: -100.0 },  1.64284868563391159e-35, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -100, x: -1.0 },  0.90442397250313899, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -100, x: -10.0 },  0.35061515251367215, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -10, m: -100, x: -100.0 },  8.19512187960476424e-09, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -100, m: -200, x: -1.0 },  0.6061497939628952629, TEST_TOL0, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -100, m: -200, x: -10.0 },  0.0063278543908877674, TEST_TOL1, "gsl_sf_hyperg_1F1_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_1F1_int_e, { n: -100, m: -200, x: -100.0 },  4.34111795007336552e-25, TEST_TOL2, "gsl_sf_hyperg_1F1_int_e");

    // 1F1

    console.log( "    ... gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0,  y: 1.5, z: 1.0 }, 2.0300784692787049755, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0,  y: 1.5, z: 10.0 },  6172.859561078406855, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0,  y: 1.5, z: 100.0 },  2.3822817898485692114e+42, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0,  y: 1.5, z: 500.0 },  5.562895351723513581e+215, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.5,  y: 2.5, z: 1.0 }, 1.8834451238277954398, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.5,  y: 2.5, z: 10.0 },  3128.7352996840916381, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 1.1, z: 1.0 },  110.17623733873889579, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 1.1, z: 10.0 },  6.146657975268385438e+09, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 1.1, z: 100.0 }, 9.331833897230312331e+55, TEST_TOL2, "gsl_sf_hyperg_1F1_e");

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 1.1,  z: 500.0 },  4.519403368795715843e+235, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 50.1, z:  2.0 },  1.5001295507968071788, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 50.1, z:  10.0 },  8.713385849265044908, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 50.1, z:  100.0 },  5.909423932273380330e+18, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 50.1, z:  500.0 },  9.740060618457198900e+165, TEST_TOL2, "gsl_sf_hyperg_1F1_e");

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: 1.0 },  5.183531067116809033e+07, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: 10.0 },  1.6032649110096979462e+28, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: 100.0 },  1.1045151213192280064e+110, TEST_TOL2, "gsl_sf_hyperg_1F1_e");

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 50.1, z: 1.0 },  7.222953133216603757, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 50.1, z: 10.0 },  1.0998696410887171538e+08, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 50.1, z: 100.0 },  7.235304862322283251e+63, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0, y: 1.5, z: -1.0 }, 0.5380795069127684191, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0, y: 1.5, z: -10.0 },  0.05303758099290164485, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0, y: 1.5, z: -100.0 }, 0.005025384718759852803, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0, y: 1.5, z: -500.0 }, 0.0010010030151059555322, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0, y: 1.1, z: -500.0 }, 0.00020036137599690208265, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 1.1, z: -1.0 },     0.07227645648935938168, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 1.1, z: -10.0 },    0.0003192415409695588126, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    //
    // sensitive to the pair_ratio hack in hyperg_1F1.c
    // TEST_SF_RLX(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -100, &r),  -8.293425316123158950e-16, 50.0*TEST_SNGL);
    //
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 10.0, y: 1.1, z: -500.0 },  -3.400379216707701408e-23, TEST_TOL2, "gsl_sf_hyperg_1F1_e");

    //TEST_SF_RLX(s, gsl_sf_hyperg_1F1_e, (50, 1.1, -90, &r),   -7.843129411802921440e-22, TEST_SQRT_TOL0);
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 50.0, y: 1.1, z: -100.0 },   4.632883869540640460e-24, TEST_SQRT_TOL0, "gsl_sf_hyperg_1F1_e");

    // FIXME:
    // tolerance is poor, but is consistent within reported error
    //
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 50.0, y: 1.1, z: -110.0 }, 5.642684651305310023e-26, 0.03, "gsl_sf_hyperg_1F1_e");

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: -1.0 },    0.0811637344096042096, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: -10.0 },   0.00025945610092231574387, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: -50.0 },   2.4284830988994084452e-13, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: -90.0 },   2.4468224638378426461e-22, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: -99.0 },   1.0507096272617608461e-23, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: -100.0 },  1.8315497474210138602e-24, TEST_TOL2, "gsl_sf_hyperg_1F1_e");

    // FIXME:
    // Reported error is too small.
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: -101.0 }, -2.3916306291344452490e-24, 0.04, "gsl_sf_hyperg_1F1_e");

    // FIXME:
    // Reported error is too small.
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 1.1, z: -110.0 }, -4.517581986037732280e-26, TEST_TOL0, "gsl_sf_hyperg_1F1_e");

    // FIXME:
    // Result is terrible, but reported error is very large, so consistent.
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100.0, y: 10.1, z: -220.0 }, -4.296130300021696573e-64, TEST_TOL1, "gsl_sf_hyperg_1F1_e");

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -10.0, y: -10.1, z: 10.0 }, 10959.603204633058116, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -10.0, y: -10.1, z: 1000.0 }, 2.0942691895502242831e+23, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -10.0, y: -100.1, z: 10.0 },  2.6012036337980078062, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -1000.0, y: -1000.1, z: 10.0 },  22004.341698908631636, TEST_TOL3, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -1000.0, y: -1000.1, z: 200.0 }, 7.066514294896245043e+86, TEST_TOL3, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -8.1, y: -10.1, z: -10.0 },    0.00018469685276347199258, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -8.1, y: -1000.1, z: -10.0 },  0.9218280185080036020, TEST_TOL0, "gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -10,  y: -5.1, z: 1 },  16.936141866089601635, TEST_TOL2, "gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -10,  y: -5.1, z: 10 },  771534.0349543820541, TEST_TOL2, "gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -10,  y: -5.1, z: 100 },  2.2733956505084964469e+17, TEST_TOL2, "gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100, y: -50.1, z: -1 },  0.13854540373629275583, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100, y: -50.1, z: -10 },  -9.142260314353376284e+19, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100, y: -50.1, z: -100 },  -1.7437371339223929259e+87, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100, y: -50.1, z: 1 },  7.516831748170351173, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100, y: -50.1, z: 10 },  1.0551632286359671976e+11, TEST_SQRT_TOL0, "gsl_sf_hyperg_1F1_e" );

    // These come out way off. On the other hand, the error estimates
    // are also very large; so much so that the answers are consistent
    // within the reported error. Something will need to be done about
    // this eventually
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100, y: -50.1, z: 50  },  -7.564755600940346649e+36, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100, y: -50.1, z: 100 },  4.218776962675977e+55, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -10.5,  y: -8.1,   z: 0.1   },  1.1387201443786421724, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -10.5,  y: -11.1,  z: 1.0   },  2.5682766147138452362, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.5, y: -80.1,  z: 10.0  },  355145.4517305220603, TEST_TOL3, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.5, y: -102.1, z: 10.0  },  18678.558725244365016, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.5, y: -500.1, z: 10.0  },  7.342209011101454, TEST_TOL0, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.5, y: -500.1, z: 100.0 },  1.2077443075367177662e+8, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -500.5, y: -80.1,  z: 2.0   },  774057.8541325341699, TEST_TOL4, "gsl_sf_hyperg_1F1_e" );

    // UNIMPL
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100, y: -10.1, z: 1 },  -2.1213846338338567395e+12, TEST_TOL0, "gsl_sf_hyperg_1F1_e" );
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100, y: -10.1, z: 10 },  -6.624849346145112398e+39, TEST_TOL0, "gsl_sf_hyperg_1F1_e" );
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100, y: -10.1, z: 100 },  -1.2413466759089171904e+129, TEST_TOL0, "gsl_sf_hyperg_1F1_e" );

    // UNIMPL
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100, y: -10.1, z: -1 },  34456.29405305551691, TEST_TOL0, "gsl_sf_hyperg_1F1_e" );
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100, y: -10.1, z: -10 },  -7.809224251467710833e+07, TEST_TOL0, "gsl_sf_hyperg_1F1_e" );
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 100, y: -10.1, z: -100 },   -5.214065452753988395e-07, TEST_TOL0, "gsl_sf_hyperg_1F1_e" );

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 1.1,  z: 1.0 },  0.21519810496314438414, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 1.1,  z: 10.0 },  8.196123715597869948, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 1.1,  z: 100.0 },  -1.4612966715976530293e+20, TEST_TOL1, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 20.1, z:  1.0 },  0.0021267655527278456412, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 20.1, z:  10.0 },   2.0908665169032186979e-11, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 20.1, z:  100.0 },  -0.04159447537001340412, TEST_TOL2, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 1.1,  z: -1.0 },  2.1214770215694685282e+07, TEST_TOL3, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 1.1,  z: -10.0 },  1.0258848879387572642e+24, TEST_TOL3, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 1.1,  z: -100.0 },  1.1811367147091759910e+67, TEST_TOL3, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 50.1, z:  -1.0 },  6.965259317271427390, TEST_TOL3, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 50.1, z:  -10.0 },  1.0690052487716998389e+07, TEST_TOL3, "gsl_sf_hyperg_1F1_e");
    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -100.0, y: 50.1, z:  -100.0 },  6.889644435777096248e+36, TEST_TOL3, "gsl_sf_hyperg_1F1_e");

    // Bug report from Fernando Pilotto

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -2.05, y: 1.0, z: 5.05 }, 3.79393389516785e+00, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );

    // Bug reports from Ivan Liu

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -26.0, y:  2.0, z: 100.0 }, 1.444786781107436954e+19, TEST_TOL3, "gsl_sf_hyperg_1F1_e");

    // FIXME
    // This one is computed with a huge error, there is loss of
    // precision but the error estimate flags the problem (assuming the
    // user looks at it).  We should probably trap any return with
    // err>|val| and signal loss of precision */
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -26.1, y: 2.0, z: 100.0 }, 1.341557199575986995e+19, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );

    // Bug report H.Moseby

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.2, y: 1.1e-15, z: 1.5 }, 8254503159672429.02, TEST_TOL3, "gsl_sf_hyperg_1F1_e");

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0, y: 1000000.5, z: 0.8e6 + 0.5 }, 4.999922505099443804e+00, TEST_TOL3, "gsl_sf_hyperg_1F1_e");

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.0, y: 1000000.5, z: 1001000.5 }, 3480.3699557431856166, TEST_TOL4, "gsl_sf_hyperg_1F1_e");

    // FIXME -- FIX THESE NEXT RELEASE
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.1, y: 1000000.5, z: 1001000.5 }, 7304.6126942641350122, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );
    // TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 0.9, y: 1000000.5, z: 1001000.5 }, 1645.4879293475410982, TEST_TOL3, "gsl_sf_hyperg_1F1_e" );

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: -1.1,y:  1000000.5, z:  1001000.5 }, -5.30066488697455e-04, TEST_TOL3, "gsl_sf_hyperg_1F1_e");

    TEST_SF_3D( s, gsl_sf_hyperg_1F1_e, { x: 1.5, y: 1000000.5,  z: 0.8e6 + 0.5 }, 11.18001288977894650469927615, TEST_TOL4, "gsl_sf_hyperg_1F1_e");

    // U for integer parameters

    console.log( "    ... gsl_sf_hyperg_U_int_e" );
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 1, x: 0.0001 },  8.634088070212725330, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 1, x: 0.01 },  4.078511443456425847, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 1, x: 0.5 },  0.9229106324837304688, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 1, x: 2.0 },  0.3613286168882225847, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 1, x: 100.0 },  0.009901942286733018406, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 1, x: 1000.0 },  0.0009990019940238807150, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 8, x: 0.01 },  7.272361203006010000e+16, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 8, x: 1.0 },  1957.0, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 8, x: 5.0 },  1.042496, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 8, x: 8.0 },  0.3207168579101562500, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 8, x: 50.0 },  0.022660399001600000000, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 8, x: 100.0 },  0.010631236727200000000, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 8, x: 1000.0 },  0.0010060301203607207200, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 20, x: 1.0 },  1.7403456103284421000e+16, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 20, x: 20.0 },  0.22597813610531052969, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 50, x: 1.0 },  3.374452117521520758e+61, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 50, x: 50.0 },  0.15394136814987651785, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 100, x: 0.1 },  1.0418325171990852858e+253, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 100, x: 1.0 },  2.5624945006073464385e+154, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 100, x: 50.0 },  3.0978624160896431391e+07, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 100, x: 100.0 },  0.11323192555773717475, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 100, x: 200.0 },  0.009715680951406713589, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 100, x: 1000.0 },  0.0011085142546061528661, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: 1000, x: 2000.0 },  0.0009970168547036318206, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -1, x: 1.0 },  0.29817368116159703717, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -1, x: 10.0 },  0.07816669698940409380, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -10, x: 1.0 },  0.08271753756946041959, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -10, x: 5.0 },  0.06127757419425055261, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -10, x: 10.0 },  0.04656199948873187212, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -10, x: 20.0 },  0.031606421847946077709, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -100, x: 0.01 },  0.009900000099999796950, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -100, x: 1.0 },  0.009802970197050404429, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -100, x: 10.0 },  0.009001648897173103447, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -100, x: 20.0 },  0.008253126487166557546, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -100, x: 50.0 },  0.006607993916432051008, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -100, x: 90.0 },  0.005222713769726871937, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -100, x: 110.0 },  0.004727658137692606210, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
// ----  TEST_SF_RLX(s, gsl_sf_hyperg_U_int_e, (1, -1000, 1, &r),   0.0009980029970019970050, TEST_SQRT_TOL0);
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 1, m: -1000, x: 1010.0 },  0.0004971408839859245170, TEST_TOL4, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 8, m: 1, x: 0.001 },  0.0007505359326875706975, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 8, m: 1, x: 0.5 },  6.449509938973479986e-06, TEST_TOL3, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 8, m: 1, x: 8.0 },  6.190694573035761284e-10, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 8, m: 1, x: 20.0 },  3.647213845460374016e-12, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 8, m: 8, x: 1.0 },  0.12289755012652317578, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 8, m: 8, x: 10.0 },  5.687710359507564272e-09, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 8, m: 8, x: 20.0 },  2.8175404594901039724e-11, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 100, m: 100, x: 0.01 },  1.0099979491941914867e+196, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 100, m: 100, x: 0.1 },  1.0090713562719862833e+97, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 100, m: 100, x: 1.0 },  0.009998990209084729106, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: 100, m: 100, x: 20.0 },  1.3239363905866130603e-131, TEST_TOL2, "gsl_sf_hyperg_U_int_e");

    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 1, x: 0.01 },  3.274012540759009536e+06, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 1, x: 1.0 },  1.5202710000000000000e+06, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 1, x: 10.0 },  1.0154880000000000000e+08, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 1, x: 100.0 },  3.284529863685482880e+19, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 10, x: 1.0 },  1.1043089864100000000e+11, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 100, x: 1.0 },  1.3991152402448957897e+20, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 100, x: 10.0 },  5.364469916567136000e+19, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 100, x: 100.0 },  3.909797568000000000e+12, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -10, m: 100, x: 500.0 },  8.082625576697984130e+25, TEST_TOL0, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: 1, x: 0.01 },  1.6973422555823855798e+64, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: 1, x: 1.0 },  7.086160198304780325e+63, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: 1, x: 10.0 },  5.332862895528712200e+65, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: 10, x: 1.0 },  -7.106713471565790573e+71, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: 100, x: 1.0 },  2.4661377199407186476e+104, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: 10, x: 10.0 },  5.687538583671241287e+68, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: 100, x: 10.0 },  1.7880761664553373445e+102, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -90, m: 1, x: 0.01 },  4.185245354032917715e+137, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -90, m: 1, x: 0.1 },  2.4234043408007841358e+137, TEST_TOL3, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -90, m: 1, x: 10.0 },  -1.8987677149221888807e+139, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -90, m: 10, x: 10.0 },  -5.682999988842066677e+143, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -90, m: 100, x: 10.0 },  2.3410029853990624280e+189, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -90, m: 1000, x: 10.0 },  1.9799451517572225316e+271, TEST_TOL3, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: -1, x: 10.0 },  -9.083195466262584149e+64, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: -10, x: 10.0 },  -1.4418257327071634407e+62, TEST_TOL1, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: -100, x: 0.01 },  3.0838993811468983931e+93, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -50, m: -100, x: 10.0 },  4.014552630378340665e+95, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -100, m: -100, x: 10.0 },  2.0556466922347982030e+162, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -100, m: -200, x: 10.0 },  1.1778399522973555582e+219, TEST_TOL2, "gsl_sf_hyperg_U_int_e");
    TEST_SF_IID( s, gsl_sf_hyperg_U_int_e, { n: -100, m: -200, x: 100.0 },  9.861313408898201873e+235, TEST_TOL3, "gsl_sf_hyperg_U_int_e");

    // U

    console.log( "    ... gsl_sf_hyperg_U_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 0.0001, z: 0.0001 }, 1.0000576350699863577, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 0.0001, z: 1.0 }, 0.9999403679233247536, TEST_TOL0, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 0.0001, z: 100.0 }, 0.9995385992657260887, TEST_TOL0, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 1.0, z: 0.0001 }, 1.0009210608660065989, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 1.0, z: 1.0 }, 0.9999999925484179084, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 10.0, z: 1.0 }, 13.567851006281412726, TEST_TOL3, "gsl_sf_hyperg_U_e");
// ----  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 10, 5, &r), 1.0006265020064596364, TEST_SQRT_TOL0);
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 10.0, z: 10.0 }, 0.9999244381454633265, TEST_TOL0, "gsl_sf_hyperg_U_e");
// ----  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 1, &r),  2.5890615708804247881e+150, TEST_SQRT_TOL0);
// ----  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 10, &r),  2.3127845417739661466e+55, TEST_SQRT_TOL0);
// ----  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 50, &r), 6402.818715083582554, TEST_SQRT_TOL0);
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 100.0, z: 98.0 }, 0.9998517867411840044, TEST_TOL2, "gsl_sf_hyperg_U_e");
// ----  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 300, &r),  2.5389557274938010716e+213, TEST_SQRT_TOL0);
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 1000.0, z: 999.0 }, 0.9997195294193261604, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.0001, y: 1000.0, z: 1100.0 },  0.9995342990014584713, TEST_TOL1, "gsl_sf_hyperg_U_e");
// ----  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.5, 1000, 300, &r), 1.1977955438214207486e+217, TEST_SQRT_TOL0);
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.5, y: 1000.0, z: 800.0 }, 9.103916020464797207e+08, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.5, y: 1000.0, z: 998.0 }, 0.21970269691801966806, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 0.5, y: 0.5, z: 1.0 }, 0.7578721561413121060, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 1.0, y: 0.0001, z: 0.0001 }, 0.9992361337764090785, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 1.0, y: 0.0001, z: 1.0 }, 0.4036664068111504538, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 1.0, y: 0.0001, z: 100.0 }, 0.009805780851264329587, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 1.0, y: 1.2, z: 2.0 }, 0.3835044780075602550, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 1.0, y: -0.0001, z: 1.0 }, 0.4036388693605999482, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 8.0, y: 10.5, z: 1.0 },  27.981926466707438538, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 8.0, y: 10.5, z: 10.0 },  2.4370135607662056809e-8, TEST_TOL0, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 8.0, y: 10.5, z: 100.0 },  1.1226567526311488330e-16, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 10.0, y: -2.5, z: 10.0 },  6.734690720346560349e-14, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 10.0, y: 2.5, z: 10.0 },  6.787780794037971638e-13, TEST_TOL0, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: 10.0, y: 2.5, z: 50.0 },  2.4098720076596087125e-18, TEST_TOL0, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 1.1, z: 1.0 },  -3.990841457734147e+6, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 1.1, z: 10.0 },  1.307472052129343e+8, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 1.1, z: 50.0 },  3.661978424114088e+16, TEST_TOL0, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 1.1, z: 90.0 },  8.09469542130868e+19, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 1.1, z: 99.0 },  2.546328328942063e+20, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 1.1, z: 100.0 },  2.870463201832814e+20, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 1.1, z: 200.0 },  8.05143453769373e+23, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 10.1, z: 0.1 },  -3.043016255306515e+20, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 10.1, z: 1.0 },  -3.194745265896115e+12, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 10.1, z: 4.0 },  -6.764203430361954e+07, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 10.1, z: 10.0 },  -2.067399425480545e+09, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 10.1, z: 50.0 },  4.661837330822824e+14, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 100.4, z: 10.0 },  -6.805460513724838e+66, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 100.4, z: 50.0 },  -2.081052558162805e+18, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 100.4, z: 80.0 },  2.034113191014443e+14, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 100.4, z: 100.0 },  6.85047268436107e+13, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -10.5, y: 100.4, z: 200.0 },  1.430815706105649e+20, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -19.5, y: 82.1, z: 10.0 },  5.464313196201917432e+60, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -50.5, y: 100.1, z: 10.0 },  -5.5740216266953e+126, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -50.5, y: 100.1, z: 40.0 },  5.937463786613894e+91, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -50.5, y: 100.1, z: 50.0 },  -1.631898534447233e+89, TEST_TOL1, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -50.5, y: 100.1, z: 70.0 },  3.249026971618851e+84, TEST_TOL2, "gsl_sf_hyperg_U_e");
    TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -50.5, y: 100.1, z: 100.0 },  1.003401902126641e+85, TEST_TOL1, "gsl_sf_hyperg_U_e");

    // Bug report from Stefan Gerlach

    // FIXME
    // TEST_SF_3D( s, gsl_sf_hyperg_U_e, { x: -2.0, y: 4.0, z: 1.0 },  11.0, TEST_TOL0, "gsl_sf_hyperg_U_e" );

    // 2F1

    console.log( "    ... gsl_sf_hyperg_2F1_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.0, dx: 1.0,  y: 1.0, dy: 0.5 }, 2.0, TEST_TOL0, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: 8.0,  y: 1.0, dy: 0.5 }, 12451584.0, TEST_TOL0, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: -8.0, y: 1.0, dy: 0.5 }, 0.13671875, TEST_TOL0, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: -8.1, y: 1.0, dy: 0.5 }, 0.14147385378899930422, TEST_TOL4, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: -8.0, y: 1.0, dy: -0.5 }, 4945.136718750000000, TEST_TOL0, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: -8.0, y: -5.5, dy: 0.5 },  -906.6363636363636364, TEST_TOL0, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: -8.0, y: -5.5, dy: -0.5 }, 24565.363636363636364, TEST_TOL0, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: 8.0,  y: 1.0, dy: -0.5 }, -0.006476312098196747669, TEST_TOL2, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: 8.0,  y: 5.0, dy: 0.5 }, 4205.714285714285714, TEST_TOL0, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 8.0, dx: 8.0,  y: 5.0, dy: -0.5 }, 0.0028489656290296436616, TEST_TOL2, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: 1.0, dy: 0.99 }, 1.2363536673577259280e+38 , TEST_TOL2, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: -1.5, dy: 0.99 }, 3.796186436458346579e+46, TEST_TOL2, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: -1.5, dy: -0.99 }, 0.14733409946001025146, TEST_TOL1, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: -8.5, dy: 0.99 }, -1.1301780432998743440e+65, TEST_TOL2, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: -8.5, dy: -0.99 }, -8.856462606575344483, TEST_TOL1, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: -21.5, dy: 0.99 }, 2.0712920991876073253e+95, TEST_TOL3, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: -21.5, dy: -0.99 }, -74.30517015382249216, TEST_TOL2, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: -100.5, dy: 0.99 },  -3.186778061428268980e+262, TEST_TOL3, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 9.0, dx: 9.0,  y: -100.5, dy: -0.99 },  2.4454358338375677520, TEST_TOL1, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 25.0, dx:  25.0, y: 1.0, dy: -0.5 }, -2.9995530823639545027e-06, TEST_SQRT_TOL0, "gsl_sf_hyperg_2F1_e");

    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: 1.0-1.0 / 64.0 }, 3.17175539044729373926, TEST_TOL3, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: 1.0-1.0 / 128.0 }, 3.59937243502024563424, TEST_TOL2, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: 1.0-1.0 / 256.0 }, 4.03259299524392504369, TEST_TOL1, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: 1.0-1.0 / 1024.0 }, 4.90784159359675398250, TEST_TOL1, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: 1.0-1.0 / 65536.0 }, 7.552266033399683914, TEST_TOL1, "gsl_sf_hyperg_2F1_e");
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: 1.0-1.0 / 16777216.0 }, 11.08235454026043830363, TEST_TOL1, "gsl_sf_hyperg_2F1_e");

    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: -1.0+1.0/1024.0 }, 0.762910940909954974527, TEST_TOL0, "gsl_sf_hyperg_2F1_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: -1.0+1.0/65536.0 }, 0.762762124908845424449, TEST_TOL0, "gsl_sf_hyperg_2F1_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 2.0, dy: -1.0+1.0/1048576.0 }, 0.762759911089064738044, TEST_TOL0, "gsl_sf_hyperg_2F1_e" );

    // added special handling with x == 1.0 , Richard J. Mathar, 2008-01-09

    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: 0.5, y: 3.0, dy: 1.0 }, 1.6976527263135502482014268 , TEST_TOL2, "gsl_sf_hyperg_2F1_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 1.5, dx: -4.2,y:  3.0, dy: 1.0 }, 0.15583601560025710649555254 , TEST_TOL2, "gsl_sf_hyperg_2F1_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: -7.4,dx:  0.7,y:  -1.5, dy: 1.0 }, -0.34478866959246584996859 , TEST_TOL2, "gsl_sf_hyperg_2F1_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_e, { x: 0.1, dx: -2.7,y:  -1.5, dy: 1.0 }, 1.059766766063610122925 , TEST_TOL2, "gsl_sf_hyperg_2F1_e" );

    // 2F1 conj

    console.log( "    ... gsl_sf_hyperg_2F1_conj_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_e, { x:  1.0, dx:  1.0, y: 1.0, dy:  0.5 }, 3.352857095662929028, TEST_TOL0, "gsl_sf_hyperg_2F1_conj_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_e, { x:  8.0, dx:  8.0, y: 1.0, dy:  0.5 }, 1.7078067538891293983e+09, TEST_TOL0, "gsl_sf_hyperg_2F1_conj_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_e, { x:  8.0, dx:  8.0, y: 5.0, dy:  0.5 }, 285767.15696901140627, TEST_TOL1, "gsl_sf_hyperg_2F1_conj_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_e, { x:  8.0, dx:  8.0, y: 1.0, dy: -0.5 }, 0.007248196261471276276, TEST_TOL3, "gsl_sf_hyperg_2F1_conj_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_e, { x:  8.0, dx:  8.0, y: 5.0, dy: -0.5 }, 0.00023301916814505902809, TEST_TOL3, "gsl_sf_hyperg_2F1_conj_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_e, { x: 25.0, dx: 25.0, y: 1.0, dy: -0.5 }, 5.1696944096e-06, TEST_SQRT_TOL0, "gsl_sf_hyperg_2F1_conj_e" );

    // updated correct values, testing enabled, Richard J. Mathar, 2008-01-09

    console.log( "    ... gsl_sf_hyperg_2F0_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_2F0_e, { x: 0.01, y: 1.0,  z: -0.02 }, 0.99980388665511730901180717   , TEST_TOL0, "gsl_sf_hyperg_2F0_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_2F0_e, { x: 0.1,  y: 0.5,  z: -0.02 }, 0.99901595171179281891589794   , TEST_TOL0, "gsl_sf_hyperg_2F0_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_2F0_e, { x: 1.0,  y:  1.0, z: -0.02 }, 0.98075549650574351826538049000    , TEST_TOL0, "gsl_sf_hyperg_2F0_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_2F0_e, { x: 8.0,  y:  8.0, z: -0.02 }, 0.32990592849626965538692141   , TEST_TOL0, "gsl_sf_hyperg_2F0_e" );
    TEST_SF_3D( s, gsl_sf_hyperg_2F0_e, { x: 50.0, y: 50.0, z: -0.02 }, 0.2688995263772964415245902e-12 , TEST_TOL0, "gsl_sf_hyperg_2F0_e" );

    // 2F1 renorm

    console.log( "    ... gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 1.0, dx:  1.0, y:    1.0, dy:   0.5 }, 2.0, TEST_TOL0, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 8.0, dx:  8.0, y:    1.0, dy:   0.5 }, 12451584.0, TEST_TOL0, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 8.0, dx: -8.0, y:    1.0, dy:   0.5 }, 0.13671875, TEST_TOL0, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 8.0, dx: -8.0, y:    1.0, dy:  -0.5 }, 4945.13671875, TEST_TOL0, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 8.0, dx: -8.0, y:   -5.5, dy:   0.5 }, -83081.19167659493609, TEST_TOL2, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 8.0, dx: -8.0, y:   -5.5, dy:  -0.5 }, 2.2510895952730178518e+06, TEST_TOL2, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 8.0, dx:  8.0, y:    5.0, dy:   0.5 }, 175.2380952380952381, TEST_TOL1, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 9.0, dx:  9.0, y:   -1.5, dy:  0.99 }, 1.6063266334913066551e+46, TEST_TOL2, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 9.0, dx:  9.0, y:   -1.5, dy: -0.99 }, 0.06234327316254516616, TEST_TOL2, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 5.0, dx:  5.0, y:   -1.0, dy:   0.5 }, 4949760.0, TEST_TOL1, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 5.0, dx:  5.0, y:  -10.0, dy:   0.5 }, 139408493229637632000.0, TEST_TOL2, "gsl_sf_hyperg_2F1_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_renorm_e, { x: 5.0, dx:  5.0, y: -100.0, dy:   0.5 }, 3.0200107544594411315e+206, TEST_TOL3, "gsl_sf_hyperg_2F1_renorm_e" );
  
    // 2F1 conj renorm

    console.log( "    ... gsl_sf_hyperg_2F1_conj_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_renorm_e, { x: 9.0, dx: 9.0, y:   -1.5, dy:  0.99 }, 5.912269095984229412e+49,   TEST_TOL2, "gsl_sf_hyperg_2F1_conj_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_renorm_e, { x: 9.0, dx: 9.0, y:   -1.5, dy: -0.99 }, 0.10834020229476124874,     TEST_TOL2, "gsl_sf_hyperg_2F1_conj_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_renorm_e, { x: 5.0, dx: 5.0, y:   -1.0, dy:   0.5 }, 1.4885106335357933625e+08,  TEST_TOL2, "gsl_sf_hyperg_2F1_conj_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_renorm_e, { x: 5.0, dx: 5.0, y:  -10.0, dy:   0.5 }, 7.968479361426355095e+21,   TEST_TOL2, "gsl_sf_hyperg_2F1_conj_renorm_e" );
    TEST_SF_4D( s, gsl_sf_hyperg_2F1_conj_renorm_e, { x: 5.0, dx: 5.0, y: -100.0, dy:   0.5 }, 3.1113180227052313057e+208, TEST_TOL3, "gsl_sf_hyperg_2F1_conj_renorm_e" );

    return s;

} // test_hyperg


// ****************************************************************************

export function test_coulomb( )
{

    var status  = { Integer: 0 };
    var s       = { Integer: 0 };
    var s0      = 0;
    var k_G     = 0;
    var Fe      = { Double: 0.0 };
    var Ge      = { Double: 0.0 };
    var lam_min = 0.0;
    var lam_F   = 0.0;
    var eta     = 0.0;
    var x       = 0.0;
    var F  = { val: 0.0, err: 0.0 }; // Result;
    var Fp = { val: 0.0, err: 0.0 }; // Result;
    var G  = { val: 0.0, err: 0.0 }; // Result;
    var Gp = { val: 0.0, err: 0.0 }; // Result;

    const WKB_TOL = 1.0e+04 * TEST_SQRT_TOL0;

    console.log( "Test Coulomb Functions ..." );

    console.log( "    ... gsl_sf_hydrogenicR_1_e" );
    TEST_SF_DD( s, gsl_sf_hydrogenicR_1_e, { x: 3.0, y:  2.0 },  0.025759948256148471036, TEST_TOL0, "gsl_sf_hydrogenicR_1_e" );
    TEST_SF_DD( s, gsl_sf_hydrogenicR_1_e, { x: 3.0, y: 10.0 }, 9.724727052062819704e-13, TEST_TOL1, "gsl_sf_hydrogenicR_1_e" );
    status.Integer = status.Integer + s.Integer;

    console.log( "    ... gsl_sf_hydrogenicR_e" );
    TEST_SF_IIDD( s, gsl_sf_hydrogenicR_e, { i1: 4, i2: 1, x: 3.0, y: 0.0 },  0.0,  TEST_TOL0, "gsl_sf_hydrogenicR_e" );
    TEST_SF_IIDD( s, gsl_sf_hydrogenicR_e, { i1: 4, i2: 0, x: 3.0, y: 2.0 }, -0.03623182256981820062,  TEST_TOL2, "gsl_sf_hydrogenicR_e" );
    TEST_SF_IIDD( s, gsl_sf_hydrogenicR_e, { i1: 4, i2: 1, x: 3.0, y: 2.0 }, -0.028065049083129581005, TEST_TOL2, "gsl_sf_hydrogenicR_e" );
    TEST_SF_IIDD( s, gsl_sf_hydrogenicR_e, { i1: 4, i2: 2, x: 3.0, y: 2.0 },  0.14583027278668431009,  TEST_TOL0, "gsl_sf_hydrogenicR_e" );
    status.Integer = status.Integer + s.Integer;

    TEST_SF_IIDD( s, gsl_sf_hydrogenicR_e, { i1: 100, i2:  0, x: 3.0, y: 2.0 }, -0.00007938950980052281367, TEST_TOL3, "gsl_sf_hydrogenicR_e" );
    TEST_SF_IIDD( s, gsl_sf_hydrogenicR_e, { i1: 100, i2: 10, x: 3.0, y: 2.0 },  7.112823375353605977e-12,  TEST_TOL2, "gsl_sf_hydrogenicR_e" );
    TEST_SF_IIDD( s, gsl_sf_hydrogenicR_e, { i1: 100, i2: 90, x: 3.0, y: 2.0 },  5.845231751418131548e-245, TEST_TOL2, "gsl_sf_hydrogenicR_e" );
    status.Integer = status.Integer + s.Integer;

    console.log( "    ... gsl_sf_coulomb_wave_FG_e" );
    lam_F = 0.0;
    k_G   = 0;
    eta   = 1.0;
    x     = 5.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  0.6849374120059439677, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp, -0.7236423862556063963, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G, -0.8984143590920205487, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -0.5108047585190350106, TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(1.0, 5.0, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 10.0;
    k_G   = 2;
    eta   = 1.0;
    x     = 5.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  0.0006423773354915823698, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  0.0013299570958719702545, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  33.27615734455096130,     TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -45.49180102261540580,     TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s,"  gsl_sf_coulomb_wave_FG_e(1.0, 5.0, lam_F=10, lam_G=8)" );
    status.Integer = status.Integer + s0;

    lam_F = 4.0;
    k_G   = 2;
    eta   = 50.0;
    x     = 120.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  0.0735194711823798495, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  0.6368149124126783325, TEST_TOL3 );
    s0 = s0 + s.Integer;
    //
    //s += test_sf_check_result(message_buff,  G,  , TEST_TOL5);
    //s += test_sf_check_result(message_buff, Gp, , TEST_TOL5);
    //
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(50.0, 120.0, lam_F=4, lam_G=2)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = -1000.0;
    x     = 1.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  9.68222518991341e-02, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  5.12063396274631e+00, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  1.13936784379472e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -4.30243486522438e+00, TEST_TOL3) ;
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(-1000.0, 1.0, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_min = 0.0;
    eta     = -50.0;
    x       = 5.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  1.52236975714236e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  2.03091041166137e+00, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  4.41680690236251e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -6.76485374766869e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(-50.0, 5.0, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_min = 0.0;
    eta     = -50.0;
    x       = 1000.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F, -0.2267212182760888523, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp, -0.9961306810018401525, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G, -0.9497684438900352186, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp,  0.2377656295411961399, TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(-50.0, 1000.0, lam_F=0, lam_G=0)");
    status.Integer = status.Integer + s0;

    lam_F = 10.0;
    k_G   = 0;
    eta   = -50.0;
    x     = 5.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F, -3.681143602184922e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  1.338467510317215e+00, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  3.315883246109351e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp,  1.510888628136180e+00, TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(-50.0, 5.0, lam_F=10, lam_G=10)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = -4.0;
    x     = 5.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  4.078627230056172e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  1.098212336357310e+00, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  6.743270353832442e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -6.361104272804447e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(-4.0, 5.0, lam_F=0, lam_G=0" );
    status.Integer = status.Integer + s0;

    lam_F = 3.0;
    k_G   = 0;
    eta   = -4.0;
    x     = 5.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F, -2.568630935581323e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  1.143229422014827e+00, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  7.879899223927996e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp,  3.859905878106713e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(-4.0, 5.0, lam_F=3, lam_G=3" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = 1.0;
    x     = 2.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  6.61781613832681e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  4.81557455709949e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  1.27577878476828e+00, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -5.82728813097184e-01, TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(1.0, 2.0, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = 1.0;
    x     = 0.5;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  0.08315404535022023302, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  0.22693874616222787568, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  3.1060069279548875140,  TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -3.549156038719924236,   TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(1.0, 0.5, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.5;
    k_G   = 0;
    eta   = 1.0;
    x     = 0.5;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  0.04049078073829290935, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  0.14194939168094778795, TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  4.720553853049677897,   TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -8.148033852319180005,   TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(1.0, 0.5, lam_F=0.5, lam_G=0.5)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.1;
    k_G   = 0;
    eta   = 1.0;
    x     = 0.5;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  0.07365466672379703418, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  0.21147121807178518647, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  3.306705446241024890, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -4.082931670935696644, TEST_TOL5 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(1.0, 0.5, lam_F=0.1, lam_G=0.1)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = 8.0;
    x     = 1.05;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  9.882706082810274357e-09, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  4.005167028235547770e-08, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  1.333127992006686320e+07, TEST_SQRT_TOL0 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -4.715914530842402330e+07, TEST_SQRT_TOL0 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(8.0, 1.05, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.1;
    k_G   = 0;
    eta   = 8.0;
    x     = 1.05;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  9.611416736061987761e-09, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  3.909628126126824140e-08, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  1.365928464219262581e+07, 4.0*TEST_SQRT_TOL0 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -4.848117385783386850e+07, 4.0*TEST_SQRT_TOL0 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(8.0, 1.05, lam_F=0.1, lam_G=0.1)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = 50.0;
    x     = 0.1;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  2.807788027954216071e-67, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  9.677600748751576606e-66, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  5.579810686998358766e+64, TEST_SQRT_TOL0 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -1.638329512756321424e+66, TEST_SQRT_TOL0 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(50.0, 0.1, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = 10.0;
    x     = 5.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  1.7207454091787930614e-06, 10.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  3.0975994706405458046e-06, 10.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  167637.56609459967623, 10.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -279370.76655361803075, 10.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(10.0, 5.0, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = 25.0;
    x     = 10.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  1.5451274501076114315e-16, 5.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  3.1390869393378630928e-16, 5.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  1.6177129008336318136e+15, 5.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -3.1854062013149740860e+15, 5.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(25.0, 10.0, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    k_G   = 0;
    eta   = 1.0;
    x     = 9.2;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F, -0.25632012319757955655, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  0.91518792286724220370, TEST_TOL5 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  1.03120585918973466110, TEST_SQRT_TOL0 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp,  0.21946326717491250193, TEST_SQRT_TOL0 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(1.0, 9.2, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    eta   = 10.0;
    x     = 10.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  0.0016262711250135878249, WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  0.0017060476320792806014, WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  307.87321661090837987, WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -291.92772380826822871, WKB_TOL );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(10.0, 10.0, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    lam_F = 0.0;
    eta   = 100.0;
    x     = 1.0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  8.999367996930662705e-126, 10.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  1.292746745757069321e-124, 10.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  3.936654148133683610e+123, 10.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp, -5.456942268061526371e+124, 10.0 * WKB_TOL );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(100.0, 1.0, lam_F=0, lam_G=0)" );
    status.Integer = status.Integer + s0;

    // compute F_1(eta=0,x=3.25), F'_1 and G_1(eta=0,x=3.25), G'_1

    lam_F = 1.0;
    eta   = 0.0;
    x     = 3.25;
    k_G   = 0;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  Math.sin( x ) / x - Math.cos( x ), TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  -Math.sin( x ) / (x * x) + Math.cos( x ) / x + Math.sin( x ), TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  Math.cos( x ) / x + Math.sin( x ), TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp,  -Math.cos( x ) / (x * x) - Math.sin( x ) / x + Math.cos( x ), TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(3.25, 0.0, lam_F=1, lam_G=1)" );
    status.Integer = status.Integer + s0;

    // compute F_1(eta=0,x=3.25), F'_1 and G_0(eta=0,x=3.25), G'_0

    lam_F = 1.0;
    eta   = 0.0;
    x     = 3.25;
    k_G   = 1;
    gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_G, F, Fp, G, Gp, Fe, Ge );
    s.Integer  = 0;
    s0 = 0;
    test_sf_check_result( s,  F,  Math.sin( x ) / x - Math.cos( x ), TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Fp,  -Math.sin( x ) / (x * x) + Math.cos( x ) / x + Math.sin( x ), TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s,  G,  Math.cos( x ), TEST_TOL3 );
    s0 = s0 + s.Integer;
    test_sf_check_result( s, Gp,  -Math.sin( x ), TEST_TOL3 );
    s0 = s0 + s.Integer;
    gsl_test( s, "  gsl_sf_coulomb_wave_FG_e(3.25, 0.0, lam_F=1, lam_G=0)" );
    status.Integer = status.Integer + s0;

    return status;

 } // test_coulomb


// END SF.Test;

// ----------------------------------------------------------------------------
// EOF SF-Test.mjs
