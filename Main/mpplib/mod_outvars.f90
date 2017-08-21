!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
module mod_outvars

  use mod_realkinds

  public

  real(rk8) , dimension(:,:) , pointer :: xlon_out => null()
  real(rk8) , dimension(:,:) , pointer :: xlat_out => null()
  real(rk8) , dimension(:,:) , pointer :: topo_out => null()
  real(rk8) , dimension(:,:) , pointer :: mask_out => null()
  real(rk8) , dimension(:,:) , pointer :: ps_out => null()

  real(rk8) , dimension(:,:) , pointer :: sub_xlon_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_xlat_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_topo_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_mask_out => null()

  real(rk8) , dimension(:,:) , pointer :: sub_ps_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: atm_u_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_v_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_t_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_omega_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_qv_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_qc_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_tke_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_kth_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_kzm_out => null()
  real(rk8) , dimension(:,:) , pointer :: atm_tgb_out => null()
  real(rk8) , dimension(:,:) , pointer :: atm_tpr_out => null()
  real(rk8) , dimension(:,:) , pointer :: atm_tsw_out => null()

  real(rk8) , dimension(:,:) , pointer :: srf_uvdrag_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_tg_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_tlef_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_evp_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_scv_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_sena_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_flw_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_fsw_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_fld_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_sina_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_tpr_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_prcv_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_zpbl_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_aldirs_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_aldifs_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_sund_out => null()
  real(rk8) , dimension(:,:) , pointer :: srf_seaice_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: srf_u10m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: srf_v10m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: srf_t2m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: srf_q2m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: srf_smw_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: srf_runoff_out => null()

  real(rk8) , dimension(:,:) , pointer :: sts_tgmax_out => null()
  real(rk8) , dimension(:,:) , pointer :: sts_tgmin_out => null()
  real(rk8) , dimension(:,:) , pointer :: sts_pcpmax_out => null()
  real(rk8) , dimension(:,:) , pointer :: sts_pcpavg_out => null()
  real(rk8) , dimension(:,:) , pointer :: sts_sund_out => null()
  real(rk8) , dimension(:,:) , pointer :: sts_psmin_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: sts_t2max_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: sts_t2min_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: sts_t2avg_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: sts_w10max_out => null()

  real(rk8) , dimension(:,:) , pointer :: sub_uvdrag_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_tg_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_tlef_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_evp_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_scv_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_sena_out => null()
  real(rk8) , dimension(:,:) , pointer :: sub_tlake_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: sub_u10m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: sub_v10m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: sub_t2m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: sub_q2m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: sub_smw_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: sub_runoff_out => null()

  real(rk8) , dimension(:,:) , pointer :: rad_frsa_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_frla_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_clrst_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_clrss_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_clrls_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_clrlt_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_solin_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_sabtp_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_totcf_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_totcl_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_totci_out => null()
  real(rk8) , dimension(:,:) , pointer :: rad_firtp_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: rad_cld_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: rad_clwp_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: rad_qrs_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: rad_qrl_out => null()

  real(rk8) , dimension(:,:) , pointer :: lak_tg_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_tpr_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_scv_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_sena_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_sina_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_fsw_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_flw_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_fld_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_evp_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_aldirs_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_aldifs_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_aveice_out => null()
  real(rk8) , dimension(:,:) , pointer :: lak_hsnow_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: lak_tlake_out => null()

  real(rk8) , dimension(:,:) , pointer :: opt_acstoarf_out => null()
  real(rk8) , dimension(:,:) , pointer :: opt_acstsrrf_out => null()
  real(rk8) , dimension(:,:) , pointer :: opt_acstalrf_out => null()
  real(rk8) , dimension(:,:) , pointer :: opt_acssrlrf_out => null()
  real(rk8) , dimension(:,:) , pointer :: opt_aod_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: opt_aext8_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: opt_assa8_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: opt_agfu8_out => null()

  real(rk8) , dimension(:,:) , pointer :: che_wdrflx_out => null()
  real(rk8) , dimension(:,:) , pointer :: che_wdcflx_out => null()
  real(rk8) , dimension(:,:) , pointer :: che_ddflx_out => null()
  real(rk8) , dimension(:,:) , pointer :: che_emflx_out => null()
  real(rk8) , dimension(:,:) , pointer :: che_ddvel_out => null()
  real(rk8) , dimension(:,:) , pointer :: che_burden_out => null()
  real(rk8) , dimension(:,:) , pointer :: che_pblten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_mixrat_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_cheten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_advhten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_advvten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_difhten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_cuten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_tuten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_raiten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_wasten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_bdyten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_sedten_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: che_emten_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: slab_qflx_out => null()

end module mod_outvars
