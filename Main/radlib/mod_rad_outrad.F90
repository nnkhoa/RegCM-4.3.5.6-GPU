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

module mod_rad_outrad

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_rad_common
  use mod_runparams
  use mod_outvars

  private

  public :: allocate_mod_rad_outrad , radout

  integer(ik4) :: npr

  contains

  subroutine allocate_mod_rad_outrad
    implicit none
    npr = (jci2-jci1+1)*(ici2-ici1+1)
  end subroutine allocate_mod_rad_outrad
!
  subroutine radout(lout,solin,sabtp,frsa,clrst,clrss,qrs,firtp,         &
                    frla,clrlt,clrls,qrl,slwd,sols,soll,solsd,solld,     &
                    totcf,totcl,totci,cld,clwp,abv,sol,aeradfo,aeradfos, &
                    aerlwfo,aerlwfos,tauxar3d,tauasc3d,gtota3d,deltaz)
!
! copy radiation output quantities to model buffer
!
! change units of the radiative fluxes from cgs to mks
!
! compute the total radiative heat flux at the surface for
! the surface temperature computation
!
    implicit none
!
!     input/output arguments
!
! solin  - instantaneous incident solar
! sabtp  - total column absorbed solar flux
! frsa   - surface absorbed solar flux
! clrst  - clear sky total column abs solar flux
! clrss  - clear sky surface absorbed solar flux
! qrs    - solar heating rate
! firtp  - net up flux top of model (up-dwn flx)
! frla   - longwave cooling of surface (up-dwn flx)
! clrlt  - clr sky net up flx top of model (up-dwn f
! clrls  - clr sky lw cooling of srf (up-dwn flx)
! qrl    - longwave cooling rate
! slwd   - surface longwave down flux
! cld    - cloud fractional cover
! clwp   - cloud liquid water path
! soll   - Downward solar rad onto surface (lw direct)
! solld  - Downward solar rad onto surface (lw diffuse)
! sols   - Downward solar rad onto surface (sw direct)
! solsd  - Downward solar rad onto surface (sw diffuse)
!
    logical , intent(in) :: lout ! Preapre data for outfile
    real(rk8) , pointer , dimension(:) :: clrls , clrlt ,  &
                clrss , clrst , firtp , frla , frsa ,      &
                sabtp , slwd , solin , soll , solld ,      &
                sols , solsd , totcf , totcl , totci , abv , sol
    real(rk8) , pointer , dimension(:,:) :: cld , clwp , qrl , qrs , deltaz
    real(rk8) , pointer , dimension(:,:,:) :: tauxar3d , tauasc3d , gtota3d
    real(rk8) , pointer , dimension(:) :: aeradfo , aeradfos
    real(rk8) , pointer , dimension(:) :: aerlwfo , aerlwfos
    intent (in) cld , clrls , clrlt , clrss , clrst ,             &
                clwp , firtp , frla , frsa , qrl , qrs , sabtp ,  &
                slwd , solin , soll , solld , sols , solsd ,      &
                totcf , totcl , totci , aeradfo , aeradfos,       &
                aerlwfo , aerlwfos , deltaz
!
    integer(ik4) :: i , j , k , n
    !
    ! total heating rate in deg/s
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          heatrt(j,i,k) = qrs(n,k) + qrl(n,k)
          n = n + 1
        end do
      end do
    end do
!
!   surface absorbed solar flux in watts/m2
!   net up longwave flux at the surface
!
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        srfabswflx(j,i) = frsa(n)
        abveg(j,i) = abv(n)
        solar(j,i) = sol(n)
        srflwflxup(j,i) = frla(n)
        srflwflxdw(j,i) = slwd(n)  ! BATS Output
        n = n + 1
      end do
    end do
!
!   for coupling with bats
!
!   for now assume abveg (solar absorbed by vegetation) is equal
!   to frsa (solar absorbed by surface). possible problems are
!   over sparsely vegetated areas in which vegetation and ground
!   albedo are significantly different
!
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        totsol(j,i) = soll(n) + sols(n) + solsd(n) + solld(n)
        soldir(j,i) = sols(n)
        soldif(j,i) = solsd(n)
#ifdef CLM
        solswdir(j,i) = sols(n)
        sollwdir(j,i) = soll(n)
        solswdif(j,i) = solsd(n)
        sollwdif(j,i) = solld(n)
#endif
        n = n + 1
      end do
    end do

    if ( ktau == 0 ) return

    if ( ifchem .and. iaerosol == 1 ) then
      call copy4d_div(tauxar3d,opt_aext8_out,8,deltaz)
      call copy4d_div(tauasc3d,opt_assa8_out,8,deltaz)
      call copy4d_div(gtota3d,opt_agfu8_out,8,deltaz)
      call copy2d_integrate_from3(tauxar3d,opt_aod_out,8)
      call copy2d_add(aeradfo,opt_acstoarf_out)
      call copy2d_add(aeradfos,opt_acstsrrf_out)
      call copy2d_add(aerlwfo,opt_acstalrf_out)
      call copy2d_add(aerlwfos,opt_acssrlrf_out)
    end if
!
    if ( ifrad ) then
      if ( lout ) then
        call copy3d(cld,rad_cld_out)
        call copy3d(clwp,rad_clwp_out)
        call copy3d(qrs,rad_qrs_out)
        call copy3d(qrl,rad_qrl_out)

        call copy2d(frsa,rad_frsa_out)
        call copy2d(frla,rad_frla_out)
        call copy2d(clrst,rad_clrst_out)
        call copy2d(clrss,rad_clrss_out)
        call copy2d(clrlt,rad_clrlt_out)
        call copy2d(clrls,rad_clrls_out)
        call copy2d(solin,rad_solin_out)
        call copy2d(sabtp,rad_sabtp_out)
        call copy2d(totcf,rad_totcf_out)
        call copy2d(totcl,rad_totcl_out)
        call copy2d(totci,rad_totci_out)
        call copy2d(firtp,rad_firtp_out)
      end if
    end if

  end subroutine radout

  subroutine copy2d(a,b)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:) :: b
    integer(ik4) :: i , j , n
    if ( associated(b) ) then
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          b(j,i) = a(n)
          n = n + 1
        end do
      end do
    end if
  end subroutine copy2d

  subroutine copy2d_integrate_from3(a,b,l)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      b(:,:) = d_zero
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i) = b(j,i) + a(n,k,l)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy2d_integrate_from3

  subroutine copy3d(a,b)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:,:) :: b
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy3d

  subroutine copy4d(a,b,l)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k,l)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d

  subroutine copy4d_mult(a,b,l,c)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(in) , dimension(:,:) :: c
    real(rk8) , pointer , intent(out) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k,l) * c(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d_mult

  subroutine copy4d_div(a,b,l,c)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(in) , dimension(:,:) :: c
    real(rk8) , pointer , intent(out) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k,l) / c(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d_div

  subroutine copy2d_add(a,b)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:) :: b
    integer(ik4) :: i , j , n
    if ( associated(b) ) then
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          b(j,i) = b(j,i) + a(n)
          n = n + 1
        end do
      end do
    end if
  end subroutine copy2d_add

  subroutine copy3d_add(a,b)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:,:) :: b
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = b(j,i,k) + a(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy3d_add

  subroutine copy4d_add(a,b,l)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(out) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = b(j,i,k) + a(n,k,l)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d_add

end module mod_rad_outrad
