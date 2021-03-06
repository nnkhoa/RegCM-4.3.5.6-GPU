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

module mod_pbl_interface

  use mod_realkinds
  use mod_service
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_mppparam
  use mod_atm_interface
  use mod_pbl_common
  use mod_pbl_holtbl
  use mod_pbl_uwtcm
  use mod_runparams , only : iqc , iqv , dt , rdt , ichem , hsigma , dsigma

  public

  contains

  subroutine init_pbl(atm2,atms,aten,holtten,uwten,adf,heatrt,chiten, &
                      remdrd,cchifxuw,psdot,sfs,mddom,ldmsk,chtrdpv)
    implicit none
    type (atmstate) , intent(in) :: atm2 , aten , holtten , uwten
    type (slice) , intent(in) :: atms
    type (diffx) , intent(in) :: adf
    type (domain) , intent(in) :: mddom
    type (surfstate) , intent(in) :: sfs
    real(rk8) , pointer , dimension(:,:,:) :: heatrt
    real(rk8) , pointer , dimension(:,:,:,:) :: chiten
    real(rk8) , pointer , dimension(:,:,:) :: cchifxuw
    real(rk8) , pointer , dimension(:,:,:) :: remdrd
    integer(ik4) , pointer , dimension(:,:) :: ldmsk
    real(rk8) , pointer , dimension(:,:) :: psdot
    real(rk8) , pointer , dimension(:,:) :: chtrdpv

    tkemin = 1.0D-8

    call assignpnt(aten%u,uten)
    call assignpnt(aten%v,vten)
    call assignpnt(aten%t,tten)
    call assignpnt(aten%tke,tketen)
    call assignpnt(aten%qx,qxten)
    call assignpnt(uwten%u,uuwten)
    call assignpnt(uwten%v,vuwten)
    call assignpnt(uwten%t,tuwten)
    call assignpnt(uwten%tke,tkeuwten)
    call assignpnt(uwten%qx,qxuwten)
    call assignpnt(atm2%tke,tkests)
    call assignpnt(atms%ubx3d,uxatm)
    call assignpnt(atms%vbx3d,vxatm)
    call assignpnt(atms%ubd3d,udatm)
    call assignpnt(atms%vbd3d,vdatm)
    call assignpnt(atms%tb3d,tatm)
    call assignpnt(atms%qxb3d,qxatm)
    call assignpnt(atms%chib3d,chmx)
    call assignpnt(atms%thx3d,thxatm)
    call assignpnt(atms%za,za)
    call assignpnt(atms%zq,zq)
    call assignpnt(atms%dzq,dzq)
    call assignpnt(atms%rhox2d,rhox2d)
    call assignpnt(adf%difft,difft)
    call assignpnt(adf%diffqx,diffqx)
    call assignpnt(heatrt,radheatrt)
    call assignpnt(holtten%qx,diagqx)
    call assignpnt(chiten,chten)
    call assignpnt(remdrd,drmr)
    call assignpnt(sfs%psb,sfcps)
    call assignpnt(psdot,sfcpd)
    call assignpnt(sfs%tgb,tg)
    call assignpnt(sfs%qfx,qfx)
    call assignpnt(sfs%hfx,hfx)
    call assignpnt(sfs%uvdrag,uvdrag)
    call assignpnt(mddom%coriol,coriolis)
    call assignpnt(mddom%msfx,mapfcx)
    call assignpnt(ldmsk,landmsk)
    call assignpnt(chtrdpv,depvel)
    call assignpnt(cchifxuw,chifxuw)
  end subroutine init_pbl

  subroutine get_data_from_tcm(tcmstate,tcmtend,aten,atm1,atm2,bRegridWinds)
    use mod_runparams , only : ibltyp
    implicit none
    type(atmstate) , intent(inout) :: tcmtend , aten , atm1 , atm2
    type(tcm_state) , intent(inout) :: tcmstate
    logical , intent(in) :: bRegridWinds

    ! Don't update the model variables if we are the diagnostic mode
    ! (Holtslag running, and UW updating tke)
    if ( ibltyp /= 99 ) then
      !
      ! Put the t and qv tendencies in to difft and diffq for
      ! application of the sponge boundary conditions (see mod_tendency)
      !
      diffqx(jci1:jci2,ici1:ici2,:,iqv) =  &
                  diffqx(jci1:jci2,ici1:ici2,:,iqv) +  &
                  tcmtend%qx(jci1:jci2,ici1:ici2,:,iqv)
      difft(jci1:jci2,ici1:ici2,:) =  &
                  difft(jci1:jci2,ici1:ici2,:) +  &
                  tcmtend%t(jci1:jci2,ici1:ici2,:)

      ! Put the cloud water tendency in aten
      aten%qx(jci1:jci2,ici1:ici2,:,iqc) =  &
                 aten%qx(jci1:jci2,ici1:ici2,:,iqc) +   &
                 tcmtend%qx(jci1:jci2,ici1:ici2,:,iqc)

      ! Put the tracer tendencies in chiuwten
      ! TODO: may want to calcuate rmdr here following holtbl
      if ( ichem == 1 ) then
        chten(jci1:jci2,ici1:ici2,:,:) = &
              chten(jci1:jci2,ici1:ici2,:,:)+chiuwten(jci1:jci2,ici1:ici2,:,:)
      end if

      if ( .not. bRegridWinds ) then
        !
        ! If the TCM calculations were done on the dot grid, then
        ! the u and v tendencies need not to be regridded to the dot grid
        !
        aten%u(jci1:jci2,ici1:ici2,:) = &
                 aten%u(jci1:jci2,ici1:ici2,:) +    &
                 tcmtend%u(jci1:jci2,ici1:ici2,:)
        aten%v(jci1:jci2,ici1:ici2,:) = &
                 aten%v(jci1:jci2,ici1:ici2,:) +    &
                 tcmtend%v(jci1:jci2,ici1:ici2,:)
      end if

      zpbl(jci1:jci2,ici1:ici2) = &
                tcmstate%zpbl(jci1:jci2,ici1:ici2)
    end if

!   !
!   ! Interpolate kzm and kth from the interfaces to the midpoints
!   ! if the diffusivities are to be used in the holtslag model
!   !
!   if ( ibltyp == 99 ) then
!     do k = 1 , kz
!       tcmstate%kzm(:,:,k) = sqrt(tcmstate%kzm(:,:,k)*tcmstate%kzm(:,:,k+1))
!       tcmstate%kth(:,:,k) = sqrt(tcmstate%kth(:,:,k)*tcmstate%kth(:,:,k+1))
!     end do
!   end if
!   !
!   ! Shift kth and kzm for output if the UW model is running
!   !
!   if ( ibltyp == 2 ) then
!     do k = 1 , kz
!       tcmstate%kzm(:,:,k) = tcmstate%kzm(:,:,k+1)
!       tcmstate%kth(:,:,k) = tcmstate%kth(:,:,k+1)
!     end do
!   end if

    aten%tke(jci1:jci2,ici1:ici2,:) = &
                tcmtend%tke(jci1:jci2,ici1:ici2,:)
    !
    ! Set the surface tke (diagnosed)
    !
    atm1%tke(jci1:jci2,ici1:ici2,kzp1) =  &
               tcmstate%srftke(jci1:jci2,ici1:ici2)
    atm2%tke(jci1:jci2,ici1:ici2,kzp1) =  &
               tcmstate%srftke(jci1:jci2,ici1:ici2)

!   tcmtend%qx = 0.0d0
!   tcmtend%t = 0.0d0
!   tcmtend%u = 0.0d0
!   tcmtend%v = 0.0d0
!   tcmtend%tke = 0.0d0


  end subroutine get_data_from_tcm

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     This subroutine computes the horizontal flux-divergence terms   c
!     for tke.  Second-order difference is used.                      c
!                                                                     c
!     dxx    : is the horizontal distance.                            c
!     j      : is the j'th slice of f anf ften.                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine hadvtke(tcmstate,atm,twt,dxx)
    implicit none
    real(rk8) , intent(in) :: dxx
    type(atmstate) , intent(in) :: atm
    type(tcm_state) , intent(inout) :: tcmstate
    real(rk8) , pointer , dimension(:,:) :: twt
    integer(ik4) :: i , k , j
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'hadvtke'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2

        !
        ! Interpoalte the winds to the full sigma levels
        ! while the advection term is calculated
        !
          tcmstate%advtke(j,i,k) = tcmstate%advtke(j,i,k)  &
           -(((atm%u(j+1,i+1,k-1)+atm%u(j+1,i,k-1))*twt(k,2)  &
             +(atm%u(j+1,i+1,k)  +atm%u(j+1,i,k))*twt(k,1)) &
             *( atm%tke(j,i,k)+atm%tke(j+1,i,k))  &
            -((atm%u(j,i+1,k-1)+atm%u(j,i,k-1))*twt(k,2)  &
             +(atm%u(j,i+1,k)  +atm%u(j,i,k))*twt(k,1)) &
             *( atm%tke(j,i,k)+atm%tke(j-1,i,k))  &
            +((atm%v(j,i+1,k-1)+atm%v(j+1,i+1,k-1))*twt(k,2)  &
             +(atm%v(j,i+1,k)  +atm%v(j+1,i+1,k))*twt(k,1)) &
             *( atm%tke(j,i,k)+atm%tke(j,i+1,k))  &
            -((atm%v(j,i,k-1)+atm%v(j+1,i,k-1))*twt(k,2)  &
             +(atm%v(j,i,k)  +atm%v(j+1,i,k))*twt(k,1)) &
             *( atm%tke(j,i,k)+atm%tke(j,i-1,k))) &
             /(dxx*mapfcx(j,i)*mapfcx(j,i))
!TAO Debug:
!         tcmstate%advtke(j,i,k)= tcmstate%advtke(j,i,k) + 1d-8*atm%tke(j+1,i,k)
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine hadvtke

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the vertical flux-divergence terms.    c
!                                                                     c
!     j      : jth slice of variable fa.                              c
!                                                                     c
!     ind = 1 : Bretherton's vertical advection method                c
!           2 : Alternate vertical advection method (unknown origin)  c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine vadvtke(tcmstate,qdot,ind)
    implicit none
    integer(ik4) , intent(in) :: ind
    type(tcm_state) :: tcmstate
    real(rk8) , dimension(:,:,:) , pointer , intent(in) :: qdot
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'vadvtke'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
        
    !
    ! Use Bretherton's method for tke advection
    !
    if ( ind == 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            dotqdot(j,i,k) = (qdot(j,i,k)+qdot(j,i,k+1))
            ftmp(j,i,k) = 0.5D0*(tcmstate%tkeps(j,i,k)+tcmstate%tkeps(j,i,k+1))
          end do
        end do
      end do

      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tcmstate%advtke(j,i,k) = tcmstate%advtke(j,i,k) - &
                         (dotqdot(j,i,k)*ftmp(j,i,k) -  &
                         dotqdot(j,i,k-1)*ftmp(j,i,k-1))/ &
                         (dsigma(k)+dsigma(k-1))
          end do
        end do
      end do
    !
    ! Use an alternative method (this came from where?)
    !
    else
      do i = ici1 , ici2
        do j = jci1 , jci2
          tcmstate%advtke(j,i,1) = tcmstate%advtke(j,i,1)-  &
                         qdot(j,i,2)*tcmstate%tkeps(j,i,2)/dsigma(1)
        end do
      end do
      do k = 2 , kzm1
        do i = ici1 , ici2
          do j = jci1 , jci2
            tcmstate%advtke(j,i,k) = tcmstate%advtke(j,i,k) &
                         -(qdot(j,i,k+1)*tcmstate%tkeps(j,i,k+1)   &
                         - qdot(j,i,k)*tcmstate%tkeps(j,i,k))/dsigma(k)
          end do
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          tcmstate%advtke(j,i,kz) = tcmstate%advtke(j,i,kz)+  &
                         qdot(j,i,kz)*tcmstate%tkeps(j,i,kz)/dsigma(kz)
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine vadvtke

  subroutine set_tke_bc(atm1,atm2)
    implicit none
    type(atmstate) , intent(inout) :: atm1 , atm2
    !
    ! Set tke boundary conditions
    !
    if ( ma%has_bdyleft ) then
      atm1%tke(jce1,:,:) = tkemin ! East boundary
      atm2%tke(jce1,:,:) = tkemin ! East boundary
    end if
    if ( ma%has_bdyright ) then
      atm1%tke(jce2,:,:) = tkemin ! West boundary
      atm2%tke(jce2,:,:) = tkemin ! West boundary
    end if
    if ( ma%has_bdytop ) then
      atm1%tke(:,ice2,:) = tkemin  ! South boundary
      atm2%tke(:,ice2,:) = tkemin  ! South boundary
    end if
    if ( ma%has_bdybottom ) then
      atm1%tke(:,ice1,:) = tkemin  ! North boundary
      atm2%tke(:,ice1,:) = tkemin  ! North boundary
    end if
    !
    ! End set tke boundary conditions
    !
  end subroutine set_tke_bc

  subroutine check_conserve_qt(rcmqxten,tcmtend,tcmstate,kmax)
    implicit none
    type(atmstate) , intent(in) :: tcmtend
    type(tcm_state) :: tcmstate
    real(rk8) , dimension(:,:,:,:) , intent(in) :: rcmqxten
    integer(ik4) , intent(in) :: kmax
    real(rk8) , dimension(kmax) :: rho1d , rhobydpdz1d
    real(rk8) :: qwtcm , qwrcm , qwanom , dtops , xps , ps2 , dza
    integer(ik4) :: i , j , k

    do k = 1 , kzm1
      do j = jci1 , jci2
        do i = ici1 , ici2
          xps = (hsigma(k)*sfcps(j,i)+ptop)
          ps2 = (hsigma(k+1)*sfcps(j,i)+ptop)
          dza = za(j,i,k) - za(j,i,k+1)
          rhobydpdz1d(k) = d_1000*(ps2-xps)/(egrav*dza)
        end do
        rhobydpdz1d(kz) = rhox2d(j,i)
        dtops = dt/sfcps(j,i)

        rho1d = d_1000*(hsigma*sfcps(j,i) + ptop) / &
                    ( rgas * tatm(j,i,:) *  &
                    (d_one + ep1* qxatm(j,i,:,iqv) - qxatm(j,i,:,iqc)) )


        qwtcm = sum((tcmtend%qx(j,i,:,iqv) + tcmtend%qx(j,i,:,iqc))  &
                      *rho1d*dzq(j,i,:))*dtops
        qwrcm = sum((rcmqxten(j,i,:,iqv) + rcmqxten(j,i,:,iqc))   &
                      *rho1d*dzq(j,i,:))*dtops
!                     *rhobydpdz1d*dzq(j,i,:))*dtops

!       qwanom = qwtcm - qwrcm
!       qwanom = qwrcm-qfx(j,i)*dt
        qwanom = qwtcm - dt*qfx(j,i)
        tcmstate%kzm(j,i,1) = qwanom
      end do
    end do
  end subroutine check_conserve_qt
!
  ! Set the net surface flux (wet/dry dep + emission) to send to the UW TCM
  !  This routine is called in mod_tendency.F90 just prior to the call of
  !  uwtcm()
  subroutine set_tracer_surface_fluxes()
  implicit none 
  integer(ik4) :: itr

    !Set the variable chifxuw to be the net flux for each tracer
    ! (dummy declaration below -- declared and allocated in mod_pbl_common.F90)
    tracerfluxloop:  &
    do itr = 1,ntr
      !chifxuw(:,:,ntr) = ...
    end do tracerfluxloop

  end subroutine set_tracer_surface_fluxes

end module mod_pbl_interface
