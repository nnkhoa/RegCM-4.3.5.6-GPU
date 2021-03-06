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

module mod_slice
!
! Fill 3D spaces for calculations
!
  use mod_constants
  use mod_runparams
  use mod_atm_interface
  use mod_che_interface
  use mod_pbl_interface
  use mod_rad_interface
!
  private
!
  public :: mkslice
!
  contains 
!
  subroutine mkslice
    implicit none
!
    real(rk8) :: cell , pres , satvp
    integer(ik4) :: i , j , k , n
!
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%tb3d(j,i,k) = atm2%t(j,i,k)/sfs%psb(j,i)
        end do
      end do
    end do
    do n = 1 , nqx
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atms%qxb3d(j,i,k,n) = atm2%qx(j,i,k,n)/sfs%psb(j,i)
            if ( atms%qxb3d(j,i,k,n) < minqx ) atms%qxb3d(j,i,k,n) = minqx
          end do
        end do
      end do
    end do
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              atms%chib3d(j,i,k,n) = chib(j,i,k,n)/sfs%psb(j,i)
            end do
          end do
        end do
      end do
    end if
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%ubx3d(j,i,k) = d_rfour* & 
              (atm2%u(j,i,k)+atm2%u(j,i+1,k)+ &
               atm2%u(j+1,i,k)+atm2%u(j+1,i+1,k))/sfs%psb(j,i)
          atms%vbx3d(j,i,k) = d_rfour* &
              (atm2%v(j,i,k)+atm2%v(j,i+1,k)+ &
               atm2%v(j+1,i,k)+atm2%v(j+1,i+1,k))/sfs%psb(j,i)
        end do
      end do
    end do
 
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          atms%ubd3d(j,i,k) = atm2%u(j,i,k)/psdot(j,i)
          atms%vbd3d(j,i,k) = atm2%v(j,i,k)/psdot(j,i)
        end do
      end do
    end do

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%pb3d(j,i,k) = hsigma(k)*sfs%psb(j,i)+ptop
          atms%thx3d(j,i,k) = atms%tb3d(j,i,k) * &
                  ((sfs%psb(j,i)+ptop)/atms%pb3d(j,i,k))**rovcp
                   
        end do
      end do
    end do
    !
    ! Find troposphere hgt.
    !
    iloop: &
    do i = ici1 , ici2
      jloop: &
      do j = jci1 , jci2
        kloop: &
        do k = kz , 1 , -1
          if ( atms%pb3d(j,i,k)*d_10 < ptrop(j,i)*d_r100 ) then
            ktrop(j,i) = kz-k-1
            cycle jloop
          end if
        end do kloop
      end do jloop
    end do iloop
    !
    ! compute the height at full (za) and half (zq) sigma levels:
    !
    do i = ice1 , ice2
      do j = jce1 , jce2
        atms%zq(j,i,kzp1) = d_zero
      end do
    end do
    do k = kz , 1 , -1
      do i = ice1 , ice2
        do j = jce1 , jce2
          cell = ptop/sfs%psb(j,i)
          atms%zq(j,i,k) = atms%zq(j,i,k+1) + rovg*atms%tb3d(j,i,k) *  &
                    dlog((sigma(k+1)+cell)/(sigma(k)+cell))
        end do
      end do
    end do
!
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%za(j,i,k) = d_half*(atms%zq(j,i,k)+atms%zq(j,i,k+1))
          atms%dzq(j,i,k) = atms%zq(j,i,k) - atms%zq(j,i,k+1)
        end do
      end do
    end do
   
!-----Calculate the relative humidity and air density

    do i = ice1 , ice2
      do j = jce1 , jce2
        pres = (sfs%psb(j,i)+ptop)*d_1000
        atms%rhox2d(j,i) = pres/(rgas*atms%tb3d(j,i,kz))
      end do
    end do

    atms%rhb3d(:,:,:) = d_zero

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          pres = atms%pb3d(j,i,k)*d_1000
          atms%rhob3d(j,i,k) = pres/(rgas*atms%tb3d(j,i,k)) !air density
          if ( atms%tb3d(j,i,k) > tzero ) then
            satvp = svp1*d_1000*dexp(svp2*(atms%tb3d(j,i,k)-tzero) / &
                                          (atms%tb3d(j,i,k)-svp3))
          else
            satvp = svp4*d_1000*dexp(svp5-svp6/atms%tb3d(j,i,k))
          end if
          atms%qsb3d(j,i,k) = (ep2*satvp)/(pres-satvp)
          if ( atms%qsb3d(j,i,k) > d_zero ) then
            atms%rhb3d(j,i,k) = atms%qxb3d(j,i,k,iqv)/atms%qsb3d(j,i,k)
          end if
        end do
      end do
    end do
   
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      do k = 1 , kzp1
        do i = ice1 , ice2
          do j = jce1 , jce2
            atms%tkeb3d(j,i,k) = atm2%tke(j,i,k)
          end do
        end do
      end do
    end if

  end subroutine mkslice
!
end module mod_slice
