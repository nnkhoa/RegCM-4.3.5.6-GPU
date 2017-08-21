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

module mod_che_mppio
!
  use mod_intkinds
  use mod_realkinds
  use mod_mppparam
  use mod_runparams
  use mod_dynparam
  use mod_memutil
  use mod_mpmessage
  use mod_che_param
  use mod_che_common
  use mod_che_species
!
  public

  real(rk8) , pointer , dimension(:,:,:,:) :: remlsc_io , remcvc_io
  real(rk8) , pointer , dimension(:,:,:) :: remdrd_io
  real(rk8) , pointer , dimension(:,:) :: ssw2da_io , sdeltk2d_io ,   &
                                         sdelqk2d_io , sfracv2d_io , &
                                         sfracb2d_io , sfracs2d_io , &
                                         svegfrac2d_io
  real(rk8) , pointer , dimension(:,:,:,:) :: chia_io , chib_io
  real(rk8) , pointer , dimension(:,:,:,:) :: chemall_io
  real(rk8) , pointer , dimension(:,:,:,:) :: taucldsp_io
!
  contains 
    !
    ! This routines allocate all the arrays contained in the module
    !
    subroutine allocate_mod_che_mppio
      implicit none

      if ( ichem == 1 ) then
        if ( myid == iocpu ) then
          call getmem4d(remlsc_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:remlsc_io')
          call getmem4d(remcvc_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:remcvc_io')
          call getmem3d(remdrd_io,jcross1,jcross2,icross1,icross2, &
                        1,ntr,'che_mppio:remdrd_io')

          call getmem2d(ssw2da_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:ssw2da_io')
          call getmem2d(sdelqk2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sdelqk2d_io')
          call getmem2d(sdeltk2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sdeltk2d_io')
          call getmem2d(sfracb2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracb2d_io')
          call getmem2d(sfracs2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracs2d_io')
          call getmem2d(sfracv2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracv2d_io')
          call getmem2d(svegfrac2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:svegfrac2d_io')

          call getmem4d(chia_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:chia_io')
          call getmem4d(chib_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:chib_io')
          if ( igaschem == 1 .and. ichsolver > 0 ) then
            call getmem4d(chemall_io,jcross1,jcross2,icross1,icross2, &
                          1,kz,1,totsp,'che_mppio:chemall_io')
            call getmem4d(taucldsp_io,jcross1,jcross2,icross1,icross2, &
                          0,kz,1,nspi,'che_mppio:chemall_io')
          end if
        end if
      end if
    end subroutine allocate_mod_che_mppio
!
end module mod_che_mppio
