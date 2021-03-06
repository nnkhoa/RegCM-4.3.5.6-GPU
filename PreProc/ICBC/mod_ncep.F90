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

module mod_ncep

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil
  use mod_message
  use mod_memutil
  use mod_nchelper

  private

  integer(ik4) , parameter :: klev = 13 , jlat = 73 , ilon = 144

  real(rk8) , dimension(ilon,jlat) :: psvar
  real(rk8) , dimension(jlat) :: glat
  real(rk8) , dimension(ilon) :: glon
  real(rk8) , dimension(klev) :: sigma1 , sigmar

  real(rk8) , dimension(ilon,jlat,klev) :: wvar

  real(rk8) , target , dimension(ilon,jlat,klev*3) :: b2
  real(rk8) , target , dimension(ilon,jlat,klev*2) :: d2
  real(rk8) , pointer , dimension(:,:,:) :: b3
  real(rk8) , pointer , dimension(:,:,:) :: d3
  integer(2) , pointer , dimension(:,:,:) :: work
  
  real(rk8) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rk8) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rk8) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rk8) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

  integer(ik4) :: year , month , day , hour

  public :: getncep , headernc

  contains

  subroutine getncep(idate,itype)
  implicit none
!
  integer(ik4) , intent(in) :: itype
  type(rcm_time_and_date) , intent(in) :: idate
!
!     D      BEGIN LOOP OVER NTIMES
!
!
  call split_idate(idate,year,month,day,hour)
!
  if ( itype == 1 ) then
    call cdc6hour(idate,globidate1)
  else if ( itype == 2 ) then
    call cdc6hour2(idate,globidate1)
  end if

  write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
  call bilinx2(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
  call bilinx2(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
  call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!HH:  CHANGE THE VERTICAL ORDER.
  call top2btm(t3,jx,iy,klev)
  call top2btm(q3,jx,iy,klev)
  call top2btm(h3,jx,iy,klev)
  call top2btm(u3,jx,iy,klev)
  call top2btm(v3,jx,iy,klev)
!HH:OVER
!
!     ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
  call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,klev)
 
  call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
  if(i_band == 1) then
     call p1p2_band(b3pd,ps4,jx,iy)
  else
     call p1p2(b3pd,ps4,jx,iy)
  endif
 
!
!     F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
  call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)
 
  call readsst(ts4,idate)

!     F2  DETERMINE P* AND HEIGHT.
!
!     F3  INTERPOLATE U, V, T, AND Q.
  call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
  call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
!
  call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
 
  call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
  call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
!
!     F4  DETERMINE H
  call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
  end subroutine getncep

  subroutine cdc6hour(idate,idate0)
  use netcdf
  implicit none
!
  type(rcm_time_and_date) , intent (in) :: idate , idate0
!
  integer(ik4) :: i , ilev , inet , it , j , kkrec , m , k , nlev , istatus
  character(len=21) :: inname
  character(len=256) :: pathaddname
  character(len=5) , dimension(7) :: varname
  real(rk8) :: xadd , xscale
  integer(ik4) , dimension(10) , save :: icount , istart
  integer(ik4) , dimension(7) , save :: inet7 , ivar7
  real(rk8) , dimension(7) , save :: xoff , xscl
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers.
!
  data varname/'air' , 'hgt' , 'rhum' , 'uwnd' , 'vwnd' , 'omega' , 'pres'/
!
!     Below in the ncopen call is the file name of the netCDF file.
!     You may want to add code to read in the file name and the
!     variable name.
!     OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
  xadd = 0.0D0
  xscale = 1.0D0
!fix  do kkrec = 1,7
  nlev = 0
  do kkrec = 1 , 5
    if ( dattyp == 'NNRP1' ) then
      if ( kkrec == 1 .or. kkrec == 2 .or. kkrec == 4 .or. kkrec == 5 ) nlev = klev
      if ( kkrec == 6 ) nlev = 12
      if ( kkrec == 3 ) nlev = 8
      if ( kkrec == 7 ) nlev = 0
    else if ( dattyp == 'NNRP2' ) then
      if ( kkrec <= 6 ) nlev = klev
      if ( kkrec == 7 ) nlev = 0
    else
    end if
    if ( idate == idate0 .or. (lfdoyear(idate) .and. lmidnight(idate))) then
      if ( kkrec == 1 ) then
        write (inname,99001) year , 'air.' , year
      else if ( kkrec == 2 ) then
        write (inname,99001) year , 'hgt.' , year
      else if ( kkrec == 3 ) then
        write (inname,99002) year , 'rhum.' , year
      else if ( kkrec == 4 ) then
        write (inname,99002) year , 'uwnd.' , year
      else if ( kkrec == 5 ) then
        write (inname,99002) year , 'vwnd.' , year
      else if ( kkrec == 6 ) then
        write (inname,99003) year , 'omega.' , year
      else if ( kkrec == 7 ) then
        write (inname,99004) year , 'pres.sfc.' , year
      else
      end if
 
      if ( dattyp == 'NNRP1' ) then
        pathaddname = trim(inpglob)//'/NNRP1/'//inname
      else if ( dattyp == 'NNRP2' ) then
        pathaddname = trim(inpglob)//'/NNRP2/'//inname
      else
      end if
      istatus = nf90_open(pathaddname,nf90_nowrite,inet7(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(pathaddname))
      istatus = nf90_inq_varid(inet7(kkrec),varname(kkrec),ivar7(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__, &
           'Variable '//varname(kkrec)//' error in file'//trim(pathaddname))
      istatus = nf90_get_att(inet7(kkrec),ivar7(kkrec), &
                            'scale_factor',xscl(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    ':scale_factor in file'//trim(pathaddname))
      istatus = nf90_get_att(inet7(kkrec),ivar7(kkrec), &
                             'add_offset',xoff(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    ':add_offset in file'//trim(pathaddname))
      write (stdout,*) inet7(kkrec) , pathaddname , xscl(kkrec) , xoff(kkrec)
    end if
 
    it = (day-1)*4 + hour/6 + 1
    if ( month == 2 ) it = it + 31*4
    if ( month == 3 ) it = it + 59*4
    if ( month == 4 ) it = it + 90*4
    if ( month == 5 ) it = it + 120*4
    if ( month == 6 ) it = it + 151*4
    if ( month == 7 ) it = it + 181*4
    if ( month == 8 ) it = it + 212*4
    if ( month == 9 ) it = it + 243*4
    if ( month == 10 ) it = it + 273*4
    if ( month == 11 ) it = it + 304*4
    if ( month == 12 ) it = it + 334*4
    if ( mod(day,4) == 0 .and. month > 2 ) it = it + 4
    if ( mod(day,100) == 0 .and. month > 2 ) it = it - 4
    if ( mod(day,400) == 0 .and. month > 2 ) it = it + 4
!bxq_
    do m = 1 , 4
      istart(m) = 1
    end do
    do m = 5 , 10
      istart(m) = 0
      icount(m) = 0
    end do
    icount(1) = ilon
    icount(2) = jlat
    icount(4) = 1460
    if ( mod(day,4) == 0 ) icount(4) = 1464
    if ( mod(day,100) == 0 ) icount(4) = 1460
    if ( mod(day,400) == 0 ) icount(4) = 1464
    istart(4) = it
    icount(4) = 1
    inet = inet7(kkrec)
    if ( nlev > 0 ) then
      icount(3) = nlev
      istatus = nf90_get_var(inet,ivar7(kkrec),work,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    'read error in file'//trim(pathaddname))
      xscale = xscl(kkrec)
      xadd = xoff(kkrec)
      do ilev = 1 , nlev
        if ( kkrec == 1 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              tvar(i,jlat+1-j,14-ilev) = real(dble(work(i,j,ilev))*xscale+xadd)
            end do
          end do
        else if ( kkrec == 2 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              hvar(i,jlat+1-j,14-ilev) = real(dble(work(i,j,ilev))*xscale+xadd)
            end do
          end do
        else if ( kkrec == 3 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              rhvar(i,jlat+1-j,14-ilev) = real(dmin1((dble(work(i,j,ilev))* &
                            xscale+xadd)*0.01D0,1.D0))
            end do
          end do
        else if ( kkrec == 4 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              uvar(i,jlat+1-j,14-ilev) = real(dble(work(i,j,ilev))*xscale+xadd)
            end do
          end do
        else if ( kkrec == 5 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              vvar(i,jlat+1-j,14-ilev) = real(dble(work(i,j,ilev))*xscale+xadd)
            end do
          end do
        else if ( kkrec == 6 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              wvar(i,jlat+1-j,14-ilev) = real(dble(work(i,j,ilev))*xscale+xadd)
            end do
          end do
        else
        end if
      end do
    else if ( nlev == 0 ) then
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar7(kkrec),work,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    'read error in file'//trim(pathaddname))
      if ( kkrec == 7 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            psvar(i,jlat+1-j) = real(dble(work(i,j,1))*xscale + xadd)
          end do
        end do
      end if
    else
    end if
    if ( dattyp == 'NNRP1' ) then
!         It's a pity that we have to nudge the values by the following
!         way
      do k = 5 , 1 , -1
        do j = 1 , jlat
          do i = 1 , ilon
            rhvar(i,j,k) = rhvar(i,j,k+1)
          end do
        end do
      end do
 
      do j = 1 , jlat
        do i = 1 , ilon
          wvar(i,j,1) = 0.0
        end do
      end do
    end if
  end do


99001 format (i4,'/',a4,i4,'.nc')
99002 format (i4,'/',a5,i4,'.nc')
99003 format (i4,'/',a6,i4,'.nc')
99004 format (i4,'/',a9,i4,'.nc')
!
  end subroutine cdc6hour

  subroutine cdc6hour2(idate,idate0)
  use netcdf
  use mod_grid
  implicit none
!
  type(rcm_time_and_date) , intent(in) :: idate , idate0
!
  integer(ik4) :: i , ii , ilev , inet , it , j , jj , kkrec , m , nlev , istatus
  character(len=24) :: inname
  character(len=256) :: pathaddname
  character(len=5) , dimension(7) :: varname
  integer(ik4) :: iii , jjj
  real(rk8) :: xadd , xscale
  integer(ik4) , dimension(10) :: icount , istart
  integer(ik4) , dimension(7) :: inet7 , ivar7
  real(rk8) , dimension(7) :: xoff , xscl
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers.
!
!     DATA ARRAY AND WORK ARRAY
!
  data varname/'air' , 'hgt' , 'rhum' , 'uwnd' , 'vwnd' , 'omega' , 'pres'/
!
  xadd = 0.0D0
  xscale = 1.0D0

  if ( idate == idate0 ) then
    i0 = idnint(lon0/2.5D0) + 1
    if ( i0 <= 0 ) i0 = i0 + ilon
    if ( i0 > ilon ) i0 = i0 - ilon
    i1 = idnint(lon1/2.5D0) + 1
    if ( i1 <= 0 ) i1 = i1 + ilon
    if ( i1 > ilon ) i1 = i1 - ilon
    j0 = idnint(lat0/2.5D0) + 36
    j1 = idnint(lat1/2.5D0) + 1

    iii = i1 - i0
    jjj = j1 - j0
  end if
!
  nlev = 0
  do kkrec = 1 , 5
    if ( kkrec <= 6 ) nlev = klev
    if ( kkrec == 7 ) nlev = 0
    if ( idate == idate0 .or. (lfdoyear(idate) .and. lmidnight(idate))) then
      if ( kkrec == 1 ) then
        write (inname,99001) year , 'air.WIN.' , year
      else if ( kkrec == 2 ) then
        write (inname,99001) year , 'hgt.WIN.' , year
      else if ( kkrec == 3 ) then
        write (inname,99002) year , 'rhum.WIN.' , year
      else if ( kkrec == 4 ) then
        write (inname,99002) year , 'uwnd.WIN.' , year
      else if ( kkrec == 5 ) then
        write (inname,99002) year , 'vwnd.WIN.' , year
      else if ( kkrec == 6 ) then
        write (inname,99003) year , 'omega.WIN.' , year
      else if ( kkrec == 7 ) then
        write (inname,99004) year , 'pres.sfc.WIN.' , year
      else
      end if
 
      pathaddname = trim(inpglob)//'/NNRP2/'//inname
      istatus = nf90_open(pathaddname,nf90_nowrite,inet7(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(pathaddname))
      istatus = nf90_inq_varid(inet7(kkrec),varname(kkrec),ivar7(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    'error in file'//trim(pathaddname))
      istatus = nf90_get_att(inet7(kkrec),ivar7(kkrec), &
                             'scale_factor',xscl(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    ':scale_factor in file'//trim(pathaddname))
      istatus = nf90_get_att(inet7(kkrec),ivar7(kkrec), &
                             'add_offset',xoff(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    ':add_offset in file'//trim(pathaddname))
      write (stdout,*) inet7(kkrec) , pathaddname , xscl(kkrec) , xoff(kkrec)
    end if
 
    it = (day-1)*4 + hour/6 + 1
    if ( month == 2 ) it = it + 31*4
    if ( month == 3 ) it = it + 59*4
    if ( month == 4 ) it = it + 90*4
    if ( month == 5 ) it = it + 120*4
    if ( month == 6 ) it = it + 151*4
    if ( month == 7 ) it = it + 181*4
    if ( month == 8 ) it = it + 212*4
    if ( month == 9 ) it = it + 243*4
    if ( month == 10 ) it = it + 273*4
    if ( month == 11 ) it = it + 304*4
    if ( month == 12 ) it = it + 334*4
    if ( mod(day,4) == 0 .and. month > 2 ) it = it + 4
    if ( mod(day,100) == 0 .and. month > 2 ) it = it - 4
    if ( mod(day,400) == 0 .and. month > 2 ) it = it + 4
!bxq_
    do m = 1 , 4
      istart(m) = 1
    end do
    do m = 5 , 10
      istart(m) = 0
      icount(m) = 0
    end do
    icount(1) = iii
    icount(2) = jjj
    icount(4) = 1460
    if ( mod(day,4) == 0 ) icount(4) = 1464
    if ( mod(day,100) == 0 ) icount(4) = 1460
    if ( mod(day,400) == 0 ) icount(4) = 1464
    istart(4) = it
    icount(4) = 1
    inet = inet7(kkrec)
    if ( nlev > 0 ) then
      icount(3) = nlev + 1
      istatus = nf90_get_var(inet,ivar7(kkrec),work,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    'read error in file'//trim(pathaddname))
      xscale = xscl(kkrec)
      xadd = xoff(kkrec)
      do ilev = 1 , nlev
        if ( kkrec == 1 ) then
          do j = 1 , jjj
            jj = j0 + j
            if ( i0 > i1 ) then
              do ii = i0 , ilon
                i = ii - i0 + 1
                tvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
              do ii = 1 , i1
                i = ii + (ilon-i0) + 1
                tvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            else
              do ii = i0 , i1
                i = ii - i0 + 1
                tvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            end if
          end do
 
        else if ( kkrec == 2 ) then
          do j = 1 , jjj
            jj = j0 + j
            if ( i0 > i1 ) then
              do ii = i0 , ilon
                i = ii - i0 + 1
                hvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
              do ii = 1 , i1
                i = ii + (ilon-i0) + 1
                hvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            else
              do ii = i0 , i1
                i = ii - i0 + 1
                hvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            end if
          end do
        else if ( kkrec == 3 ) then
          do j = 1 , jjj
            jj = j0 + j
            if ( i0 > i1 ) then
              do ii = i0 , ilon
                i = ii - i0 + 1
                rhvar(ii,jj,ilev) = real(dmin1((dble(work(i,j,ilev+1))*xscale+ &
                                    xadd)*0.01D0,1.D0))
              end do
              do ii = 1 , i1
                i = ii + (ilon-i0) + 1
                rhvar(ii,jj,ilev) = real(dmin1((dble(work(i,j,ilev+1))*xscale+ &
                                    xadd)*0.01D0,1.D0))
              end do
            else
              do ii = i0 , i1
                i = ii - i0 + 1
                rhvar(ii,jj,ilev) = real(dmin1((dble(work(i,j,ilev+1))*xscale+ &
                                    xadd)*0.01D0,1.D0))
              end do
            end if
          end do
        else if ( kkrec == 4 ) then
          do j = 1 , jjj
            jj = j0 + j
            if ( i0 > i1 ) then
              do ii = i0 , ilon
                i = ii - i0 + 1
                uvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
              do ii = 1 , i1
                i = ii + (ilon-i0) + 1
                uvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            else
              do ii = i0 , i1
                i = ii - i0 + 1
                uvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            end if
          end do
        else if ( kkrec == 5 ) then
          do j = 1 , jjj
            jj = j0 + j
            if ( i0 > i1 ) then
              do ii = i0 , ilon
                i = ii - i0 + 1
                vvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
              do ii = 1 , i1
                i = ii + (ilon-i0) + 1
                vvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            else
              do ii = i0 , i1
                i = ii - i0 + 1
                vvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            end if
          end do
        else if ( kkrec == 6 ) then
          do j = 1 , jjj
            jj = j0 + j
            if ( i0 > i1 ) then
              do ii = i0 , ilon
                i = ii - i0 + 1
                wvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
              do ii = 1 , i1
                i = ii + (ilon-i0) + 1
                wvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            else
              do ii = i0 , i1
                i = ii - i0 + 1
                wvar(ii,jj,ilev) = real(dble(work(i,j,ilev+1))*xscale + xadd)
              end do
            end if
          end do
        end if
      end do
    else if ( nlev == 0 ) then
      icount(3) = nlev
      istatus = nf90_get_var(inet,ivar7(kkrec),work,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    'read error in file'//trim(pathaddname))
      if ( kkrec == 7 ) then
        do j = 1 , jjj
          jj = j0 + j
          if ( i0 > i1 ) then
            do ii = i0 , ilon
              i = ii - i0 + 1
              psvar(ii,jj) = real(dble(work(i,j,1))*xscale + xadd)
            end do
            do ii = 1 , i1
              i = ii + (ilon-i0) + 1
              psvar(ii,jj) = real(dble(work(i,j,1))*xscale + xadd)
            end do
          else
            do ii = i0 , i1
              i = ii - i0 + 1
              psvar(ii,jj) = real(dble(work(i,j,1))*xscale + xadd)
            end do
          end if
        end do
      end if
    else
    end if
  end do


99001 format (i4,'/',a8,i4,'.nc')
99002 format (i4,'/',a9,i4,'.nc')
99003 format (i4,'/',a10,i4,'.nc')
99004 format (i4,'/',a13,i4,'.nc')
!
  end subroutine cdc6hour2

  subroutine headernc
  implicit none
!
  integer(ik4) :: i , j , k , kr
!
!     X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X
!     X X A1
!
!     ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!     jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!     klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!
  sigmar(1) = .07
  sigmar(2) = .1
  sigmar(3) = .15
  sigmar(4) = .2
  sigmar(5) = .25
  sigmar(6) = .3
  sigmar(7) = .4
  sigmar(8) = .5
  sigmar(9) = .6
  sigmar(10) = .7
  sigmar(11) = .85
  sigmar(12) = .925
  sigmar(13) = 1.0
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
  do i = 1 , ilon
    glon(i) = float(i-1)*2.5
  end do
  do j = 1 , jlat
    glat(j) = -90.0 + float(j-1)*2.5
  end do
!HH:OVER
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
  do k = 1 , klev
    kr = klev - k + 1
    sigma1(k) = sigmar(kr)
  end do
 
  call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_ncep:work')
  call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ncep:b3')
  call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ncep:d3')

!     Set up pointers

  u3 => d3(:,:,1:klev)
  v3 => d3(:,:,klev+1:2*klev)
  t3 => b3(:,:,1:klev)
  h3 => b3(:,:,klev+1:2*klev)
  q3 => b3(:,:,2*klev+1:3*klev)
  uvar => d2(:,:,1:klev)
  vvar => d2(:,:,klev+1:2*klev)
  tvar => b2(:,:,1:klev)
  hvar => b2(:,:,klev+1:2*klev)
  rhvar => b2(:,:,2*klev+1:3*klev)

  end subroutine headernc
!
end module mod_ncep
