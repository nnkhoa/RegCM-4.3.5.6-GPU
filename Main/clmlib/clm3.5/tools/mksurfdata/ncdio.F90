module ncdio

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdioMod
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_sys_mod    , only : shr_sys_flush
  use mkvarsur       , only : spval
!
! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  save
  public :: check_ret   ! checks return status of netcdf calls
  public :: check_var   ! determine if variable is on netcdf file
  public :: check_dim   ! validity check on dimension
!
! !REVISION HISTORY:
!
!EOP
!
! !PRIVATE METHODS:
!
  interface ncd_ioglobal
     module procedure ncd_ioglobal_int_var
     module procedure ncd_ioglobal_real_var
     module procedure ncd_ioglobal_int_1d
     module procedure ncd_ioglobal_real_1d
     module procedure ncd_ioglobal_int_2d
     module procedure ncd_ioglobal_real_2d
     module procedure ncd_ioglobal_int_3d
     module procedure ncd_ioglobal_real_3d
  end interface
  logical  :: masterproc = .true. ! always use 1 proc
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dim
!
! !INTERFACE:
  subroutine check_dim(ncid, dimname, value)
!
! !DESCRIPTION:
! Validity check on dimension
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer, intent(in) :: value
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: dimid, dimlen    ! temporaries
!-----------------------------------------------------------------------

    call check_ret(nf_inq_dimid (ncid, trim(dimname), dimid), 'check_dim')
    call check_ret(nf_inq_dimlen (ncid, dimid, dimlen), 'check_dim')
    if (dimlen /= value) then
       write (6,*) 'CHECK_DIM error: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for variable ',trim(dimname)
       call abort()
    end if

  end subroutine check_dim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var
!
! !INTERFACE:
  subroutine check_var(ncid, varname, varid, readvar)
!
! !DESCRIPTION:
! Check if variable is on netcdf file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: varname
    integer, intent(out)         :: varid
    logical, intent(out)         :: readvar
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ret     ! return value
!-----------------------------------------------------------------------

    readvar = .true.
    if (masterproc) then
       ret = nf_inq_varid (ncid, varname, varid)
       if (ret/=NF_NOERR) then
          write(6,*)'CHECK_VAR: variable ',trim(varname),' is not on initial dataset'
          call shr_sys_flush(6)
          readvar = .false.
       end if
    end if
  end subroutine check_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, calling)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ret
    character(len=*) :: calling
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then
       write(6,*)'netcdf error from ',trim(calling)
       write(6,*)'error ',ret,' = ',nf_strerror(ret)
       call abort()
    end if

  end subroutine check_ret

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar
!
! !INTERFACE:
  subroutine ncd_defvar(ncid, varname, xtype, &
       dim1name, dim2name, dim3name, dim4name, dim5name, &
       long_name, units, cell_method, missing_value, fill_value, &
       imissing_value, ifill_value)
!
! !DESCRIPTION:
!  Define a netcdf variable
!
! !ARGUMENTS:
    implicit none
    integer         , intent(in)  :: ncid                    ! input unit
    character(len=*), intent(in)  :: varname                 ! variable name
    integer         , intent(in)  :: xtype                   ! external type
    character(len=*), intent(in), optional :: dim1name       ! dimension name
    character(len=*), intent(in), optional :: dim2name       ! dimension name
    character(len=*), intent(in), optional :: dim3name       ! dimension name
    character(len=*), intent(in), optional :: dim4name       ! dimension name
    character(len=*), intent(in), optional :: dim5name       ! dimension name
    character(len=*), intent(in), optional :: long_name      ! attribute
    character(len=*), intent(in), optional :: units          ! attribute
    character(len=*), intent(in), optional :: cell_method    ! attribute
    real(r8)        , intent(in), optional :: missing_value  ! attribute for real
    real(r8)        , intent(in), optional :: fill_value     ! attribute for real
    integer         , intent(in), optional :: imissing_value ! attribute for int
    integer         , intent(in), optional :: ifill_value    ! attribute for int
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n              ! indices
    integer :: ndims          ! dimension counter
    integer :: dimid(5)       ! dimension ids
    integer :: varid          ! variable id
    integer :: itmp           ! temporary
    character(len=256) :: str ! temporary
    character(len=32) :: subname='NCD_DEFVAR_REAL' ! subroutine name
!-----------------------------------------------------------------------

    if (.not. masterproc) return

    ! Determine dimension ids for variable

    dimid(:) = 0

    if (present(dim1name)) then
       call check_ret(nf_inq_dimid(ncid, dim1name, dimid(1)), subname)
    end if
    if (present(dim2name)) then
       call check_ret(nf_inq_dimid(ncid, dim2name, dimid(2)), subname)
    end if
    if (present(dim3name)) then
       call check_ret(nf_inq_dimid(ncid, dim3name, dimid(3)), subname)
    end if
    if (present(dim4name)) then
       call check_ret(nf_inq_dimid(ncid, dim4name, dimid(4)), subname)
    end if
    if (present(dim5name)) then
       call check_ret(nf_inq_dimid(ncid, dim5name, dimid(5)), subname)
    end if

    ! Define variable

    if (present(dim1name)) then
       ndims = 0
       do n = 1, size(dimid)
          if (dimid(n) /= 0) ndims = ndims + 1
       end do
       call check_ret(nf_def_var(ncid, trim(varname), xtype, ndims, dimid(1:ndims), varid), subname)
    else
       call check_ret(nf_def_var(ncid, varname, xtype, 0, 0, varid), subname)
    end if
    if (present(long_name)) then
       call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
    end if
    if (present(units)) then
       call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call check_ret(nf_put_att_text(ncid, varid, 'cell_method', len_trim(str), trim(str)), subname)
    end if
    if (present(fill_value)) then
       call check_ret(nf_put_att_double(ncid, varid, '_FillValue', xtype, 1, fill_value), subname)
    end if
    if (present(missing_value)) then
       call check_ret(nf_put_att_double(ncid, varid, 'missing_value', xtype, 1, missing_value), subname)
    end if
    if (present(ifill_value)) then
       call check_ret(nf_put_att_int(ncid, varid, '_FillValue', xtype, 1, ifill_value), subname)
    end if
    if (present(imissing_value)) then
       call check_ret(nf_put_att_int(ncid, varid, 'missing_value', xtype, 1, imissing_value), subname)
    end if

  end subroutine ncd_defvar

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_var(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! I/O of integer variable
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), intent(in)  :: varname            ! variable name
    integer         , intent(in)  :: data               ! local decomposition data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: dimid(1)                       ! dimension id
    integer :: start(1), count(1)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_INT_VAR' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = nt; count(1) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_int(ncid, varid, data), subname)
       end if
#ifdef SPMD
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call abort()
       end if
       if (varpresent) then
          call mpi_bcast(data, 1, MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call abort()
          end if
       end if
#endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_var(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! I/O of real variable
!

! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), intent(in)  :: varname            ! variable name
    real(r8)        , intent(in)  :: data               ! local decomposition data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: dimid(1)                       ! dimension id
    integer :: start(1), count(1)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_VAR' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = nt; count(1) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_double(ncid, varid, data), subname)
       end if
#ifdef SPMD
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call abort()
       end if
       if (varpresent) then
          call mpi_bcast(data, 1, MPI_REAL8, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call abort()
          end if
       end if
#endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_1d(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! Master I/O for 1d integer data
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:)          ! local decomposition data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(2), ndims                ! dimension ids
    integer :: start(2), count(2)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_INT_1D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data)
             start(2) = nt; count(2) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_int(ncid, varid, data), subname)
       end if
#ifdef SPMD
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call abort()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call abort()
          end if
       end if
#endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_1d(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! Master I/O for 1d real data
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data(:)          ! local decomposition input data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: dimid(2), ndims                ! dimension ids
    integer :: start(2), count(2)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_1D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data)
             start(2) = nt; count(2) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_double(ncid, varid, data), subname)
       end if
#ifdef SPMD
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call abort()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call abort()
          end if
       end if
#endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_2d(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global 2d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:,:)        ! local decomposition input data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_2D_INT_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_int(ncid, varid, data), subname)
       end if
#ifdef SPMD
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call abort()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call abort()
          end if
       end if
#endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_2d(varname, data, long_name, units, flag, &
                                  ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global 2d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data(:,:)        ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_2D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_double(ncid, varid, data), subname)
       end if
#ifdef SPMD
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call abort()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call abort()
          end if
       end if
#endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_3d(varname, data, long_name, units, flag, &
                                 ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global 3d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:,:,:)      ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: dimid(4), ndims          ! dimension ids
    integer :: start(4), count(4)       ! output bounds
    integer :: ier                      ! error code
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_3D_INT_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = 1;  count(3) = size(data, dim=3)
             start(4) = nt; count(4) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_int(ncid, varid, data), subname)
       end if
#ifdef SPMD
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call abort()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call abort()
          end if
       end if
#endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_3d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_3d(varname, data, long_name, units, flag, &
                                  ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global 3d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data(:,:,:)      ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: ier                      ! error code
    integer :: dimid(4), ndims          ! dimension ids
    integer :: start(4), count(4)       ! output bounds
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_3D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = 1;  count(3) = size(data, dim=3)
             start(4) = nt; count(4) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_double(ncid, varid, data), subname)
       end if
#ifdef SPMD
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call abort()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call abort()
          end if
       end if
#endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_3d

end module ncdio
