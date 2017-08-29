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

module mod_regcm_interface
  use omp_lib
  use mod_memutil
  use mod_service
  use mod_che_interface
  use mod_atm_interface
  use mod_pbl_interface
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_ncio
  use mod_ncout
  use mod_output
  use mod_split
  use mod_bdycod
  use mod_init
  use mod_header
  use mod_params
  use mod_tendency
  use mod_tstep
  use mod_service
  use mod_mppio
  use mod_cloud_s1
  use mod_sun
#ifdef CLM
  use perf_mod
  use mod_mtrxclm
  use spmdMod, only: mpicom
  use clm_varsur , only : numdays
#endif
#ifdef ESMFCPL
  use mod_update, only: ocn_put => RCM_PutExportData
  use mod_update, only: ocn_get => RCM_GetImportData
#endif
  implicit none
  include 'mpif.h'
!
  private
  public :: RCM_initialize
  public :: RCM_run
  public :: RCM_finalize

  real(rk8) :: dtinc
  real(rk8) :: extime

  real(rk8) :: tend_time
  real(rk8) :: split_time
  real(rk8) :: solar_time
  real(rk8) :: new_bound_time
  real(rk8) :: fill_bound_time
  integer(ik4) :: no_loop

  data extime /d_zero/
  contains
 
  subroutine RCM_initialize(mpiCommunicator)
    implicit none
    integer, intent(in), optional :: mpiCommunicator
!
    integer(ik4) :: ierr
    character(256) :: namelistfile, prgname
! 
!**********************************************************************
!
!   MPI Initialization
!
!**********************************************************************
!
    if (present(mpiCommunicator)) then
      mycomm = mpiCommunicator
    else
      mycomm = MPI_COMM_WORLD
    end if
    call mpi_comm_rank(mycomm, myid, ierr)
    call mpi_comm_size(mycomm, nproc, ierr)
#ifndef MPI_SERIAL
    call mpi_comm_set_errhandler(mycomm, mpi_errors_return, ierr)
#endif
!
    call whoami(myid)
    call setup_mesg(myid)
!
#ifdef DEBUG 
    call activate_debug()
#endif
!
!**********************************************************************
!
!   Read input global namelist
!
!**********************************************************************
!
    if ( myid == iocpu ) then
      call getarg(0, prgname)
      call getarg(1, namelistfile)
      call initparam(namelistfile, ierr)
      if ( ierr /= 0 ) then
        write ( 6, * ) 'Parameter initialization not completed'
        write ( 6, * ) 'Usage : '
        write ( 6, * ) '          ', trim(prgname), ' regcm.in'
        write ( 6, * ) ' '
        write ( 6, * ) 'Check argument and namelist syntax'
        stop
      end if
    end if
!
    call broadcast_params

    call memory_init
!
    call header(myid,nproc)
    call set_nproc
    call setup_model_indexes
!
#ifdef DEBUG 
    call start_debug()
#endif 
!
!**********************************************************************
!
!   Parameter Setup
!
!**********************************************************************
!
    call param
    dtinc = dt
!
!**********************************************************************
!
!   Read initial data and boundary conditions
!
!**********************************************************************
!
    ! Update solar constant from TSI dataset
    solcon = solar_irradiance( )
    scon = solcon*d_1000
    !
    ! Calculate solar declination angle at startup
    !
    if ( myid == italk ) then
      write (stdout,*) 'Calculate solar declination angle at ',toint10(idatex)
    end if
#ifdef CLM
    numdays = dayspy
    call solar_clm(idatex,calday,declin,xyear)
#else
    call solar1
#endif
    call init_bdy
!
!**********************************************************************
!
!   Initialize data (from IC or restart)
!
!**********************************************************************
!
    call init
!
!**********************************************************************
!
!   Initialize split explicit scheme
!
!**********************************************************************
!
    call spinit
!
!**********************************************************************
!
!   Setup the output files
!
!**********************************************************************
!
    call init_output_streams(do_parallel_netcdf_out)
    call output
!
!**********************************************************************
!
!   Set the boundary conditions for this timestep
!
!**********************************************************************
!
    call bdyval(xbctime)
!
!**********************************************************************
!
!   Calculate Zenital Angle
!
!**********************************************************************
!
#ifdef CLM
    call zenit_clm(coszrs)
#else
    call zenitm(coszrs)
#endif
!
!**********************************************************************
!
!   Clean up and logging
!
!**********************************************************************
!
#ifdef DEBUG
    call time_print(6,'inizialization phase')
    call time_reset()
#endif
  end subroutine RCM_initialize
!
!=======================================================================
!                                                                      !
!     This routine runs RegCM model from specified starting (TimeStr)  !
!     to ending (TimeEnd) time-steps.                                  !
!                                                                      !
!=======================================================================
!
  subroutine RCM_run(timestr, timeend)
    implicit none
    real(rk8) , intent(in) :: timestr   ! starting time-step
    real(rk8) , intent(in) :: timeend   ! ending   time-step
    real(rk8) :: start_sub_time
    real(rk8) :: finish_sub_time
    real(rk8) :: start_loop_time
    real(rk8) :: end_loop_time
    integer(ik4) :: no_iteration, i
    character(len=32) :: appdat
!
#ifdef DEBUG
    ! if ( enable_newmicro ) call grid_nc_create('qqxp',cross,zqxn,qqxp)
    ! call grid_nc_create('qxatm',cross,atm1%qx,nc_4d)
#endif
    if (myid == italk) then
      write(stdout,*) 'Execution Time = ', extime
      write(stdout,*) 'Start Time = ', timestr
      write(stdout,*) 'End Time = ', timeend
    endif
    
    no_iteration = (timeend - timestr)/dtinc

    !$OMP PARALLEL DO  
    do i = 1, no_iteration
    !do while ( extime >= timestr .and. extime < timeend)
      if ( extime < timestr .or. extime >= timeend)
        exit
      end if  
      call cpu_time(start_loop_time)
      no_loop = no_loop + 1
#ifdef DEBUG
      ! call grid_nc_write(nc_4d)
#endif
      if ( mod(ktau,kday) == 0 ) then
        solcon = solar_irradiance( )
        scon = solcon*d_1000
      end if
      !
      ! Refined start
      !
      if ( .not. ifrest ) then
        if ( rfstrt ) then
          if ( (ktau == 0) .or. dtinc /= deltmx ) then
            call tstep(extime,dtinc)
            if ( myid == italk ) then
              write(stdout, 99001) extime , dtinc , dt , dt2 , &
                                   dtsec , ktau , xyear
            end if
          end if
        end if
      end if
      !
      ! Get information from ocean model
      !
#ifdef ESMFCPL
      if ( iocncpl == 1 ) then
        if (ktau > ntcpl) then 
          call ocn_get(myid)
        end if
      end if
#endif
      !
      ! Compute tendencies
      !
      !write(stderr,*) 'Computing tendencies'
      call cpu_time(start_sub_time)
      call tend
      call cpu_time(finish_sub_time)
      tend_time = tend_time + (finish_sub_time - start_sub_time)
      !write(stderr,*) 'Tendencies Computation Time: ', &
      !        (finish_sub_time - start_sub_time)
      !
      ! Split modes
      !
      !write(stderr,*) 'Split Mode Start'
      call cpu_time(start_sub_time)
      call splitf
      call cpu_time(finish_sub_time)
      split_time = split_time + (finish_sub_time - start_sub_time)
      !write(stderr,*) 'Split Mode Time: ', &
      !        (finish_sub_time - start_sub_time)
      !
      ! Boundary code (do not execute at the end of run)
      !
      if ( ktau /= mtau ) then
        if ( nbdytime == 0 ) then
          !
          ! recalculate solar declination angle if reading bdy
          !
          if ( myid == italk ) then
            write (stdout,*) &
              'Calculate solar declination angle at ',toint10(idatex)
          end if
          !write(stderr,*) 'Solar declination angle Start'
          call cpu_time(start_sub_time)
#ifdef CLM
          if ( myid == italk ) then
            write(stderr,*) 'solar_clm'
          end if
          call solar_clm(idatex,calday,declin,xyear)
#else
          if ( myid == italk ) then
            write(stderr,*) 'solar1'
          end if
          call solar1
#endif
          call cpu_time(finish_sub_time)
          solar_time = solar_time + (finish_sub_time - start_sub_time)
          !write(stdout,*) 'Solar declination angle Time: ', &
          !    (finish_sub_time - start_sub_time)
          !
          ! Read in new boundary conditions
          !
          !write(stdout,*) 'Read new Boundary Start'
          call cpu_time(start_sub_time)
          call bdyin
          call cpu_time(finish_sub_time)
          new_bound_time = new_bound_time + (finish_sub_time - start_sub_time)
          !write(stdout,*) 'Read new Boundary Time: ', &
          !    (finish_sub_time - start_sub_time)

        end if
        !
        ! fill up the boundary values for xxb and xxa variables:
        !
        !write(stdout,*) 'Fill Boundary Start'
        call cpu_time(start_sub_time)
        call bdyval(xbctime)
        call cpu_time(finish_sub_time)
        !write(stdout,*) 'Fill Boundary Time: ', &
        !      (finish_sub_time - start_sub_time)
        fill_bound_time = fill_bound_time + (finish_sub_time - start_sub_time)

      end if
      !
      ! Write output for this timestep if requested
      !
      call output
      !
      ! Send information to ocean model
      !
#ifdef ESMFCPL
      if ( iocncpl == 1 ) then
        call ocn_put(myid)
      end if
#endif
      !
      ! Increment execution time
      !
      extime = extime + dtinc
      if ( debug_level > 3 ) then
        if ( myid == italk ) then
          appdat = tochar(idatex)
          write(6,'(a,a,f12.2)') 'Simulation time: ', appdat, extime
        end if
      end if
      call cpu_time(end_loop_time)
      if ( myid == italk) then
        write(stdout,*) 'Iteration time: ',  end_loop_time - start_loop_time
      end if
    end do
    !$OMP END PARALLEL DO
  
#ifdef DEBUG
    call stop_debug()
    ! if ( enable_newmicro ) call grid_nc_destroy(qqxp)
    ! call grid_nc_destroy(nc_4d)
    call time_print(6,'evolution phase')
#endif
!
99001 format (6x,'large domain: extime = ',f7.1,' dtinc = ',f7.1,       &
        & ' dt = ',f7.1,' dt2 = ',f7.1,' dtsec = ',f6.1,' ktau = ', &
        & i7,' in year ',i4)

  end subroutine RCM_run

  subroutine RCM_finalize
    implicit none
    character(len=32) :: appdat
!
    if ( myid == italk ) then
      appdat = tochar(idate2)
      write(stdout,*) 'Restart file for next run is written at time ',appdat
    end if
!
    call close_icbc
    if ( ichem == 1 ) call close_chbc
    call dispose_output_streams
!
#ifdef CLM
    call t_prf('timing_all',mpicom)
    call t_finalizef()
#endif
!
    call memory_destroy
    call finaltime(myid)
!
    if ( myid == italk ) then
      write(stdout,*) 'Number of loop: ', no_loop
      write(stdout,*) 'Tendencies time: ', tend_time
      write(stdout,*) 'Split time: ', split_time
      write(stdout,*) 'Solar time: ', solar_time
      write(stdout,*) 'New Boundary time: ', new_bound_time
      write(stdout,*) 'Fill Boundary time: ', fill_bound_time
      write(stdout,*) 'RegCM V4 simulation successfully reached end'

    end if
!
  end subroutine RCM_finalize
!
end module mod_regcm_interface
