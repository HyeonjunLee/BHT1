!-----------------------------------------------------------
! File: bht.f90
! Description: Transient heat conduction in a composite cylindrical shell
! with non-uniform grid, firing events, optional chrome coating,
! and input parameters read from "input.txt".
!-----------------------------------------------------------
program transient_heat_conduction_input
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  ! Simulation parameter variables (read from file)
  real(dp) :: T_outside, inner_diam, steel_thick, d_chrome_val, V
  real(dp) :: T_air_initial, T_fire, total_time, dt
  integer :: include_coating, N_coat, N_steel, num_firing, i
  integer, allocatable :: firing_events(:)
  
  ! Variables for simulation results
  real(dp), allocatable :: r(:)
  real(dp), allocatable :: T_record(:,:)  ! (num_snapshots, grid nodes)
  real(dp), allocatable :: times(:)
  real(dp), allocatable :: T_air_record(:)
  real(dp) :: r1, inner_d, r_chrome_out

  !-----------------------------------------------------------------
  ! Read input parameters from "input.txt"
  !-----------------------------------------------------------------
  open(unit=10, file="input.txt", status="old", action="read")
    read(10,*) T_outside         ! Ambient temperature in °C
    read(10,*) inner_diam        ! Inner diameter (m)
    read(10,*) steel_thick       ! Steel thickness (m)
    read(10,*) d_chrome_val      ! Chrome coating thickness (m)
    read(10,*) include_coating   ! 1: include, 0: exclude
    read(10,*) V                 ! Internal air flow velocity (m/s)
    read(10,*) T_air_initial     ! Initial internal air temperature (°C)
    read(10,*) T_fire            ! Temperature after firing (°C)
    read(10,*) total_time        ! Total simulation time (s)
    read(10,*) dt                ! Time step (s)
    read(10,*) N_coat            ! Number of grid points in coating region
    read(10,*) N_steel           ! Number of grid points in steel region
    read(10,*) num_firing        ! Number of firing events
    allocate(firing_events(num_firing))
    do i = 1, num_firing
       read(10,*) firing_events(i)
    end do
  close(10)

  ! Call the simulation subroutine (parameters for the shell region are passed in mm)
  call transient_heat_conduction_nonuniform_with_firing( T_outside, inner_diam*1000.0_dp, &
       steel_thick*1000.0_dp, d_chrome_va

