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
       steel_thick*1000.0_dp, d_chrome_val*1000.0_dp, V, 1.2_dp, 1000.0_dp, T_air_initial, &
       total_time, dt, N_coat, N_steel, firing_events, T_fire, (include_coating == 1), &
       r, T_record, times, T_air_record, r1, inner_d, r_chrome_out )

  ! Check for NaN values in the results
  if (any(isnan(T_air_record)) .or. any(isnan(T_record))) then
     print *, 'Error: NaN values detected in the simulation results.'
     stop
  end if

  print '(A, F8.2)', 'Final internal air temperature (°C): ', T_air_record(size(T_air_record))
  print *, 'Final temperature profile in shell:'
  do i = 1, size(r)
     print '(F8.2, F8.2)', r(i), T_record(size(T_record,1), i)
  end do

contains

  subroutine transient_heat_conduction_nonuniform_with_firing( T_outside, inner_diameter_mm, steel_thickness_mm, d_chrome_mm, &
         V, rho_air, c_p_air, T_air_initial, total_time, dt, N_coat, N_steel, firing_events, T_fire, &
         include_coating, r, T_record, times, T_air_record, r1, inner_d, r_chrome_out )
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    ! Input parameters
    real(dp), intent(in) :: T_outside
    real(dp), intent(in) :: inner_diameter_mm, steel_thickness_mm, d_chrome_mm
    real(dp), intent(in) :: V, rho_air, c_p_air, T_air_initial, total_time, dt, T_fire
    integer, intent(in) :: N_coat, N_steel
    integer, intent(in) :: firing_events(:)
    logical, intent(in) :: include_coating
    ! Output arrays
    real(dp), allocatable, intent(out) :: r(:)
    real(dp), allocatable, intent(out) :: T_record(:,:)  ! (num_snapshots, N)
    real(dp), allocatable, intent(out) :: times(:)
    real(dp), allocatable, intent(out) :: T_air_record(:)
    ! Output geometry
    real(dp), intent(out) :: r1, inner_d, r_chrome_out

    ! Local variables
    integer :: N, nstep, num_steps, num_snapshots, i, idx
    real(dp) :: time, dr_min, dt_max, h_current, F_conv, F_cond, dTdt, dT_air
    real(dp) :: A_inner, m_air, r_face0, V0, k0, cp0, rho0, k1, k_eff
    real(dp) :: r_im, r_ip, V_i, k_i, cp_i, rho_i, F_left, F_right
    real(dp) :: k_left, k_right  ! 추가된 변수 선언
    real(dp) :: T_air           ! 내부 공기 온도
    integer :: num_firing_local
    real(dp), allocatable :: T(:), T_new(:)
    real(dp), allocatable :: temp_snapshots(:,:)
    real(dp), allocatable :: time_snapshots(:)
    real(dp), allocatable :: T_air_snapshots(:)
    real(dp) :: r2_local

    ! Set number of snapshots to 50
    num_snapshots = 50

    ! Geometry conversion: inner_d in m
    inner_d = inner_diameter_mm / 1000.0_dp
    r1 = inner_d / 2.0_dp
    if (include_coating) then
      r_chrome_out = r1 + d_chrome_mm / 1000.0_dp
    else
      r_chrome_out = r1
    end if
    r2_local = r_chrome_out + steel_thickness_mm / 1000.0_dp

    ! Generate non-uniform grid over shell region
    allocate(r(N_coat + N_steel - 1))
    call generate_nonuniform_grid(r, r1, r_chrome_out, r2_local, N_coat, N_steel)
    N = size(r)
    dr_min = minval( r(2:N) - r(1:N-1) )
    dt_max = dr_min**2 / (2.0_dp * 1.0_dp)   ! (여기서는 대략 상수 1.0 사용)
    if (dt > dt_max) then
       print *, 'Warning: dt too large. dt_max = ', dt_max
    end if

    num_steps = int(total_time / dt)
    allocate(time_snapshots(num_snapshots))
    allocate(temp_snapshots(num_snapshots, N))
    allocate(T_air_snapshots(num_snapshots))

    allocate(T(N))
    allocate(T_new(N))
    T = T_outside
    T_air = T_air_initial

    A_inner = 2.0_dp * 3.141592653589793_dp * r1
    m_air = rho_air * 3.141592653589793_dp * r1**2

    idx = 0
    do nstep = 0, num_steps
      time = nstep * dt

      do i = 1, size(firing_events)
         if (abs(time - real(firing_events(i), dp)) < dt/2.0_dp) then
            T_air = T_fire
         end if
      end do

      if (mod(nstep, max(1, num_steps/num_snapshots)) == 0) then
         idx = idx + 1
         time_snapshots(idx) = time
         temp_snapshots(idx, :) = T
         T_air_snapshots(idx) = T_air
      end if

      ! --- Node 1 update (inner boundary) ---
      r_face0 = 0.5_dp * (r(1) + r(2))
      V0 = 3.141592653589793_dp * (r_face0**2 - r1**2)
      if (include_coating .and. (r(1) <= r_chrome_out)) then
         k0  = k_chrome(T(1))
         cp0 = cp_chrome(T(1))
         rho0= 7190.0_dp
      else
         k0  = k_steel(T(1))
         cp0 = cp_steel(T(1))
         rho0= 7850.0_dp
      end if
      if (include_coating .and. (r(2) <= r_chrome_out)) then
         k1 = k_chrome(T(2))
      else
         k1 = k_steel(T(2))
      end if
      k_eff = 0.5_dp * (k0 + k1)
      F_cond = - k_eff * (T(2) - T(1)) / (r(2) - r(1)) * (2.0_dp * 3.141592653589793_dp * r_face0)
      h_current = dittus_boelter_h(V, inner_d, T_air + 273.15_dp, 'cooling')
      F_conv = h_current * A_inner * (T_air - T(1))
      dTdt = (F_conv + F_cond) / (rho0 * cp0 * V0)
      T_new(1) = T(1) + dt * dTdt

      ! --- Internal nodes update: i = 2 to N-1 ---
      do i = 2, N - 1
         r_im = 0.5_dp * (r(i) + r(i-1))
         r_ip = 0.5_dp * (r(i+1) + r(i))
         V_i = 3.141592653589793_dp * (r_ip**2 - r_im**2)
         if (include_coating .and. (r(i) <= r_chrome_out)) then
            k_i  = k_chrome(T(i))
            cp_i = cp_chrome(T(i))
            rho_i = 7190.0_dp
         else
            k_i  = k_steel(T(i))
            cp_i = cp_steel(T(i))
            rho_i = 7850.0_dp
         end if
         if (include_coating .and. (r(i-1) <= r_chrome_out)) then
            k_left = k_chrome(0.5_dp*(T(i)+T(i-1)))
         else
            k_left = k_steel(0.5_dp*(T(i)+T(i-1)))
         end if
         F_left = - k_left * (T(i) - T(i-1)) / (r(i) - r(i-1)) * (2.0_dp * 3.141592653589793_dp * r_im)
         if (include_coating .and. (r(i+1) <= r_chrome_out)) then
            k_right = k_chrome(0.5_dp*(T(i)+T(i+1)))
         else
            k_right = k_steel(0.5_dp*(T(i)+T(i+1)))
         end if
         F_right = - k_right * (T(i+1) - T(i)) / (r(i+1) - r(i)) * (2.0_dp * 3.141592653589793_dp * r_ip)
         dTdt = (F_left - F_right) / (rho_i * cp_i * V_i)
         T_new(i) = T(i) + dt * dTdt
      end do

      T_new(N) = T_outside

      dT_air = - (h_current * A_inner * (T_air - T(1))) / (rho_air * c_p_air * (3.141592653589793_dp * r1**2))
      T_air = T_air + dt * dT_air

      T = T_new

      ! Print variables at each step
      print '(A, I8)', 'Step: ', nstep
      print '(A, F8.2)', 'Time: ', time
      print '(A, F8.2)', 'T_air: ', T_air
      print '(A, F8.2)', 'T: ', T
      print '(A, F8.2)', 'T_new: ', T_new
      print '(A, F8.2)', 'dT_air: ', dT_air
      print '(A, F8.2)', 'F_conv: ', F_conv
      print '(A, F8.2)', 'F_cond: ', F_cond
      print '(A, F8.2)', 'dTdt: ', dTdt

      ! Check for NaN values in the temperature array
      if (any(isnan(T))) then
         print *, 'Error: NaN values detected at time step ', nstep
         stop
      end if

    end do

    allocate(T_record(idx, N))
    T_record = temp_snapshots(1:idx, :)
    allocate(times(idx))
    times = time_snapshots(1:idx)
    allocate(T_air_record(idx))
    T_air_record = T_air_snapshots(1:idx)
  end subroutine transient_heat_conduction_nonuniform_with_firing

  pure function k_chrome(T) result(k)
    real(dp), intent(in) :: T
    real(dp) :: k
    k = 93.0_dp - 0.008163_dp * (T - 20.0_dp)
  end function k_chrome

  pure function k_steel(T) result(k)
    real(dp), intent(in) :: T
    real(dp) :: k
    k = 80.0_dp - 0.010204_dp * (T - 20.0_dp)
  end function k_steel

  pure function cp_chrome(T) result(cp)
    real(dp), intent(in) :: T
    real(dp) :: cp
    cp = 450.0_dp + 0.05102_dp * (T - 20.0_dp)
  end function cp_chrome

  pure function cp_steel(T) result(cp)
    real(dp), intent(in) :: T
    real(dp) :: cp
    cp = 450.0_dp + 0.07143_dp * (T - 20.0_dp)
  end function cp_steel

  pure function dittus_boelter_h(V, D, T_air, mode) result(h)
    real(dp), intent(in) :: V, D, T_air
    character(len=*), intent(in) :: mode
    real(dp) :: h, Re, Pr, Nu, n
    real(dp), parameter :: rho = 1.2_dp, mu = 1.8e-5_dp, k_air = 0.026_dp, cp_air = 1000.0_dp
    Re = (rho * V * D) / mu
    Pr = (cp_air * mu) / k_air
    if (mode == 'cooling') then
      n = 0.3_dp
    else
      n = 0.4_dp
    end if
    Nu = 0.023_dp * Re**0.8_dp * Pr**n
    h = (Nu * k_air) / D
  end function dittus_boelter_h

  subroutine generate_nonuniform_grid(r, r1, r_chrome_out, r2, N_coat, N_steel)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(out) :: r(:)
    real(dp), intent(in) :: r1, r_chrome_out, r2
    integer, intent(in) :: N_coat, N_steel
    integer :: i, total_points
    real(dp) :: dr_coat, dr_steel
    total_points = N_coat + N_steel - 1
    dr_coat = (r_chrome_out - r1) / real(N_coat - 1, dp)
    do i = 1, N_coat
       r(i) = r1 + (i - 1) * dr_coat
    end do
    dr_steel = (r2 - r_chrome_out) / real(N_steel - 1, dp)
    do i = 2, N_steel
       r(N_coat + i - 1) = r_chrome_out + (i - 1) * dr_steel
    end do
  end subroutine generate_nonuniform_grid

end program transient_heat_conduction_input

