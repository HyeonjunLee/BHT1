!-----------------------------------------------------------
! File: transient_heat_conduction.f90
! Description: Transient heat conduction in a composite cylindrical shell
!   with non-uniform grid, firing events (reset internal air temperature),
!   and optional chrome coating.
!-----------------------------------------------------------
program transient_heat_conduction
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  ! Simulation parameters
  real(dp), parameter :: T_outside    = 25.0_dp       ! Ambient temperature (°C)
  real(dp), parameter :: inner_diam   = 0.155_dp      ! Inner diameter (m)
  real(dp), parameter :: steel_thick  = 0.0775_dp     ! Steel thickness (m)
  real(dp), parameter :: d_chrome_val = 0.0001_dp     ! Chrome coating thickness (m)
  logical, parameter :: include_coating = .true.       ! Option: .true. -> include coating, .false. -> steel-only
  real(dp), parameter :: V            = 2.0_dp        ! Internal air flow velocity (m/s)
  real(dp), parameter :: T_air_initial= 2500.0_dp     ! Initial internal air temperature (°C)
  real(dp), parameter :: T_fire       = 2500.0_dp     ! Temperature after firing (°C)
  real(dp), parameter :: total_time   = 300.0_dp      ! Total simulation time (s)
  real(dp), parameter :: dt           = 0.001_dp      ! Time step (s)
  integer, parameter :: N_coat        = 10            ! Grid points in coating region
  integer, parameter :: N_steel       = 140           ! Grid points in steel region

  ! Derived geometry
  real(dp) :: r1, r_chrome_out, r2   ! r1: inner radius; r_chrome_out: interface; r2: outer radius
  real(dp) :: inner_d                  ! inner diameter in m (same as inner_diam)
  integer :: N_total
  ! Grid and temperature arrays
  real(dp), allocatable :: r(:)       ! Non-uniform radial grid (m)
  real(dp), allocatable :: T(:)       ! Temperature at grid nodes (°C)
  real(dp), allocatable :: T_new(:)
  
  ! Time integration variables
  integer :: n, nsteps, i
  real(dp) :: time, dTdt, dT_air, h_current
  real(dp) :: T_air                   ! Internal air temperature (°C)
  real(dp) :: A_inner, m_air          ! Inner surface area and mass of air (per unit length)
  
  ! Variables for finite volume update (node 1: inner boundary)
  real(dp) :: r_face0, V0, k0, cp0, rho0, k1, k_eff, F_cond, F_conv, T_new1
  ! Variables for internal nodes update
  real(dp) :: r_im, r_ip, V_i, k_i, cp_i, rho_i, k_left, k_right, F_left, F_right

  ! Constants
  real(dp), parameter :: pi = 3.141592653589793_dp
  real(dp), parameter :: rho_air_const = 1.2_dp      ! kg/m³ for air
  real(dp), parameter :: cp_air_const  = 1000.0_dp     ! J/(kg·K) for air

  ! For convenience, set inner_d equal to inner_diam
  inner_d = inner_diam

  ! Compute geometry
  r1 = inner_d / 2.0_dp
  if (include_coating) then
    r_chrome_out = r1 + d_chrome_val
  else
    r_chrome_out = r1
  end if
  r2 = r_chrome_out + steel_thick

  ! Generate non-uniform grid:
  ! In coating region: from r1 to r_chrome_out, N_coat points (uniform)
  ! In steel region: from r_chrome_out to r2, N_steel points (uniform), avoiding duplicate interface
  N_total = N_coat + N_steel - 1
  allocate(r(N_total))
  call generate_nonuniform_grid(r, r1, r_chrome_out, r2, N_coat, N_steel)

  ! Allocate temperature arrays
  allocate(T(N_total))
  allocate(T_new(N_total))
  
  ! Set initial conditions: shell temperature = T_outside everywhere
  T = T_outside

  ! Set initial internal air temperature
  T_air = T_air_initial

  ! Pre-calculate inner surface area and air mass (per unit length)
  A_inner = 2.0_dp * pi * r1
  m_air = rho_air_const * pi * r1**2

  ! Total number of time steps
  nsteps = int(total_time / dt)

  ! Main time integration loop
  do n = 1, nsteps
     time = n * dt

     ! Firing events: if current time is within dt/2 of any event (60, 120, 180 s), reset T_air
     if ( (abs(time - 60.0_dp) < dt/2.0_dp) .or. &
          (abs(time - 120.0_dp) < dt/2.0_dp) .or. &
          (abs(time - 180.0_dp) < dt/2.0_dp) ) then
        T_air = T_fire
     end if

     ! Compute convection coefficient using Dittus–Boelter (inner diameter used)
     h_current = dittus_boelter_h(V, inner_d, T_air + 273.15_dp, 'cooling')

     ! --- Node 1 update (inner boundary) ---
     ! Control volume for node 1: from r = r1 to face between node 1 and node 2
     r_face0 = 0.5_dp * (r(1) + r(2))
     V0 = pi * (r_face0**2 - r1**2)
     ! Get material properties at node 1:
     if (include_coating .and. (r(1) <= r_chrome_out)) then
        k0  = k_chrome(T(1))
        cp0 = cp_chrome(T(1))
        rho0= 7190.0_dp
     else
        k0  = k_steel(T(1))
        cp0 = cp_steel(T(1))
        rho0= 7850.0_dp
     end if
     ! For node 2:
     if (include_coating .and. (r(2) <= r_chrome_out)) then
        k1 = k_chrome(T(2))
     else
        k1 = k_steel(T(2))
     end if
     k_eff = 0.5_dp * (k0 + k1)
     F_cond = - k_eff * (T(2) - T(1)) / (r(2) - r(1)) * (2.0_dp * pi * r_face0)
     F_conv = h_current * A_inner * (T_air - T(1))
     dTdt = (F_conv + F_cond) / (rho0 * cp0 * V0)
     T_new1 = T(1) + dt * dTdt

     ! --- Internal nodes update: for i = 2 to N_total-1 ---
     do i = 2, N_total - 1
        r_im = 0.5_dp * (r(i) + r(i-1))
        r_ip = 0.5_dp * (r(i+1) + r(i))
        V_i = pi * (r_ip**2 - r_im**2)
        ! Get properties at node i:
        if (include_coating .and. (r(i) <= r_chrome_out)) then
           k_i  = k_chrome(T(i))
           cp_i = cp_chrome(T(i))
           rho_i = 7190.0_dp
        else
           k_i  = k_steel(T(i))
           cp_i = cp_steel(T(i))
           rho_i = 7850.0_dp
        end if
        ! Left interface: between node i-1 and i
        if (include_coating .and. (r(i-1) <= r_chrome_out)) then
           k_left = k_chrome(0.5_dp*(T(i) + T(i-1)))
        else
           k_left = k_steel(0.5_dp*(T(i) + T(i-1)))
        end if
        F_left = - k_left * (T(i) - T(i-1)) / (r(i) - r(i-1)) * (2.0_dp * pi * r_im)
        ! Right interface: between node i and i+1
        if (include_coating .and. (r(i+1) <= r_chrome_out)) then
           k_right = k_chrome(0.5_dp*(T(i) + T(i+1)))
        else
           k_right = k_steel(0.5_dp*(T(i) + T(i+1)))
        end if
        F_right = - k_right * (T(i+1) - T(i)) / (r(i+1) - r(i)) * (2.0_dp * pi * r_ip)
        dTdt = (F_left - F_right) / (rho_i * cp_i * V_i)
        T_new(i) = T(i) + dt * dTdt
     end do

     ! --- Outer boundary: node N_total is held at T_outside
     T_new(N_total) = T_outside

     ! --- Update internal air temperature (lumped model)
     dT_air = - (h_current * A_inner * (T_air - T(1))) / (rho_air_const * cp_air_const * (pi * r1**2))
     T_air = T_air + dt * dT_air

     ! --- Update T array for next time step ---
     T(1) = T_new1
     do i = 2, N_total - 1
        T(i) = T_new(i)
     end do
     T(N_total) = T_outside

  end do  ! End of time loop

  ! 출력: 최종 내부 공기 온도와 셸의 온도 분포
  print *, 'Final internal air temperature (°C): ', T_air
  print *, 'Final temperature profile in shell:'
  do i = 1, N_total
     print *, r(i), T(i)
  end do

contains

  !-----------------------------------------------------------
  ! Material property functions (T in °C)
  !-----------------------------------------------------------
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

  !-----------------------------------------------------------
  ! Dittus–Boelter correlation: compute convective heat transfer coefficient
  ! V: velocity (m/s), D: characteristic diameter (m), T_air: air temperature in Kelvin.
  !-----------------------------------------------------------
  pure function dittus_boelter_h(V, D, T_air, mode) result(h)
    real(dp), intent(in) :: V, D, T_air
    character(len=*), intent(in) :: mode
    real(dp) :: h, Re, Pr, Nu, n
    ! Use constant air properties
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

  !-----------------------------------------------------------
  ! Generate a non-uniform grid from r1 to r2.
  ! In the coating region [r1, r_chrome_out]: use N_coat uniformly spaced points.
  ! In the steel region [r_chrome_out, r2]: use N_steel uniformly spaced points.
  ! The interface (r_chrome_out) is not duplicated.
  !-----------------------------------------------------------
  subroutine generate_nonuniform_grid(r, r1, r_chrome_out, r2, N_coat, N_steel)
    real(dp), intent(out) :: r(:)
    real(dp), intent(in) :: r1, r_chrome_out, r2
    integer, intent(in) :: N_coat, N_steel
    integer :: i
    real(dp) :: dr_coat, dr_steel
    ! Coating region
    dr_coat = (r_chrome_out - r1) / real(N_coat - 1, dp)
    do i = 1, N_coat
       r(i) = r1 + (i - 1) * dr_coat
    end do
    ! Steel region (starting from r_chrome_out; skip duplicate)
    dr_steel = (r2 - r_chrome_out) / real(N_steel - 1, dp)
    do i = 2, N_steel
       r(N_coat + i - 1) = r_chrome_out + (i - 1) * dr_steel
    end do
  end subroutine generate_nonuniform_grid

end program transient_heat_conduction

