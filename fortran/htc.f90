program heat_transfer_cylinder
    implicit none
    ! Constants
    real, parameter :: pi = 3.141592653589793
    real, parameter :: k_chrome = 93.0    ! Thermal conductivity of chrome (W/mK)
    real, parameter :: k_steel = 50.2     ! Thermal conductivity of steel (W/mK)
    real, parameter :: h_air = 1000.0     ! Convective heat transfer coefficient of air (W/m^2K)
    real, parameter :: T_air_initial = 2500.0 ! Initial temperature of air inside the barrel (C)
    real, parameter :: T_ambient = 25.0   ! Ambient temperature outside the barrel (C)
    real, parameter :: r_inner = 0.0775   ! Inner radius of the barrel (m)
    real, parameter :: r_chrome = 0.0785  ! Outer radius of the chrome layer (m)
    real, parameter :: r_outer = 0.155    ! Outer radius of the steel layer (m)
    real, parameter :: dr = 0.001         ! Radial step size (m)
    real, parameter :: dt = 0.01          ! Time step size (s)
    real, parameter :: t_end = 10.0       ! End time (s)
    integer, parameter :: n_r = int((r_outer - r_inner) / dr) + 1
    integer, parameter :: n_t = int(t_end / dt) + 1

    ! Arrays
    real :: T(n_r), T_new(n_r)
    real :: r(n_r)
    integer :: i, j

    ! Initialize radial positions and temperature distribution
    do i = 1, n_r
        r(i) = r_inner + (i - 1) * dr
        if (r(i) <= r_chrome) then
            T(i) = T_air_initial
        else
            T(i) = T_ambient
        end if
    end do

    ! Time-stepping loop
    do j = 1, n_t
        ! Update temperature distribution
        do i = 2, n_r-1
            if (r(i) <= r_chrome) then
                T_new(i) = T(i) + dt * k_chrome * ( &
                    (T(i+1) - T(i)) / (r(i+1) - r(i)) - &
                    (T(i) - T(i-1)) / (r(i) - r(i-1)) ) / (r(i) * dr)
            else
                T_new(i) = T(i) + dt * k_steel * ( &
                    (T(i+1) - T(i)) / (r(i+1) - r(i)) - &
                    (T(i) - T(i-1)) / (r(i) - r(i-1)) ) / (r(i) * dr)
            end if
        end do

        ! Boundary conditions
        T_new(1) = T(1) + dt * h_air * (T_air_initial - T(1)) / (r_inner * dr)
        T_new(n_r) = T(n_r) + dt * h_air * (T_ambient - T(n_r)) / (r_outer * dr)

        ! Update temperature array
        T = T_new
    end do

    ! Output final temperature distribution
    open(unit=1, file='temperature_distribution.txt', status='unknown')
    do i = 1, n_r
        write(1,*) r(i), T(i)
    end do
    close(1)

end program heat_transfer_cylinder