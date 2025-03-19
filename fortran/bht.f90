program barrel_heat_transfer
    use hdf5
    implicit none

    ! 물리 상수 및 초기 조건
    real(8), parameter :: pi = 3.141592653589793
    real(8), parameter :: barrel_diameter = 0.155       ! 포신 지름 (m)
    real(8), parameter :: barrel_thickness = 0.0775     ! 포신 두께 (m)
    real(8), parameter :: coating_thickness = 0.001     ! 크롬 코팅 두께 (m)
    real(8), parameter :: initial_air_temp = 2500.0     ! 내부 공기 초기 온도 (°C)
    real(8), parameter :: ambient_temp = 25.0           ! 외부 공기 온도 (°C)
    real(8), parameter :: sim_time = 60.0               ! 시뮬레이션 시간 (s)
    real(8), parameter :: dt = 0.01                     ! 시간 간격 (s)
    real(8), parameter :: dx = 0.0001                   ! 공간 격자 크기 (m)

    ! 열전도율 및 밀도, 비열 (철과 크롬의 값)
    real(8), parameter :: k_steel = 50.2                ! 철의 열전도율 (W/m·K)
    real(8), parameter :: rho_steel = 7850.0            ! 철의 밀도 (kg/m³)
    real(8), parameter :: cp_steel = 470.0              ! 철의 비열 (J/kg·K)
    real(8), parameter :: k_chrome = 93.7               ! 크롬의 열전도율 (W/m·K)
    real(8), parameter :: rho_chrome = 7190.0           ! 크롬의 밀도 (kg/m³)
    real(8), parameter :: cp_chrome = 460.0             ! 크롬의 비열 (J/kg·K)

    ! 대류 열전달 계수
    real(8), parameter :: h_inner = 100.0               ! 내부 대류 열전달 계수 (W/m²·K)
    real(8), parameter :: h_outer = 10.0                ! 외부 대류 열전달 계수 (W/m²·K)

    ! 격자 설정
    integer, parameter :: n_nodes = int((barrel_thickness + coating_thickness) / dx) + 1
    real(8) :: temp(n_nodes), temp_new(n_nodes)
    real(8) :: r(n_nodes)
    integer :: i
    real(8) :: time
    real(8) :: start_time, end_time

    ! HDF5 관련 변수
    integer(HID_T) :: hdf5_file, hdf5_space, hdf5_dataset
    integer :: hdf5_status
    integer(8) :: dims(1)

    ! 파일 존재 여부 확인 및 삭제
    logical :: file_exists  ! Declare as logical
    inquire(file="temperature_distribution.h5", exist=file_exists)
    if (file_exists) then
        call execute_command_line("rm -f temperature_distribution.h5")
    end if

    ! 초기화
    do i = 1, n_nodes
        r(i) = (i - 1) * dx
        if (r(i) <= coating_thickness) then
            temp(i) = initial_air_temp
        else
            temp(i) = ambient_temp
        end if
        print *, "r(", i, ") = ", r(i), ", temp(", i, ") = ", temp(i)
    end do

    if (any(r < 0.0)) then
        print *, "Error: Negative radius values detected."
        stop
    end if

    ! 시작 시간 기록
    call cpu_time(start_time)

    ! 시간 루프
    time = 0.0
    do while (time < sim_time)
        ! 열 전달 계산
        do i = 2, n_nodes - 1
            if (r(i) <= coating_thickness) then
                temp_new(i) = temp(i) + dt * k_chrome / (rho_chrome * cp_chrome) * &
                                            ((temp(i+1) - 2.0 * temp(i) + temp(i-1)) / dx**2)
            else
                temp_new(i) = temp(i) + dt * k_steel / (rho_steel * cp_steel) * &
                                            ((temp(i+1) - 2.0 * temp(i) + temp(i-1)) / dx**2)
            end if
        end do

        ! 경계 조건
        temp_new(1) = temp(1) + dt * h_inner / (rho_chrome * cp_chrome) * &
                                    (initial_air_temp - temp(1))
        temp_new(n_nodes) = temp(n_nodes) + dt * h_outer / (rho_steel * cp_steel) * &
                                                (ambient_temp - temp(n_nodes))

        ! 온도 갱신
        temp = temp_new
        time = time + dt
    end do

    ! 종료 시간 기록
    call cpu_time(end_time)

    ! HDF5 파일 생성
    call h5fcreate_f("temperature_distribution.h5", H5F_ACC_TRUNC_F, hdf5_file, hdf5_status)
    if (hdf5_status /= 0) then
        print *, "Error creating HDF5 file."
        stop
    end if

    ! 데이터 공간 생성
    dims(1) = n_nodes
    call h5screate_simple_f(1, dims, hdf5_space, hdf5_status)
    if (hdf5_status /= 0) then
        print *, "Error creating HDF5 dataspace."
        stop
    end if

    ! 데이터셋 생성 및 쓰기
    call h5dcreate_f(hdf5_file, "radius", H5T_NATIVE_DOUBLE, hdf5_space, hdf5_dataset, hdf5_status)
    if (hdf5_status /= 0) then
        print *, "Error creating dataset 'radius'."
        stop
    end if

    call h5dwrite_f(hdf5_dataset, H5T_NATIVE_DOUBLE, r, dims, hdf5_status)
    if (hdf5_status /= 0) then
        print *, "Error writing dataset 'radius'."
        stop
    end if

    call h5dclose_f(hdf5_dataset, hdf5_status)

    call h5dcreate_f(hdf5_file, "temperature", H5T_NATIVE_DOUBLE, hdf5_space, hdf5_dataset, hdf5_status)
    if (hdf5_status /= 0) then
        print *, "Error creating dataset 'temperature'."
        stop
    end if

    call h5dwrite_f(hdf5_dataset, H5T_NATIVE_DOUBLE, temp, dims, hdf5_status)
    if (hdf5_status /= 0) then
        print *, "Error writing dataset 'temperature'."
        stop
    end if

    call h5dclose_f(hdf5_dataset, hdf5_status)

    ! 데이터 공간 및 파일 닫기
    call h5sclose_f(hdf5_space, hdf5_status)
    call h5fclose_f(hdf5_file, hdf5_status)

    print *, "Simulation complete. Results have been saved to temperature_distribution.h5."
    print *, "Simulation runtime: ", end_time - start_time, " seconds."

end program barrel_heat_transfer