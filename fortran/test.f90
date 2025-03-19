PROGRAM HeatConduction1D
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = 10       ! Number of nodes
    REAL, PARAMETER :: L = 1.0         ! Length of the rod (m)
    REAL, PARAMETER :: k = 200.0       ! Thermal conductivity (W/m·K)
    REAL, PARAMETER :: A = 0.01        ! Cross-sectional area (m^2)
    REAL, PARAMETER :: q_gen = 1000.0  ! Heat generation per unit volume (W/m^3)
    REAL, PARAMETER :: T_left = 100.0  ! Left boundary temperature (°C)
    REAL, PARAMETER :: T_right = 50.0  ! Right boundary temperature (°C)

    ! Variables
    REAL :: dx, rhs(n), T(n)
    REAL, DIMENSION(n, n) :: K_matrix
    INTEGER :: i, j

    ! Initialize
    dx = L / (n - 1)
    K_matrix = 0.0
    rhs = 0.0
    T = 0.0

    ! Assemble the stiffness matrix and RHS vector
    DO i = 1, n
        IF (i == 1) THEN
            K_matrix(i, i) = 1.0
            rhs(i) = T_left
        ELSE IF (i == n) THEN
            K_matrix(i, i) = 1.0
            rhs(i) = T_right
        ELSE
            K_matrix(i, i - 1) = -k * A / dx
            K_matrix(i, i) = 2.0 * k * A / dx
            K_matrix(i, i + 1) = -k * A / dx
            rhs(i) = q_gen * A * dx
        END IF
    END DO

    ! Solve the system of equations K_matrix * T = rhs
    CALL SolveLinearSystem(K_matrix, rhs, T, n)

    ! Output results
    PRINT *, "Node", "Temperature (°C)"
    DO i = 1, n
        PRINT *, i, T(i)
    END DO

CONTAINS

    SUBROUTINE SolveLinearSystem(A, b, x, n)
        REAL, INTENT(IN) :: A(n, n)
        REAL, INTENT(IN) :: b(n)
        REAL, INTENT(OUT) :: x(n)
        INTEGER, INTENT(IN) :: n
        REAL :: tempA(n, n), tempb(n)
        INTEGER :: i, j, k

        tempA = A
        tempb = b

        ! Gaussian elimination
        DO k = 1, n - 1
            DO i = k + 1, n
                tempb(i) = tempb(i) - tempA(i, k) / tempA(k, k) * tempb(k)
                DO j = k + 1, n
                    tempA(i, j) = tempA(i, j) - tempA(i, k) / tempA(k, k) * tempA(k, j)
                END DO
                tempA(i, k) = 0.0
            END DO
        END DO

        ! Back substitution
        x(n) = tempb(n) / tempA(n, n)
        DO i = n - 1, 1, -1
            x(i) = (tempb(i) - SUM(tempA(i, i + 1:n) * x(i + 1:n))) / tempA(i, i)
        END DO
    END SUBROUTINE SolveLinearSystem

END PROGRAM HeatConduction1D