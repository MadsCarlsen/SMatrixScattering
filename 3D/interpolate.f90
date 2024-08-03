module interpolate
use utils 
implicit none 
private 

integer, parameter :: wp = dp 

public cubic_spline, eval_cubic_spline

contains 

! Solves the tridiaonal matrix system Ax = r using the Thomas algorithm
subroutine tridiag_elimination(a, b, c, u, r, N)
    integer, intent(in) :: N  ! Size of the matrix 
    complex(wp), intent(in) :: a(N), b(N), c(N), r(N)  ! Lower, main, upper diagonal and right hand side 
    complex(wp), intent(out) :: u(N)  ! Solution vector
    complex(wp) :: beta, gam(N)
    integer :: i

    ! Initial setup 
    beta = b(1)
    u(1) = r(1) / beta 

    ! Forward substitution to clear the lower diagonal 
    do i=2, N 
        gam(i-1) = c(i-1) / beta 
        beta = b(i) - a(i) * gam(i-1)
        ! Here we should in principle check if beta is 0...
        u(i) = (r(i) - a(i)*u(i-1)) / beta 
    end do 

    ! Now backwards substitution to clear the upper diagonal. Last element is already done. 
    do i=N-1, 1, -1
        u(i) = u(i) - gam(i)*u(i+1)
    end do 
end subroutine tridiag_elimination

! Subroutine to obtain the second derivatives of the a spline 
subroutine cubic_spline(x, y, spline_res, N)
    integer, intent(in) :: N
    complex(wp), intent(in) :: y(N)
    real(wp), intent(in) :: x(N)
    complex(wp), intent(out) :: spline_res(N)
    complex(wp) :: a(N-2), b(N-2), c(N-2), d(N-2)  ! Minus two, since endpoints are fixed 
    integer :: i, j

    ! Endpoints set to make natural cubic spline 
    spline_res(1) = 0._wp 
    spline_res(N) = 0._wp

    ! Prepare tridiag matrix arrays 
    b(1) = (x(3)-x(1)) / 3._wp 
    d(1) = (y(3)-y(2)) / (x(3)-x(2)) - (y(2)-y(1)) / (x(2)-x(1))
    do i=2, N-2 
        j = i+1 
        c(i-1) = (x(i+1)-x(i)) / 6._wp
        b(i) = (x(j+1)-x(j-1)) / 3._wp
        a(i) = (x(j+1)-x(j)) / 6._wp
        d(i) = (y(j+1)-y(j)) / (x(j+1)-x(j)) - (y(j)-y(j-1)) / (x(j)-x(j-1))
    end do 

    ! Solve the tridiagonal matrix system to find the second derivatives
    call tridiag_elimination(a, b, c, spline_res(2:N-1), d, N-2)
end subroutine cubic_spline 

! Subroutine to evaluate a cubic spline in random x point 
function eval_cubic_spline(x, x_list, y_list, y_deriv_list, N) result(y)
    integer, intent(in) :: N 
    real(wp), intent(in) :: x
    real(wp), intent(in) :: x_list(N)
    complex(wp), intent(in) :: y_list(N), y_deriv_list(N)
    real(wp) :: h, A, B, C, D
    complex(wp) :: y 
    integer :: k_low, k_high, k

    k_low = 1 
    k_high = N 
    ! First determine the indices that bracket x using binary search
    do while (k_high - k_low > 1) 
        k = (k_low + k_high) / 2  
        if (x_list(k) > x) then 
            k_high = k
        else 
            k_low = k 
        end if 
    end do 

    ! Now k_low and k_high bracket the value x. Evaluate the spline polynomial in this interval 
    h = x_list(k_high) - x_list(k_low)
    A = (x_list(k_high) - x) / h
    B = (x - x_list(k_low)) / h
    C = 1._wp / 6._wp * (A**3 - A) * h**2
    D = 1._wp / 6._wp * (B**3 - B) * h**2
    y = A*y_list(k_low) + B*y_list(k_high) + C*y_deriv_list(k_low) + D*y_deriv_list(k_high)
end function eval_cubic_spline 

end module interpolate