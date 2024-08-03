module utils 
implicit none 
private

! Different kinds of floating point numbers
integer, parameter :: sp = kind(1.0)
integer, parameter :: dp = kind(1.0d0) 

! Some constants 
real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp  ! Do I need the dp here? 

! Define the working precision for this module! 
integer, parameter :: wp = dp  

! Public constants and kinds 
public sp, dp, pi 

! Interface for overloading 
interface trapz 
    module procedure trapz_real, trapz_complex
end interface trapz 

! Public routines
public linspace, trapz, savetxt


contains 
    ! Just a standard linspace subroutine, nothing fancy 
    subroutine linspace(arr, start, end, N)
        integer, intent(in) :: N 
        real(wp), intent(in) :: start, end 
        real(wp), intent(out) :: arr(N)
        real(wp) :: step 
        integer :: i 

        step = (end - start) / (N-1)

        do i = 1, N
            arr(i) = start + (i-1)*step 
        end do 
    end subroutine linspace

    ! A simple trapezoidal rule integration function 
    function trapz_real(func_vec, step_size)
        real(wp), intent(in) :: step_size  ! The step size of the grid used for function evaluation 
        real(wp), intent(in) :: func_vec(:)  ! The function to integrate evaluated at evenly spaced points x
        real(wp) :: trapz_real
        integer :: N  ! To hold the size of the array
        
        N = size(func_vec)
        trapz_real = 0.5_dp * step_size * (func_vec(1) + func_vec(N)) + step_size * sum(func_vec(2:N-1))
    end function trapz_real

    ! A simple trapezoidal rule integration function 
    function trapz_complex(func_vec, step_size)
        real(wp), intent(in) :: step_size  ! The step size of the grid used for function evaluation 
        complex(wp), intent(in) :: func_vec(:)  ! The function to integrate evaluated at evenly spaced points x
        complex(wp) :: trapz_complex
        integer :: N  ! To hold the size of the array
        
        N = size(func_vec)
        trapz_complex = 0.5_dp * step_size * (func_vec(1) + func_vec(N)) + step_size * sum(func_vec(2:N-1))
    end function trapz_complex

    ! Simple subroutine to save a 1D real array to a data file (txt)
    subroutine savetxt(file_name, array, include_Nr_of_elements, index)
        character(*), intent(in) :: file_name
        real(wp), intent(in) :: array(:)
        logical, optional, intent(in) :: include_Nr_of_elements 
        integer, optional, intent(in) :: index
        logical :: include_num
        integer :: i, index_to_use

        ! Set default value if not provided
        if (present(include_Nr_of_elements)) then 
            include_num = include_Nr_of_elements
        else 
            include_num = .false. 
        end if

        if (present(index)) then 
            index_to_use = index + 6  ! To avoid 5 and 6?
        else 
            index_to_use = 1
        end if

        ! Save the array to file 
        open(index_to_use, file=file_name)
        if (include_num) then 
            write(index_to_use,*) size(array)
        end if
        do i=1, size(array)
            write(index_to_use,*) array(i)
        end do 
    end subroutine savetxt
end module utils 