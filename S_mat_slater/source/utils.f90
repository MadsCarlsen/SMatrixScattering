module utils 
implicit none 
private

! Different kinds of floating point numbers
integer, parameter :: sp = kind(1.0)
integer, parameter :: dp = kind(1.0d0) 

! Some constants 
real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp   

! Define the working precision for this module! 
integer, parameter :: wp = dp  

! Public constants and kinds 
public sp, dp, pi 

! Interface for overloading 
interface trapz 
    module procedure trapz_real, trapz_complex
end interface trapz 

! Public routines
public linspace, trapz, savetxt, init_slater, save_array

! Define Slater state type 
type, public :: SlaterState
    integer :: Ns  ! Number of basis functions in the state 
    real(wp), allocatable :: c(:) ! Expansions coefficients (Including the normalization constant)
    real(wp), allocatable :: alpha(:)  ! Exponents 
    integer, allocatable :: n(:) ! Radial exponents 
    integer :: l  ! Angular momentum and magnetic quantum number
    integer :: m 
    real(wp) :: E  ! Energy of the state 
    complex(wp) :: init_coeff 
    integer :: index
end type 


contains 
    subroutine init_slater(state, Ns, c_arr, alpha_arr, n_arr, l, m, E)
        type(SlaterState), intent(out) :: state
        integer, intent(in) :: Ns
        integer, intent(in) :: n_arr(Ns)
        real(wp), intent(in) :: alpha_arr(Ns), c_arr(Ns)
        integer, intent(in) :: l, m
        real(wp), intent(in) :: E
        integer :: i 
         
        ! Allocate the arrays
        allocate(state%c(Ns))
        allocate(state%alpha(Ns))
        allocate(state%n(Ns))

        ! Assign the values 
        state%Ns = Ns
        state%c = c_arr
        state%alpha = alpha_arr
        state%n = n_arr
        state%l = l
        state%m = m
        state%E = E

        ! Calculate the norms 
        !do i=1, Ns
        !    norm_i = sqrt((2._wp*alpha_arr(i))**(2*n_arr(i)+1) / gamma(2._wp*n_arr(i)+1._wp))
        !    state%c(i) = c_arr(i) * norm_i
        !end do
    end subroutine init_slater

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

    subroutine save_array(file_name, array, unit)
        character(*), intent(in) :: file_name
        real(wp), intent(in) :: array(:,:)
        integer, optional, intent(in) :: unit 
        integer :: index 

        if (present(unit)) then 
            index = unit + 6  ! To avoid 5 and 6?
        else 
            index = 1
        end if

        open(file=file_name, unit=index, form='unformatted')
        write(1) array
        close(1)
    end subroutine save_array
end module utils 