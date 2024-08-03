program main 
use GaussLegendre
use utils 
use T_matrix
implicit none 

! --- LOAD PARAMETERS FOR SIMULATION --- !

! Define working precision 
integer, parameter :: wp = dp 

! Parameters for simulation (here standard values)
integer :: max_m = 40  ! Maximum magnetic quantum number to calculate T for 
integer :: N_quad = 100  ! Nr of quad points to use in T-matrix calculation
logical :: elastic = .true.  ! Toggles whether or not we perform elastic or multi channel calculation 
logical :: k_single = .true.  ! Toggels wheter to perform calculation for a single k or a range of k's
real(wp) :: k = 0.1_wp  ! k-vector to calculate T-matrix for (if k_single = .true.)
real(wp) :: k_min = 0.1_wp  ! Lower range of k values to calculate for (if k_single = .false.)
real(wp) :: k_max = 2._wp  ! Upper range of k values to calculate for (if k_single = .false.)
integer :: N_k = 50  ! Nr. of k values to calculate for (if k_single = .false.)
logical :: hydrogen_bool = .false.
integer :: Nr = 3000  ! Nr of radial points to use in calculation of potential matrix elements
real(wp) :: r_max = 20._wp  ! Maximum radius in potential matrix element calculations 


! Variables needed in simulation 
complex(wp), allocatable :: T_res(:,:), T_list(:), T_res_channel2(:,:)
complex(wp), allocatable :: T_list_channel2(:,:)   ! Used to store the channel 2 matrix elements if non-elastic scattering 
real(wp), allocatable :: k_quad(:), w_quad(:)
real(wp), allocatable :: k_list(:)
integer :: i 

! Array calculation variables 
integer :: ierr, array_cal_index
character(len=10) :: command_line_input 
character(len=100) :: save_file1 
character(len=100) :: save_file2

! Load parameters from the settings file 
namelist /T_MAT_SETTINGS/ max_m, N_quad, elastic, k_single, k, k_min, k_max, N_k, r_max, hydrogen_bool, Nr
open(file='settings_3D.nml', unit=1)
read(nml=T_MAT_SETTINGS, unit=1)
close(1)

! Allocated based on settings 
if (k_single) then 
    N_k = 1 
    allocate(k_list(N_k))
    k_list(1) = k
else 
    allocate(k_list(N_k))
    call linspace(k_list, k_min, k_max, N_k)
end if 

allocate(T_res(max_m+1, N_k))
allocate(k_quad(N_quad))
allocate(w_quad(N_quad))

! Get index for array calculations 
if (command_argument_count() /= 1) then
    write(*,*) 'Provide config file name!'
else
    ! Get the command line argument
    call get_command_argument(1, command_line_input, ierr)
    read(command_line_input,*,iostat=ierr) array_cal_index
end if 

! Make check to see if array-calculation index is in range???

! Make T-matrix code load the potential data 
call load_data(r_max, Nr, hydrogen_bool)


! --- SIMULATION --- !

! Calculate Gaussian quadrature points and weights and rescale them to [0, inf[
call gauss_leg(N_quad, k_quad, w_quad, 0._wp, pi/2._wp)
w_quad = w_quad / (cos(k_quad)**2)
k_quad = tan(k_quad)

! Calculate T-matrix element for the specific k value given by the array-calculation index 
write(*,*) 'Performing calculation for k = ', k_list(array_cal_index)
allocate(T_list(max_m+1))
call calculate_several_m_elastic(T_list, max_m, k_list(array_cal_index), N_quad, k_quad, w_quad, hydrogen_bool)

! Save result to file with index name 
write(save_file1, "(A,A,A)") 'data/', trim(command_line_input), '_real.dat'
write(save_file2, "(A,A,A)") 'data/', trim(command_line_input), '_imag.dat' 

open(file=save_file1, unit=1, form='unformatted')
open(file=save_file2, unit=2, form='unformatted')
write(1) real(T_list)
write(2) aimag(T_list)
close(1)
close(2)

end program main 