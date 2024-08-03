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

! Make T-matrix code load the potential data 
call load_data(r_max, Nr, hydrogen_bool)


! --- SIMULATION --- !

! Calculate Gaussian quadrature points and weights and rescale them to [0, inf[
call gauss_leg(N_quad, k_quad, w_quad, 0._wp, pi/2._wp)
w_quad = w_quad / (cos(k_quad)**2)
k_quad = tan(k_quad)

!write(*,*) w_quad

! Call T-matrix code based on settings
if (elastic) then 
    ! Elastic call (This is not done in parallel, but the different m's are...)
    write(*,*) 'Performing elastic scattering calculation...'
    allocate(T_list(max_m+1))
    do i=1, N_k
        write(*,*) i, '/', N_k
        call calculate_several_m_elastic(T_list, max_m, k_list(i), N_quad, k_quad, w_quad, hydrogen_bool)
        T_res(:,i) = T_list
    end do 

    ! Now save the data 
    open(file='T00_real.dat', unit=1, form='unformatted')
    open(file='T00_imag.dat', unit=2, form='unformatted')
    write(1) real(T_res)
    write(2) aimag(T_res)
    
    open(file='k_list.dat', unit=3)
    do i=1, N_k
        write(3,*) k_list(i)
    end do
    close(1)
    close(2)
    close(3)
else  
    ! Non-elastic (two channel) call 
    write(*,*) 'Performing two channel scattering calculation...'
    allocate(T_list_channel2(max_m+1, 2))
    allocate(T_res_channel2(max_m+1, N_k))
    do i=1, N_k
        call calculate_several_m_two_channel(T_list_channel2, max_m, k_list(i), N_quad, k_quad, w_quad)
        T_res(:,i) = T_list_channel2(:,1)
        T_res_channel2(:,i) = T_list_channel2(:,2) 
    end do 

    ! Now save the data 
    open(file='data/T00_real.dat', unit=1, form='unformatted')
    open(file='data/T00_imag.dat', unit=2, form='unformatted')
    open(file='data/T01_real.dat', unit=3, form='unformatted')
    open(file='data/T01_imag.dat', unit=4, form='unformatted')
    write(1) real(T_res)
    write(2) aimag(T_res)
    write(3) real(T_res_channel2)
    write(4) aimag(T_res_channel2)
    
    open(file='data/k_list.dat', unit=7)
    do i=1, N_k
        write(7,*) k_list(i)
    end do
    close(1)
    close(2)
    close(3)
    close(4)
    close(7)

end if 


end program main 