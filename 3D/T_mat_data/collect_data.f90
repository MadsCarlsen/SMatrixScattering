program collect_data 
implicit none 

integer, parameter :: wp = kind(1d0)
character(len=50) :: target_folder
character(len=100) :: file_path_real, file_path_imag 
integer :: ierr, i 

! Define parameters needed to load the config file 
integer :: max_m   ! Maximum magnetic quantum number to calculate T for 
integer :: N_quad   ! Nr of quad points to use in T-matrix calculation
logical :: elastic   ! Toggles whether or not we perform elastic or multi channel calculation 
logical :: k_single  ! Toggels wheter to perform calculation for a single k or a range of k's
real(wp) :: k   ! k-vector to calculate T-matrix for (if k_single = .true.)
real(wp) :: k_min   ! Lower range of k values to calculate for (if k_single = .false.)
real(wp) :: k_max   ! Upper range of k values to calculate for (if k_single = .false.)
integer :: N_k  ! Nr. of k values to calculate for (if k_single = .false.)
logical :: hydrogen_bool 
integer :: Nr   ! Nr of radial points to use in calculation of potential matrix elements
real(wp) :: r_max   ! Maximum radius in potential matrix element calculations 

! Define matrix to contain all the T-matrix elements 
real(wp), allocatable :: T_real_list(:), T_imag_list(:)
real(wp), allocatable :: T_real_mat(:,:), T_imag_mat(:,:)

! Define namelists 
namelist /T_MAT_SETTINGS/ max_m, N_quad, elastic, k_single, k, k_min, k_max, N_k, r_max, hydrogen_bool, Nr

! Get target folder from command line argument 
if (command_argument_count() /= 1) then
    write(*,*) 'Provide config file name!'
else
    ! Get the command line argument
    call get_command_argument(1, target_folder, ierr)
end if 

! Load settings from file - only interested in Nk and max_m here! 
open(file='settings_3D.nml', unit=1)
read(nml=T_MAT_SETTINGS, unit=1)
close(1)

! Allocate the different arrays 
allocate(T_real_mat(max_m+1, N_k))
allocate(T_imag_mat(max_m+1, N_k))
allocate(T_real_list(max_m+1))
allocate(T_imag_list(max_m+1))

! Now load each file one after one and save to the matrices 
do i=1, 800
    ! First get the name of the file 
    write(file_path_real, '(A,i0,A)') 'data/', i, '_real.dat' 
    write(file_path_imag, '(A,i0,A)') 'data/', i, '_imag.dat' 
    
    ! Load and add to arrays 
    open(1, file=file_path_real, form='unformatted')
    open(2, file=file_path_imag, form='unformatted')
    read(1) T_real_list
    read(2) T_imag_list 
    T_real_mat(:,i) = T_real_list 
    T_imag_mat(:,i) = T_imag_list
    close(1)
    close(2)
    !write(*,*) T_real_list
end do 

! Now save the full arrays 
open(file='T00_real.dat', unit=1, form='unformatted')
open(file='T00_imag.dat', unit=2, form='unformatted')
write(1) T_real_mat
write(2) T_imag_mat
close(1)
close(2)

end program collect_data 
