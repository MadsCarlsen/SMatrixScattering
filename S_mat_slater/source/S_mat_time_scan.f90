program main 
use utils 
use Wigner3j
use born_elements
use GaussLegendre
implicit none

integer, parameter :: wp = dp 

! Struct to contain the WP settings and type to make it easy to pass around / implement new WP types
type :: WP_parameters 
    integer :: wp_type = 0  ! Toggle for different WP types implemented below 
    real(wp) :: k0z, sigma_z, sigma_perp, z0, tz, t_perp, t_bunch, b, norm_factor, g_factor, omega
end type 

! Parameters determining what to calculate 
logical :: cal_tot_signal
! Parameters for calculating of the S-matrix and general integration    
integer :: Nk
real(wp) :: k_low, k_high 
integer :: N_theta_int, N_theta_final, N_phi_int, N_phi_final
real(wp) :: theta_int_min, theta_int_max, theta_final_min, theta_final_max
real(wp) :: phi_int_min = 0._wp, phi_int_max = 2._wp*pi
real(wp) :: phi_final_min, phi_final_max 

! Parameters for the target 
complex(wp), allocatable :: state_coeffs(:) 
type(SlaterState), allocatable :: state_array(:)
integer :: N_states
real(wp) :: ta0  ! Initial time to add to the atomic system cohrent superposition 

! Parameters for the WP 
real(wp) :: k0z, sigma_time, sigma_perp_factor, omega, g_factor, b  ! These are loaded 
integer :: WP_type, max_N  ! This is loaded 
real(wp) :: t_bunch, z0, tz, t_perp, sigma_space_z, sigma_space_perp, sigma_z, sigma_perp 
type(WP_parameters) :: Wp_p

! Arrays for integration and results 
real(wp), allocatable :: theta_quad(:), theta_final(:), k_list(:), phi_quad(:), phi_final(:)
real(wp), allocatable :: theta_w(:), phi_w(:)
integer :: i,j 
real(wp), allocatable :: res(:,:)
real(wp) :: coeff_sum
 
! Command line input to control the final channel (for array calculation)
integer :: ierr, array_cal_index
character(len=10) :: command_line_input 
character(len=100) :: save_file

! Constants 
real(wp), parameter :: sqrt2 = sqrt(2._wp)

! TEST
integer :: Nt
real(wp), allocatable :: delay_times(:)

! Define namelists 
namelist /CALCULATE_SETTINGS/ cal_tot_signal
namelist /S_MAT_SETTINGS/ Nk, N_theta_int, N_theta_final, theta_int_min &
                        , theta_int_max, theta_final_min, theta_final_max, k_low, k_high, &
                        phi_final_min, phi_final_max, N_phi_int, N_phi_final
namelist /TARGET_SETTINGS/ b, ta0
namelist /WAVEPACKET_SETTINGS/ WP_type, k0z, sigma_time, sigma_perp, omega, g_factor, max_N 

! Load data 
open(file='settings.nml', unit=1)
read(nml=CALCULATE_SETTINGS, unit=1)
close(1)
open(file='settings.nml', unit=1)
read(nml=S_MAT_SETTINGS, unit=1)
close(1)
open(file='settings.nml', unit=1)
read(nml=TARGET_SETTINGS, unit=1)
close(1)
open(file='settings.nml', unit=1)
read(nml=WAVEPACKET_SETTINGS, unit=1)
close(1)

! Get command line input 
if (command_argument_count() /= 1) then
    write(*,*) 'Provide final channel integer as command line input!'
    stop 
else
    ! Get the command line argument
    call get_command_argument(1, command_line_input, ierr)
    read(command_line_input,*,iostat=ierr) array_cal_index
end if 

! Setup arrays for integration
allocate(res(N_theta_final, N_phi_final))
allocate(k_list(Nk))
allocate(theta_final(N_theta_final))
allocate(phi_final(N_phi_final))
allocate(theta_quad(N_theta_int))
allocate(phi_quad(N_phi_int))
allocate(theta_w(N_theta_int))
allocate(phi_w(N_phi_int))
call linspace(k_list, k_low, k_high, Nk)
if (N_theta_final == 1) then 
    theta_final = theta_final_min
else   
    call linspace(theta_final, theta_final_min, theta_final_max, N_theta_final)
end if 
if (N_phi_final == 1) then 
    phi_final = phi_final_min  
else
    call linspace(phi_final, pi*phi_final_min, pi*phi_final_max, N_phi_final)
end if 

! Calculate theta and phi integration quadrature points  
!call linspace(theta_quad, theta_int_min, theta_int_max, N_theta_int)
!call linspace(phi_quad, phi_int_min, phi_int_max, N_phi_int)
call gauss_leg(N_theta_int, theta_quad, theta_w, theta_int_min, theta_int_max)
call gauss_leg(N_phi_int, phi_quad, phi_w, phi_int_min, phi_int_max)



! Calculate parameters 
sigma_space_z = sigma_time * k0z 
sigma_z = 1._wp / (2._wp * sigma_space_z)   
t_bunch = k0z**2 / (g_factor * omega**2)

! Update wavepacket struct with parameters 
Wp_p = WP_parameters(WP_type, k0z, sigma_z, sigma_perp, z0, tz, t_perp, t_bunch, b, 1._wp, g_factor, omega)


! SETUP TARGET AND INITIAL STATE
N_states = 5  ! 1s, 2s, 2p-1, 2p0, 2p1
allocate(state_array(N_states))
allocate(state_coeffs(N_states))

! Initial state expansion coeffs (1s, 2s, 2p-1, 2p0, 2p1) 
state_coeffs = [1._wp/sqrt(2._wp), 0._wp, 1._wp/sqrt(4._wp), 0._wp, -1._wp/sqrt(4._wp)]

! Setup states in terms of slater orbitals 
call init_slater(state_array(1), 1, [2._wp], [1._wp], [1], 0, 0, -0.5_wp)
call init_slater(state_array(2), 2, [1._wp/sqrt(2._wp), -1._wp/(2._wp*sqrt(2._wp))], [0.5_wp, 0.5_wp], [1, 2], 0, 0, -0.125_wp)
call init_slater(state_array(3), 1, [1._wp/(2._wp*sqrt(6._wp))], [0.5_wp], [2], 1, -1, -0.125_wp)
state_array(4) = state_array(3)
state_array(5) = state_array(3)
state_array(4)%m = 0 
state_array(5)%m = 1

! Normalize the initial state in terms of expansion coefficients 
coeff_sum = 0._wp
do i=1, size(state_coeffs)
    coeff_sum = coeff_sum + abs(state_coeffs(i))**2 
end do 
if (coeff_sum /= 1._wp) then 
    write(*,*) 'Normalizing the initial state...'
    state_coeffs = state_coeffs / sqrt(coeff_sum)
end if 

! Now set the values in the wp slater_state structs 
do i=1, N_states 
    state_array(i)%init_coeff = state_coeffs(i)
    state_array(i)%index = i
end do 


! SIMULATIONS 
tz = 0._wp 
t_perp = -t_bunch  
z0 = k0z * t_bunch 
!ta0 = 0._wp  ! Time propagation to add to the atomic system coherent superposition 

! Set the WP settings 
WP_p%t_perp = t_perp 
WP_p%z0 = z0 
WP_p%tz = tz 
call normalize_WP(WP_p)  ! Normalize WP on the grid (neeeded if E-gauss, since not analytically normalized)
write(*,*) 'WP settings: ', Wp_p%k0z, Wp_p%sigma_z, Wp_p%sigma_perp, Wp_p%z0, Wp_p%tz, Wp_p%t_perp, Wp_p%t_bunch, Wp_p%b
write(*,*) 'Calculations done for WP type: ', WP_type

! Delay times 
Nt = 4
allocate(delay_times(Nt))
call linspace(delay_times, 0._wp, 2._wp*pi/(-0.125_wp + 0.5_wp), Nt)  ! Scan over one period 

! Now get the angular distributions!
i = 3 
write(*,*) 'Calculating channel ', i, '/'
write(*,*) 'Performing calculations for time ', array_cal_index 

! Add the atomic system initial time 
do j=1, N_states
    state_array(j)%init_coeff = state_coeffs(j) * exp(cmplx(0._wp, -state_array(j)%E*delay_times(array_cal_index), wp))
end do 
res = angular_dist_channel(k_list, state_array(i), state_array, Wp_p, theta_final, phi_final) 

write(save_file, "(A,I0,A,I0,A)") 'data/angular_dist_channel_', i, '_time_', array_cal_index, '.txt'
!call save_array(save_file, res)
open(2, file=save_file)
do i = 1, size(theta_final)
    write(2,*) res(i,1)
end do
 
! Also save the dimensions of the arrays in a .txt file (if this is channel 1 calculation)
if (array_cal_index == 1) then
    open(1, file='data/dim.txt')
    write(1,*) N_theta_final
    write(1,*) N_phi_final 
    write(1,*) theta_final_min
    write(1,*) theta_final_max
    write(1,*) phi_final_min
    write(1,*) phi_final_max
    write(1,*) Nt 
    write(1,*) delay_times(1)
    write(1,*) delay_times(Nt)
    close(1)
end if 


contains 

! Function to calculate the angular distribution, for a given final channel 
function angular_dist_channel(k_list, state_f, state_arr, WP_p, theta_final, phi_final) result(res)
    real(wp), intent(in) :: k_list(:), theta_final(:), phi_final(:)
    type(SlaterState), intent(in) :: state_f, state_arr(:) 
    type(WP_parameters), intent(in) :: WP_p
    real(wp) :: res(size(theta_final), size(phi_final)) 
    complex(wp) :: S_list(size(k_list))
    complex(wp) :: psi_f
    integer :: i, j, m, state_index
    logical :: elastic 
    real(wp), allocatable :: wigner_0(:), wigner_m(:)
    type(SlaterState) :: state_i
    real(wp) :: li, lf, mi, mf 
     
    !!$OMP PARALLEL DO PRIVATE(i, S_list, j, m, state_index, psi_f, elastic, state_i, wigner_0, wigner_m, li, lf, mi, mf)
    do i=1, size(theta_final)
        !write(*,*) i, '/', size(theta_final)
        do j=1, size(phi_final)
            !write(*,*) j, '/', size(phi_final)
            S_list = 0._wp  ! Reset the S-matrix list 
            do m=1, size(k_list)
                psi_f = eval_WP(k_list(m), theta_final(i), WP_p)
                
                ! Loop over all the initial state channels 
                do state_index = 1, size(state_arr) 
                    if (state_index == state_f%index) then
                        elastic = .true.
                    else 
                        elastic = .false.
                    end if 
                    
                    ! Check if the initial state has any components in the channel before calculation... 
                    if (abs(state_arr(state_index)%init_coeff) /= 0._wp) then 
                        
                        state_i = state_arr(state_index)

                        ! Calculate Wigner-coeffs here? To perform checks on wether or not it is total 0 (is it ever?, else move to S-mat)
                        lf = state_f%l
                        li = state_i%l
                        mf = state_f%m
                        mi = state_i%m
                        
                        call wigner_3j_right(li, lf, mf-mi, mi, -mf, wigner_m)
                        call wigner_3j_right(li, lf, 0._wp, 0._wp, 0._wp, wigner_0)

                        S_list(m) = S_list(m) + calculate_S_element(k_list(m), theta_final(i), phi_final(j), & 
                                                state_f, state_i, WP_p, wigner_0, wigner_m, elastic) 
                                                
                        ! Deallocate the wigner arrays so they are ready for new iteration? 
                        deallocate(wigner_0)
                        deallocate(wigner_m)
                    end if 
                end do 
                
                !Add the final WP state to get the full S-matrix element 
                if (cal_tot_signal) then 
                    S_list(m) = S_list(m) + psi_f*state_f%init_coeff 
                end if 

            end do
            
            ! Now calculate the k-integral 
            res(i,j) = trapz(k_list**2 * abs(S_list)**2, k_list(2)-k_list(1)) 
        end do 
    end do 
    !!$OMP END PARALLEL DO
end function angular_dist_channel


! Function to calculate the scattering part of the S-matrix element 
function calculate_S_element(kf, theta_f, phi_f, state_f, state_i, WP_p, wigner_0, wigner_m, elastic) result (S_val)
    real(wp), intent(in) :: kf, theta_f, phi_f  ! Final momentum values 
    type(SlaterState), intent(in) :: state_f, state_i  ! Final and initial WP states
    type(WP_parameters), intent(in) :: WP_p  ! WP parameters
    logical, intent(in) :: elastic 
    complex(wp) :: S_val
    real(wp) :: ki, kx_f, ky_f, kz_f, delta_x, delta_y, delta_z, x_temp_i, delta
    real(wp), intent(in) :: wigner_0(:), wigner_m(:)
    complex(wp) :: T_elem, phi_integrand(N_phi_int), theta_integrand(N_theta_int), psi_i(N_theta_int)
    complex(wp) :: T_elem_ana
    integer :: i, j 

    ! First determine the corresponding initial momentum 
    if (elastic) then 
        ki = kf 
    else 
        ki = sqrt(kf**2 + 2._wp*(state_f%E - state_i%E))    ! Sign here 
    end if
    
    ! Evaluate the WP in the final angles (ASSUMES ONLY THETA DEPENDENCE)
    do i=1, N_theta_int 
        psi_i(i) = eval_WP(ki, theta_quad(i), WP_p)
    end do 

    ! Perform double integral over theta and phi. Do phi first, since WP's are independent of phi 
    kx_f = kf*sin(theta_f)*cos(phi_f)
    ky_f = kf*sin(theta_f)*sin(phi_f)
    kz_f = kf*cos(theta_f)

    do i=1, N_theta_int 
        delta_z = ki*cos(theta_quad(i)) - kz_f
        x_temp_i = ki*sin(theta_quad(i)) 

        do j=1, N_phi_int 
            delta_x = x_temp_i * cos(phi_quad(j)) - kx_f
            delta_y = x_temp_i * sin(phi_quad(j)) - ky_f  ! Using x_temp_i since theta part is the same for x/y
            delta = sqrt(delta_x**2 + delta_y**2 + delta_z**2)
            !write(*,*) delta, elastic 
            T_elem = calculate_born_element(delta, delta_z/delta, atan2(delta_y, delta_x), state_f, state_i, &
                                            wigner_0, wigner_m, elastic)

            !T_elem_ana = calculate_born_analytic(delta, delta_z/delta, atan2(delta_y,delta_x), state_f, state_i)
            
            !write(*,*) T_elem 
            !write(*,*) T_elem_ana 
            !write(*,*)

            phi_integrand(j) = exp(cmplx(0._wp, WP_p%b*ki*sin(theta_quad(i))*cos(phi_quad(j)),wp)) * T_elem 
        end do 
        
        ! Perform the phi integral 
        theta_integrand(i) = sin(theta_quad(i)) * psi_i(i) * sum(phi_integrand * phi_w) ! trapz(phi_integrand, phi_quad(2)-phi_quad(1))
    end do 
    
    ! Last, perform the theta integral and multiply with the state coefficient  
    S_val = sum(theta_integrand * theta_w) * cmplx(0._wp, -2._wp*pi*ki, wp) * &  ! trapz(theta_integrand, theta_quad(2)-theta_quad(1)) * 
            exp(cmplx(0._wp, -WP_p%b*kf*sin(theta_f)*cos(phi_f),wp)) * state_i%init_coeff 
                  !write(*,*) S_mat

    !if (kf > 11.9992660550459) then 
    !    write(*,*) kf, theta_f, phi_f, delta, T_elem 
    !    write(*,*) S_val 
    !    write(*,*) ki, theta_quad(1), sigma_z  
    !    write(*,*) psi_i(1)
    !    ! Kill program
    !    stop 
    !end if
end function calculate_S_element


! Wrapper function to evaluate the WP, gets the correct WP type 
function eval_WP(k, theta, WP_params) result(res)
    real(wp) :: k, theta
    type(WP_parameters), intent(in) :: WP_params
    complex(wp) :: res 

    select case (WP_params%wp_type)
    case (0)
        res = Gaussian_wp(k, theta, WP_params)
    case (1)
        res = E_Gauss_wp(k, theta, WP_params)
    case(2)
        res = PINEM_wp(k, theta, WP_params)
    end select 
end function eval_WP 


! Gaussian WP 
function Gaussian_wp(k, theta, WP_p) result(res)
    real(wp) :: k, theta 
    type(WP_parameters), intent(in) :: WP_p
    complex(wp) :: res 
    real(wp) :: kz, k_perp
    
    kz = k * cos(theta)
    k_perp = k * sin(theta) 

    !write(*,*) WP_p%tz 

    res = 1._wp/(2._wp*pi*WP_p%sigma_z**2)**(1._wp/4._wp) * &
            1._wp/sqrt(2._wp*pi*WP_p%sigma_perp**2)
    res = res * exp(-(kz-WP_p%k0z)**2/(4._wp*WP_p%sigma_z**2)) * exp(-k_perp**2/(4._wp*WP_p%sigma_perp**2))
    res = res * exp(cmplx(0._wp, -kz**2/2._wp * WP_p%tz, wp)) * &
                exp(cmplx(0._wp, -k_perp**2/2._wp * WP_p%t_perp, wp)) * & 
                exp(cmplx(0._wp, kz*WP_p%z0, wp))
end function Gaussian_wp

! E-Gaussian WP 
function E_Gauss_wp(k, theta, WP_p) result(res)
    real(wp) :: k, theta 
    type(WP_parameters), intent(in) :: WP_p
    complex(wp) :: res 
    real(wp) :: kz, k_perp

    kz = k * cos(theta)
    k_perp = k * sin(theta)

    res = WP_p%norm_factor * exp(-(k-WP_p%k0z)**2/(2._wp*WP_p%sigma_z**2) - sin(theta)**2/(2._wp*WP_p%sigma_perp**2))
    res = res * exp(cmplx(0._wp, -kz**2/2._wp * WP_p%tz, wp)) * &
                exp(cmplx(0._wp, -k_perp**2/2._wp * WP_p%t_perp, wp)) * & 
                exp(cmplx(0._wp, kz*WP_p%z0, wp))
end function E_Gauss_wp

function Gaussian_one_dim(k, k0, sigma) result(res)
    real(wp) :: k, k0, sigma
    complex(wp) :: res 
    res = 1._wp/(2._wp*pi*sigma**2)**(1._wp/4._wp) * &
            exp(-(k-k0)**2/(4._wp*sigma**2))
end function Gaussian_one_dim

function PINEM_wp(k, theta, Wp_p) result(res)
    real(wp) :: k, theta
    type(WP_parameters), intent(in) :: WP_p
    complex(wp) :: res 
    real(wp) :: delta, kz, k_perp
    real(wp) :: phi_N = 0._wp 
    integer :: N 

    delta = WP_p%omega/Wp_p%k0z
    kz = k*cos(theta)
    k_perp = k*sin(theta)
    
    ! Build the PINEM part in x direction
    res = bessel_jn(0, WP_p%g_factor) * exp(cmplx(0._wp, phi_N, wp)) * &
          Gaussian_one_dim(kz, Wp_p%k0z, Wp_p%sigma_z)
    do N=1, max_N 
        res = res + bessel_jn(N, Wp_p%g_factor) * exp(cmplx(0._wp, phi_N, wp)) * &
              Gaussian_one_dim(kz, Wp_p%k0z + N*delta, Wp_p%sigma_z)
        res = res + bessel_jn(-N, Wp_p%g_factor) * exp(cmplx(0._wp, phi_N, wp)) * &
              Gaussian_one_dim(kz, Wp_p%k0z - N*delta, Wp_p%sigma_z)
    end do 

    ! Then add the Gaussian in y direction 
    res = res * 1._wp/sqrt(2._wp*pi*Wp_p%sigma_perp**2) * exp(-k_perp**2/(4._wp*Wp_p%sigma_perp**2))

    ! Now add phases 
    res = res * exp(cmplx(0._wp, -kz**2/2 * Wp_p%tz, wp)) * &
                exp(cmplx(0._wp, -k_perp**2/2 * Wp_p%t_perp, wp)) * & 
                exp(cmplx(0._wp, kz*Wp_p%z0, wp))
end function PINEM_wp

! Function to determine the normalization factor of a WP on the grid. Assumes no phi dependence.
subroutine normalize_WP(WP_p)
    type(WP_parameters), intent(inout) :: WP_p
    real(wp) :: norm
    real(wp) :: k_integrand(size(k_list)), theta_integrand(size(theta_quad))
    integer :: i,j 

    WP_p%norm_factor = 1._wp  ! Start by setting this to one (in case WP depends on it)
    norm = 2._wp * pi  ! From phi integral

    write(*,*) 'Performing numerical normalization of WP on grid...'
    
    ! Perform double integral over k/theta. 
    do i=1, size(k_list)
        do j=1, size(theta_quad)
            theta_integrand(j) = sin(theta_quad(j)) * abs(eval_WP(k_list(i), theta_quad(j), WP_p))**2
        end do 
        k_integrand(i) = k_list(i)**2 * sum(theta_integrand * theta_w)
    end do 
    norm = norm * trapz(k_integrand, k_list(2)-k_list(1))
    WP_p%norm_factor = 1._wp / sqrt(norm)
    write(*,*) 'Done! WP norm on grid was: ', norm
end subroutine normalize_WP

end program main