program main 
use utils 
use interpolate
implicit none 

integer, parameter :: wp = dp 

! Parameters determining what to calculate 
logical :: cal_PINEM, cal_Gauss, cal_tot_elastic, cal_scatt_elastic, cal_scatt_inelastic, cal_interference 

! Parameters for calculating of the S-matrix and general integration    
integer :: Nk, Nk_T_mat
real(wp) :: k_low, k_high 
integer :: N_theta_int, N_theta_final, N_phi_int, N_phi_final
real(wp) :: theta_int_min, theta_int_max, theta_final_min, theta_final_max
real(wp) :: phi_int_min = 0._wp, phi_int_max = 2._wp*pi
real(wp) :: phi_final_min, phi_final_max 
logical :: interpolate_bool 

! Parameters for the target 
real(wp) :: E0, E1  ! These are loaded
integer :: max_l 
real(wp) :: E10, period 
character(100) :: run_id 
character(:), allocatable :: run_id_trim

! Parameters for the WP 
real(wp) :: k0z, sigma_time, sigma_perp_factor, omega, g_factor, b  ! These are loaded 
integer :: max_N  ! This is loaded 
real(wp) :: t_bunch, z0, tz, t_perp, sigma_space_z, sigma_space_perp, sigma_z, sigma_perp 

! Arrays for integration and results 
real(wp), allocatable :: theta_int(:), theta_final(:), k_list(:), phi_int(:), phi_final(:)
complex(wp), allocatable :: T00(:,:), T01(:,:)
complex(wp) :: T_vec(2)
integer :: i,j 
real(wp) :: ki, kf 
real(wp) :: mass = 1._wp
real(wp), allocatable :: res(:,:)
 
! Spline arrays 
integer, parameter :: N_spline = 12
real(wp) :: x_spline(N_spline) 
complex(wp) :: y_spline(N_spline), spline_coeffs(N_spline)

! TEST 
integer, parameter :: Nt = 100
integer :: ti
real(wp) :: t_scan_list(Nt), time_scan_res(Nt)
character(200) :: save_file
real(wp) :: delta_k


! Define namelists 
namelist /CALCULATE_SETTINGS/ cal_PINEM, cal_Gauss, cal_tot_elastic, cal_scatt_elastic, cal_scatt_inelastic, cal_interference

namelist /S_MAT_SETTINGS/ Nk, N_theta_int, N_theta_final, theta_int_min &
                        , theta_int_max, theta_final_min, theta_final_max, k_low, k_high, &
                        phi_final_min, phi_final_max, N_phi_int, N_phi_final, interpolate_bool 

namelist /TARGET_SETTINGS/ E0, E1, b, run_id 

namelist /WAVEPACKET_SETTINGS/ k0z, sigma_time, sigma_perp_factor, omega, g_factor, max_N 

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


! Calculate parameters 
E10 = E1 - E0 
period = 2._wp * pi / E10 
sigma_time = sigma_time * period 
sigma_space_z = sigma_time * k0z 
sigma_space_perp = sigma_space_z * sigma_perp_factor 
sigma_z = 1.0505169575771722e-05_wp !1._wp / (2._wp * sigma_space_z)
sigma_perp = 0.13555319735133256_wp !1._wp / (2._wp * sigma_space_perp)  ! Should this be different when perp? 
!omega = 0.02278167626453977_wp!0.06509050361297077_wp !0.02278167626453977_wp !omega * E10 
t_bunch = k0z**2 / (g_factor * omega**2)
delta_k = omega / k0z 

write(*,*) k0z, k_low, k_high 
write(*,*) sigma_z, sigma_perp, omega, t_bunch

! Load the dims 
run_id_trim = trim(run_id)
!run_id = trim(run_id)
write(*,*) 'Reading T-mat data from : ', run_id_trim
call read_dims(max_l, Nk_T_mat, run_id_trim)

! Check if splined call, else Nk must match Nk_T_mat
if (.not. interpolate_bool) then 
    if (Nk /= Nk_T_mat) then 
        write(*,*) 'Nk does not match the T-matrix dimensions!'
        stop 
    end if 
end if

! Setup arrays for integration

! Hacks for time integration 
N_theta_final = 1 

allocate(res(N_theta_final, N_phi_final))
allocate(k_list(Nk))
allocate(theta_int(N_theta_int))
allocate(theta_final(N_theta_final))
allocate(phi_int(N_phi_int))
allocate(phi_final(N_phi_final))
allocate(T00(max_l+1, Nk_T_mat))
allocate(T01(max_l+1, Nk_T_mat))
call linspace(theta_int, theta_int_min, theta_int_max, N_theta_int)
!call linspace(theta_final, theta_final_min, theta_final_max, N_theta_final)
theta_final = acos(1 - delta_k / k0z)
write(*,*) 'theta_final = ', theta_final

call linspace(k_list, k_low, k_high, Nk)
call linspace(phi_int, phi_int_min, phi_int_max, N_phi_int)
if (N_phi_final == 1) then 
    phi_final = 0._wp 
else
    call linspace(phi_final, pi*phi_final_min, pi*phi_final_max, N_phi_final)
end if 


! LOAD THE T-MATRIX ELEMENTS - only reads T01 if inealstic==.true. 
! NB! SHOULD CALCULATE THE T-MATRIX ELEMENTS FOR KF INSTEAD OF KI!!!! 
call read_T_mat(max_l, Nk_T_mat, run_id_trim, T00, T01)

! Interpolate if chosen - will be used to evaluate the T-matrix at arbitary angles
! NB: This picks one k value, only works if the T-matrix does not vary much in k (which it does not)
if (interpolate_bool) then 
    call linspace(x_spline, cos(2._wp*theta_int_max), 1._wp, N_spline)
    do i=1, N_spline 
        y_spline(i) = eval_T(T00(:,1), x_spline(i))
    end do
    call cubic_spline(x_spline, y_spline, spline_coeffs, N_spline)

    ! Fix T00 and T01 to be arrays of correct size to be used in funcs.. Stupid hack, but they are not used in this case...
    deallocate(T00)
    deallocate(T01)
    allocate(T00(1, Nk))
    allocate(T01(1, Nk))
    T00 = 0._wp 
    T01 = 0._wp
end if 

! CALCULATE ANGULAR DISTRIBUTIONS
tz = 0._wp 
t_perp = -t_bunch !* 1.9_wp*pi   
z0 = k0z * t_bunch !* 1.9_wp*pi

! MANUAL CALCULATIONS CAN GO HERE  
!allocate(P_test(N_theta_final))
!write(*,*) 'Calculating PINEM elastic scattering...'
!P_test = angular_dist_single(k_list, T00, theta_final, theta_int, phi_int, tz, t_perp, z0, b, omega, .true.)
!call savetxt('data/P00_PINEM.txt', P_test)


call linspace(t_scan_list, 1._wp, 3._wp*pi, Nt)

do ti=1, Nt 

tz = 0._wp 
t_perp = -t_bunch * t_scan_list(ti)   
z0 = k0z * t_bunch * t_scan_list(ti)
    

write(*,*) 'Calculating time scan ', ti, '/', Nt

! BOOLEAN TOGGLE FOR DIFFERENT CALCULATIONS 
if (cal_PINEM) then 
    if (cal_tot_elastic) then 
        !write(*,*) 'Calculating PINEM elastic, total...'
        res = angular_dist_elastic(k_list, T00, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, .true.)
        !write(save_file, '(A,I0,A)') 'data_time_scan/PINEM_tot_elastic_', ti, '.dat'
        !call save_array(save_file, res)
        time_scan_res(ti) = res(1,1)
    end if 

    if (cal_scatt_elastic) then 
        !write(*,*) 'Calculating PINEM elastic, scattered only...'
        res = angular_dist_scattered(k_list, T00, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, &
        .true., .true.)
        !write(save_file, '(A,I0,A)') 'data_time_scan/PINEM_scatt_elastic_', ti, '.dat'
        !call save_array(save_file, res)
        time_scan_res(ti) = res(1,1)
    end if 

    if (cal_scatt_inelastic) then 
        !write(*,*) 'Calculating PINEM inelastic...'
        res = angular_dist_scattered(k_list, T01, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, &
        .true., .false.)
        !write(save_file, '(A,I0,A)') 'data_time_scan/PINEM_scatt_inelastic_', ti, '.dat'
        !call save_array(save_file, res)
        time_scan_res(ti) = res(1,1)
    end if
    
    if (cal_interference) then 
        !write(*,*) 'Calculating PINEM interference...'
        res = angular_dist_interference(k_list, T00, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, .true.) 
        !write(save_file, '(A,I0,A)') 'data_time_scan/PINEM_interference_', ti, '.dat'
        !call save_array(save_file, res)
        time_scan_res(ti) = res(1,1)
    end if 
end if 

if (cal_Gauss) then 
    if (cal_tot_elastic) then 
        !write(*,*) 'Calculating Gauss elastic, total...'
        res = angular_dist_elastic(k_list, T00, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, .false.)
        !write(save_file, '(A,I0,A)') 'data_time_scan/Gauss_tot_elastic_', ti, '.dat'
        !call save_array(save_file, res)
        time_scan_res(ti) = res(1,1)
    end if 

    if (cal_scatt_elastic) then 
        !write(*,*) 'Calculating Gauss elastic, scattered only...'
        res = angular_dist_scattered(k_list, T00, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, &
        .false., .true.)
        !write(save_file, '(A,I0,A)') 'data_time_scan/Gauss_scatt_elastic_', ti, '.dat'
        !call save_array(save_file, res)
        time_scan_res(ti) = res(1,1)
    end if 

    if (cal_scatt_inelastic) then 
        !write(*,*) 'Calculating Gauss inelastic...'
        res = angular_dist_scattered(k_list, T01, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, &
        .false., .false.)
        !write(save_file, '(A,I0,A)') 'data_time_scan/Gauss_scatt_inelastic_', ti, '.dat'
        !call save_array(save_file, res)
        time_scan_res(ti) = res(1,1)
    end if
    
    if (cal_interference) then 
        !write(*,*) 'Calculating Gauss interference...'
        res = angular_dist_interference(k_list, T00, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, .false.) 
        !write(save_file, '(A,I0,A)') 'data_time_scan/Gauss_interference_', ti, '.dat'
        !call save_array(save_file, res)
        time_scan_res(ti) = res(1,1)
    end if 
end if 

end do 

! Also save the dimensions of the arrays in a .txt file
call savetxt('data_time_scan/times.txt', t_scan_list*t_bunch)
call savetxt('data_time_scan/res.txt', time_scan_res)


contains 

! Functtion to calculate the angular distribution of the scattered WP 
function angular_dist_elastic(k_list, T_data, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, PINEM_bool) 
    real(wp) :: k_list(:), theta_final(:), phi_final(:), theta_int(:), phi_int(:), tz, t_perp, z0, b, omega
    logical :: PINEM_bool
    complex(wp) :: T_data(:,:)
    real(wp) :: angular_dist_elastic(size(theta_final), size(phi_final)) 
    complex(wp) :: S_list(size(k_list))
    complex(wp) :: wp_arr(1)
    integer :: i,j,m

    !do i=1, size(theta_final)
        !write(*,*) i, '/', size(theta_final)
        !do j=1, size(phi_final)
            !write(*,*) j, '/', size(phi_final)
            !$OMP PARALLEL DO PRIVATE(wp_arr)
            do m=1, size(k_list)
                if (PINEM_bool) then 
                    wp_arr = PINEM_wp(k_list(m), [theta_final(1)], k0z, sigma_z, sigma_perp, tz, &
                                        t_perp, z0, omega, max_N, g_factor)
                else 
                    wp_arr = Gaussian_wp(k_list(m), [theta_final(1)], sigma_z, sigma_perp, k0z, tz, t_perp, z0)
                end if     
                S_list(m) = calculate_S_mat(k_list(m), k_list(m), T_data(:,m), theta_final(1), phi_final(1), theta_int, &
                            phi_int, tz, t_perp, z0, b, omega, PINEM_bool)  + wp_arr(1)  
    
                !write(*,*) S_list(m)
            end do
            !$OMP END PARALLEL DO

            ! Now calculate the k-integral 
            angular_dist_elastic(1,1) = trapz(k_list**2 * abs(S_list)**2, k_list(2)-k_list(1)) 
        !end do 
    !end do 

end function angular_dist_elastic 

! Functtion to calculate the angular distribution of the scattered WP 
function angular_dist_scattered(k_list, T_data, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, &
                                PINEM_bool, elastic_bool) 
    real(wp) :: k_list(:), theta_final(:), phi_final(:), theta_int(:), phi_int(:), tz, t_perp, z0, b, omega
    logical :: PINEM_bool, elastic_bool
    complex(wp) :: T_data(:,:)
    real(wp) :: angular_dist_scattered(size(theta_final), size(phi_final)) 
    complex(wp) :: S_list(size(k_list))
    real(wp) :: ki_list(size(k_list))
    integer :: i,j,m

    ! Find the correct initial momentum to be used in S-matrix calculations 
    if (elastic_bool) then 
        ki_list = k_list 
    else 
        ki_list = sqrt(k_list**2 + 2*mass*(E1-E0))
    end if

    !$OMP PARALLEL DO PRIVATE(S_list, j, m)
    do i=1, size(theta_final)
        do j=1, size(phi_final)
            do m=1, size(k_list)  
                S_list(m) = calculate_S_mat(k_list(m), ki_list(m), T_data(:,m), theta_final(i), phi_final(j), theta_int, &
                            phi_int, tz, t_perp, z0, b, omega, PINEM_bool)        
            end do
            ! Now calculate the k-integral 
            angular_dist_scattered(i,j) = trapz(k_list**2 * abs(S_list)**2, k_list(2)-k_list(1)) 
        end do 
    end do 
    !$OMP END PARALLEL DO
end function angular_dist_scattered

! Functtion to calculate the angular distribution of the scattered WP  
function angular_dist_interference(k_list, T_data, theta_final, phi_final, theta_int, phi_int, tz, t_perp, z0, b, omega, PINEM_bool) 
    real(wp) :: k_list(:), theta_final(:), phi_final(:), theta_int(:), phi_int(:), tz, t_perp, z0, b, omega
    logical :: PINEM_bool
    complex(wp) :: T_data(:,:)
    real(wp) :: angular_dist_interference(size(theta_final), size(phi_final)) 
    complex(wp) :: S_list(size(k_list))
    complex(wp) :: wp_arr(1)
    integer :: i,j,m

    !$OMP PARALLEL DO PRIVATE(S_list, wp_arr, j, m)
    do i=1, size(theta_final)
        do j=1, size(phi_final)
            do m=1, size(k_list)
                if (PINEM_bool) then 
                    wp_arr = PINEM_wp(k_list(m), [theta_final(i)], k0z, sigma_z, sigma_perp, tz, &
                                        t_perp, z0, omega, max_N, g_factor)
                else 
                    wp_arr = Gaussian_wp(k_list(m), [theta_final(i)], sigma_z, sigma_perp, k0z, tz, t_perp, z0)
                end if     
                S_list(m) = calculate_S_mat(k_list(m), k_list(m), T_data(:,m), theta_final(i), phi_final(j), theta_int, &
                            phi_int, tz, t_perp, z0, b, omega, PINEM_bool) * cmplx(0._wp, 2._wp) * conjg(wp_arr(1)) 
            end do
            ! Now calculate the k-integral 
            angular_dist_interference(i,j) = trapz(k_list**2 * aimag(S_list), k_list(2)-k_list(1)) 
        end do 
    end do 
    !$OMP END PARALLEL DO
end function angular_dist_interference 

! Function to calculate the S-matrix element at a given final momentum 
function calculate_S_mat(kf, ki, Tm_list, theta_f, phi_f, theta_int_list, phi_int_list, tz, &
                        t_perp, z0, b, omega, PINEM_bool) result(S_mat)
    logical :: PINEM_bool
    real(wp) :: kf, ki, theta_f, tz, t_perp, z0, omega, theta_int_list(:), b, phi_f, phi_int_list(:)
    complex(wp) :: Tm_list(:)
    complex(wp) :: S_mat
    complex(wp) :: T_list(size(theta_int_list)), theta_integrand(size(theta_int_list)), phi_integrand(size(phi_int_list))
    complex(wp) :: a_wp(size(theta_int_list))
    real(wp) :: cos_fac, sin_fac, cos_angle 
    integer :: i,j
    complex(wp) :: T_element
 
    if (PINEM_bool) then 
        a_wp = PINEM_wp(ki, theta_int_list, k0z, sigma_z, sigma_perp, tz, t_perp, z0, omega, max_N, g_factor)
    else 
        a_wp = Gaussian_wp(ki, theta_int_list, sigma_z, sigma_perp, k0z, tz, t_perp, z0)
    end if
    ! Perform double integral over theta and phi. Do phi first, since WP's are independent of phi 
    do i=1, size(theta_int_list)
        cos_fac = cos(theta_f) * cos(theta_int_list(i))
        sin_fac = sin(theta_f) * sin(theta_int_list(i))
        
        do j=1, size(phi_int_list)
            cos_angle = sin_fac * cos(phi_f-phi_int_list(j)) + cos_fac  
            
            if (interpolate_bool) then 
                T_element = eval_cubic_spline(cos_angle, x_spline, y_spline, spline_coeffs, N_spline)
            else 
                T_element = eval_T(Tm_list, cos_angle)
            end if 
            phi_integrand(j) = exp(cmplx(0._wp, b*ki*sin(theta_int_list(i))*cos(phi_int_list(j)),wp)) * T_element 
        end do 
        
        ! Perform the phi integral 
        theta_integrand(i) = sin(theta_int_list(i)) * a_wp(i) * trapz(phi_integrand, phi_int_list(2)-phi_int_list(1))        
    end do 
    
    ! Last, perform the theta integral 
    S_mat = trapz(theta_integrand, theta_int_list(2)-theta_int_list(1)) * cmplx(0._wp, -2._wp*pi*ki, wp) * & 
                  exp(cmplx(0._wp, -b*kf*sin(theta_f)*cos(phi_f),wp)) 
    !write(*,*) S_mat
end function calculate_S_mat

! Function to evaluate the T-matrix element for a final theta
function eval_T(Tm_list, x) result(T)
    complex(wp) :: Tm_list(:) 
    real(wp) :: x 
    complex(wp) :: T
    integer :: i
    real(wp) :: poly, poly1, polyl 
    
    ! Starting value for recursive evaluation of Legendre poly (We are now assuming at least 2 l in Tm_list...)
    poly = 1._wp 
    poly1 = x
    T = poly * Tm_list(1) + 3._wp*poly1*Tm_list(2)

    ! Now recursively evaluate the rest 
    do i=2, size(Tm_list)-1
        polyl = (x*(2._wp*i-1._wp)*poly1 - (i-1._wp)*poly)/i
        T = T + (2._wp*i+1._wp) * polyl * Tm_list(i+1)
        poly = poly1 
        poly1 = polyl  
    end do 
end function eval_T 

! Gaussian WP 
function Gaussian_wp(k, theta_list, sigma_z, sigma_perp, k0z, tz, t_perp, z0) result(res)
    real(wp) :: k, sigma_z, sigma_perp, k0z, tz, t_perp, z0, theta_list(:)
    complex(wp) :: res(size(theta_list))
    real(wp) :: kz(size(theta_list)), k_perp(size(theta_list))

    kz = k * cos(theta_list)
    k_perp = k * sin(theta_list) 

    res = 1._wp/(2._wp*pi*sigma_z**2)**(1._wp/4._wp) * &
            1._wp/sqrt(2._wp*pi*sigma_perp**2)
    res = res * exp(-(kz-k0z)**2/(4._wp*sigma_z**2)) * exp(-k_perp**2/(4._wp*sigma_perp**2))
    res = res * exp(cmplx(0._wp, -kz**2/2 * tz, wp)) * &
                exp(cmplx(0._wp, -k_perp**2/2 * t_perp, wp)) * & 
                exp(cmplx(0._wp, kz*z0, wp))
end function Gaussian_wp

function Gaussian_one_dim(k, k0, sigma) result(res)
    real(wp) :: k(:), k0, sigma
    complex(wp) :: res(size(k))
    res = 1._wp/(2._wp*pi*sigma**2)**(1._wp/4._wp) * &
            exp(-(k-k0)**2/(4._wp*sigma**2))
end function Gaussian_one_dim

function PINEM_wp(k, theta_list, k0z, sigma_z, sigma_perp, tz, t_perp, z0, omega, max_N, g_factor) result(res)
    real(wp) :: k, k0z, sigma_z, sigma_perp, tz, t_perp, z0, omega, g_factor, theta_list(:)
    integer :: max_N 
    complex(wp) :: res(size(theta_list))
    real(wp) :: delta, kz(size(theta_list)), k_perp(size(theta_list))
    real(wp) :: phi_N = 0._wp 
    integer :: N 

    delta = omega/k0z
    kz = k*cos(theta_list)
    k_perp = k*sin(theta_list)
    
    ! Build the PINEM part in x direction
    res = bessel_jn(0, g_factor) * exp(cmplx(0._wp, phi_N, wp)) * Gaussian_one_dim(kz, k0z, sigma_z)
    do N=1, max_N 
        res = res + bessel_jn(N, g_factor) * exp(cmplx(0._wp, phi_N, wp)) * Gaussian_one_dim(kz, k0z + N*delta, sigma_z)
        res = res + bessel_jn(-N, g_factor) * exp(cmplx(0._wp, phi_N, wp)) * Gaussian_one_dim(kz, k0z - N*delta, sigma_z)
    end do 

    ! Then add the Gaussian in y direction 
    res = res * 1._wp/sqrt(2._wp*pi*sigma_perp**2) * exp(-k_perp**2/(4._wp*sigma_perp**2))

    ! Now add phases 
    res = res * exp(cmplx(0._wp, -kz**2/2 * tz, wp)) * &
                exp(cmplx(0._wp, -k_perp**2/2 * t_perp, wp)) * & 
                exp(cmplx(0._wp, kz*z0, wp))
end function PINEM_wp

! Subroutine to load k and N dims.
subroutine read_dims(max_l, Nk, run_id)
    integer, intent(out) :: max_l, Nk
    character(*), intent(in) :: run_id
    open(1, file='T_mat_data/' // run_id // '/dims.txt')
    read(1,*) max_l
    read(1,*) Nk
    close(1)
end subroutine read_dims

subroutine read_T_mat(max_l, Nk, run_id, T00, T01)
    integer, intent(in) :: max_l, Nk
    character(*), intent(in) :: run_id
    complex(wp), intent(out) :: T00(max_l+1, Nk), T01(max_l+1, Nk)

    real(wp) :: T00_data_r(max_l+1, Nk)
    real(wp) :: T00_data_i(max_l+1, Nk)
    real(wp) :: T01_data_r(max_l+1, Nk)
    real(wp) :: T01_data_i(max_l+1, Nk)

    open(1, file='T_mat_data/' // run_id // '/T00_real.dat', form='unformatted')
    open(2, file='T_mat_data/' // run_id // '/T00_imag.dat', form='unformatted')
    read(1) T00_data_r 
    read(2) T00_data_i 
    T00 = cmplx(T00_data_r, T00_data_i, wp)

    if (cal_scatt_inelastic) then 
        open(3, file='T_mat_data/' // run_id // '/T01_real.dat', form='unformatted')
        open(4, file='T_mat_data/' // run_id // '/T01_imag.dat', form='unformatted')
        read(3) T01_data_r 
        read(4) T01_data_i
        T01 = cmplx(T01_data_r, T01_data_i, wp)
    else 
        T01 = 0._wp 
    end if
end subroutine read_T_mat

end program main 
