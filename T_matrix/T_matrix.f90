module T_matrix 
use utils
use spherical_bessel
use GaussLegendre
use LAPACK_routines
implicit none 
private 

! Define working precision
integer, parameter :: wp = dp 

! Variables used to load data file with potential 
integer, save :: N_data  ! Used to hold the nr. of data points in the files 
real(wp), save, allocatable :: r_data(:), V00_data(:), V01_data(:), V11_data(:)  ! Used to hold the data points
real(wp), save :: dr  ! Used to hold the step size of the radial grid 
real(wp), save :: E0, E1  ! Used to store the energy of the two level system (TODO LOAD THIS IN SOMEHOW!)
real(wp), save :: sqrt2 = sqrt(2._wp)
real(wp), save :: sqrt_pi = sqrt(pi)
real(wp), save :: pi_fac = 1._wp / (2._wp * pi**2)
real(wp), save :: log2 = log(2._wp)

real(wp), save, allocatable :: test(:)

! Define public routines 
public calculate_T_matrix_elastic, calculate_several_m_elastic, load_data, calculate_several_m_two_channel

contains 

! Subroutine to load data (TODO some warnings?)
subroutine load_data(r_max, Nr, hydrogen_bool)
    integer, intent(in) :: Nr 
    real(wp), intent(in) :: r_max 
    logical, intent(in) :: hydrogen_bool 

    ! Set the number of data points in the numerical integration 
    N_data = Nr

    ! Allocate arrays 
    allocate(r_data(Nr))
    allocate(V00_data(Nr))
    allocate(V01_data(Nr))
    allocate(V11_data(Nr))
    call linspace(r_data, 0._wp, r_max, Nr)

    if (hydrogen_bool) then 
        ! Make hydrogen dataset - premult the V-lists with r^2 here to avoid divergences 
        V00_data = -(r_data**2+ r_data) * exp(-2._wp * r_data)
        V01_data = 4._wp/(9._wp * sqrt(2._wp)) * exp(-3._wp/2._wp * r_data) * (r_data + 2._wp/3._wp) * r_data**2 
        V11_data = -exp(-r_data) * (6._wp + 2._wp*r_data + r_data**2 + 8._wp/r_data)/8._wp * r_data**2
        E0 = -0.5_wp
        E1 = -0.125_wp

    else 
        ! Load data or something here 
        V00_data = -0.48_wp 
        V11_data = -0.92_wp 
        V01_data = 0.1_wp 
        E0 = -0.5015219985837196_wp 
        E1 = -0.12955528032839306_wp 
    end if 

    write(*,*) "Data loaded!"
    write(*,*) E0, E1
    dr = r_data(2) - r_data(1)  ! Get the radial step size
end subroutine load_data

! Function to calculate the hypergeometric series by summing the series. Should work for abs(x) <= 0.5
function hypgeo_sum(a,b,c,x) result(res)
    real(wp), intent(in) :: a,b,c,x
    real(wp) :: eps = 1e-10_wp  ! Tolerance to determine convergence 
    real(wp) :: factor, res, aa, bb, cc
    integer :: i 
    res = 1._wp 
    factor = 1._wp 
    aa = a 
    bb = b 
    cc = c 
    do i=1, 1000
        factor = factor * (aa*bb/cc) * x/i
        res = res + factor 
        aa = aa + 1._wp
        bb = bb + 1._wp
        cc = cc + 1._wp
        if (abs(factor)<eps) then 
            exit
        end if 
    end do 
end function hypgeo_sum 

! Function to evaluate the Legendre function of the second kind using the hypergeometric series
function Ql_hyp(l,x)
    real(wp), intent(in) :: x 
    integer, intent(in) :: l 
    real(wp) Ql_hyp
    Ql_hyp = sqrt_pi * hypgeo_sum(0.5_wp*l+1._wp, 0.5_wp*l+0.5_wp, l+3._wp/2._wp, 1._wp/x**2) * &
            exp(log_gamma(l+1._wp) - (l+1._wp)*log2 - (l+1._wp)*log(x) - log_gamma(l+3._wp/2._wp))
end function Ql_hyp

! Function to calculate potential matrix elements for hydrogen using Laplace transform 
function hyd_V_mel_1s1s(l, kf, ki) result(res)
    integer, intent(in) :: l 
    real(wp), intent(in) :: kf, ki 
    real(wp) :: Ql, Ql1, Ql_deriv, x, res, int1, int2 
    
    x = (4._wp + kf**2 + ki**2)/(2._wp*kf*ki)

    if (x .ge. sqrt2) then 
        ! Use the Legendre expression 
        Ql = Ql_hyp(l, x)
        Ql1 = Ql_hyp(l+1, x)
        Ql_deriv = (l+1._wp)*(Ql1 - x*Ql)/(x**2 - 1._wp)  ! Recursion formula for derviative 
        int1 = 1._wp/(2._wp*kf*ki) * Ql 
        int2 = 1._wp/(2._wp*kf*ki) * Ql_deriv * 2._wp/(kf*ki)
        res = -(int1 - int2) * pi_fac 

        ! Perform sign test?! 
        !write(*,*) '1', res
        !res = pi_fac * trapz(sph_jn_list(l, ki*r_data) * &
        !sph_jn_list(l, kf*r_data) * V00_data, dr)
        !write(*,*) '2', res

    else 
        ! Not in easy convergence region of hypgeo, just do trapz integral 
        res = pi_fac * trapz(sph_jn_list(l, ki*r_data) * &
                sph_jn_list(l, kf*r_data) * V00_data, dr)  ! Remember r^2 is included in V!
    end if 
end function hyd_V_mel_1s1s

! Function to calculate the elastic T-matrix element for a momentum value, given the Guassian quadrature
! This is specific for hydrogen 1s -> 1s scattering, using the Legendre expression to eval matrix elements
function calculate_T_matrix_elastic_hydrogen(k0, Nk, m, k_quad, w_quad) result(T)
    real(wp), intent(in) :: k0  ! Momentum value for which to calculate the T-matrix
    integer, intent(in) :: Nk  ! Number of momentum points to use in the Guass quad 
    integer, intent(in) :: m  ! Magnetic quantum number
    real(wp), intent(in) :: k_quad(Nk), w_quad(Nk)  
    complex(wp) :: T 
    complex(wp) :: T_list(Nk+1)
    real(wp) :: Vm_list(Nk+1)
    real(wp) :: k_list(Nk+1)
    complex(wp) :: Mat(Nk+1, Nk+1)
    real(wp) :: D_sum 
    complex(wp) :: D 
    integer :: i,j 
    real(wp) :: mass = 1._wp  ! CHANGE THIS LATER? 
    
    k_list(1:Nk) = k_quad 
    k_list(Nk+1) = k0

    ! First calculate the Vm_list 
    do i=1, Nk+1
        Vm_list(i) = hyd_V_mel_1s1s(m, k_list(i), k0)
    end do

    ! Then calculate the matrix used for getting the T_list 
    do j=1, Nk+1
        do i=1, Nk+1
            ! Calculate the half on shell part 
            if (j == Nk+1) then 
                D_sum = sum(w_quad / (k0**2 - k_quad**2))
                D = cmplx(-8._wp*pi * mass*k0**2*D_sum, -4._wp*pi**2 * mass, wp)
            else 
                D = cmplx(8._wp*pi * mass * w_quad(j)*k_quad(j)**2 / (k0**2 - k_quad(j)**2), 0._wp, wp)
            end if 
            Mat(i,j) = -hyd_V_mel_1s1s(m, k_list(i), k_list(j)) * D

            if (i == j) then 
                Mat(i,j) = Mat(i,j) + 1._wp
            end if 
        end do
    end do

    ! Now take the inverse and find the T-matrix elements 
    call complex_inverse(Mat, Nk+1)
    T_list = matmul(Mat, Vm_list)
    T = T_list(Nk+1)
end function calculate_T_matrix_elastic_hydrogen


! Function to calculate the elastic T-matrix element for a momentum value, given the Guassian quadrature
function calculate_T_matrix_elastic(k0, Nk, m, k_quad, w_quad) result(T)
    real(wp), intent(in) :: k0  ! Momentum value for which to calculate the T-matrix
    integer, intent(in) :: Nk  ! Number of momentum points to use in the Guass quad 
    integer, intent(in) :: m  ! Magnetic quantum number
    real(wp), intent(in) :: k_quad(Nk), w_quad(Nk)  
    complex(wp) :: T 
    complex(wp) :: T_list(Nk+1)
    real(wp) :: Vm_list(Nk+1)
    real(wp) :: k_list(Nk+1)
    complex(wp) :: Mat(Nk+1, Nk+1)
    real(wp) :: D_sum 
    complex(wp) :: D 
    integer :: i,j 
    real(wp) :: mass = 1._wp  ! CHANGE THIS LATER? 
    
    k_list(1:Nk) = k_quad 
    k_list(Nk+1) = k0

    ! First calculate the Vm_list 
    do i=1, Nk+1
        Vm_list(i) = trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * sph_jn_list(m, k0*r_data) * V00_data, dr)
    end do
    Vm_list = Vm_list / (2._wp * pi**2)

    ! Then calculate the matrix used for getting the T_list 
    do j=1, Nk+1
        do i=1, Nk+1
            ! Calculate the half on shell part 
            if (j == Nk+1) then 
                D_sum = sum(w_quad / (k0**2 - k_quad**2))
                D = cmplx(-8._wp*pi * mass*k0**2*D_sum, -4._wp*pi**2 * mass, wp)
            else 
                D = cmplx(8._wp*pi * mass * w_quad(j)*k_quad(j)**2 / (k0**2 - k_quad(j)**2), 0._wp, wp)
            end if 
            Mat(i,j) = - 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
            sph_jn_list(m, k_list(j)*r_data) * V00_data, dr) * D

            if (i == j) then 
                Mat(i,j) = Mat(i,j) + 1._wp
            end if 

            write(*,*) 'End!'
        end do
    end do

    ! Now take the inverse and find the T-matrix elements 
    call complex_inverse(Mat, Nk+1)
    T_list = matmul(Mat, Vm_list)
    T = T_list(Nk+1)
end function calculate_T_matrix_elastic


! Subtoruine to calculate elastic T-matrix for several m values 
subroutine calculate_several_m_elastic(T_list, max_m, k0, Nk, k_quad, w_quad, hydrogen_bool)
    integer, intent(in) :: max_m, Nk 
    complex(wp), intent(inout) :: T_list(max_m+1)
    real(wp), intent(in) :: k0 
    real(wp), intent(in) :: k_quad(Nk), w_quad(Nk)
    logical, intent(in) :: hydrogen_bool
    integer :: i 

    !$OMP PARALLEL DO
    do i=1, max_m+1
        print *, "Calculating T-matrix for m = ", i-1
        if (hydrogen_bool) then 
            T_list(i) = calculate_T_matrix_elastic_hydrogen(k0, Nk, i-1, k_quad, w_quad)
            !write(*,*), i, T_list(i)
        else 
            T_list(i) = calculate_T_matrix_elastic(k0, Nk, i-1, k_quad, w_quad)
        end if 
    end do 
    !$OMP END PARALLEL DO
end subroutine calculate_several_m_elastic


! Function to calculate the two channel T-matrix elements for a momentum value.
! NB this function only calculates the T-matrix for 0->0 and 0->1 channels (not 1->1)
function calculate_T_matrix_two_channel(k0, N_quad, m, k_quad, w_quad) result(T)
    real(wp), intent(in) :: k0  ! Momentum value for which to calculate the T-matrix
    integer, intent(in) :: N_quad  ! Number of momentum points to use in the Guass quad 
    integer, intent(in) :: m  ! Magnetic quantum number
    real(wp), intent(in) :: k_quad(N_quad), w_quad(N_quad) 
    real(wp) :: k1_square, k1  ! Critical momentum value for second channel  
    complex(wp) :: T(2)
    complex(wp) :: T_list(2*N_quad+4)
    real(wp) :: Vm_list(2*N_quad+4)
    real(wp) :: k_list(2*N_quad+4)
    complex(wp) :: Mat(2*N_quad+4, 2*N_quad+4)
    complex(wp) :: D 
    integer :: i,j
    real(wp) :: mass = 1._wp  ! CHANGE THIS LATER? 
    logical :: channel_open = .false.

    ! Determine whether the second channel is open and find the critical momentum value 
    k1_square = k0**2 + 2._wp * mass * (E0-E1)
    if (k1_square < 0._wp) then 
        ! The channel is closed 
        !write(*,*) 'Channel is closed'
        k1 = 0._wp  ! Set this value to something arbitrary, it is not going to be used? 
        channel_open = .false.
    else 
        !write(*,*) 'Channel is open'
        k1 = sqrt(k1_square)
        channel_open = .true.
    end if 
    
    ! Add momenta values to collective array 
    k_list(1:N_quad) = k_quad 
    k_list(N_quad+1) = k0
    k_list(N_quad+2) = k1 
    k_list(N_quad+3 : 2*N_quad+2) = k_quad
    k_list(2*N_quad+3) = k0 
    k_list(2*N_quad+4) = k1 

    ! Calculate the cyldrincal expansions of the potentials 
    do i=1, N_quad 
        Vm_list(i) = trapz(r_data**2 * sph_jn_list(m, k_quad(i)*r_data) * sph_jn_list(m, k0*r_data) * V00_data, dr)
    end do
    Vm_list(N_quad+1) = trapz(r_data**2 * sph_jn_list(m, k0*r_data) * sph_jn_list(m, k0*r_data) * V00_data, dr)
    
    if (channel_open) then  ! CHECK IF CHANNEL IS EVEN OPEN? (Should I have this?)
        Vm_list(N_quad+2) = trapz(r_data**2 * sph_jn_list(m, k1*r_data) * sph_jn_list(m, k0*r_data) * V00_data, dr)

        do i=1, N_quad
            Vm_list(N_quad+2 + i) = trapz(r_data**2 * sph_jn_list(m, k_quad(i)*r_data) * sph_jn_list(m, k1*r_data) * V01_data, dr)
        end do
        Vm_list(2*N_quad+3) = trapz(r_data**2 * sph_jn_list(m, k0*r_data) * sph_jn_list(m, k1*r_data) * V01_data, dr)
        Vm_list(2*N_quad+4) = trapz(r_data**2 * sph_jn_list(m, k1*r_data) * sph_jn_list(m, k1*r_data) * V01_data, dr)
    else 
        Vm_list(N_quad+2 : 2*N_quad+4) = 0._wp
    end if 

    Vm_list = Vm_list / (2._wp * pi**2)

    ! Calculate the matrix used for getting the T_list 
    do j=1, 2*N_quad+4
        do i=1, 2*N_quad+4
            ! Determine the channel that kj belongs to 
            ! First the three cases where kj is in the first channel 
            D = 0._wp 
            if (j <= N_quad) then 
                ! Channel 0 - not critical momentum 
                D = 8._wp*mass*pi * w_quad(j) * k_list(j)**2 / (k0**2 - k_list(j)**2)
                ! Now deternine the correct potential matrix element to calculate 
                if (i <= N_quad+2) then 
                    D = D * 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
                    sph_jn_list(m, k_list(j)*r_data) * V00_data, dr)
                else 
                    D = D * 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
                    sph_jn_list(m, k_list(j)*r_data) * V01_data, dr)
                end if 

            else if (j == N_quad+1) then 
                ! Channel 0 - ciritical momentum k0 
                D = -4._wp*mass*pi * cmplx(2._wp*k0**2 * sum(w_quad / (k0**2 - k_quad**2)), pi*k0, wp) 
                ! Now deternine the correct potential matrix element to calculate (Same as above)
                if (i <= N_quad+2) then 
                    D = D * 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
                    sph_jn_list(m, k_list(j)*r_data) * V00_data, dr)
                else 
                    D = D * 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
                    sph_jn_list(m, k_list(j)*r_data) * V01_data, dr)
                end if 

            else if (j == N_quad+2) then 
                ! Channel 0 - critical momentum k1 
                D = 0._wp
                

            ! Now the three cases where kj is in the second channel 
            else if (j <= 2*N_quad+2) then   
                ! Channel 1 - not critical momentum 
                if (channel_open) then  ! SHOULD I CHECK IF CHANNEL IS OPEN HERE???
                    D = 8._wp*mass*pi * w_quad(j-N_quad-2) * k_list(j)**2 / (k1**2 - k_list(j)**2)
                    ! Now deternine the correct potential matrix element to calculate (Same as above)
                    if (i <= N_quad+2) then 
                        D = D * 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
                        sph_jn_list(m, k_list(j)*r_data) * V01_data, dr)
                    else 
                        D = D * 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
                        sph_jn_list(m, k_list(j)*r_data) * V11_data, dr)
                    end if 
                else 
                    D = 0._wp 
                end if 

            else if (j == 2*N_quad+3) then 
                ! Channel 1 - critical momentum k0 
                D = 0._wp

            else if (j == 2*N_quad+4) then 
                ! Channel 1 - critical momentum k1 
                if (channel_open) then 
                    D = -4._wp*pi*mass * cmplx(2._wp*k1**2 * sum(w_quad / (k1**2 - k_quad**2)), pi*k1, wp)
                    ! Now deternine the correct potential matrix element to calculate (Same as above)
                    if (i <= N_quad+2) then 
                        D = D * 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
                        sph_jn_list(m, k_list(j)*r_data) * V01_data, dr)
                    else 
                        D = D * 1._wp/(2._wp*pi**2) * trapz(r_data**2 * sph_jn_list(m, k_list(i)*r_data) * &
                        sph_jn_list(m, k_list(j)*r_data) * V11_data, dr)
                    end if
                else 
                    D = 0._wp 
                end if 

            else 
                write(*,*) 'Everybody panic!'
            end if 
            
            ! Add to the matrix 
            Mat(i,j) = -1._wp * D
            if (i == j) then 
                Mat(i,j) = Mat(i,j) + 1._wp

                ! CHECK FOR NAN
                !if (Mat(i,j) /= Mat(i,j)) then 
                !    write(*,*) m, i, j, Mat(i,j)
                !end if 
            end if 
        end do
    end do

    ! Now take the inverse and find the on shell T-matrix elements 
    call complex_inverse(Mat, 2*N_quad+4)
    T_list = matmul(Mat, Vm_list)
    T(1) = T_list(N_quad+1)
    T(2) = T_list(2*N_quad+3)
    !write(*,*) T(1), T(2)
end function calculate_T_matrix_two_channel


! Subtoruine to calculate elastic T-matrix for several m values 
subroutine calculate_several_m_two_channel(T_list, max_m, k0, N_quad, k_quad, w_quad)
    integer, intent(in) :: max_m, N_quad
    complex(wp), intent(out) :: T_list(max_m+1, 2)
    real(wp), intent(in) :: k0 
    real(wp), intent(in) :: k_quad(N_quad), w_quad(N_quad)
    integer :: i 

    !$OMP PARALLEL DO
    do i=1, max_m+1
        write(*,*) "Calculating T-matrix for m = ", i-1
        T_list(i,:) = calculate_T_matrix_two_channel(k0, N_quad, i-1, k_quad, w_quad)
        !write(*,*) "Done calculating T-matrix for m = ", i-1, T_list(i,:)
    end do 
    !$OMP END PARALLEL DO
end subroutine calculate_several_m_two_channel


end module T_matrix 