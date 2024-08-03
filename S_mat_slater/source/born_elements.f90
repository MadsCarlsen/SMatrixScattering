module born_elements
use utils 
use spherical_harmonics
implicit none 
private 

! Define working precision 
integer, parameter :: wp = dp 

public :: calculate_born_element, calculate_born_analytic

contains 

! Main function collects all the elements of caculating the Born elements - input are two slater states, initial and final momentum 
function calculate_born_element(delta, cos_delta_theta, delta_phi, slater_f, slater_i, wigner_0, wigner_m, elastic) result(res)
    type(SlaterState), intent(in) :: slater_f, slater_i
    real(wp), intent(in) :: wigner_0(:), wigner_m(:)   ! Arrays of neede Wigner-3j symbols 
    real(wp), intent(in) :: delta, cos_delta_theta, delta_phi
    logical, intent(in) :: elastic
    complex(wp) :: res 
    complex(wp) :: ylm_fac
    integer :: lf, li, mf, mi, l, mu_i, mu_f, n_sum, w_count, w_count0, Nl  ! Intial and final angular momentum and magnetic quantum numbers
    real(wp) :: alpha_sum, R_sum
     
    complex(wp) :: R_test

    ! Get quantum numbers 
    lf = slater_f%l
    li = slater_i%l
    mf = slater_f%m
    mi = slater_i%m

    ! Should do check on elastic case and also prevent divergence! 

    ! Calculate the series of Wigner-3j symbols needed (done recursively so get all l values needed for the sum)
    ! These should be found outside, if several calls are needed and only delta varies 
    !call wigner3j_j1(li, lf, 0, 0, 0, wigner_000)
    !call wigner3j_j1(li, lf, mf-mi, mi, mf, wigner_mmm)

    w_count = 1  ! Counter to get the correct wigner symbols 
    if (abs(mf-mi) > abs(lf-li)) then 
        ! Here we need to shift the wigner_0 index, since we are not starting at lowest possible l
        w_count0 = abs(mf-mi) - abs(lf-li) + 1
    else 
        w_count0 = 1 
    end if 

    res = 0._wp 
    do l = max(abs(lf-li), abs(mf-mi)), lf+li 
        ! Calculate spherical harmonics factor 
        ylm_fac = sph_harm(cos_delta_theta, delta_phi, l, mf - mi)
        ylm_fac = conjg(ylm_fac) * cmplx(0._wp, 1._wp, wp)**(l) * (-1)**mf * &
              sqrt((2._wp*li+1._wp)*(2._wp*l+1._wp)*(2._wp*lf+1._wp)/(4._wp*pi))

        ! Calculate sum over radial integrals of the Slater functions 
        R_sum = 0._wp
        do mu_f=1, slater_f%Ns 
            do mu_i=1, slater_i%Ns
                alpha_sum = slater_f%alpha(mu_f) + slater_i%alpha(mu_i)
                n_sum = slater_f%n(mu_f) + slater_i%n(mu_i)
                R_sum = R_sum + slater_f%c(mu_f) * slater_i%c(mu_i) * &
                        radial_integral(delta, n_sum, alpha_sum, l)
            end do 
        end do 
    
        ! Add to the l sum 
        res = res + ylm_fac * wigner_0(w_count0) * wigner_m(w_count) * R_sum
        w_count = w_count + 1
        w_count0 = w_count0 + 1 
    end do
    
    ! Mult with remaining factors and return 
    res = 2._wp/(pi*delta**2) * res 

    ! If elasic scattering need to subtract the core term
    if (elastic) then 
        res = res - 1._wp / (2._wp*pi**2 * delta**2)
    end if 
end function calculate_born_element

function calculate_born_analytic(delta, cos_delta_theta, delta_phi, slater_f, slater_i) result(res)
    type(SlaterState), intent(in) :: slater_f, slater_i 
    real(wp), intent(in) :: delta, cos_delta_theta, delta_phi
    complex(wp) :: res 
    integer :: lf, li, mf, mi   ! Intial and final angular momentum and magnetic quantum numbers

    ! Get quantum numbers 
    lf = slater_f%l
    li = slater_i%l
    mf = slater_f%m
    mi = slater_i%m

    ! For p-1 -> 2s
    !res = 1._wp/(2._wp*sqrt(12._wp)) * (8._wp*delta/(delta**2+1._wp)**3 + & 
    !         4._wp*delta*(delta**2-5._wp)/(delta**2+1._wp)**4)
    !res = res * -1._wp/(2*pi**2 * delta**2) * (-1._wp)**(0) * cmplx(0._wp, 1._wp, wp) * & 
    !         sqrt(4._wp*pi) * conjg(sph_harm(cos_delta_theta, delta_phi, 1, -mi))

    ! For 2p-1 -> 2p1
    res = 4._wp*pi/(2._wp*pi**2) * sqrt(6._wp)/15._wp * conjg(sph_harm(cos_delta_theta, delta_phi, 2, 2)) * &
          sqrt(45._wp/(4._wp*pi)) * 2._wp / (delta**2 + 1._wp)**4 
end function calculate_born_analytic


! Wrapper function to evaluate the hypergeometric series. Uses linear transformation formula if abs(x) > 0.75
function eval_hypgeo(a, b, c, x) result(res)
    real(wp), intent(in) :: a,b,c,x
    real(wp) :: res
    real(wp) :: gamma_fac1, gamma_fac2 
    res = hypgeo_sum(a,b,c,x)

    
    !if (abs(x) < 0.75_wp) then 
    !    res = hypgeo_sum(a,b,c,x)
    !else 
    !    if (c == 0 .or. (c-a-b) == 0 .or. (c-a)==0 .or. (c-b)==0) then 
    !        write(*,*) c, c-a-b, c-a, c-b
    !        print*, 'Hypergeometric series diverges'
    !        print*, x
    !        return
    !    end if
    !    ! Here we have to use a linear transformation formula to hopefully get a convering series 
    !    gamma_fac1 = exp(log_gamma(c) + log_gamma(c-a-b) - log_gamma(c-a) - log_gamma(c-b))
    !    gamma_fac2 = exp(log_gamma(c) + log_gamma(a+b-c) - log_gamma(a) - log_gamma(b))
    !    res = gamma_fac1 * hypgeo_sum(a, b, a+b-c+1._wp, 1._wp-x) + &
    !          (1._wp-x)**(c-a-b) * gamma_fac2 * hypgeo_sum(c-a, c-b, c-a-b+1._wp, 1._wp-x) 
    !end if
end function eval_hypgeo


! Function to calculate the hypergeometric series by summing the series. Should work for abs(x) <= 0.5
function hypgeo_sum(a,b,c,x) result(res)
    real(wp), intent(in) :: a,b,c,x
    real(wp) :: eps = 1e-12_wp  ! Tolerance to determine convergence 
    real(wp) :: factor, res, aa, bb, cc
    integer :: i 
    
    res = 1._wp 
    factor = 1._wp 
    aa = a 
    bb = b 
    cc = c 
    do i=1, 2000
        factor = factor * (aa*bb/cc) * x/i
        res = res + factor 
        aa = aa + 1._wp
        bb = bb + 1._wp
        cc = cc + 1._wp
        if (abs(factor)<eps) then 
            !if (i>100) then 
            !    write(*,*) i, x, a, b, c, res
            !end if 
            return
        end if 
    end do 
    print*, 'Hypergeometric series did not converge'
    print*, x
    stop 
end function hypgeo_sum 

! Radial integral. Here n = n_i + n_f, alpha = alpha_i + alpha_f
function radial_integral(delta, n, alpha, l) result(res)
    real(wp) :: delta, alpha, res 
    integer :: n, l 
    real(wp) :: mu, nu, gamma_fac 
    
    mu = n + 0.5_wp 
    nu = l + 0.5_wp 
    gamma_fac = exp(log_gamma(nu+mu) - log_gamma(nu+1._wp))  ! Use log-gamma to avoid overflow 

    res = sqrt(pi)/(2**(l+1)) * delta**l * gamma_fac / sqrt((alpha**2 + delta**2)**(nu+mu)) * &
          eval_hypgeo((nu+mu)/2._wp, (1._wp-mu+nu)/2._wp, nu+1._wp, delta**2/(alpha**2+delta**2))
end function radial_integral 


function radial_integral_mod(delta, n, alpha, l) result(res)
    real(wp) :: delta, alpha, res 
    integer :: n, l 
    real(wp) :: mu, nu 

    mu = n + 0.5_wp
    nu = l + 0.5_wp
    res = (alpha**2 + delta**2)**(-0.5_wp*mu) * gamma(nu+mu) * & 
                assoc_legendre_hypgeo(alpha/sqrt(alpha**2+delta**2), mu-1._wp, -nu)

end function radial_integral_mod

function assoc_legendre_hypgeo(z, nu, mu) result(res)
    real(wp) :: z, mu, nu 
    real(wp) :: res
    real(wp) :: factor 

    write(*,*) z, mu, (z+1._wp)/(z-1._wp)
    res = ((z+1._wp)/(z-1._wp))**(mu/2._wp) / gamma(1._wp-mu) * hypgeo_sum(-nu, nu+1._wp, 1._wp - mu, (1-z)/2._wp)

    write(*,*) res 
end function assoc_legendre_hypgeo


end module born_elements
