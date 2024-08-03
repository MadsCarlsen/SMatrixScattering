module Wigner3j 

use utils 
implicit none 
private 

integer, parameter :: wp = dp 

! Module for calculating Wigner 3j symbols and Clebsch-Gordan coefficients
public :: wigner3j_j1, wigner_3j_left, wigner_3j_right 


contains 

function cal_A(j1, j2, j3, m1) result(A)
    real(wp) :: j1, j2, j3, m1 
    real(wp) :: A
    A = sqrt(j1**2 - (j2-j3)**2) * sqrt((j2+j3+1._wp)**2 -j1**2) * sqrt(j1**2 - m1**2)
end function cal_A 

function cal_B(j1, j2, j3, m1, m2, m3) result(B)
    real(wp) :: j1, j2, j3, m1, m2, m3 
    real(wp) :: B
    B = -(2._wp*j1+1._wp) * (j2*(j2+1._wp)*m1 - j3*(j3+1._wp)*m1 - j1*(j1+1._wp)*(m3-m2))
end function cal_B 


subroutine wigner_3j_right(j2, j3, m1, m2, m3, result) 
    real(wp), intent(in) :: j2, j3, m1, m2, m3 
    real(wp), intent(out), allocatable :: result(:)
    real(wp) :: j1_min, j1_max, jr, jr_0, jr_1, sum, phase_sign 
    integer :: Nj, i
    
    j1_min = max(abs(j2-j3), abs(m1))
    j1_max = j2 + j3 
    Nj = int(j1_max - j1_min) + 1

    ! Make sure m1 is within the range of non-zero elements 
    if (j1_min > j1_max) then 
        allocate(result(1))
        result(1) = 0._wp
        return
    end if

    ! Make tests on m2/m3 
    if ((m2 < -j2 .or. m2 > j2) .or. (m3 < -j3 .or. m3 > j3)) then 
        allocate(result(1))
        result(1) = 0._wp
        return
    end if 

    allocate(result(Nj))

    ! Initial value for the right side
    jr_0 = 1._wp 
    result(Nj) = jr_0

    ! Now take first step 
    if (Nj > 1) then 
        jr_1 = -cal_B(j1_max, j2, j3, m1, m2, m3) / (cal_A(j1_max, j2, j3, m1) * (j1_max+1._wp))
        result(Nj-1) = jr_1

        ! Do the rest of the recursion 
        if (Nj > 2) then 
            do i=2, Nj-1 
                jr = - ((j1_max-i+1._wp)*cal_A(j1_max-i+2._wp, j2, j3, m1)*jr_0 + cal_B(j1_max-i+1._wp, j2, j3, m1, m2, m3)*jr_1) & 
                / ((j1_max-i+2._wp)*cal_A(j1_max-i+1._wp, j2, j3, m1)) 
    
                result(Nj-i) = jr
                jr_0 = jr_1
                jr_1 = jr
            end do 
        end if
    end if 

    ! Now normalize the results 
    sum = 0._wp 
    do i=1, Nj 
        sum = sum + (2._wp*(j1_min + i -1._wp)+1._wp) * result(i)**2
    end do 
    result = result / sqrt(sum)

    ! Last determine the phase by looking at the last element
    phase_sign = sign(1._wp, result(Nj))
    if (phase_sign /= (-1._wp)**(j2-j3-m1)) then
        result = result * (-1._wp)
    end if
end subroutine wigner_3j_right 


! Recursive calculation of Wigner 3j symbols, from one side only! 
subroutine wigner_3j_left(j2, j3, m1, m2, m3, result)
    real(wp), intent(in) :: j2, j3, m1, m2, m3 
    real(wp), intent(out), allocatable :: result(:)
    real(wp) :: j1_min, j1_max, jl, jl_0, jl_1, sum, phase_sign 
    integer :: Nj, i

    j1_min = max(abs(j2-j3), abs(m1))
    j1_max = j2 + j3 
    Nj = int(j1_max - j1_min) + 1

    allocate(result(Nj))

    ! Initial value for the left side 
    jl_0 = 1._wp 
    result(Nj) = jl_0
    
    if (Nj > 1) then 
        ! Take first step 
        jl_1 = -cal_B(j1_min, j2, j3, m1, m2, m3) / (cal_A(j1_min+1._wp, j2, j3, m1) * j1_min) 
        
        result(2) = jl_1 
        write(*,*) result 
        ! Do the rest of the recursion
        if (Nj > 2) then 
            do i=2, Nj-1 
                jl = -(cal_B(j1_min+i-1._wp, j2, j3, m1, m2, m3)*jl_1 + (j1_min+i)*cal_A(j1_min+i-1._wp, j2, j3, m1)*jl_0) &
                / ((j1_min+i-1._wp)*cal_A(j1_min+i, j2, j3, m1)) 
                result(i+1) = jl
                jl_0 = jl_1
                jl_1 = jl
            end do
        end if 
    end if 

    ! Now normalize the results 
    sum = 0._wp 
    do i=1, Nj 
        sum = sum + (2._wp*(j1_min + i -1._wp)+1._wp) * result(i)**2
    end do 
    result = result / sqrt(sum)

    ! Last determine the phase by looking at the last element
    phase_sign = sign(1._wp, result(Nj))
    if (phase_sign /= (-1._wp)**(j2-j3-m1)) then
        result = result * (-1._wp)
    end if
end subroutine wigner_3j_left


! The routine below performs recursion from both sides, should be numerically stable? WIP though 
subroutine wigner3j_j1(j2, j3, m1, m2, m3, results)
    ! Calculate all the Wigner 3j for the range of j1 
    real(wp), intent(in) :: j2, j3, m1, m2, m3
    real(wp), intent(out), allocatable :: results(:)
    real(wp) :: j1_min, j1_max, j_mid, jl, jl_0, jl_1, jr_0, jr_1, jr, sum, phase_sign, lambda
    integer :: Nj, i

    ! Determine the range of j1 
    j1_min = max(abs(j2-j3), abs(m1))
    j1_max = j2 + j3 
    j_mid = (j1_max - j1_min)/2.0_wp + j1_min  ! Central j value 
    Nj = int(j1_max-j1_min)+1  ! Number of j1 values 

    write(*,*) j1_max, j1_min, j_mid, Nj

    ! Allocate the results array 
    allocate(results(Nj))

    ! Initial values for reucrsion from left/right side
    jl_0 = 1._wp 
    jr_0 = 1._wp
    results(1) = jl_0
    results(Nj) = jr_0 

    ! First recursion from each side 
    jl_1 = -cal_B(j1_min, j2, j3, m1, m2, m3) / (cal_A(j1_min+1._wp, j2, j3, m1) * j1_min) 
    jr_1 = -cal_B(j1_max, j2, j3, m1, m2, m3) / (cal_A(j1_max, j2, j3, m1) * (j1_max+1._wp)) 
    results(2) = jl_1
    results(Nj-1) = jr_1

    ! Now do recursion to the elements next to the center (if we are not already there)
    if (Nj/2 > 2) then 
        !write(*,*) Nj/2 
        do i=2, Nj/2-1 
            jl = -(cal_B(j1_min+i-1._wp, j2, j3, m1, m2, m3)*jl_1 + (j1_min+i)*cal_A(j1_min+i-1._wp, j2, j3, m1)*jl_0) &
                / ((j1_min+i-1._wp)*cal_A(j1_min+i, j2, j3, m1)) 
            jr = - ((j1_max-i+1._wp)*cal_A(j1_max-i+2._wp, j2, j3, m1)*jr_0 + cal_B(j1_max-i+1._wp, j2, j3, m1, m2, m3)*jr_1) & 
                / ((j1_max-i+2._wp)*cal_A(j1_max-i+1._wp, j2, j3, m1)) 
            
            ! These indices are dumb, maybe do recursion individually for left/right...
            
            ! If not center, save and prepare next iteration.
            !if (i < Nj/2 + 1) then
            results(i+1) = jl
            results(Nj-i) = jr
            
            !write(*,*) jl, jr 
            !write(*,*) i+1, Nj-i

            jl_0 = jl_1
            jl_1 = jl
            jr_0 = jr_1
            jr_1 = jr
            !end if 
        end do
    end if 

    ! Perform the last two steps to the center outside loop - we go past with one to match in 3 points  
    !write(*,*) jl_0, jl_1, jl, jr_0, jr_1, jr
    lambda = cal_A(j1_min, j2, j3, m1)
    !write(*,*) lambda 

    jl = -(cal_B(j_mid-1._wp, j2, j3, m1, m2, m3)*jl_1 + j_mid*cal_A(j_mid-1._wp,j2,j3,m1)*jl_0) & 
        / ((j_mid-1._wp)*cal_A(j_mid, j2, j3, m1))
    jl_0 = jl_1 
    jl_1 = jl
    jl = -(cal_B(j_mid, j2, j3, m1, m2, m3)*jl_1 + (j_mid+1._wp)*cal_A(j_mid,j2,j3,m1)*jl_0) & 
        / ((j_mid)*cal_A(j_mid+1._wp, j2, j3, m1))


    jr = - ((j_mid+1._wp)*cal_A(j_mid+2._wp,j2,j3,m1)*jr_0 + cal_B(j_mid+1._wp,j2,j3,m1,m2,m3)*jr_1) & 
        / ((j_mid+2._wp)*cal_A(j_mid+1._wp, j2, j3, m1))
    jr_0 = jr_1
    jr_1 = jr
    jr = - ((j_mid)*cal_A(j_mid+1._wp,j2,j3,m1)*jr_0 + cal_B(j_mid,j2,j3,m1,m2,m3)*jr_1) & 
        / ((j_mid+1._wp)*cal_A(j_mid, j2, j3, m1))


    ! Save the mid value 
    results(Nj/2+1) = jr_1
    !write(*,*) results
    
    ! Using three central terms we do a least square to match (this way we avoid dividing by 0 if central point is 0)
    !lambda = (jl_0*jr + jl*jr_0 + jl_1*jr_1) / (jl_0**2 + jl**2 + jl_0**2)

    ! Scale the left part with lambda 
    !do i=1, Nj/2 
    !    results(i) = results(i) * 1._wp/lambda
    !end do

    if (jr_1 /= 0._wp) then 
        lambda = jr_1/jl_1 
    else 
        lambda = jr/jl_0 
    end if 

    !write(*,*) lambda 

    ! Scale the left part by lambda 
    do i=1, Nj/2 
        results(i) = results(i) * lambda 
    end do

    ! Now normalize the results 
    sum = 0._wp 
    do i=1, Nj 
        sum = sum + (2._wp*(j1_min + i -1._wp)+1._wp) * results(i)**2
    end do 
    results = results / sqrt(sum)

    ! Last determine the phase by looking at the last element
    phase_sign = sign(1._wp, results(Nj))
    if (phase_sign /= (-1._wp)**(j2-j3-m1)) then
        results = results * (-1._wp)
    end if
end subroutine wigner3j_j1


end module 
