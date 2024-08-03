module LAPACK_routines
use utils 
implicit none 

private 

public complex_inverse 

contains 

subroutine complex_inverse(mat, N)
    integer, intent(in) :: N  ! Dimension of the square matrix 
    complex(dp), intent(inout) :: mat(N,N)
    complex(dp) :: work(N)
    integer :: IPIV(N)  ! Integer array to contain the pivot indices? 
    integer :: info 
 
    ! LU factorization of the matrix 
    call ZGETRF(N, N, mat, N, IPIV, info)
    if (info .eq. 1) then 
        write(*,*) 'Failed ZGETRF!'
    end if 

    ! Result of factorization is fed into ZGETRI to get the inverse 
    call ZGETRI(N, mat, N, IPIV, WORK, N, info)
    if (info .eq. 1) then 
        write(*,*) 'Failed ZGETRI!'
    end if
end subroutine complex_inverse 

end module LAPACK_routines