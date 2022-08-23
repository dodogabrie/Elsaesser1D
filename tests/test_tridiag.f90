! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
subroutine inv(A, Ainv, n)
   real*8, dimension(n,n), intent(in) :: A
   real*8, dimension(n,n) :: Ainv

   real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
   integer, dimension(size(A,1)) :: ipiv   ! pivot indices
   integer :: nn, info

   ! External procedures defined in LAPACK
   external DGETRF
   external DGETRI

   ! Store A in Ainv to prevent it from being overwritten by LAPACK
   Ainv = A
   nn = size(A,1)

   ! DGETRF computes an LU factorization of a general M-by-N matrix A
   ! using partial pivoting with row interchanges.
   call DGETRF(nn, n, Ainv, nn, ipiv, info)

   if (info /= 0) then
      stop 'Matrix is numerically singular!'
   end if

   ! DGETRI computes the inverse of a matrix using the LU factorization
   ! computed by DGETRF.
   call DGETRI(n, Ainv, nn, ipiv, work, nn, info)

   if (info /= 0) then
      stop 'Matrix inversion failed!'
   end if
end subroutine inv

program main
	use derivative, only: tridiag
   use IO, only: initialize
   implicit none
   real*8, dimension(:), allocatable :: a_diag, sol_p, sol_i, y, test_p, test_i
   real*8, dimension(:), allocatable :: a_lo, a_up
   real*8, dimension(:,:), allocatable :: A, Ainv
   real*8 :: e, c, d, aa, b
   integer :: i, j, n
   real :: start, finish
   n = 50
   e = 5.d0
   c = 4.d0
   d = 1.d0
   aa = 2.d0
   b = 3.d0
   call initialize()
   allocate(a_diag(n),sol_p(n), sol_i(n), y(n), A(n,n), Ainv(n,n))
   allocate(test_i(n),test_p(n))
   allocate(a_lo(n-1), a_up(n-1))
   do i = 1, n
      a_diag(i) = d
      if (i < n) then
         a_lo(i) = c
         a_up(i) = aa
      endif
   enddo

   !call random_number(a_diag)
   !call random_number(a_lo)
   !call random_number(a_lolo)
   !call random_number(a_up)
   !call random_number(a_upup)
   call random_number(y)
   do i = 1, n
      do j = 1, n
         A(i, j) = 0.d0
      enddo
   enddo
   do i=1, n
      A(i, i) = a_diag(i)
      if (i > 1) then
         A(i, i-1) = a_lo(i-1)
      endif
      if (i < n) then
         A(i, i + 1) = a_up(i)
      endif
   enddo

   call cpu_time(start)
   call tridiag(c, d, aa, y, sol_p, n)
   call cpu_time(finish)
   print '("Linear tridiagonal Time = ",f6.3," seconds.")',finish-start

   call cpu_time(start)
   call inv(A, Ainv, n)
   sol_i = matmul(Ainv, y)
   call cpu_time(finish)
   print '("Lapack inverse tridiagonal Time = ",f6.3," seconds.")',finish-start

   !write(*,*) "penta        inv "
   !do i = 1, n
   !   write(*,*) sol_p(i), sol_i(i), sol_p(i)/sol_i(i)
   !enddo
   !write(*,*) "Test reverse"
   test_i = matmul(A, sol_i)
   test_p = matmul(A, sol_p)
   !write(*,*) "init    penta        inv "

   do i = 1, n
      if (abs(y(i) - test_p(i)) > 0.000001) then
         write(*,*) i, y(i), test_p(i), test_i(i)
      endif
   enddo
end program
