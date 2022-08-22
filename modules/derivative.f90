module derivative
    real*8, dimension(:), allocatable ::  diagonal
contains
   subroutine tridiag(l, d, u, b, sol)
       ! Simple case: all elements of the diag, 
       ! sub and sup diag are equal
      use numericals, only: n
      integer :: i
      real*8, dimension(n-2) :: b, sol
      real*8 :: factor, l, d, u
      diagonal(1) = d
      do i = 2, n-2 ! Gauss reduction
        factor = l/diagonal(i-1)
        diagonal(i) = d - factor * u
        b(i) = b(i) - factor * b(i-1)
      enddo
      sol(n-2) = b(n-2)/diagonal(n-2)
      do i = n-2, 1, -1 ! Solve 
         sol(i) = (b(i) - u * sol(i+1))/diagonal(i)
      enddo
   end subroutine tridiag

   subroutine DFS1(f, der)
      use numericals
      use lele_filter, only: filter6
      ! Finite difference o(2) first derivative
      integer :: i
      real*8, dimension(n) :: f, der
      do i = 2, n-1
         der(i) = ( f(i+1) - f(i-1) )/(2 * dx)
      end do
      der(1) = ( - f(3)*0.5 + 2 * f(2) - 3 * f(1)/2 )/( dx )
      der(n) = ( 2 * f(n)/2 - 2 * f(n-1) + 0.5 * f(n-2))/( dx )
      call filter6(der)
   end subroutine DFS1

   subroutine DFS2(f, der)
      use numericals
      use lele_filter, only: filter6

      ! Finite difference o(2) second derivative 1D
      integer :: i
      real*8, dimension(n) :: f, der
      do i = 2, n-1
         der(i) = ( f(i+1) - 2 * f(i) + f(i-1) )/(dx * dx)
      end do
      der(1) = ( f(3) - 2 * f(2) + f(1) )/( dx * dx )
      der(n) = ( f(n) - 2 * f(n-1) + f(n-2))/( dx * dx )
      call filter6(der)
   end subroutine DFS2

   subroutine DFC1(f, der)
      use numericals
      ! Finite difference o(2) second derivative 1D
      integer :: i
      real*8, dimension(n) :: f, der
      !call pentadiag()
      do i = 2, n-1
         der(i) = ( f(i+1) - 2 * f(i) + f(i-1) )/(dx * dx)
      end do
      der(1) = ( f(2) - f(1) )/( dx )
      der(n) = ( f(n) - f(n-1))/( dx )
   end subroutine DFC1
end module derivative


