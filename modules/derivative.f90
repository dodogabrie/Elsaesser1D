module derivative
    real*8, dimension(:), allocatable ::  diagonal, b_mod, yy
contains
   subroutine tridiag(l, d, u, b, sol, m)
       ! Simple case: all elements of the diag, 
       ! sub and sup diag are equal
      integer :: i,m
      real*8, dimension(m) :: b, sol
      real*8 :: factor, l, d, u
      diagonal(1) = d
			b_mod(1) = b(1)
      do i = 2, m ! Gauss reduction
        factor = l/diagonal(i-1)
        diagonal(i) = d - factor * u
        b_mod(i) = b(i) - factor * b_mod(i-1)
      enddo
      sol(m) = b_mod(m)/diagonal(m)
      do i = m-1, 1, -1 ! Solve 
         sol(i) = (b_mod(i) - u * sol(i+1))/diagonal(i)
      enddo
   end subroutine tridiag

   subroutine tridiag_boundary(l, d, u, b, sol, m, bound)
       ! Simple case: all elements of the diag, 
       ! sub and sup diag are equal
      integer :: i,m
      real*8, dimension(m) :: b, sol
      real*8 :: factor, l, d, u, bound
      diagonal(1) = d
			b_mod(1) = b(1)
      factor = l/d
			diagonal(2) = d - factor * bound
			b_mod(2) = b(2) - factor * b_mod(1)
      do i = 3, m-1 ! Gauss reduction
        factor = l/diagonal(i-1)
        diagonal(i) = d - factor * u
        b_mod(i) = b(i) - factor * b_mod(i-1)
      enddo
			factor = bound / diagonal(m-1)
			diagonal(m) = d - factor * u
			b_mod(m) = b(m) - factor * b_mod(m-1)

      sol(m) = b_mod(m)/diagonal(m)
      do i = m-1, 2, -1 ! Solve 
         sol(i) = (b_mod(i) - u * sol(i+1))/diagonal(i)
      enddo
      sol(1) = (b_mod(1) - bound * sol(2))/diagonal(1)
   end subroutine tridiag_boundary


   subroutine DFS1(f, der)
      use numericals
      ! Finite difference o(2) first derivative
      integer :: i
      real*8, dimension(n) :: f, der
      do i = 2, n-1
         der(i) = ( f(i+1) - f(i-1) )/(2 * dx)
      end do
      der(1) = ( - f(3)*0.5 + 2 * f(2) - 3 * f(1)/2 )/( dx )
      der(n) = ( 3 * f(n)/2 - 2 * f(n-1) + 0.5 * f(n-2))/( dx )
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
   end subroutine DFS2

   subroutine DFC1(f, der)
      use numericals, only: n, dx
      ! Finite difference o(2) second derivative 1D
      integer :: i
      real*8, dimension(n) :: der, f
			real*8 :: lu, d, bound, inv_dx
			inv_dx = 1/dx
			bound = 3.d0
			lu = 1.d0/4.d0
			d = 1.d0
			yy(1) = inv_dx * (-17.d0/6.d0 * f(1) + 3.d0/2.d0 * f(2) + 3d0/2d0*f(3) - 1d0/6d0*f(4))
			yy(n) = - inv_dx * (-17.d0/6.d0 * f(n) + 3.d0/2.d0 * f(n-1) + 3d0/2d0*f(n-2) - 1d0/6d0*f(n-3))
			do i = 2, n-1
				yy(i) = inv_dx * 3.d0/(4.d0) * (f(i+1) - f(i-1))
			enddo
			
      call tridiag_boundary(lu, d, lu, yy, yy, n, bound)
      do i = 1, n
         der(i) = yy(i)
      end do
   end subroutine DFC1

   subroutine DFC2(f, der)
      use numericals, only: n, dx
      ! Finite difference o(2) second derivative 1D
      integer :: i
      real*8, dimension(n) :: der, f
			real*8 :: lu, d, inv_dx
			inv_dx = 1/dx
			lu = 1.d0
			d = 10.d0
			yy(1) = 1/(dx * dx) * (13.d0 * f(1) - 27.d0 * f(2) + 15.d0 *f(3) - f(4))
			yy(n) = 1/(dx * dx) * (13.d0 * f(n) - 27.d0 * f(n-1) + 15.d0 *f(n-2) - f(n-3))
			do i = 2, n-1
				yy(i) = 12.d0 * inv_dx * inv_dx  * (f(i+1) - 2 * f(i) + f(i-1))
			enddo
			
      call tridiag_boundary(lu, d, lu, yy, yy, n, 11.d0)
      do i = 1, n
         der(i) = yy(i)
      end do
   end subroutine DFC2

end module derivative


