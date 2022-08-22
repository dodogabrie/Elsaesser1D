module numericals
   implicit none
   real*8 :: dx, dt
   integer :: n, T_steps
   real*8, dimension(:), allocatable :: der1_f, der2_f, der1_g
   real*8, dimension(:), allocatable :: aux_f, aux_g, aux1_f, aux1_g
   real*8, dimension(:), allocatable :: xx
contains
   subroutine linspace()
      integer :: i
      do i = 0, n-1
         xx(i+1) = dx * i
      end do
   end subroutine linspace
end module numericals


