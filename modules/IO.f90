module IO
   implicit none
contains
   subroutine save_state(num)
      use physicals, only: f, g
      use numericals, only: n, xx
      integer :: num, i
      character(8) :: fmt
      character(5) :: x1
      fmt = '(I5.5)'
      write(x1, fmt) num
      open(1,file='data/output'//trim(x1)//'.dat',status='unknown')
      do i = 1, n
         write(1, *) xx(i), f(i), g(i)
      end do
      close(1)
   end subroutine save_state
   subroutine initialize()
      use numericals
      use physicals
      use derivative
      use lele_filter
      implicit none
      PI=4.D0*DATAN(1.D0)
      open(1,file='input.dat',status='old')
      read(1, *) n
      read(1, *) dx
      read(1, *) dt
      read(1, *) T_STEPS
      read(1, *) c
      read(1, *) nu
      L = dx * n
      allocate ( f(n) , g(n) )
      allocate ( aux_f(n) , aux1_f(n) , aux_g(n) , aux1_g(n) )
      allocate ( der1_f(n) , der1_g(n) ,der2_f(n) )
      allocate ( xx(n) )
      allocate ( alpha(n-1), beta(n-2), z(n))
      allocate ( y_aux(n-6))
      allocate ( diagonal(n-2))
      call linspace()
   end subroutine initialize
end module IO

