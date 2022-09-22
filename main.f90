! Code by Edoardo Gabrielli

subroutine RHS_l()
  ! Righten side of linear Alfèn waves
  use physicals
  use numericals
  use derivative
	use lele_filter, only: filter6
	implicit none
  integer :: i
  real*8 :: lambda_f, lambda_g, tmp_f
  call DFC1(aux_f, der1_f)
  call DFC1(aux_g, der1_g)
  do i = 1, n
    aux1_f(i) = aux_f(i) - aux_g(i)
  end do
  call DFC2(aux1_f, der2_f)
  do i = 2, n-1
    tmp_f = aux_f(i)
    aux_f(i) = c * der1_f(i) + nu * der2_f(i)
    aux_g(i) = - c * der1_g(i) - nu * der2_f(i)
  end do
  tmp_f = aux_f(1)
  lambda_f = c
  aux_f(1) = ( 1 + SIGN(1.d0, lambda_f) ) * lambda_f * 0.5d0 * der1_f(1) + nu * der2_f(1)
  lambda_g = - c 
  aux_g(1) = ( 1 + SIGN(1.d0, lambda_g) ) * lambda_g * 0.5d0 * der1_g(1)- nu * der2_f(1)

  tmp_f = aux_f(n)
  lambda_f = c
  aux_f(n) = ( 1 - SIGN(1.d0, lambda_f) ) * lambda_f * 0.5d0 * der1_f(n)+ nu * der2_f(n)
  lambda_g = - c
  aux_g(n) = ( 1 - SIGN(1.d0, lambda_g) ) * lambda_g * 0.5d0 * der1_g(n)- nu * der2_f(n)
	call filter6(aux_g)
	call filter6(aux_f)
end subroutine RHS_l


subroutine RHS_nl()
  ! Righten side of linear Alfèn waves
  use physicals
  use numericals
  use derivative
	use lele_filter, only: filter6
	implicit none
  integer :: i
  real*8 :: lambda_f, lambda_g, tmp_f
  call DFC1(aux_f, der1_f)
  call DFC1(aux_g, der1_g)
  do i = 1, n
    aux1_f(i) = aux_f(i) - aux_g(i)
  end do
  call DFC2(aux1_f, der2_f)
  do i = 2, n-1
    tmp_f = aux_f(i)
    aux_f(i) = c * der1_f(i) + aux_g(i) * der1_f(i) + nu * der2_f(i)
    aux_g(i) = - c * der1_g(i) + tmp_f * der1_g(i) - nu * der2_f(i)
  end do
  tmp_f = aux_f(1)
  lambda_f = c + aux_g(1)
  aux_f(1) = ( 1 + SIGN(1.d0, lambda_f) ) * lambda_f * 0.5d0 * der1_f(1) + nu * der2_f(1)
  lambda_g = - c + tmp_f
  aux_g(1) = ( 1 + SIGN(1.d0, lambda_g) ) * lambda_g * 0.5d0 * der1_g(1) - nu * der2_f(1)

  tmp_f = aux_f(n)
  lambda_f = c + aux_g(n)
  aux_f(n) = ( 1 - SIGN(1.d0, lambda_f) ) * lambda_f * 0.5d0 * der1_f(n) + nu * der2_f(n)
  lambda_g = - c + tmp_f
  aux_g(n) = ( 1 - SIGN(1.d0, lambda_g) ) * lambda_g * 0.5d0 * der1_g(n) - nu * der2_f(n)
	call filter6(aux_g)
	call filter6(aux_f)
end subroutine RHS_nl

subroutine RK_l(order)
! Runge Kutta of order N for time indipendent PDE
  use numericals
  use physicals
	implicit none
  integer :: order, k, i
  do i = 1, n
    aux_f(i) = f(i)
    aux_g(i) = g(i)
  end do
  do k = order, 1, -1
    call RHS_l()
    do i = 1, n
      aux_f(i) = f(i) + dt/k * aux_f(i)
      aux_g(i) = g(i) + dt/k * aux_g(i)
    enddo
  enddo
  do i = 1, n
    f(i) = aux_f(i)
    g(i) = aux_g(i)
  end do
end subroutine RK_l

subroutine RK_nl(order)
! Runge Kutta of order N for time indipendent PDE
  use numericals
  use physicals
	implicit none
  integer :: order, k, i
  do i = 1, n
    aux_f(i) = f(i)
    aux_g(i) = g(i)
  end do
  do k = order, 1, -1
    call RHS_nl()
    do i = 1, n
      aux_f(i) = f(i) + dt/k * aux_f(i)
      aux_g(i) = g(i) + dt/k * aux_g(i)
    enddo
  enddo
  do i = 1, n
    f(i) = aux_f(i)
    g(i) = aux_g(i)
  end do
end subroutine RK_nl

program main
  use physicals
  use numericals, only: xx, T_STEPS, n
  use IO
  implicit none
  integer :: order, i, state, kk
  real*8 :: mu_f, mu_g, sigma, Af, Ag, init_f, init_g
	character(4) :: nlnam, lnam
	nlnam = "nlin"
 	lnam = "line"
  order = 4
  call initialize()
  mu_f = L/5
  mu_g = 4* L/5
  sigma = L/10
	Af = 1.5d0
	Ag = 1.2d0
  kk = 1
  state = 0
!------- Non linear evolution --------------------
	write(*,*) "Starting evolution with nonlinear terms"
  do i = 1, n
    f(i) = init_f(xx(i), Af, mu_f, sigma, kk, PI)
    g(i) = init_g(xx(i), Ag, mu_g, sigma, kk, PI)
  end do
  call save_state(nlnam, state)
  do i = 1, T_STEPS
    call RK_nl(order)
    if (mod(i, 200) == 0) then
      write(*,fmt="(I10)",advance='NO') i
      state = state + 1
      call save_state(nlnam, state)
    endif
  enddo
  call save_state(nlnam, state+1)
!-------------------------------------------------
	write(*,*) ""
	write(*,*) "Starting evolution with only linear terms"
!----------- Linear evolution --------------------
  do i = 1, n
    f(i) = init_f(xx(i), Af, mu_f, sigma, kk, PI)
    g(i) = init_g(xx(i), Ag, mu_g, sigma, kk, PI)
  end do

	state = 0
  call save_state(lnam, state)
  do i = 1, T_STEPS
    call RK_l(order)
    if (mod(i, 200) == 0) then
      write(*,fmt="(I10)",advance='NO') i
      state = state + 1
      call save_state(lnam, state)
    endif
  enddo
  call save_state(lnam, state+1)
!----------------------------------------------
end program

real*8 function init_g(x, a, mu, sigma, kk, pi)
	real*8 :: x, mu, sigma, a, pi
	integer :: kk
  init_g = a * EXP( -( x - mu)*( x - mu)/(sigma * sigma) ) !* sin( kk *  2 * pi / (2 * sigma) * (x-mu))
	return
end function init_g

real*8 function init_f(x, a, mu, sigma, kk, pi)
	real*8 :: x, mu, sigma, a, pi
	integer :: kk
  init_f = a * exp( -( x - mu)*( x - mu)/(sigma * sigma) ) * sin( kk *  2 * pi / (2 * sigma) * (x-mu))
	return
end function init_f
