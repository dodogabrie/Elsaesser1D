! Code by Edoardo Gabrielli

subroutine RHS()
  ! Righten side of linear Alfèn waves
  use physicals
  use numericals
  use derivative
	use lele_filter, only: filter6
  integer :: i
  real*8 :: lambda_f, lambda_g, tmp_f
  call DFC1(aux_f, der1_f)
	call filter6(der1_f)
  call DFC1(aux_g, der1_g)
	call filter6(der1_g)
  do i = 1, n
    aux1_f(i) = aux_f(i) - aux_g(i)
  end do
  call DFC2(aux1_f, der2_f)
	call filter6(der2_f)
  do i = 2, n-1
    tmp_f = aux_f(i)
    aux_f(i) = c * der1_f(i) + aux_g(i) * der1_f(i) + nu * der2_f(i)
    aux_g(i) = - c * der1_g(i) + tmp_f * der1_g(i) - nu * der2_f(i)
  end do
  tmp_f = aux_f(1)
  lambda_f = c + aux_g(1)
  aux_f(1) = ( 1 + SIGN(1.d0, lambda_f) ) * lambda_f * 0.5d0 * der1_f(1)
  lambda_g = - c + tmp_f
  aux_g(1) = ( 1 + SIGN(1.d0, lambda_g) ) * lambda_g * 0.5d0 * der1_g(1)

  tmp_f = aux_f(n)
  lambda_f = c + aux_g(n)
  aux_f(n) = ( 1 - SIGN(1.d0, lambda_f) ) * lambda_f * 0.5d0 * der1_f(n)
  lambda_g = - c + tmp_f
  aux_g(n) = ( 1 - SIGN(1.d0, lambda_g) ) * lambda_g * 0.5d0 * der1_g(n)
end subroutine RHS

subroutine RK(order)
! Runge Kutta of order N for time indipendent PDE
  use numericals
  use physicals
  integer :: order, k
  do i = 1, n
    aux_f(i) = f(i)
    aux_g(i) = g(i)
  end do
  do k = order, 1, -1
    call RHS()
    do i = 1, n
  aux_f(i) = f(i) + dt/k * aux_f(i)
  aux_g(i) = g(i) + dt/k * aux_g(i)
    enddo
  enddo
  do i = 1, n
    f(i) = aux_f(i)
    g(i) = aux_g(i)
  end do
end subroutine RK

program main
  use physicals
  use numericals, only: xx, T_STEPS, n
  use IO
  implicit none
  integer :: order, i, state, kk
  real*8 :: mu_f, mu_g, sigma
  order = 4
  call initialize()
  mu_f = L/3
  mu_g = 2 * L/3
  sigma = L/10
  kk = 2
  state = 0
  do i = 1, n
    f(i) = 1.3 * EXP( -( xx(i) - mu_f)*( xx(i) - mu_f)/(sigma * sigma) ) * SIN( kk *  2 * PI / (4 * sigma) * (xx(i)-mu_f))
    g(i) = EXP( -( xx(i) - mu_g)*( xx(i) - mu_g)/(sigma * sigma) )! * SIN( kk *  2 * PI / L * xx(i))
  end do
  call save_state(state)
  do i = 1, T_STEPS
    call RK(order)
    if (mod(i, 50) == 0) then
      write(*,fmt="(I10)",advance='NO') i
      state = state + 1
      call save_state(state)
    endif
  enddo
  call save_state(state)
end program
