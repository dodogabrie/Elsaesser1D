program main
  use physicals
  use numericals, only: xx, T_STEPS, n, dx
  use IO
	use lele_filter
  implicit none
  integer :: i, state
  real*8 :: mu_g, sigma_g
  call initialize()
  mu_g = L
  sigma_g = 100 * dx
  state = 0
  do i = 1, n
    f(i) = EXP( -( xx(i) - mu_g)*( xx(i) - mu_g)/(sigma_g * sigma_g) ) 
    f(i) = f(i) + EXP( -( xx(i))*( xx(i))/(sigma_g * sigma_g) ) 
    g(i) = EXP( -( xx(i) - mu_g)*( xx(i) - mu_g)/(sigma_g * sigma_g) ) 
  end do
  call save_state(state)
  call filter6(f)
  call save_state(state+1)
end program
