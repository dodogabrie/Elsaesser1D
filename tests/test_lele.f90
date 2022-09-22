program main
  use physicals
  use numericals, only: xx, T_STEPS, n, dx
  use IO
	use lele_filter
  implicit none
  integer :: i, state
  real*8 :: mu_g, sigma_g
  call initialize()
  mu_g = L/2
  sigma_g =2* dx
  state = 0
  do i = 1, n
    f(i) = EXP( -( xx(i) - mu_g)*( xx(i) - mu_g)/(sigma_g * sigma_g) ) 
    g(i) = EXP( -( xx(i) - mu_g)*( xx(i) - mu_g)/(sigma_g * sigma_g) ) 
  end do
  call save_state("filt", state)
  call filter6(f)
  call save_state("filt", state+1)
end program
