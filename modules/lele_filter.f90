module lele_filter
   real*8, dimension(:), allocatable :: alpha, beta, z
   real*8, dimension(:), allocatable :: y_aux
   real*8 :: gamm, mu
   ! Module with different ways to compute derivative 1D
contains!               u_ll, u_l, u_d, u_p, u_pp  MX = Y
   subroutine pentadiag(e,     c,   d,   a,    b,    y,  sol)
      !                       gamm  mu  alp  beta
      use numericals
      real*8, dimension(n) :: d, sol, y
      real*8, dimension(n-1) :: c, a
      real*8, dimension(n-2) :: e, b
      integer :: i
      mu = d(1)
      alpha(1) = a(1)/mu
      beta(1) = b(1)/mu
      z(1) = y(1)/mu

      gamm = c(1)
      mu = d(2) - alpha(1) * gamm
      alpha(2) = (a(2)- beta(1) * gamm )/mu
      beta(2) = b(2)/mu
      z(2) = (y(2) - z(1) * gamm )/mu

      do i = 3, n-2
         gamm = c(i-1) - alpha(i-2) * e(i-2)
         mu = d(i) - beta(i-2) * e(i-2) - alpha(i-1) * gamm
         alpha(i) = (a(i) - beta(i-1) * gamm)/mu
         beta(i) = b(i)/mu
         z(i) = (y(i) - z(i-2) * e(i-2) - z(i-1) * gamm)/mu
      enddo
      gamm = c(n-2) - alpha(n-3) * e(n-3)
      mu = d(n-1) - beta(n-3) * e(n-3) - alpha(n-2) * gamm
      alpha(n-1) = (a(n-1) - beta(n-2) * gamm)/mu
      z(n-1) = (y(n-1) - z(n-3) * e(n-3) - z(n-2) * gamm)/mu

      gamm = c(n - 1) - alpha(n-2) * e(n-2)
      mu = d(n) - beta(n-2) * e(n-2) - alpha(n-1) * gamm
      z(n) = (y(n) - z(n-2) * e(n-2) - z(n-1) * gamm)/mu

      sol(n) = z(n)
      sol(n-1) = z(n-1) - alpha(n-1) * sol(n)
      do i = n - 2, 1, -1
         sol(i) = z(i) - alpha(i) * sol( i+1 ) - beta(i) * sol(i+2)
      end do
   end subroutine pentadiag

   subroutine pentadiag_simple(e,     c,   d,   a,    b,    y,  sol, m)
      !                            gamm  mu  alp  beta
      use numericals
      real*8, dimension(m) :: sol, y
      real*8 :: d, c, a, e, b
      integer :: i, m
      mu = d
      alpha(1) = a/mu
      beta(1) = b/mu
      z(1) = y(1)/mu

      gamm = c
      mu = d - alpha(1) * gamm
      alpha(2) = (a- beta(1) * gamm )/mu
      beta(2) = b/mu
      z(2) = (y(2) - z(1) * gamm )/mu

      do i = 3, m-2
         gamm = c - alpha(i-2) * e
         mu = d - beta(i-2) * e - alpha(i-1) * gamm
         alpha(i) = (a - beta(i-1) * gamm)/mu
         beta(i) = b/mu
         z(i) = (y(i) - z(i-2) * e - z(i-1) * gamm)/mu
      enddo
      gamm = c - alpha(m-3) * e
      mu = d - beta(m-3) * e - alpha(m-2) * gamm
      alpha(m-1) = (a - beta(m-2) * gamm)/mu
      z(m-1) = (y(m-1) - z(m-3) * e - z(m-2) * gamm)/mu

      gamm = c - alpha(m-2) * e
      mu = d - beta(m-2) * e - alpha(m-1) * gamm
      z(m) = (y(m) - z(m-2) * e - z(m-1) * gamm)/mu

      sol(m) = z(m)
      sol(m-1) = z(m-1) - alpha(m-1) * sol(m)
      do i = m - 2, 1, -1
         sol(i) = z(i) - alpha(i) * sol( i+1 ) - beta(i) * sol(i+2)
      end do
   end subroutine pentadiag_simple

   subroutine filter6(f)
      use numericals, only: n
      real*8, dimension(n) :: f
      real*8 :: a, b, c, d, alp, bet
      real*8 :: diag, lo, lolo, up, upup, f1, f2
			integer :: i, ii
      alp = 0.d0
      bet = 3/10.d0

      a = 0.5
      b = 0.75
      c = 3/10.d0
			d = 1.d0/20.d0
      diag = 1.d0
      lo = alp
      lolo = bet
      up = lo
      upup = lolo
      do i = 1, n-6
				ii = i+3
        y_aux(i) = a * f(ii) + b * 0.5 * (f(ii+1) + f(ii-1)) + c * 0.5 * (f(ii+2) + f(ii-2)) + d * 0.5 * (f(ii+3) + f(ii-3))
      enddo
      f1 = f(1)
      f2 = f(2)
      f(1) = 15.d0/16 * f(1) + 1.d0/16 * (4*f(2) - 6 * f(3) + 4 * f(4) - f(5))
      f(2) = 3.d0/4   * f(2) + 1.d0/16 * (  f1   + 6 * f(3) - 4 * f(4) + f(5))
      f(3) = 5.d0/8   * f(3) + 1.d0/16 * (- f1   + 4 * f2   + 4 * f(4) - f(5))
      f1 = f(n)
      f2 = f(n-1)
      f(n)   = 15.d0/16 * f(n)    + 1.d0/16*(4*f(n-1)- 6 * f(n-2) + 4 * f(n-3) - f(n-4))
      f(n-1) = 3.d0/4   * f(n-1)  + 1.d0/16*(  f1    + 6 * f(n-2) - 4 * f(n-3) + f(n-4))
      f(n-2) = 5.d0/8   * f(n-2)  + 1.d0/16*(- f1    + 4 * f2     + 4 * f(n-3) - f(n-4))
			i = 1
			ii = i + 3
			y_aux(i) = y_aux(i) - alp * f(ii - 1) - bet * f(ii - 2)
			i = 2
			ii = i + 3
			y_aux(i) = y_aux(i) - bet * f(ii - 2) 
			i = n - 6
			ii = i + 3
			y_aux(i) = y_aux(i) - alp * f(ii + 1) - bet * f(ii + 2)
			i = n - 7
			ii = i + 3
			y_aux(i) = y_aux(i) - bet * f(ii + 2)
      call pentadiag_simple(lolo, lo, diag, up, upup, y_aux, y_aux, n-6)
      do i = 1, n-6
         f(i + 3) = y_aux(i)
      enddo
   end subroutine filter6
end module lele_filter


