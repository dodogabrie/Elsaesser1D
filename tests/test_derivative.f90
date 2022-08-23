program test_DFC1
	use physicals, only: f, g, PI, L
	use numericals
	use derivative, only: DFC1, DFC2
	use IO
	integer :: i
	call initialize()

	do i = 1, n 
		f(i) = -sin(4 * PI/L * xx(i))
	end do
	call DFC2(f, g)
	call save_state(1)
end program test_DFC1
